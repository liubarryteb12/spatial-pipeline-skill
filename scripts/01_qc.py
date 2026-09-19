#!/usr/bin/env python3
"""
01_qc.py — spot 质控与过滤

**空间数据的 QC 与单细胞有一处本质不同：不能按"细胞"过滤。**
每个 spot 是 1-10 个细胞的混合物，所以：
  - `pct_counts_mt` 的分布比单细胞**窄得多**（多个细胞的线粒体信号平均掉了）
  - 单个 spot 的 UMI 少（Visium 每 spot 约 1k-10k，而 10x 单细胞是 5k-50k）
  - **过滤掉一个 spot 会在地图上留一个洞** —— 这个洞本身会影响
    空间邻域分析（邻居图会变），所以过滤阈值要谨慎

阈值是**组织特异**的：PBMC 5-10 / 淋巴结 20 / 心肌 30-50 / 肝 30+。
用 PBMC 的阈值套到心肌上会滤掉大部分 spot。

**过滤后必须检查组织是否被割裂。** 如果滤掉的 spot 集中在组织某一区域，
下游的"空间域"会把那个区域整体判成边界。这里用连通性检查发现它。
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))

import numpy as np  # noqa: E402
import scanpy as sc  # noqa: E402

from common import (df_to_records, ensure_dirs, load_config, log_info,  # noqa: E402
                    log_warn, parse_args, record_step, save_fig, set_seed,
                    spot_radius_plot_units, write_json, spatial_xy, W_DOUBLE, mm,)

HB_PREFIXES = ("HBA", "HBB", "HBD", "HBE", "HBG", "HBM", "HBQ", "HBZ")


def check_connectivity(adata, n_neighbors: int = 6) -> dict:
    """
    过滤后组织是否被割裂成多块。

    **这是空间数据特有的检查。** 单细胞数据里"滤掉一批细胞"没有几何后果，
    但空间数据里会留洞；洞足够大时组织被切成互不相邻的几块，
    而空间邻域分析、niche、空间通讯全部只在块内有效。

    用六边形邻居图（Visium 的 6 个直接邻居）的连通分量数判断。
    """
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree

    xy = spatial_xy(adata)
    # 用图像坐标建 kNN 图（Visium 网格上 k=6 就是六边形的 6 个邻居）
    k = min(n_neighbors + 1, len(xy))
    tree = cKDTree(xy)
    d, idx = tree.query(xy, k=k)
    rows = np.repeat(np.arange(len(xy)), k - 1)
    cols = idx[:, 1:].ravel()
    # 距离阈值：超过中位最近邻距离的 3 倍不算邻居（防止把割裂的两块连起来）
    med = float(np.median(d[:, 1]))
    mask = d[:, 1:].ravel() <= med * 3.0
    adj = sp.coo_matrix((np.ones(mask.sum()), (rows[mask], cols[mask])),
                        shape=(len(xy), len(xy)))
    n_comp, labels = connected_components(adj, directed=False)
    sizes = np.bincount(labels)
    sizes = np.sort(sizes)[::-1]
    return {
        "n_components": int(n_comp),
        "largest_component_size": int(sizes[0]),
        "largest_component_frac": round(float(sizes[0] / len(xy)), 4),
        "component_sizes_top5": [int(s) for s in sizes[:5]],
        "fragmented": bool(sizes[0] / len(xy) < 0.95),
    }


def run_01_qc(cfg: dict) -> dict:
    ensure_dirs(cfg)
    set_seed(cfg)
    data_dir = Path(cfg["output"]["data_dir"])
    res_dir = Path(cfg["output"]["results_dir"])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    adata = sc.read_h5ad(data_dir / "raw.h5ad")
    info = None
    try:
        import json
        with open(data_dir / "dataset_info.json", encoding="utf-8") as fh:
            info = json.load(fh)
    except Exception:  # noqa: BLE001
        pass
    organism = (info or {}).get("organism") or "Homo sapiens"
    lib = list(adata.uns["spatial"].keys())[0]
    n0 = adata.n_obs
    log_info(f"读入 {n0} spot x {adata.n_vars} 基因（{organism}）")

    # ---- 1. 指标 ------------------------------------------------------------
    # **小鼠的线粒体前缀是 mt-，人是 MT-。** 用错前缀会让 pct_counts_mt 全 0，
    # 过滤形同虚设，而日志里看不出任何异常。
    is_mouse = organism.lower().startswith("mus")
    mt_pref = ("mt-",) if is_mouse else ("MT-",)
    adata.var["mt"] = adata.var_names.str.startswith(mt_pref)
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl") if is_mouse
                                                       else ("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.startswith(HB_PREFIXES)
    n_mt = int(adata.var["mt"].sum())
    if n_mt == 0:
        raise RuntimeError(
            f"线粒体基因前缀 '{mt_pref[0]}' 一个都没匹配到 —— "
            f"物种判断可能错了（organism='{organism}'）。"
            f"前缀用错会让 pct_counts_mt 全 0，过滤形同虚设而日志无异常。")
    qc_vars = [v for v in ("mt", "ribo", "hb") if bool(adata.var[v].any())]
    sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, percent_top=None,
                               log1p=False, inplace=True)
    log_info(f"QC 指标已算（{', '.join(qc_vars)}）；线粒体基因 {n_mt} 个")

    # ---- 2. 过滤前的空间图 --------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(W_DOUBLE, mm(58)))
    xy = spatial_xy(adata)
    for ax, key, cmap in zip(axes,
                             ("total_counts", "n_genes_by_counts", "pct_counts_mt"),
                             ("viridis", "viridis", "magma")):
        s = ax.scatter(xy[:, 0], xy[:, 1], c=adata.obs[key].astype(float),
                       s=4, cmap=cmap)
        ax.set_title(key)
        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(s, ax=ax, shrink=0.8)
    fig.suptitle(f"QC metrics on tissue (n={n0})")
    save_fig(cfg, "qc_metrics_on_tissue", fig)

    # ---- 3. 过滤 ------------------------------------------------------------
    q = cfg["qc"]
    steps = []
    n_cur = adata.n_obs

    def _apply(mask, label):
        nonlocal adata, n_cur
        before = adata.n_obs
        adata = adata[mask].copy()
        steps.append({"filter": label, "before": int(before),
                      "after": int(adata.n_obs),
                      "removed": int(before - adata.n_obs)})
        n_cur = adata.n_obs
        log_info(f"  {label}: {before} -> {adata.n_obs}")

    _apply(adata.obs["total_counts"] >= float(q["min_counts"]),
           f"min_counts>={q['min_counts']}")
    if q.get("max_counts"):
        _apply(adata.obs["total_counts"] <= float(q["max_counts"]),
               f"max_counts<={q['max_counts']}")
    _apply(adata.obs["n_genes_by_counts"] >= float(q["min_genes"]),
           f"min_genes>={q['min_genes']}")
    if q.get("max_genes"):
        _apply(adata.obs["n_genes_by_counts"] <= float(q["max_genes"]),
               f"max_genes<={q['max_genes']}")
    if q.get("max_pct_mt") is not None:
        _apply(adata.obs["pct_counts_mt"] < float(q["max_pct_mt"]),
               f"max_pct_mt<{q['max_pct_mt']}")

    if adata.n_obs < int(q.get("min_spots", 100)):
        raise RuntimeError(
            f"过滤后只剩 {adata.n_obs} 个 spot < {q.get('min_spots', 100)}。\n"
            f"  阈值是**组织特异**的：PBMC 5-10 / 淋巴结 20 / 心肌 30-50 / 肝 30+。\n"
            f"  用 PBMC 的 max_pct_mt 套到心肌上会滤掉大部分 spot。\n"
            f"  过滤链: " + "; ".join(f"{s['filter']} -{s['removed']}" for s in steps))

    # ---- 4. 过滤后的连通性 --------------------------------------------------
    conn = check_connectivity(adata)
    log_info(f"组织连通性: {conn['n_components']} 块，最大块占 "
             f"{conn['largest_component_frac']:.3f}")
    if conn["fragmented"]:
        log_warn("**过滤把组织割裂了** —— 空间邻域/niche/通讯只在块内有效。"
                 "检查是不是阈值过严，或组织本身就有分离的区域")

    # ---- 5. 过滤后的空间图 --------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(W_DOUBLE, mm(72)))
    xy = spatial_xy(adata)
    s0 = axes[0].scatter(xy[:, 0], xy[:, 1], c=adata.obs["total_counts"].astype(float),
                         s=5, cmap="viridis")
    axes[0].set_title(f"After filtering (n={adata.n_obs})")
    fig.colorbar(s0, ax=axes[0], shrink=0.8)
    s1 = axes[1].scatter(xy[:, 0], xy[:, 1], c=adata.obs["pct_counts_mt"].astype(float),
                         s=5, cmap="magma")
    axes[1].set_title("pct_counts_mt")
    fig.colorbar(s1, ax=axes[1], shrink=0.8)
    for ax in axes:
        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
    save_fig(cfg, "qc_on_tissue_after", fig)

    # ---- 6. 落盘 ------------------------------------------------------------
    out = data_dir / "qc_filtered.h5ad"
    adata.write_h5ad(out)
    log_info(f"已写出 {out}（{adata.n_obs} spot x {adata.n_vars} 基因）")

    adata.obs[[c for c in ("total_counts", "n_genes_by_counts", "pct_counts_mt",
                           "pct_counts_ribo", "pct_counts_hb") if c in adata.obs.columns]] \
        .to_csv(res_dir / "qc_spots.csv")

    status = {
        "dataset_id": cfg["dataset_id"],
        "organism": organism,
        "mt_prefix": mt_pref[0],
        "n_mt_genes": n_mt,
        "n_spots_raw": int(n0),
        "n_spots_final": int(adata.n_obs),
        "filter_chain": steps,
        "thresholds": {k: q.get(k) for k in ("min_counts", "max_counts", "min_genes",
                                             "max_genes", "max_pct_mt", "min_spots")},
        "connectivity": conn,
        "threshold_note": ("空间数据的 QC 阈值是**组织特异**的："
                           "PBMC 5-10 / 淋巴结 20 / 心肌 30-50 / 肝 30+。"
                           "spot 是多个细胞的混合物，pct_counts_mt 分布比单细胞窄"),
        "status": "ok",
    }
    write_json(res_dir / "qc_status.json", status)
    return status


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        run_01_qc(cfg)
        record_step(cfg, "qc", "ok", time.time() - t0)
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "qc", "failed", time.time() - t0, message=str(e))
        raise
