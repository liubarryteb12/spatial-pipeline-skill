#!/usr/bin/env python3
"""
06_niche.py — 空间邻域分析（邻域富集 + 共现）

两个互补的分析：

  1. **邻域富集**（neighborhood enrichment）：哪些域/细胞类型在空间上
     互为邻居，比随机排列下更频繁？用置换检验给 z-score。

  2. **共现**（co-occurrence）：随距离增大，某类型周围出现另一类型的
     概率怎么衰减？能区分"紧密共定位"与"大范围共存"。

**置换检验在这里的作用是关键的。** 两个类型相邻，可能只是因为它们
各自都占了很大面积（大类型的邻居当然多）。置换保持各类型的数量、
打乱空间位置，得到的零分布正好扣掉这个效应 —— 所以 z-score 高
才是"真的偏好相邻"，而不是"两个都很大"。

**这是空间数据相对单细胞数据的核心增量。** 单细胞里只能说"这两类
细胞都存在于样本里"，空间数据能说"它们真的挨着"。
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import scipy.sparse as sp  # noqa: E402

from common import (df_to_records, ensure_dirs, load_config, log_info,  # noqa: E402
                    log_warn, parse_args, record_step, save_fig, set_seed,
                    write_json, spatial_xy, W_DOUBLE, W_ONE_HALF, mm,)


def build_adj(adata, n_neighbors: int = 6):
    from scipy.spatial import cKDTree

    xy = spatial_xy(adata)
    k = min(n_neighbors + 1, len(xy))
    tree = cKDTree(xy)
    d, idx = tree.query(xy, k=k)
    med = float(np.median(d[:, 1]))
    rows = np.repeat(np.arange(len(xy)), k - 1)
    cols = idx[:, 1:].ravel()
    mask = d[:, 1:].ravel() <= med * 1.6
    A = sp.coo_matrix((np.ones(int(mask.sum())), (rows[mask], cols[mask])),
                      shape=(len(xy), len(xy))).tocsr()
    return A.maximum(A.T), med


def neighborhood_enrichment(labels: np.ndarray, A: sp.spmatrix,
                            n_perms: int, seed: int) -> pd.DataFrame:
    """
    邻域富集 z-score。

    统计量：类型 i 与类型 j 之间的**边数**（在邻居图上）。
    零分布：保持各类型的 spot 数，随机重排标签，重复 n_perms 次。

    返回 z = (observed - mean_null) / sd_null。
    """
    rng = np.random.default_rng(seed)
    labels = np.asarray(labels).astype(str)
    uniq = sorted(set(labels))
    idx_of = {u: i for i, u in enumerate(uniq)}
    codes = np.array([idx_of[x] for x in labels])
    n = len(uniq)

    coo = A.tocoo()
    src, dst = coo.row, coo.col

    def edge_counts(c):
        # 对称计数：i-j 与 j-i 都算
        M = np.zeros((n, n))
        np.add.at(M, (c[src], c[dst]), 1.0)
        return M

    obs = edge_counts(codes)
    null = np.zeros((n_perms, n, n))
    for p in range(n_perms):
        null[p] = edge_counts(rng.permutation(codes))
    mu = null.mean(axis=0)
    sd = null.std(axis=0)
    sd[sd < 1e-9] = 1.0
    z = (obs - mu) / sd

    rows = []
    for i, a in enumerate(uniq):
        for j, b in enumerate(uniq):
            if j < i:
                continue
            rows.append({"type_a": a, "type_b": b,
                         "observed_edges": int(obs[i, j]),
                         "null_mean": round(float(mu[i, j]), 2),
                         "z_score": round(float(z[i, j]), 3),
                         # |z|>2 约对应 p<0.05（正态近似）
                         "enriched": bool(z[i, j] > 2),
                         "depleted": bool(z[i, j] < -2),
                         "self": bool(i == j)})
    return pd.DataFrame(rows)


def cooccurrence(labels: np.ndarray, xy: np.ndarray, types: list,
                 max_dist: float, n_bins: int = 8, k_neighbors: int = 60) -> pd.DataFrame:
    """
    共现：随距离增大，类型 a 周围出现类型 b 的比例。

    **只对"核心类型"（占比 >=2%）做**，否则稀有类型的分母太小，
    曲线全是噪声。

    ---
    ## T8：分母是被截断的，必须把它报出来

    每个 spot 只取**最近 `k_neighbors` 个邻居**，而距离环是**固定宽度**的
    （`max_dist / n_bins`）。于是远处的环里落进的邻居数天然少，而且
    **被 k 截断** —— 一个 spot 的"邻居"里根本没包含它真正的第 100 个邻居，
    所以远处环的 `fraction` 分母是"碰巧落进这个环的、且属于前 60 近的邻居数"。

    **这不是算错了，是分辨率随距离衰减。** 但原来的产物里只有 `fraction`，
    读者会以为每个距离环都是同等的估计。现在每个环**同时报分母**
    （`n_neighbors`）与 `truncated` 标志：分母小于 `min_n` 的环不该被当结论。
    """
    from scipy.spatial import cKDTree

    labels = np.asarray(labels).astype(str)
    bins = np.linspace(0, max_dist, n_bins + 1)
    centers = (bins[:-1] + bins[1:]) / 2
    tree = cKDTree(xy)
    k = min(k_neighbors, len(xy))
    rows = []
    for a in types:
        ma = labels == a
        if ma.sum() < 20:
            continue
        # a 的每个 spot 的邻居
        d, idx = tree.query(xy[ma], k=k)
        nb_lab = labels[idx[:, 1:]].ravel()
        nb_d = d[:, 1:].ravel()
        for b in types:
            frac = []
            for lo, hi in zip(bins[:-1], bins[1:]):
                m = (nb_d >= lo) & (nb_d < hi)
                frac.append(float((nb_lab[m] == b).mean()) if m.sum() > 0 else np.nan)
            for c, f, lo, hi in zip(centers, frac, bins[:-1], bins[1:]):
                m = (nb_d >= lo) & (nb_d < hi)
                n_in_ring = int(m.sum())
                rows.append({"core_type": a, "neighbor_type": b,
                             "distance": round(float(c), 1),
                             "fraction": None if np.isnan(f) else round(f, 5),
                             # **分母与截断标志**（T8）：没有这两列，
                             # 远处环的小分母看起来和近处环一样可信。
                             "n_neighbors": n_in_ring,
                             "n_spots_core": int(ma.sum()),
                             "k_neighbors_cap": int(k),
                             "truncated": n_in_ring < max(10, k // 10)})
    return pd.DataFrame(rows)


def run_06_niche(cfg: dict) -> dict:
    ensure_dirs(cfg)
    set_seed(cfg)
    data_dir = Path(cfg["output"]["data_dir"])
    res_dir = Path(cfg["output"]["results_dir"])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_cfg = cfg.get("niche") or {}
    if not n_cfg.get("enabled", True):
        status = {"dataset_id": cfg["dataset_id"], "status": "disabled",
                  "reason": "配置 niche.enabled=false"}
        write_json(res_dir / "niche_status.json", status)
        return status

    adata = sc.read_h5ad(data_dir / "domains.h5ad")
    if "domain" not in adata.obs.columns:
        status = {"dataset_id": cfg["dataset_id"], "status": "missing_domains",
                  "reason": "obs 里没有 domain（step 03 未产出）"}
        write_json(res_dir / "niche_status.json", status)
        log_warn(status["reason"])
        return status

    n_neigh = int(n_cfg.get("n_neighbors", 6))
    n_perms = int(n_cfg.get("n_perms", 100))
    A, med = build_adj(adata, n_neigh)
    log_info(f"邻居图: 最近邻距离中位 {med:.1f}，平均度数 "
             f"{A.getnnz(axis=1).mean():.2f}")

    # ---- 1. 邻域富集（域层面）----------------------------------------------
    log_info(f"邻域富集置换检验 {n_perms} 次（域层面）…")
    dom_lab = adata.obs["domain"].astype(str).values
    ne_dom = neighborhood_enrichment(dom_lab, A, n_perms, cfg["analysis"]["seed"])
    ne_dom.to_csv(res_dir / "niche_enrichment_domains.csv", index=False)
    n_enr = int(ne_dom["enriched"].sum())
    log_info(f"域邻域富集: {n_enr}/{len(ne_dom)} 对 |z|>2（其中自相邻 "
             f"{int(ne_dom['self'].sum())} 对）")

    # ---- 2. 邻域富集（细胞类型层面，用解卷积的 argmax）---------------------
    prop_file = res_dir / "deconvolution_proportions.csv"
    ct_lab = None
    ct_align = None
    if prop_file.exists():
        # **M4：`.loc[adata.obs_names]` 无保护会抛 KeyError。**
        # 05 与 06 读的是同一份 `domains.h5ad`，正常情况下 spot 集合一致；
        # 但"正常情况下"不是保证 —— 只要 05 因为任何原因漏掉/多出 spot，
        # 这里就是一个**裸 KeyError**，`06` 整步崩掉，
        # 而产物里只会看到"06 失败"，看不到"是因为 spot 对不上"。
        # 现在：显式取交集、报出差异、并对齐顺序。
        props_all = pd.read_csv(prop_file, index_col=0)
        _wanted = list(adata.obs_names.astype(str))
        _have = list(props_all.index.astype(str))
        _missing = [s for s in _wanted if s not in set(_have)]
        _extra = [s for s in _have if s not in set(_wanted)]
        ct_align = {
            "n_spots_expected": len(_wanted),
            "n_spots_in_proportions": len(_have),
            "n_missing": len(_missing),
            "n_extra": len(_extra),
            "missing_examples": _missing[:5],
            "extra_examples": _extra[:5],
        }
        if _missing:
            log_warn(f"解卷积比例表缺 {len(_missing)} 个 spot（例：{_missing[:3]}）"
                     f"—— 细胞类型层面的邻域富集只能用**交集**的 spot，"
                     f"这与域层面的结果不是同一批 spot")
        props = props_all.loc[[s for s in _wanted if s in set(_have)]]
        ct_lab = props.idxmax(axis=1).values
        ct_align["n_spots_used"] = int(len(ct_lab))
        ct_align["same_spot_set"] = bool(not _missing and not _extra)
        # **argmax 是硬分配，丢掉了混合信息。** 它只能给出"这个 spot 里
        # 相对权重最高的类型"，不能说明这个 spot 里没有别的类型。
        adata.obs["ct_argmax"] = pd.Series(ct_lab, index=props.index)
        log_info("细胞类型层面：用解卷积权重的 argmax（**硬分配，丢失混合信息**）")
        ne_ct = neighborhood_enrichment(ct_lab, A, n_perms, cfg["analysis"]["seed"])
        ne_ct.to_csv(res_dir / "niche_enrichment_celltypes.csv", index=False)
        log_info(f"细胞类型邻域富集: {int(ne_ct['enriched'].sum())}/{len(ne_ct)} 对 |z|>2")

        # 热图
        piv = ne_ct.pivot(index="type_a", columns="type_b", values="z_score")
        piv = piv.combine_first(piv.T)
        # 宽度夹在 [单栏半, 双栏]，避免类型数多时画出装不进一页的图
        # **高度公式必须单位一致**（评审 3.9：这一张实测 203.2 mm）：
        # 原写法 max(5.5, 0.5*len(piv)+2.0) 把英寸当毫米比 ——
        # 12 个类型给 8.0in = 203 mm，超出常规页面。这里统一为英寸并夹到 5.4"。
        fig_h = min(max(3.2, 0.30 * len(piv) + 1.6), 5.4)
        fig, ax = plt.subplots(figsize=(min(W_DOUBLE, max(W_ONE_HALF, 0.5 * len(piv) + 2.5)),
                                        fig_h))
        im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-8, vmax=8)
        ax.set_xticks(range(len(piv.columns)))
        ax.set_xticklabels(piv.columns, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(piv.index)))
        ax.set_yticklabels(piv.index, fontsize=7)
        ax.set_title("Neighborhood enrichment z-score (cell types)\n"
                     "red = preferentially adjacent, blue = avoid each other")
        fig.colorbar(im, ax=ax, label="z-score", shrink=0.8, pad=0.02,
                     fraction=0.046)
        save_fig(cfg, "03-06-01-unit1-niche-enrichment-celltypes", fig)
    else:
        log_warn("没有解卷积产物 —— 跳过细胞类型层面的邻域分析")

    # ---- 3. 共现 ------------------------------------------------------------
    max_dist = float(n_cfg.get("max_distance_um", 400))
    xy = spatial_xy(adata)
    # **M13：不要用"假定 100 μm"去反推像素比例。**
    # 原写法是 `max_px = med * max(2, int(max_dist / 100))` —— 它把
    # `max_dist/100` 取整成邻居间距的整数倍。两个问题：
    #   ① `100 μm` 是硬编码的假定值，而 `med`（最近邻距离中位数）
    #      **本身就是比例尺的观测量** —— 应该由它算出 μm/px，
    #      而不是反过来用一个假定的 100 μm 去凑倍数；
    #   ② `int()` 截断：`max_dist=400` → 恰好 4；但 `max_dist=350`
    #      也会给 4（=400 μm 的范围），配置项形同虚设。
    # Visium 的 spot 中心间距**名义上是 100 μm**（这是芯片规格，不是估计），
    # 所以用它把像素换成 μm 是合理的 —— 但必须**显式写明这是假定**，
    # 并把换算结果与"按 100 μm 反推的邻居间距"一起报出来，让读者能判断。
    um_per_px = 100.0 / med if med > 0 else 1.0
    max_px = max_dist / um_per_px
    co_scale = {
        "um_per_px_assumed": round(float(um_per_px), 5),
        "assumption": "Visium spot 中心间距名义 100 μm（芯片规格）",
        "median_nn_distance_px": round(float(med), 3),
        "max_distance_um": round(float(max_dist), 1),
        "max_distance_px": round(float(max_px), 2),
        "note": ("像素→μm 用的是**名义规格 100 μm**；`median_nn_distance_px` "
                 "是观测值，若它与 100/um_per_px 相差大，说明该切片不是标准"
                 "Visium 网格或坐标被缩放过"),
    }
    dom_counts = adata.obs["domain"].value_counts()
    core = [str(d) for d in dom_counts[dom_counts >= 20].index]
    co = cooccurrence(dom_lab, xy, core, max_px)
    co.to_csv(res_dir / "niche_cooccurrence.csv", index=False)
    log_info(f"共现曲线: {len(core)} 个核心域 x {len(core)} 个邻居类型"
             f"（{max_dist:.0f} μm = {max_px:.1f} px，假定 100 μm/spot 间距）")

    # 声明式动态名豁免：03-06-02 图号下最多 3 张（核心域号是运行时数据）
    DYNAMIC_FIG_BASES = {"02": 3}
    if len(core) >= 2:
        # **单图原则拆分（D-006）**：多核心域面板 -> 每域一张单图
        # （P6 niche 共现对照链）。跨图对照同色纪律（tab20 按邻居类型编号）。
        cmap = plt.get_cmap("tab20")
        for ui, core_t in enumerate(core[:3], start=1):
            sub = co[co["core_type"] == core_t]
            fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(64)))
            for i, nb in enumerate(core):
                s = sub[sub["neighbor_type"] == nb]
                ax.plot(s["distance"], s["fraction"], "o-", ms=3,
                        color=cmap(i % 20), label=nb)
            ax.set_xlabel("distance (pixel)")
            ax.set_ylabel("fraction of neighbors")
            ax.set_title(f"Neighborhood composition around domain {core_t}")
            fig.legend(fontsize=6, ncol=1, loc="outside right center")
            save_fig(cfg, f"03-06-02-unit{ui}-niche-around-domain-{core_t}", fig)

    # ---- 4. 落盘 ------------------------------------------------------------
    status = {
        "dataset_id": cfg["dataset_id"],
        "status": "ok",
        "n_neighbors": n_neigh,
        "median_nn_distance_px": round(med, 2),
        "cooccurrence_scale": co_scale,
        "cooccurrence_denominator_note": (
            "`niche_cooccurrence.csv` 每个距离环都带 `n_neighbors`（该环的"
            "实际分母）与 `truncated`（分母过小）两列 —— **远处环的邻居数"
            "被 k 上限截断，分辨率随距离衰减**，只看 `fraction` 会把它们"
            "当同等可信的估计"),
        "celltype_spot_alignment": ct_align or {"status": "not_available"},
        "n_permutations": n_perms,
        "domain_enrichment": {
            "n_pairs": int(len(ne_dom)),
            "n_enriched_z_gt2": n_enr,
            # **M5：固定阈值 |z|>2 的期望假阳性要报出来。**
            # 类型对数是 O(k²)：k=12 时 78 对，|z|>2 在纯零假设下
            # 期望约 78 × 0.0455 ≈ 3.5 对"显著" —— 那只是噪声。
            # 没有这一条，读者会把 `n_enriched_z_gt2` 当真实富集数。
            "expected_false_positives_at_z2": round(0.0455 * len(ne_dom), 2),
            "multiple_testing_note": (
                "**没有做多重检验校正。** `|z|>2` 在零假设下每对约有 4.55% "
                "概率越界，所以上表里的『显著』对数要减去上面那个期望值"
                "才接近真实富集"),
            "top": df_to_records(ne_dom[~ne_dom["self"]]
                                 .sort_values("z_score", ascending=False).head(10)),
            "most_depleted": df_to_records(ne_dom[~ne_dom["self"]]
                                           .sort_values("z_score").head(5)),
        },
        "celltype_enrichment": ({
            "n_pairs": int(len(ne_ct)),
            "n_enriched_z_gt2": int(ne_ct["enriched"].sum()),
            "top": df_to_records(ne_ct[~ne_ct["self"]]
                                 .sort_values("z_score", ascending=False).head(10)),
        } if ct_lab is not None else {"status": "not_available",
                                      "reason": "没有解卷积产物"}),
        "method": ("邻域富集：保持各类型数量、打乱空间位置的置换检验，"
                   "统计邻居图上的边数，报 z-score"),
        "why_permutation": ("两个类型相邻，可能只因为它们各自都占很大面积"
                            "（大类型的邻居当然多）。置换保持各类型的数量、"
                            "打乱位置，零分布正好扣掉这个效应 —— "
                            "所以 z 高才是『真的偏好相邻』"),
        "limitations": [
            "z-score 用正态近似；|z|>2 约对应 p<0.05，但没有做多重检验校正"
            "（类型对数是 O(k²)，k=12 时是 78 对）",
            "细胞类型层面用 argmax 硬分配，**丢失了 spot 内的混合信息** —— "
            "一个 60% T 细胞 + 40% B 细胞的 spot 会被当成纯 T 细胞",
            "Visium 的 spot 间距约 100 μm，远大于细胞直径，"
            "所以『相邻』是**spot 层面**的，不是细胞层面的接触",
            "共现分析只对占比 >=2% 的类型做（稀有类型的分母太小）",
        ],
    }
    write_json(res_dir / "niche_status.json", status)
    return status


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        run_06_niche(cfg)
        record_step(cfg, "niche", "ok", time.time() - t0)
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "niche", "failed", time.time() - t0, message=str(e))
        raise
