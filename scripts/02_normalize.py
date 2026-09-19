#!/usr/bin/env python3
"""
02_normalize.py — 标准化、高变基因、PCA

与单细胞流程的关键差别：

  **Visium 的 spot 是多个细胞的混合物，所以 HVG 会偏向"细胞组成差异"
  而不是"细胞状态差异"。** 这不算错 —— 空间域划分本来就想抓组成差异 ——
  但要意识到：一个基因可能只因为"这块区域的细胞构成不同"而成为 HVG，
  而不是因为它在同一类细胞里被调控了。

  另外 **必须保留原始计数层**：解卷积的 NNLS 要用计数，
  用 log 后的值会给出错的组成估计。
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))

import numpy as np  # noqa: E402
import scanpy as sc  # noqa: E402

from common import (df_to_records, ensure_dirs, load_config, log_info,  # noqa: E402
                    log_warn, parse_args, record_step, save_fig, set_seed, write_json, spatial_xy,)


def run_02_normalize(cfg: dict) -> dict:
    ensure_dirs(cfg)
    set_seed(cfg)
    data_dir = Path(cfg["output"]["data_dir"])
    res_dir = Path(cfg["output"]["results_dir"])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    adata = sc.read_h5ad(data_dir / "qc_filtered.h5ad")
    log_info(f"读入 {adata.n_obs} spot x {adata.n_vars} 基因")

    norm = cfg["norm"]
    flavor = norm.get("hvg_flavor", "seurat_v3")
    n_top = int(norm["n_top_genes"])
    target_sum = float(norm["target_sum"])

    # **原始计数单独存一层。** 解卷积的 NNLS 与任何基于计数的分析都要用它。
    adata.layers["counts"] = adata.X.copy()

    # ---- 1. HVG（顺序按口味分）---------------------------------------------
    # seurat_v3 吃原始计数，必须在 normalize 之前；另两种口味相反。
    # **顺序反了不报错，只是给出错的 HVG。**
    hvg_flavor_used = flavor
    hvg_fallback = None
    if flavor == "seurat_v3":
        try:
            import skmisc  # noqa: F401
        except ImportError:
            hvg_flavor_used = "seurat"
            hvg_fallback = {
                "requested": "seurat_v3", "used": "seurat",
                "reason": "scikit-misc (skmisc) 未安装，seurat_v3 的 loess 拟合不可用",
                "impact": "两种口味选出的 HVG 集合不同，下游 PCA/空间域/SVG 全部跟着变",
                "fix": "pip install scikit-misc",
            }
            log_warn(f"HVG 口味降级: seurat_v3 -> seurat —— {hvg_fallback['reason']}")

    if hvg_flavor_used == "seurat_v3":
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="seurat_v3",
                                    layer="counts")
        log_info("HVG (seurat_v3, 用原始计数) 先算，再 normalize")
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
    else:
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top,
                                    flavor=hvg_flavor_used)
        log_info(f"HVG ({hvg_flavor_used}, 用 log 后数据) 后算")

    n_hvg = int(adata.var["highly_variable"].sum())
    log_info(f"高变基因: {n_hvg}")
    adata.raw = adata

    sc.pl.highly_variable_genes(adata, show=False)
    fig = plt.gcf()
    fig.suptitle(f"HVG ({hvg_flavor_used}, n={n_hvg})", fontsize=10)
    save_fig(cfg, "hvg_selection", fig)

    # ---- 2. PCA -------------------------------------------------------------
    work = adata[:, adata.var["highly_variable"]].copy()
    if norm.get("scale", True):
        sc.pp.scale(work, max_value=10)
    sc.tl.pca(work, svd_solver="arpack", random_state=cfg["analysis"]["seed"])

    sc.pl.pca_variance_ratio(work, log=True, show=False)
    fig = plt.gcf()
    fig.suptitle("PCA variance ratio (elbow)", fontsize=10)
    save_fig(cfg, "pca_variance_ratio", fig)

    var_ratio = work.uns["pca"]["variance_ratio"]
    log_info(f"PC1-10 方差解释: {', '.join(f'{v:.3f}' for v in var_ratio[:10])}")

    # PCA 空间投影 —— 看主成分是否有空间结构
    # **这是空间数据特有的诊断**：如果 PC1 在组织上是随机斑点，
    # 说明主要变异是技术噪声而不是空间结构。
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    xy = spatial_xy(work)
    for ax, i in zip(axes, range(3)):
        s = ax.scatter(xy[:, 0], xy[:, 1], c=work.obsm["X_pca"][:, i],
                       s=4, cmap="RdBu_r",
                       vmin=-np.abs(work.obsm["X_pca"][:, i]).max(),
                       vmax=np.abs(work.obsm["X_pca"][:, i]).max())
        ax.set_title(f"PC{i+1} ({var_ratio[i]*100:.1f}%)", fontsize=9)
        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(s, ax=ax, shrink=0.8)
    fig.suptitle("PC scores on tissue — spatial structure means PCs capture "
                 "histology, not just noise", fontsize=10)
    save_fig(cfg, "pca_on_tissue", fig)

    # ---- 3. 落盘 ------------------------------------------------------------
    out = data_dir / "normalized.h5ad"
    work.write_h5ad(out)
    log_info(f"已写出 {out}（{work.n_obs} spot x {work.n_vars} HVG）")

    status = {
        "dataset_id": cfg["dataset_id"],
        "n_spots": int(work.n_obs),
        "n_hvg": n_hvg,
        "hvg_flavor_requested": flavor,
        "hvg_flavor_used": hvg_flavor_used,
        "hvg_fallback": hvg_fallback,
        "n_top_genes_requested": n_top,
        "target_sum": target_sum,
        "scaled_before_pca": bool(norm.get("scale", True)),
        "n_pcs": int(cfg["reduce"]["n_pcs"]),
        "pca_variance_ratio_top10": [round(float(v), 5) for v in var_ratio[:10]],
        "pca_cumvar_top10": [round(float(v), 5) for v in np.cumsum(var_ratio[:10])],
        "counts_layer_preserved": "counts" in work.layers or "counts" in adata.layers,
        "note": ("Visium 的 spot 是多个细胞的混合物，HVG 会偏向**细胞组成差异**"
                 "而非细胞状态差异 —— 这对空间域划分是想要的，"
                 "但不能解读成『同一类细胞内的调控差异』"),
        "status": "ok",
    }
    write_json(res_dir / "normalize_status.json", status)
    return status


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        run_02_normalize(cfg)
        record_step(cfg, "normalize", "ok", time.time() - t0)
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "normalize", "failed", time.time() - t0, message=str(e))
        raise
