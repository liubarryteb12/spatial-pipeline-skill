#!/usr/bin/env python3
"""
04_svg.py — 空间高变基因（Moran's I / Geary's C）

**SVG ≠ HVG，两者问的是不同问题：**
  - HVG：这个基因在**样本之间**变异大吗？（不看位置）
  - SVG：这个基因在**空间上**有结构吗？（邻居比随机更相似吗）

一个基因可以既是 HVG 又不是 SVG（在组织里均匀高表达但在样本间差异大），
也可以既是 SVG 又不是 HVG（在组织里形成清晰的梯度但整体变异不大）。

**Moran's I 的解读依赖空间权重矩阵。** 同一个基因用不同的邻居定义
会得到不同的 I 值。所以这里把权重矩阵的参数一并落盘。

**p 值用置换检验，并且必须做多重检验校正。** 20000 个基因里
一定有一堆 p<0.05，而那只是基因多。

**本工具自己实现 Moran's I，不依赖 squidpy。** 理由：squidpy 的
`spatial_autocorr` 在版本间改过参数名和返回值结构；自己实现
（20 行）比跟着上游 API 漂移更稳，而且能明确控制置换方式。
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
                    spot_radius_plot_units, write_json, spatial_xy,)


def build_weights(adata, n_neighbors: int = 6, row_standardize: bool = True):
    """
    建空间权重矩阵。

    **行标准化很重要。** 不标准化的话，度数高的 spot（组织内部的）
    对统计量的贡献大于度数低的（边缘的），而边缘 spot 恰恰是
    空间结构最有信息量的地方。
    """
    from scipy.spatial import cKDTree

    xy = spatial_xy(adata)
    k = min(n_neighbors + 1, len(xy))
    tree = cKDTree(xy)
    d, idx = tree.query(xy, k=k)
    med = float(np.median(d[:, 1]))
    thresh = med * 1.6

    rows = np.repeat(np.arange(len(xy)), k - 1)
    cols = idx[:, 1:].ravel()
    mask = d[:, 1:].ravel() <= thresh
    W = sp.coo_matrix((np.ones(int(mask.sum())), (rows[mask], cols[mask])),
                      shape=(len(xy), len(xy))).tocsr()
    W = W.maximum(W.T)
    if row_standardize:
        deg = np.asarray(W.sum(axis=1)).ravel()
        deg[deg == 0] = 1.0
        W = sp.diags(1.0 / deg) @ W
    return W, {"n_neighbors_k": int(k - 1), "median_nn_dist": round(med, 2),
               "distance_threshold": round(thresh, 2),
               "row_standardized": bool(row_standardize),
               "mean_degree": round(float(W.getnnz(axis=1).mean()), 2)}


def morans_i(X: np.ndarray, W: sp.spmatrix) -> np.ndarray:
    """
    Moran's I，按列向量化计算。

    I = (n / S0) * (x' W x) / (x' x)，其中 x 是中心化后的表达。
    行标准化后 S0 = n，公式化简为 x'Wx / x'x。

    **一次算全部基因**，不是循环 —— 循环 20000 个基因在 Python 里很慢。
    """
    n = X.shape[0]
    Xc = X - X.mean(axis=0, keepdims=True)
    denom = (Xc ** 2).sum(axis=0)
    denom[denom == 0] = 1.0
    # (W @ Xc) 再与 Xc 逐元素乘后按列求和 = x' W x
    num = (Xc * (W @ Xc)).sum(axis=0)
    S0 = float(W.sum())
    return (n / S0) * (num / denom)


def gearys_c(X: np.ndarray, W: sp.spmatrix) -> np.ndarray:
    """Geary's C。C 接近 0 表示强正空间自相关（与 Moran's I 相反）。"""
    n = X.shape[0]
    Xc = X - X.mean(axis=0, keepdims=True)
    denom = (Xc ** 2).sum(axis=0)
    denom[denom == 0] = 1.0
    # sum_ij w_ij (xi - xj)^2 的向量化：
    # = 2 * (x' diag(W.sum(1)) x - x' W x)  对对称 W 成立
    deg = np.asarray(W.sum(axis=1)).ravel()
    term = (Xc ** 2 * deg[:, None]).sum(axis=0) - (Xc * (W @ Xc)).sum(axis=0)
    S0 = float(W.sum())
    return (n - 1) / (2 * S0) * (2 * term / denom)


def permutation_pvalues(X: np.ndarray, W: sp.spmatrix, observed: np.ndarray,
                        n_perms: int, seed: int) -> np.ndarray:
    """
    置换检验 p 值。

    **置换的是表达值在 spot 上的分配**（打乱空间位置），保留表达分布。
    这检验的是"这个基因的空间排布是否比随机更结构化"。

    一次性对所有基因做：每次置换算一遍全部基因的统计量。
    比逐基因置换快得多，而且置换次数对每个基因相同（p 值可比）。
    """
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    n_ge = np.zeros(X.shape[1], dtype=np.int64)
    for _ in range(n_perms):
        perm = rng.permutation(n)
        null = morans_i(X[perm], W)
        n_ge += (null >= observed)
    # +1 校正：避免 p=0（有限次置换下不可能有比 1/(n+1) 更小的 p）
    return (n_ge + 1) / (n_perms + 1)


def run_04_svg(cfg: dict) -> dict:
    ensure_dirs(cfg)
    set_seed(cfg)
    data_dir = Path(cfg["output"]["data_dir"])
    res_dir = Path(cfg["output"]["results_dir"])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    s = cfg.get("svg") or {}
    if not s.get("enabled", True):
        status = {"dataset_id": cfg["dataset_id"], "status": "disabled",
                  "reason": "配置 svg.enabled=false"}
        write_json(res_dir / "svg_status.json", status)
        return status

    adata = sc.read_h5ad(data_dir / "domains.h5ad")
    log_info(f"读入 {adata.n_obs} spot x {adata.n_vars} HVG")

    # ---- 1. 权重矩阵 --------------------------------------------------------
    W, winfo = build_weights(adata, int(s.get("n_neighbors", 6)))
    log_info(f"空间权重矩阵: 平均度数 {winfo['mean_degree']}，"
             f"行标准化 {winfo['row_standardized']}")

    # ---- 2. 用全基因集还是 HVG？--------------------------------------------
    # **用全基因集**（adata.raw）。HVG 是按"样本间变异"选的，会漏掉
    # 空间结构强但整体变异不大的基因 —— 而那正是 SVG 要找的。
    # 但全基因集有 3 万多个，置换检验会很慢，所以先按方差粗筛。
    if adata.raw is not None:
        Xall = adata.raw.X
        Xall = Xall.toarray() if sp.issparse(Xall) else np.asarray(Xall)
        genes_all = list(adata.raw.var_names)
        var = Xall.var(axis=0)
        # 保留方差 > 0 的基因（全零基因的 I 无定义）
        keep = np.flatnonzero(var > 0)
        # 上限：按方差取前 N 个，控制置换检验耗时
        max_genes = int(s.get("max_genes", 4000))
        if len(keep) > max_genes:
            keep = keep[np.argsort(var[keep])[::-1][:max_genes]]
            log_info(f"全基因集 {len(genes_all)} 个，按方差取前 {max_genes} 个做置换检验")
        X = np.ascontiguousarray(Xall[:, keep].astype(np.float64))
        genes = [genes_all[i] for i in keep]
    else:
        X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
        X = np.ascontiguousarray(X.astype(np.float64))
        genes = list(adata.var_names)
        log_warn("adata.raw 为空 —— 只用 HVG 做 SVG。"
                 "**HVG 是按样本间变异选的，会漏掉空间结构强但整体变异不大的基因**")

    log_info(f"SVG 候选基因: {X.shape[1]} 个")

    # ---- 3. 统计量 ----------------------------------------------------------
    method = s.get("method", "moran")
    if method == "moran":
        stat = morans_i(X, W)
        stat_name = "morans_I"
        # Moran's I 的期望值（无空间自相关时）≈ -1/(n-1)
        expected = -1.0 / (X.shape[0] - 1)
    elif method == "geary":
        stat = gearys_c(X, W)
        stat_name = "gearys_C"
        expected = 1.0
    else:
        raise ValueError(f"不支持的 svg.method: {method}（moran / geary）")
    log_info(f"{stat_name} 计算完成（期望值 {expected:.5f}）")

    # ---- 4. 置换 p 值 -------------------------------------------------------
    n_perms = int(s.get("n_perms", 100))
    log_info(f"置换检验 {n_perms} 次（打乱表达值在 spot 上的分配）…")
    pvals = permutation_pvalues(X, W, stat, n_perms, cfg["analysis"]["seed"])

    # ---- 5. 多重检验校正 ----------------------------------------------------
    from statsmodels.stats.multitest import multipletests
    p_adj = multipletests(pvals, method="fdr_bh")[1]

    res = pd.DataFrame({"gene": genes, stat_name: stat,
                        "expected": expected, "p_value": pvals,
                        "p_adj_bh": p_adj})
    # Moran's I：越大越有空间结构；Geary's C：越小越有
    res["spatial_score"] = (res[stat_name] - expected if method == "moran"
                            else expected - res[stat_name])
    res = res.sort_values("spatial_score", ascending=False).reset_index(drop=True)
    res.to_csv(res_dir / "svg_results.csv", index=False)

    n_sig = int((res["p_adj_bh"] < 0.05).sum())
    log_info(f"SVG: {n_sig}/{len(res)} 个基因 BH 校正后 p<0.05")
    top = res.head(10)
    log_info("  top10: " + ", ".join(
        f"{r.gene}({getattr(r, stat_name):.3f})" for r in top.itertuples()))

    # ---- 6. 出图 ------------------------------------------------------------
    n_top = int(s.get("n_top", 30))
    top_genes = res.head(n_top)["gene"].tolist()
    ncol = 6
    nrow = int(np.ceil(len(top_genes) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.5 * ncol, 2.6 * nrow))
    axes = np.atleast_1d(axes).ravel()
    sf = float(adata.uns["spatial"][list(adata.uns["spatial"])[0]]
               ["scalefactors"]["tissue_hires_scalef"])
    xy = spatial_xy(adata, sf)
    gi = {g: i for i, g in enumerate(genes)}
    for ax, g in zip(axes, top_genes):
        v = X[:, gi[g]]
        sc_ = ax.scatter(xy[:, 0], xy[:, 1], c=v, s=3, cmap="viridis")
        row = res[res["gene"] == g].iloc[0]
        ax.set_title(f"{g}\nI={row[stat_name]:.3f}", fontsize=7)
        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
    for ax in axes[len(top_genes):]:
        ax.axis("off")
    fig.suptitle(f"Top {len(top_genes)} spatially variable genes "
                 f"({stat_name}, BH p<0.05: {n_sig})", fontsize=11)
    save_fig(cfg, "svg_top_genes", fig)

    # 统计量分布
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    ax.hist(res[stat_name], bins=60, color="#2C7FB8", alpha=0.85)
    ax.axvline(expected, color="#B2182B", ls="--", lw=1.2,
               label=f"expected (no autocorr) = {expected:.4f}")
    ax.set_xlabel(stat_name); ax.set_ylabel("number of genes")
    ax.set_title(f"{stat_name} distribution across {len(res)} genes", fontsize=10)
    ax.legend(fontsize=8)
    save_fig(cfg, "svg_stat_distribution", fig)

    status = {
        "dataset_id": cfg["dataset_id"],
        "status": "ok",
        "method": method,
        "stat_name": stat_name,
        "expected_under_no_autocorrelation": round(float(expected), 5),
        "n_genes_tested": int(len(res)),
        "n_significant_bh": n_sig,
        "n_permutations": n_perms,
        "weights": winfo,
        "top_genes": df_to_records(res.head(20)),
        "gene_set_note": ("用的是**全基因集**（adata.raw）按方差粗筛后的结果。"
                          "HVG 是按样本间变异选的，会漏掉空间结构强但"
                          "整体变异不大的基因 —— 而那正是 SVG 要找的"),
        "method_note": ("**SVG ≠ HVG。** HVG 问『样本间变异大吗』（不看位置），"
                        "SVG 问『邻居比随机更相似吗』。两者可以完全不同"),
        "limitations": [
            f"{stat_name} 的数值**依赖空间权重矩阵**：换邻居数或标准化方式会得到不同值",
            "置换检验打乱的是表达值在 spot 上的分配，保留了表达分布；"
            "这检验的是空间排布，不是表达差异",
            f"置换 {n_perms} 次，最小可能 p 值是 1/{n_perms + 1} —— "
            f"更小的 p 值无法区分",
            "Visium 的 spot 含 1-10 个细胞，所以 SVG 可能是『细胞组成在空间上变化』"
            "而非『同一类细胞内该基因被调控』",
        ],
    }
    write_json(res_dir / "svg_status.json", status)
    return status


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        run_04_svg(cfg)
        record_step(cfg, "svg", "ok", time.time() - t0)
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "svg", "failed", time.time() - t0, message=str(e))
        raise
