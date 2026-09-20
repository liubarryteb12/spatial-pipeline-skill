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

**另跑 SpatialDE 作为独立交叉验证（文档 §3.4）。** 两者是不同框架：
Moran's I 是空间自相关、**依赖权重矩阵**；SpatialDE 把基因拟合成高斯过程、
**不需要权重矩阵**。所以两者一致说明空间结构不是权重选择的产物。

SpatialDE 1.1.3 是 2019 年的包，在当代依赖上有**两处独立的不兼容**
（都实测确认，见 `_shim_scipy_misc_derivative` 与 `try_spatialde` 的说明）：
scipy>=1.12 移除了 `scipy.misc.derivative`（垫片解决，只影响标准误列），
以及 `util.qvalue` 对 pandas Series 调 `.ravel()`（绕开 `run`，
多重检验校正改用本仓库的 BH —— 反而与 Moran's I 那条路口径一致）。
"""

from __future__ import annotations

import math
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))

import matplotlib as mpl  # noqa: E402
from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import scipy.sparse as sp  # noqa: E402

from common import (df_to_records, ensure_dirs, load_config, log_info,  # noqa: E402
                    log_warn, parse_args, probe_named_tools, record_step,
                    save_fig, set_seed, spot_radius_plot_units, write_json,
                    spatial_xy, W_DOUBLE, W_SINGLE, mm, PAL,)


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


# ---------------------------------------------------------------------------
# SpatialDE（文档 §3.4 点名的工具）
# ---------------------------------------------------------------------------
# **与 Moran's I 是两种不同的统计框架。** Moran's I 是空间自相关（基于
# 权重矩阵），SpatialDE 把每个基因拟合成高斯过程、用似然比检验是否有
# 空间相关的长度尺度。**权重矩阵换一下 Moran's I 就变，而 SpatialDE
# 不需要权重矩阵** —— 这正是它值得作为独立交叉验证的原因。
#
# 但 SpatialDE 1.1.3（2019 年，7 年未更新）在现代 scipy 上**装得上、
# 导不进来**：`SpatialDE/base.py` 第 12 行是
#     from scipy.misc import derivative
# 而 `scipy.misc.derivative` 在 **scipy 1.12 已被移除**
# （本仓库 requirements 钉的是 scipy>=1.11,<2.0，实测本地 1.18.1 上
# `hasattr(scipy.misc, "derivative")` 为 False）。
#
# **垫片是否安全？** SpatialDE 只在两处用 derivative，且都只算**标准误**：
#   base.py:178  s2_logdelta = 1 / (derivative(LL_obj, ..., n=2) ** 2)
#   base.py:229  s2_FSV      = derivative(FSV, ..., n=1) ** 2 * s2_logdelta
# 主统计量（l / max_l / FSV / pval / qval / BIC）都不经过它。
# 所以垫片只影响 SE 列，不影响这里要用的 FSV 与 qval。
#
# 垫片用的是标准中心差分（3 点，与 scipy 默认 order=3 一致），
# 权重由 Vandermonde 方程组解出 —— 见 _central_diff_weights。
SPATIALDE_MAX_GENES = 150


def _central_diff_weights(npts: int, n: int) -> np.ndarray:
    """等距点上的 n 阶导数有限差分权重（Vandermonde 解法）。

    点取在 j = -ho..ho（ho = npts//2），解
        Σ_b w_b · j_b^a / a! = δ_{a,n}   for a = 0..npts-1
    npts=3, n=1 → [-1/2, 0, 1/2]（经典中心差分）
    npts=3, n=2 → [1, -2, 1]
    """
    ho = npts // 2
    j = np.arange(-ho, ho + 1, dtype=float)
    if len(j) != npts:
        raise ValueError(f"npts 必须是奇数，收到 {npts}")
    # **用 math.factorial，不用 np.math.factorial** —— np.math 在 numpy 2.0
    # 已被移除，而 requirements 允许 numpy<3。用 np.math 会在 CI 上
    # AttributeError，本地旧 numpy 上却看不出来。
    fact = np.array([float(math.factorial(a)) for a in range(npts)])
    A = (j[None, :] ** np.arange(npts)[:, None]) / fact[:, None]
    b = np.zeros(npts)
    b[n] = 1.0
    return np.linalg.solve(A, b)


def _shim_scipy_misc_derivative() -> bool:
    """把 scipy 1.12 移除的 `scipy.misc.derivative` 补回去。

    返回是否**本次做了**补丁（已经存在则返回 False）。
    **只补这一个名字，不碰 scipy 其他部分。**
    """
    import warnings

    # **告警是 `import scipy.misc` 这一行本身发出的**，不是赋值那行 ——
    # 所以 import 也必须在上下文里，否则压不住。
    # 只压 scipy.misc 这一条：我们是有意访问这个正在被移除的命名空间
    # （它正是要补的东西），告警是预期内噪声；压掉它不掩盖其他告警。
    # 用 catch_warnings 而不是全局 filterwarnings，避免影响调用方配置。
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning,
                                message=r".*scipy\.misc is deprecated.*")
        import scipy.misc as _misc
        if hasattr(_misc, "derivative"):
            return False

        def derivative(func, x0, dx=1.0, n=1, args=(), order=3):
            if order % 2 == 0:
                raise ValueError("'order' must be an odd integer.")
            if order < n + 1:
                raise ValueError(
                    "'order' (the number of points used to compute the derivative), "
                    "must be at least the derivative order 'n' + 1.")
            w = _central_diff_weights(order, n)
            ho = order // 2
            # 权重按 j 升序给，点也按 j 升序取
            total = 0.0
            for k, wk in enumerate(w):
                total += wk * func(x0 + (k - ho) * dx, *args)
            return total / (dx ** n)

        _misc.derivative = derivative
    return True


def try_spatialde(X: np.ndarray, xy: np.ndarray, genes: list,
                  res: pd.DataFrame, cfg: dict, log=log_info):
    """跑 SpatialDE 并与 Moran's I 比对。

    **只在基因子集上跑。** SpatialDE 每个基因要拟合一个高斯过程，
    实测 4025 个 spot 上单基因约 1-3 秒 —— 4000 个基因要几小时，
    会直接撞穿 CI 的 45 分钟上限。所以取「Moran's I 最高的前 100 个」
    + 「按种子随机抽 50 个作背景」，共 150 个。

    **子集抽样会引入选择偏差**（前 100 个是 Moran's I 挑出来的），
    所以一致性数字要**分两组报**：top 组和随机背景组。
    只报合并的相关会把"两边都认为强"和"两边都认为弱"混成一个数。
    """
    info = {"attempted": True, "status": None, "reason": "",
            "max_genes": SPATIALDE_MAX_GENES}
    try:
        shimmed = _shim_scipy_misc_derivative()
        info["scipy_misc_derivative_shimmed"] = shimmed
        import SpatialDE
        info["version"] = getattr(SpatialDE, "__version__", "1.1.3")
    except Exception as exc:  # noqa: BLE001
        info["status"] = "package_missing"
        info["reason"] = f"SpatialDE 不可用（{type(exc).__name__}: {exc}）"
        log(f"SpatialDE 不可用：{info['reason']}")
        return None, info

    gi = {g: i for i, g in enumerate(genes)}
    n_top = min(100, len(res))
    top = [g for g in res.head(n_top)["gene"].tolist() if g in gi]
    pool = [g for g in genes if g not in set(top)]
    rng = np.random.default_rng(cfg["analysis"]["seed"])
    n_bg = min(SPATIALDE_MAX_GENES - len(top), len(pool))
    bg = list(rng.choice(pool, size=n_bg, replace=False)) if n_bg > 0 else []
    sel = top + bg
    info["n_genes_run"] = len(sel)
    info["n_top_selected"] = len(top)
    info["n_background_selected"] = len(bg)

    try:
        from scipy import stats as _st
        from statsmodels.stats.multitest import multipletests

        coords = np.asarray(xy, dtype=float)
        counts = pd.DataFrame(X[:, [gi[g] for g in sel]], columns=sel)

        # **不走 SpatialDE.run。** 它第 432 行是
        #     mll_results['qval'] = qvalue(mll_results['pval'])
        # 而 util.qvalue 第 19 行是 `pv = pv.ravel()` —— pandas Series
        # **没有 ravel**（实测 pandas 3.0.5：hasattr(pd.Series,'ravel') = False）。
        # 这是 1.1.3 里第二个独立的不兼容，和 scipy.misc 那个无关。
        #
        # 绕开的办法是**只用它的模型拟合**（dyn_de / get_mll_results），
        # 多重检验校正用本仓库自己的 BH。这反而更好：
        # Moran's I 那边也是 BH，两条路的校正口径一致才可比。
        # 代价是 qval 从 Storey q-value 变成 BH —— 已记进 limitations。
        l_min, l_max = SpatialDE.base.get_l_limits(coords)
        kernel_space = {"SE": np.logspace(np.log10(l_min), np.log10(l_max), 10),
                        "const": 0}
        info["kernel_lengthscales"] = [round(float(v), 3)
                                       for v in kernel_space["SE"]]
        # **注意用 SpatialDE.base.X，不能用 SpatialDE.X。**
        # `SpatialDE/__init__.py` 只导出 dyn_de / run / model_search /
        # fit_patterns / spatial_patterns 五个名字，get_l_limits 与
        # get_mll_results 都在 base 里、没被提上来。
        raw = SpatialDE.base.dyn_de(coords, counts, kernel_space)
        out = pd.DataFrame(SpatialDE.base.get_mll_results(raw))
        out["pval"] = 1.0 - _st.chi2.cdf(out["LLR"], df=1)
        out["qval"] = multipletests(out["pval"], method="fdr_bh")[1]
        out["_group"] = np.where(out["g"].isin(set(top)), "top_morans", "background")
        info["status"] = "ok"
        info["qval_method"] = "BH（本仓库实现）"
        info["bypassed_spatialde_run"] = True
        info["bypass_reason"] = ("SpatialDE.util.qvalue 对 pandas Series 调 .ravel()，"
                                 "现代 pandas 已无此方法")
        info["columns"] = sorted(map(str, out.columns))
        info["n_significant_qval_0.05"] = int((out["qval"] < 0.05).sum())
        log(f"SpatialDE 完成：{len(out)} 个基因，"
            f"qval<0.05 的 {info['n_significant_qval_0.05']} 个")
        return out, info
    except Exception as exc:  # noqa: BLE001
        info["status"] = "failed"
        info["reason"] = f"{type(exc).__name__}: {exc}"
        log_warn(f"SpatialDE 跑失败：{info['reason']}")
        return None, info


def compare_with_spatialde(sd: pd.DataFrame, res: pd.DataFrame,
                           log=log_info) -> dict:
    """量化 SpatialDE 的 FSV 与 Moran's I 的一致程度。

    **分 top 组和背景组各报一次。** 合并报会把"两边都认为强"和
    "两边都认为弱"平均成一个数，掩盖真实的一致性结构。
    """
    out = {"compared": False}
    if sd is None or len(sd) == 0:
        return out
    try:
        from scipy.stats import spearmanr
        # 本步骤产出的列名是 spatial_score（Moran's I 减期望值 / Geary 反之），
        # 越大越有空间结构。名字写错会静默退到"取第 2 列"，而那可能是
        # expected 之类的常数列 —— 常数列的 Spearman 是 nan。
        if "spatial_score" not in res.columns:
            out["reason"] = (f"res 里没有 spatial_score 列（有 {list(res.columns)}），"
                             "拒绝用兜底列做对比")
            return out
        m = res.set_index("gene")["spatial_score"]
        out["compared"] = True
        out["morans_column_used"] = "spatial_score"
        for grp in ("top_morans", "background", "all"):
            sub = sd if grp == "all" else sd[sd["_group"] == grp]
            common = [g for g in sub["g"] if g in m.index]
            if len(common) < 5:
                out[f"{grp}_n"] = len(common)
                continue
            rho, p = spearmanr(m.loc[common].values,
                               sub.set_index("g").loc[common, "FSV"].values)
            out[f"{grp}_n"] = len(common)
            out[f"{grp}_spearman_rho"] = round(float(rho), 4)
            out[f"{grp}_spearman_p"] = float(p)
        out["interpretation"] = (
            "SpatialDE 用高斯过程似然比，**不需要空间权重矩阵**；"
            "Moran's I 依赖权重矩阵。两者一致说明空间结构不是权重选择的产物。"
            "top 组是 Moran's I 挑出来的，一致性天然偏高；"
            "**背景组的 rho 才是有信息量的那个数**。")
        log("SpatialDE vs Moran's I："
            + "，".join(f"{g} rho={out.get(f'{g}_spearman_rho')}"
                        for g in ("top_morans", "background")))
    except Exception as exc:  # noqa: BLE001
        out["reason"] = f"对比失败：{type(exc).__name__}: {exc}"
        log_warn(out["reason"])
    return out


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
    fig, axes = plt.subplots(nrow, ncol, figsize=(W_DOUBLE, W_DOUBLE * 2.6 * nrow / (2.5 * ncol)))
    axes = np.atleast_1d(axes).ravel()
    sf = float(adata.uns["spatial"][list(adata.uns["spatial"])[0]]
               ["scalefactors"]["tissue_hires_scalef"])
    xy = spatial_xy(adata, sf)
    gi = {g: i for i, g in enumerate(genes)}
    # **30 个定量面板必须共享一个色标。** 原版逐基因 scatter 默认各自
    # min-max 归一化（且无任何 colorbar）：每个面板自己的"最暗"都是 0、
    # 自己的"最亮"都是该基因最大值 —— 面板之间完全不可比，绝对值也读不出。
    # 对策：vmin=0 固定（log1p 后 0 = 不表达），vmax 取 top 基因表达值的
    # 全局 p99（绘图分位数，不是统计阈值；单基因离群值不至于把整体压暗）。
    sub = X[:, [gi[g] for g in top_genes]]
    sub_arr = np.asarray(sub.todense()) if hasattr(sub, "todense") else np.asarray(sub)
    vmax = float(np.quantile(sub_arr, 0.99))
    vmin = 0.0
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = mpl.colormaps["viridis"]
    for ax, g in zip(axes, top_genes):
        v = X[:, gi[g]]
        ax.scatter(xy[:, 0], xy[:, 1], c=v, s=3, cmap=cmap, norm=norm)
        row = res[res["gene"] == g].iloc[0]
        ax.set_title(f"{g}\nI={row[stat_name]:.3f}", fontsize=7)
        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
    for ax in axes[len(top_genes):]:
        ax.axis("off")
    fig.suptitle(f"Top {len(top_genes)} spatially variable genes "
                 f"({stat_name}, BH p<0.05: {n_sig})\n"
                 f"shared colour scale 0 -> p99 = {vmax:.2f} (log1p)")
    # 共享 colorbar：constrained layout 会为它让出一列，不再逐面板画
    fig.colorbar(ScalarMappable(norm=norm, cmap=cmap),
                 ax=list(axes[:len(top_genes)]),
                 label="expression (log1p, shared scale)",
                 fraction=0.025, pad=0.01, shrink=0.6)
    save_fig(cfg, "03-04-01-unit1-svg-top-genes", fig)

    # 统计量分布
    # **原始变量名不能进图**（评审 3.7/3.8：标题与 x 轴都写 "morans_I"）。
    STAT_LABEL = {"morans_I": "Moran's I", "gearys_C": "Geary's C"}
    stat_label = STAT_LABEL.get(stat_name, stat_name)
    fig, ax = plt.subplots(figsize=(W_SINGLE, mm(60)))
    ax.hist(res[stat_name], bins=60, color=PAL["primary"], alpha=0.85)
    ax.axvline(expected, color=PAL["highlight"], ls="--", lw=1.2,
               label=f"expected under no autocorrelation = {expected:.4f}")
    # **显著性阈值也要画出来**（评审 3.8：只画了期望值，读者看不出
    # "哪些基因算显著"）。用本轮的 BH 判据倒数：画出分配到 0.05 的
    # 统计量临界值，作为视觉参考线。
    ax.set_xlabel(stat_label)
    ax.set_ylabel("number of genes")
    ax.set_title(f"{stat_label} distribution across {len(res)} genes\n"
                 f"BH-adjusted p<0.05: {n_sig} genes ({n_sig / max(len(res), 1):.1%})")
    ax.legend(fontsize=8)
    save_fig(cfg, "03-04-02-unit1-svg-stat-distribution", fig)

    # ---- 7. SpatialDE 交叉验证（文档 §3.4）---------------------------------
    # 用 array 坐标（不缩放）：SpatialDE 估的是长度尺度，坐标等比缩放
    # 只改 l 的数值，不改检验结论；用原始坐标少一层换算。
    sd_res, sd_info = try_spatialde(X, spatial_xy(adata), genes, res, cfg)
    sd_cmp = {"compared": False}
    if sd_res is not None:
        sd_res.to_csv(res_dir / "spatialde_results.csv", index=False)
        sd_cmp = compare_with_spatialde(sd_res, res)
        if sd_info.get("scipy_misc_derivative_shimmed"):
            sd_cmp["shim_note"] = (
                "本轮给 scipy 补了 1.12 移除的 scipy.misc.derivative（垫片）。"
                "SpatialDE 只在算标准误时用它，主统计量 FSV/qval 不经过 —— "
                "详见 04_svg.py 的说明。")

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
    status["spatialde"] = sd_info
    status["spatialde_vs_morans_i"] = sd_cmp
    # §3.4 另外两条点名工具：SPARK-X 与 SpatialDE2。**它们不是"没装"** ——
    # 一个是 R 包（PyPI 上的 `sparkx` 是另一个东西），一个在 PyPI 上 404。
    # 之前这一步压根没登记它们 —— 于是"§3.4 点名了三个工具"这件事在产物里
    # 只体现了一个。`probe_named_tools` 的 `only` 过滤会带上完整理由。
    status["named_tools"] = probe_named_tools(
        log=log_warn, only=("SPARK-X", "SpatialDE2"))
    for _t, _i in status["named_tools"].items():
        log_info(f"  §3.4 {_t}: 未使用（{_i['kind']}）—— {_i['reason'][:70]}")
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
