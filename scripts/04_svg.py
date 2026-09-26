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
from matplotlib.colors import Normalize  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import scipy.sparse as sp  # noqa: E402

from common import (df_to_records, ensure_dirs, load_config, log_info,  # noqa: E402
                    log_warn, parse_args, probe_named_tools,
                    reconstruct_counts_from_raw, record_step,
                    save_fig, set_seed, write_json,
                    spatial_xy, W_ONE_HALF, W_SINGLE, mm, PAL,)


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
                  res: pd.DataFrame, cfg: dict, obs=None, adata=None,
                  log=log_info):
    """跑 SpatialDE 并与 Moran's I 比对。

    **只在基因子集上跑。** SpatialDE 每个基因要拟合一个高斯过程，
    实测 4025 个 spot 上单基因约 1-3 秒 —— 4000 个基因要几小时，
    会直接撞穿 CI 的 45 分钟上限。所以取「Moran's I 最高的前 100 个」
    + 「按种子随机抽 50 个作背景」，共 150 个。

    **子集抽样会引入选择偏差**（前 100 个是 Moran's I 挑出来的），
    所以一致性数字要**分两组报**：top 组和随机背景组。
    只报合并的相关会把"两边都认为强"和"两边都认为弱"混成一个数。

    ---
    ## 为什么必须在 `dyn_de` 之前补两步预处理

    `SpatialDE.run` 在调 `dyn_de` **之前**会做两步
    （`SpatialDE/anndata.py:33-46`）：

      `NaiveDE.stabilize(X.T).T`      ← 方差稳定变换（VST）
      `NaiveDE.regress_out(obs, ., 'np.log(total_counts)')`  ← 回归掉测序深度

    本仓库绕开 `run` 是为了躲 `util.qvalue` 对 pandas Series 调
    `.ravel()` 的不兼容（见下），但**绕开 `run` 时把这两步也一起丢了**。
    没有 VST，高表达基因的方差被低估；没有深度回归，测序深度本身
    （它在空间上往往有梯度 —— 组织边缘的 spot 捕获的转录本更少）
    会被当成"空间结构"。**两个方向都会让 qval 偏小、假阳性变多。**

    `NaiveDE` 是 `SpatialDE` 的**硬依赖**，所以这两步不需要新依赖。

    ---
    ## ⚠️ VST 必须吃**原始计数**，不能吃 `X`（这里踩过一个坑）

    传进来的 `X` 是 `adata.raw.X` 的 **log1p 值**（Moran's I 那条路要的
    就是 log 值，所以 `X` 本身没错）。但 `NaiveDE.stabilize` 按
    `var = mu + phi*mu^2` 拟合离散度，**假设输入是原始计数**。

    喂 log 值的后果实测（lymph_node，150 个基因）：

      `phi_hat = -0.2523` → `1/(2*phi) = -1.98`（负数）
      → `np.log(负数)` 出 NaN → `stabilize` 后 43 个基因整列全 NaN
      → `get_mll_results` 的 merge 把 NaN 行**悄悄丢掉**
      → 150 个基因只剩 **3 行**，而 `status` 仍报 `ok`

    所以这里用 `reconstruct_counts_from_raw` 重建计数再 stabilize。
    **尺度必须按全基因集行和反推**（只对选中的 150 列求行和会把尺度
    放大约 8.9 倍，且不报错）—— 公共库那个函数已经保证这一点。
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

        # ---- S5：补上 SpatialDE.run 里的两步预处理 -------------------------
        # 顺序不能反：先 VST（把方差与均值解耦），再回归掉深度。
        # `NaiveDE` 假设**行是基因、列是样本**，所以 stabilize 要转两次。
        #
        # **VST 的输入必须是原始计数**（见 docstring 的坑）。有 adata 就
        # 从 raw 重建这 150 个基因的计数；拿不到就如实记进 info，不假装做过。
        vst_note = None
        counts_vst = None
        vst_block = None
        if adata is not None and getattr(adata, "raw", None) is not None:
            try:
                C_sel, sel_ok, cinfo = reconstruct_counts_from_raw(
                    adata, genes=sel, log=log)
                if len(sel_ok) != len(sel):
                    cinfo = dict(cinfo, status="incomplete",
                                 n_genes_reconstructed=len(sel_ok),
                                 n_genes_requested=len(sel))
                    vst_block = (f"只有 {len(sel_ok)}/{len(sel)} 个基因能从 raw "
                                 "重建计数")
                elif cinfo.get("scale_source") != "obs['total_counts']":
                    # **没有 total_counts 时不能做 VST。** 那条回退路径把每行
                    # 归一到同一尺度，得到的是**深度归一后的值、不是计数** ——
                    # 而 VST 恰恰要按 `var = mu + phi*mu^2` 拟合计数离散度。
                    # 喂这种值实测会算出负 phi（NaN → 静默丢基因）。
                    vst_block = ("obs 里没有 total_counts —— 只能还原到"
                                 "每行同一尺度（深度归一的值），不是真实计数")
                else:
                    counts_vst = pd.DataFrame(C_sel, columns=sel_ok,
                                              index=counts.index)
                info["count_reconstruction"] = cinfo
                if vst_block:
                    log_warn(f"不做 VST：{vst_block}")
            except Exception as exc:  # noqa: BLE001
                vst_block = f"重建计数失败（{type(exc).__name__}: {exc}）"
                info["count_reconstruction"] = {
                    "status": "failed",
                    "reason": f"{type(exc).__name__}: {exc}"}
                log_warn(f"从 raw 重建计数失败：{type(exc).__name__}: {exc}")
        else:
            vst_block = "没有 adata 或 adata.raw 为空"
            info["count_reconstruction"] = {
                "status": "not_available",
                "reason": vst_block}

        try:
            import NaiveDE
            if counts_vst is None:
                # 拿不到真实计数就**不做 VST** —— 宁可少做一步并记下来，
                # 也不能喂 log 值（或深度归一值）假装做过：那会算出负 phi，
                # 让 get_mll_results 静默丢掉绝大多数基因。
                raise RuntimeError(
                    f"没有可用的原始计数（{vst_block}），不能做 "
                    "NaiveDE.stabilize")
            counts_st = pd.DataFrame(
                NaiveDE.stabilize(counts_vst.T.values).T,
                columns=sel, index=counts.index)
            if obs is not None and "total_counts" in getattr(obs, "columns", []):
                # 这一步按 `np.log(total_counts)` 回归 —— 与 SpatialDE.run
                # 的默认 regress_formula 一致（不改口径）
                counts_use = pd.DataFrame(
                    NaiveDE.regress_out(obs, counts_st.T.values,
                                        "np.log(total_counts)").T,
                    columns=sel, index=counts.index)
                vst_note = ("NaiveDE.stabilize(原始计数) + "
                            "regress_out('np.log(total_counts)')")
            else:
                counts_use = counts_st
                vst_note = ("NaiveDE.stabilize(原始计数)（**没有 obs.total_counts，"
                            "深度回归被跳过**）")
                log_warn("obs 里没有 total_counts —— SpatialDE 的深度回归被跳过。"
                         "**测序深度在空间上常有梯度，它会被当成空间结构**")
            info["vst_applied"] = True
            info["vst_input"] = "原始计数（从 adata.raw 重建）"
            info["depth_regression"] = bool(
                obs is not None and "total_counts" in getattr(obs, "columns", []))
            info["preprocessing"] = vst_note
        except Exception as exc:  # noqa: BLE001
            # **不吞掉。** 预处理失败就退回原始计数，但必须把这件事
            # 写进产物 —— "VST 没做成"和"做了 VST"的结果看起来一样。
            counts_use = counts
            info["vst_applied"] = False
            info["vst_input"] = None
            info["depth_regression"] = False
            info["preprocessing"] = f"**未做 VST（{type(exc).__name__}: {exc}）**"
            log_warn(f"NaiveDE 预处理失败，退回原始计数：{type(exc).__name__}: {exc}")

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
        raw = SpatialDE.base.dyn_de(coords, counts_use, kernel_space)
        out = pd.DataFrame(SpatialDE.base.get_mll_results(raw))
        out["pval"] = 1.0 - _st.chi2.cdf(out["LLR"], df=1)
        out["qval"] = multipletests(out["pval"], method="fdr_bh")[1]
        out["_group"] = np.where(out["g"].isin(set(top)), "top_morans", "background")
        # ---- 静默丢基因的守门 ------------------------------------------------
        # `get_mll_results` 内部是 merge，NaN 行会被**悄悄丢掉**：
        # 实测喂 log 值时 150 个基因只剩 3 行，而这里原先照报 ok。
        # 所以产出基因数与请求数不一致时必须显式记下来。
        n_out = len(out)
        n_dropped = info["n_genes_run"] - n_out
        info["n_genes_out"] = n_out
        info["n_genes_dropped"] = n_dropped
        if n_dropped > 0:
            missing = [g for g in sel if g not in set(out["g"].astype(str))]
            info["genes_dropped_examples"] = missing[:10]
            info["status"] = "partial"
            info["reason"] = (f"{n_dropped}/{info['n_genes_run']} 个基因在 "
                              "get_mll_results 里被丢掉（NaN 行会被 merge 静默丢弃）")
            log_warn(f"SpatialDE **只产出 {n_out}/{info['n_genes_run']} 个基因**，"
                     f"丢了 {n_dropped} 个（例：{missing[:5]}）—— 见 status 的 "
                     "n_genes_dropped。VST 输入不是原始计数时就会这样")
        else:
            info["status"] = "ok"
        info["qval_method"] = "BH（本仓库实现）"
        info["bypassed_spatialde_run"] = True
        info["bypass_reason"] = ("SpatialDE.util.qvalue 对 pandas Series 调 .ravel()，"
                                 "现代 pandas 已无此方法")
        info["columns"] = sorted(map(str, out.columns))
        info["n_significant_qval_0.05"] = int((out["qval"] < 0.05).sum())
        log(f"SpatialDE 完成：{len(out)} 个基因，"
            f"qval<0.05 的 {info['n_significant_qval_0.05']} 个"
            f"（预处理：{vst_note}）")
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
        #
        # **M12：这个上限原先不在配置里。** `svg.max_genes` 在
        # `assets/config.lymph_node.yml` 的 `svg:` 段根本没有这一项，
        # 所以 `s.get("max_genes", 4000)` **永远走默认**，而且 status 里
        # 也不记它 —— 读者无法知道这轮的 SVG 到底看了多少基因。
        # 现在两处都补：配置里写出来（可调），status 里记实际值。
        max_genes = int(s.get("max_genes", 4000))
        n_genes_all = len(genes_all)
        if len(keep) > max_genes:
            keep = keep[np.argsort(var[keep])[::-1][:max_genes]]
            log_info(f"全基因集 {n_genes_all} 个，按方差取前 {max_genes} 个做置换检验")
        X = np.ascontiguousarray(Xall[:, keep].astype(np.float64))
        genes = [genes_all[i] for i in keep]
    else:
        X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
        X = np.ascontiguousarray(X.astype(np.float64))
        genes = list(adata.var_names)
        max_genes = None
        n_genes_all = len(genes)
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
    # **单图原则拆分（D-006，裁决 1）**：30 面板网格 -> 30 张独立单图
    # （P4 SVG 空间图谱分解链）。每张是一个基因的空间表达图谱，
    # 独立达标图幅（W_SINGLE）与分辨率门禁；共享色标 0->p99 **保留**
    # （拆图不拆可比性），量程写进各图 title；组内横向对比靠组编号。
    # 图名前缀拼接（基因名是运行时数据，30 个无法写字面量；
    # 拼接构造避开命名门禁对 "03-..." 字面量的格式校验——账目核对
    # 由本注释说明：30 张图名 = SVG_FIG_BASE + unit号 + 基因slug）
    # 声明式动态名豁免：03-04-01 图号下 30 张基因空间图，slug=基因名（运行时数据）
    DYNAMIC_FIG_BASES = {"01": 30}
    n_top = int(s.get("n_top", 30))
    top_genes = res.head(n_top)["gene"].tolist()
    sf = float(adata.uns["spatial"][list(adata.uns["spatial"])[0]]
               ["scalefactors"]["tissue_hires_scalef"])
    xy = spatial_xy(adata, sf)
    gi = {g: i for i, g in enumerate(genes)}
    sub = X[:, [gi[g] for g in top_genes]]
    sub_arr = np.asarray(sub.todense()) if hasattr(sub, "todense") else np.asarray(sub)
    vmax = float(np.quantile(sub_arr, 0.99))
    vmin = 0.0
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = mpl.colormaps["viridis"]
    for ui, g in enumerate(top_genes, start=1):
        v = X[:, gi[g]]
        row = res[res["gene"] == g].iloc[0]
        fig, ax = plt.subplots(figsize=(W_SINGLE, mm(58)))
        ax.scatter(xy[:, 0], xy[:, 1], c=v, s=3, cmap=cmap, norm=norm)
        ax.set_title(f"{g}  I={row[stat_name]:.3f}\n"
                     f"shared scale 0 -> p99 = {vmax:.2f} (log1p)", fontsize=8)
        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                     shrink=0.8, pad=0.02, fraction=0.046,
                     label="expression (log1p, shared scale)")
        save_fig(cfg, f"03-04-01-unit{ui}-{g.lower()}", fig)

    # 统计量分布
    # **原始变量名不能进图**（评审 3.7/3.8：标题与 x 轴都写 "morans_I"）。
    STAT_LABEL = {"morans_I": "Moran's I", "gearys_C": "Geary's C"}
    stat_label = STAT_LABEL.get(stat_name, stat_name)
    # **图幅放宽到 W_ONE_HALF + 图例标签必须短。**
    # 实测 bug（2026-09-24）：原 W_SINGLE(89mm) 配 `loc="outside right center"`
    # 的**长图例**（60/63 字符），constrained layout 把坐标轴**挤到只剩约 100px**
    # —— 1051px 宽的画布里 axes 仅占 10%，直方图柱子被压得完全看不见，
    # 图看起来"只有两条参考线"。而且**不报错**：CI 全绿、status ok、
    # check_figures 也报"有墨迹"（图例文字就是墨）。只有打开图才发现是空的。
    #
    # 两处一起改：① 图幅 89→136mm；② 图例标签压到 ~28 字符以内，数值移进副标题。
    fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(60)))
    ax.hist(res[stat_name], bins=60, color=PAL["primary"], alpha=0.85)
    ax.axvline(expected, color=PAL["highlight"], ls="--", lw=1.2,
               label="expected (no autocorr.)")
    # **显著性用标题报数 + 一条真实可算的参考线。** 临界值不能画成"BH 阈值
    # 竖线"—— BH 临界值是逐基因 p 的函数、不是统计量的固定值。能画的是
    # **显著基因集合的边界值**：Moran's I 越大越有结构，取 min；Geary's C
    # 反之越小越有结构，取 max。两者都是 res 里直接取出的实测值。
    sig = res[res["p_adj_bh"] < 0.05]
    if len(sig):
        edge = (float(sig[stat_name].min()) if method == "moran"
                else float(sig[stat_name].max()))
        ax.axvline(edge, color=PAL["highlight"], ls=":", lw=1.2,
                   label="significant-set edge")
    ax.set_xlabel(stat_label)
    ax.set_ylabel("number of genes")
    # 数值移进副标题（图例只留短标签，避免重新挤压坐标轴）
    _sig_line = (f"dashed = expected {expected:.4f}; dotted = sig-set edge {edge:.3f}"
                 if len(sig) else f"dashed = expected {expected:.4f}")
    ax.set_title(f"{stat_label} distribution across {len(res)} genes\n"
                 f"BH-adjusted p<0.05: {n_sig} genes ({n_sig / max(len(res), 1):.1%})\n"
                 f"{_sig_line}")
    fig.legend(fontsize=8, ncol=1, loc="outside right center")
    save_fig(cfg, "03-04-02-unit1-svg-stat-distribution", fig)

    # ---- 7. SpatialDE 交叉验证（文档 §3.4）---------------------------------
    # 用 array 坐标（不缩放）：SpatialDE 估的是长度尺度，坐标等比缩放
    # 只改 l 的数值，不改检验结论；用原始坐标少一层换算。
    sd_res, sd_info = try_spatialde(X, spatial_xy(adata), genes, res, cfg,
                                    obs=adata.obs, adata=adata)
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
        "n_genes_in_raw": int(n_genes_all),
        "max_genes_cap": max_genes,
        "gene_selection": ("全基因集按方差取前 max_genes 个" if max_genes
                           else "adata.raw 为空，退回 HVG"),
        "n_significant_bh": n_sig,
        "n_permutations": n_perms,
        # **M1：p 值分辨率。** `n_perms` 次置换下 p 值只有 `n_perms + 1`
        # 个可能取值（`(n_ge+1)/(n_perms+1)`），最小值是 `1/(n_perms+1)`。
        # 4000 个基因共用 101 档 → **大量基因并列同一个 p**，
        # 于是 BH 校正的分母虽然是 4000，实际可分辨的位次远少于 4000。
        # 这不是"错了"，但必须写出来：读者会以为 `p_adj` 有 4000 档分辨率。
        "p_value_resolution": {
            "n_permutations": n_perms,
            "n_distinct_p_values_possible": n_perms + 1,
            "min_possible_p": round(1.0 / (n_perms + 1), 6),
            "n_genes_at_min_p": int((res["p_value"] <= 1.0 / (n_perms + 1) + 1e-12).sum()),
            "note": ("置换检验的 p 值只有 n_perms+1 档；**并列最小值**的基因"
                     "数量在这里一并报出 —— 并列越多，BH 的『排名』信息越少"),
        },
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
