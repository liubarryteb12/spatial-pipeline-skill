#!/usr/bin/env python3
"""
08_spatial_trajectory.py — 空间拟时序（空间感知 vs 朴素的定量对比）

**为什么空间数据要单独做轨迹。**
把单细胞的拟时序工具（Monocle / Slingshot / DPT）直接套到空间数据上，
它们只看表达相似度、**不知道 spot 在哪**。结果是一条在空间上破碎的
排序：相邻的 spot 可以拿到相差很远的拟时序值。

本脚本把这件事**量化**出来，而不是声称"我们的方法更好"：
  - `expr_pseudotime`   —— 直接在原始 PCA 嵌入上跑 DPT（= 传统工具的做法）
  - `spatial_pseudotime`—— 先在**空间平滑过的**嵌入上跑 DPT
  两者的空间自相关（Moran's I）一比，改善有多少就是多少。

**这不是发育轨迹。** 空间数据里没有时间信息，"空间拟时序"报的是
沿组织主导轴的**空间排序**。淋巴结的滤泡↔T 细胞区是空间分区，
不是分化先后。产物里必须写清楚，否则一个排序会被读成发育过程。

**没有真正用到的具名工具**（StPedf / SpaceFlow / ISORT / Stereopy-TGPI /
stLearn）都在 status 的 `named_tools_not_used` 里写明原因。
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402

from common import (ensure_dirs, load_config, log_info,  # noqa: E402
                    log_warn, parse_args, probe_named_tools, record_step, save_fig,
                    set_seed,
                    spatial_xy, write_json, W_ONE_HALF, mm,)

# 复用 03 的空间平滑与 04 的 Moran's I —— 不重复实现。
# 目录名以数字开头，不能直接 import，所以按文件路径加载。
def _load(modname: str, fnname: str):
    """按文件路径加载兄弟脚本里的函数。"""
    import importlib.util

    p = Path(__file__).resolve().parent / f"{modname}.py"
    spec = importlib.util.spec_from_file_location(modname, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return getattr(mod, fnname)


N_GENES_ALONG = 300


def run_08_spatial_trajectory(cfg: dict) -> dict:
    ensure_dirs(cfg)
    set_seed(cfg)
    data_dir = Path(cfg["output"]["data_dir"])
    res_dir = Path(cfg["output"]["results_dir"])
    seed = cfg["analysis"]["seed"]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc
    from scipy.stats import spearmanr

    st = cfg.get("spatial_trajectory") or {}
    if not st.get("enabled", True):
        status = {"dataset_id": cfg["dataset_id"], "status": "disabled",
                  "reason": "配置 spatial_trajectory.enabled=false"}
        write_json(res_dir / "spatial_trajectory_status.json", status)
        log_info("空间拟时序已按配置关闭")
        return status

    src = data_dir / "domains.h5ad"
    if not src.exists():
        status = {"dataset_id": cfg["dataset_id"], "status": "not_available",
                  "reason": f"{src.name} 不存在（step 03 未产出）"}
        write_json(res_dir / "spatial_trajectory_status.json", status)
        log_warn(status["reason"])
        return status

    adata = sc.read_h5ad(src)
    n_obs = int(adata.n_obs)
    if n_obs < 50:
        status = {"dataset_id": cfg["dataset_id"], "status": "not_applicable",
                  "reason": f"只有 {n_obs} 个 spot，轨迹分析无从谈起"}
        write_json(res_dir / "spatial_trajectory_status.json", status)
        log_warn(status["reason"])
        return status
    if "X_pca" not in adata.obsm:
        status = {"dataset_id": cfg["dataset_id"], "status": "missing_pca",
                  "reason": "obsm 里没有 X_pca（step 02 未产出）"}
        write_json(res_dir / "spatial_trajectory_status.json", status)
        log_warn(status["reason"])
        return status

    # ---- 1. 空间权重与平滑 ------------------------------------------------
    build_weights = _load("04_svg", "build_weights")
    morans_i = _load("04_svg", "morans_i")
    smooth_embeddings = _load("03_spatial_domains", "smooth_embeddings")
    spatial_neighbor_graph = _load("03_spatial_domains", "spatial_neighbor_graph")

    W, w_info = build_weights(adata, n_neighbors=int(st.get("n_neighbors", 6)),
                              row_standardize=True)
    adj, adj_info = spatial_neighbor_graph(adata, n_neighbors=int(st.get("n_neighbors", 6)))
    log_info(f"空间权重：{w_info}")

    alpha = float(st.get("smoothing_alpha", 0.5))
    xy = spatial_xy(adata)

    # ---- 2. 选根 -----------------------------------------------------------
    # 判据：转录组复杂度（表达基因数）最高的 spot。
    # **这是启发式，不是生物学判据。** 空间数据里没有时间轴，"根"只是
    # 排序的起点；换一个根，整条轴会反向（见 limitations）。
    #
    # **M3：复杂度必须在全基因集上数。** `adata` 是 `domains.h5ad`，
    # 它的 `.X` 是**只含 2000 个 HVG 的 scaled 矩阵**（02_normalize 写盘时
    # 就裁过）—— 于是在 `.X` 上数 `> 0` 得到的是"有多少个 HVG 的
    # z-score 为正"，一个**恒在 2000 附近**的量，几乎不区分 spot。
    # 复杂度要问的是"这个 spot 测到了多少种基因"，那只有 `adata.raw` 有。
    if adata.raw is not None:
        _Xr = adata.raw.X
        _Xr = _Xr.toarray() if hasattr(_Xr, "toarray") else np.asarray(_Xr)
        complexity = np.asarray((_Xr > 0).sum(axis=1)).ravel().astype(float)
        complexity_source = "adata.raw（全基因集，>0 计数）"
    else:
        complexity = np.asarray((adata.X > 0).sum(axis=1)).ravel().astype(float)
        complexity_source = ("adata.X（**只有 HVG 的 scaled 矩阵** —— "
                             "数出来的是『z-score 为正的 HVG 个数』，"
                             "不是『检测到的基因数』）")
        log_warn(f"adata.raw 为空，复杂度只能退回 HVG 矩阵：{complexity_source}")
    root_idx = int(np.argmax(complexity))
    root_record = {
        "method": "auto_max_detected_genes",
        "complexity_source": complexity_source,
        "complexity_range": [round(float(np.min(complexity)), 2),
                             round(float(np.max(complexity)), 2)],
        "root_spot": str(adata.obs_names[root_idx]),
        "root_xy": [round(float(xy[root_idx, 0]), 2), round(float(xy[root_idx, 1]), 2)],
        "reason": ("取表达基因数最多的 spot 作排序起点。**这是启发式**："
                   "空间数据没有时间轴，根只是排序的起点；"
                   "换根会让整条轴反向，方向本身不携带生物学结论"),
    }
    cfg_root = st.get("root_spot")
    if cfg_root is not None:
        names = list(adata.obs_names.astype(str))
        if str(cfg_root) not in names:
            status = {"dataset_id": cfg["dataset_id"], "status": "bad_root",
                      "reason": f"spatial_trajectory.root_spot='{cfg_root}' 不在 spot 列表里"}
            write_json(res_dir / "spatial_trajectory_status.json", status)
            log_warn(status["reason"])
            return status
        root_idx = names.index(str(cfg_root))
        root_record = {"method": "configured", "root_spot": str(cfg_root),
                       "root_xy": [round(float(xy[root_idx, 0]), 2),
                                   round(float(xy[root_idx, 1]), 2)],
                       "reason": "来自配置 spatial_trajectory.root_spot"}

    # ---- 3. 两条拟时序 -----------------------------------------------------
    def dpt_on(rep: np.ndarray, tag: str):
        """在给定嵌入上跑扩散图 + DPT。返回拟时序数组。"""
        tmp = sc.AnnData(X=np.zeros((n_obs, 1), dtype=np.float32))
        tmp.obs_names = adata.obs_names
        tmp.obsm["X_emb"] = np.asarray(rep, dtype=float)
        sc.pp.neighbors(tmp, n_neighbors=int(cfg["analysis"].get("n_neighbors", 15)),
                        use_rep="X_emb", random_state=seed)
        sc.tl.diffmap(tmp, random_state=seed)
        tmp.uns["iroot"] = root_idx
        sc.tl.dpt(tmp)
        pt = tmp.obs["dpt_pseudotime"].astype(float).values
        if not np.isfinite(pt).any():
            raise ValueError(f"{tag}: DPT 全为 nan")
        return np.nan_to_num(pt, nan=float(np.nanmax(pt)) if np.isfinite(pt).any() else 0.0)

    try:
        expr_pt = dpt_on(adata.obsm["X_pca"], "expr")
        log_info(f"朴素拟时序（原始 PCA 嵌入）：范围 {expr_pt.min():.3f}-{expr_pt.max():.3f}")
    except Exception as e:  # noqa: BLE001
        status = {"dataset_id": cfg["dataset_id"], "status": "failed",
                  "reason": f"朴素 DPT 失败: {type(e).__name__}: {e}"}
        write_json(res_dir / "spatial_trajectory_status.json", status)
        log_warn(status["reason"])
        return status

    try:
        smooth = smooth_embeddings(adata, adj, alpha=alpha, key="X_pca",
                                   out_key="X_pca_spatial")
        spatial_pt = dpt_on(smooth, "spatial")
        log_info(f"空间感知拟时序（平滑 alpha={alpha}）："
                 f"范围 {spatial_pt.min():.3f}-{spatial_pt.max():.3f}")
    except Exception as e:  # noqa: BLE001
        status = {"dataset_id": cfg["dataset_id"], "status": "failed",
                  "reason": f"空间感知 DPT 失败: {type(e).__name__}: {e}"}
        write_json(res_dir / "spatial_trajectory_status.json", status)
        log_warn(status["reason"])
        return status

    # ---- 4. 核心量化：空间自相关 ------------------------------------------
    # 这就是文档 §3.6 说的"传统工具直接套用会给出噪声大的结果"的可检验形式。
    I_expr = float(morans_i(expr_pt.reshape(-1, 1), W)[0])
    I_spatial = float(morans_i(spatial_pt.reshape(-1, 1), W)[0])
    expected = -1.0 / (n_obs - 1)

    # 邻居间拟时序差（越小越平滑）—— 与 Moran's I 互为佐证
    def neighbor_gap(v):
        deg = np.asarray(adj.sum(axis=1)).ravel()
        deg[deg == 0] = 1.0
        nb_mean = (sp.diags(1.0 / deg) @ adj @ v.reshape(-1, 1)).ravel()
        return float(np.mean(np.abs(v - nb_mean)))

    gap_expr = neighbor_gap(expr_pt)
    gap_spatial = neighbor_gap(spatial_pt)
    rho_two = float(spearmanr(expr_pt, spatial_pt).correlation)

    log_info(f"Moran's I（拟时序的空间自相关）：朴素 {I_expr:.4f} → 空间感知 {I_spatial:.4f}"
             f"（期望值 {expected:.4f}）")
    log_info(f"邻居间拟时序平均绝对差：朴素 {gap_expr:.4f} → 空间感知 {gap_spatial:.4f}")

    # ---- 5. 与空间轴的相关（梯度方向）--------------------------------------
    grad_rows = []
    for axis_name, axis_vals in (("x", xy[:, 0]), ("y", xy[:, 1])):
        for tag, v in (("expr", expr_pt), ("spatial", spatial_pt)):
            r = float(spearmanr(v, axis_vals).correlation)
            grad_rows.append({"pseudotime": tag, "spatial_axis": axis_name,
                              "spearman_rho": round(r, 4)})
    pd.DataFrame(grad_rows).to_csv(res_dir / "spatial_trajectory_gradient.csv", index=False)

    # ---- 6. 沿空间拟时序变化的基因 -----------------------------------------
    gene_rows = []
    genes_along_error = None
    try:
        use = adata.raw.to_adata() if adata.raw is not None else adata
        X = use.X
        X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
        X = X.astype(np.float32)
        expressed = (X > 0).sum(axis=0) >= max(10, int(0.01 * n_obs))
        idx = np.flatnonzero(expressed)
        Xs = X[:, idx]
        # 抽列一次算完，避免在 36601 基因上重复扫
        rho_pt = np.nan_to_num(np.array([spearmanr(Xs[:, k], spatial_pt).correlation
                                         for k in range(Xs.shape[1])]))
        I_g = morans_i(Xs, W)          # 每个基因自己的空间自相关
        order = np.argsort(-np.abs(rho_pt))[:N_GENES_ALONG]
        names_g = [use.var_names[i] for i in idx]
        gene_rows = [{
            "gene": names_g[i],
            "rho_with_spatial_pseudotime": round(float(rho_pt[i]), 4),
            "morans_I": round(float(I_g[i]), 4),
            "direction": "increases" if rho_pt[i] > 0 else "decreases",
        } for i in order]
        pd.DataFrame(gene_rows).to_csv(res_dir / "spatial_trajectory_genes.csv", index=False)
        log_info(f"沿空间拟时序变化的基因：表达基因 {len(idx)} 个，报前 {len(gene_rows)} 个")
    except Exception as e:  # noqa: BLE001
        # **S3：不再静默吞掉。** 整块原先只 `log_warn`，而
        # `spatial_trajectory_genes.csv` 在 main_analysis 里**没有消费者** ——
        # "基因没算出来"和"算出来了"在验收层看起来一样。
        # 现在：① 失败原因进 status（`genes_along_pseudotime`）；
        # ② 落一个**只有表头**的 CSV，让消费方能区分"文件不存在"与"空结果"。
        genes_along_error = f"{type(e).__name__}: {e}"
        log_warn(f"沿空间拟时序的基因分析失败: {genes_along_error}")
        pd.DataFrame(columns=["gene", "rho_with_spatial_pseudotime",
                              "morans_I", "direction"]).to_csv(
            res_dir / "spatial_trajectory_genes.csv", index=False)

    # ---- 7. 每 spot 落盘 ---------------------------------------------------
    out = pd.DataFrame({
        "spot": adata.obs_names.astype(str),
        "x": xy[:, 0].round(2), "y": xy[:, 1].round(2),
        "expr_pseudotime": expr_pt, "spatial_pseudotime": spatial_pt,
    })
    if "domain" in adata.obs.columns:
        out["domain"] = adata.obs["domain"].astype(str).values
    out.to_csv(res_dir / "spatial_pseudotime.csv", index=False)

    # **单图原则拆分（D-006）**：两个组合图 -> 5 张独立单图（P8 空间拟时序
    # 方法对照链）。unit1/2/3 = A/|A-B| 对照结构（unit3 依赖前两张，
    # 组内阅读顺序 1->2->3）；unit4/5 = 分布对照 + 相关散点（证据链）。
    # Moran's I 各写进各图标题。
    fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(60)))
    s0 = ax.scatter(xy[:, 0], xy[:, 1], c=expr_pt, s=7, cmap="viridis")
    ax.set_title(f"Expression-only pseudotime\nMoran's I = {I_expr:.3f}")
    fig.colorbar(s0, ax=ax, label="pseudotime",
                 fraction=0.046, pad=0.02, shrink=0.8)
    ax.set_xlabel("x (fullres px)"); ax.set_ylabel("y (fullres px)")
    ax.set_aspect("equal")
    save_fig(cfg, "03-08-01-unit1-pseudotime-expr-only", fig)

    fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(60)))
    s1 = ax.scatter(xy[:, 0], xy[:, 1], c=spatial_pt, s=7, cmap="viridis")
    ax.set_title(f"Spatially smoothed pseudotime (alpha={alpha})\nMoran's I = {I_spatial:.3f}")
    fig.colorbar(s1, ax=ax, label="pseudotime",
                 fraction=0.046, pad=0.02, shrink=0.8)
    ax.set_xlabel("x (fullres px)"); ax.set_ylabel("y (fullres px)")
    ax.set_aspect("equal")
    save_fig(cfg, "03-08-01-unit2-pseudotime-smoothed", fig)

    fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(60)))
    s2 = ax.scatter(xy[:, 0], xy[:, 1],
                    c=np.abs(expr_pt - spatial_pt), s=7, cmap="magma")
    ax.set_title("|difference| expression-only vs smoothed")
    fig.colorbar(s2, ax=ax, label="|delta pseudotime|",
                 fraction=0.046, pad=0.02, shrink=0.8)
    ax.set_xlabel("x (fullres px)"); ax.set_ylabel("y (fullres px)")
    ax.set_aspect("equal")
    save_fig(cfg, "03-08-01-unit3-pseudotime-difference", fig)

    fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(64)))
    ax.hist(expr_pt, bins=40, alpha=0.65, label="expression-only", color="#B2182B")
    ax.hist(spatial_pt, bins=40, alpha=0.65, label="spatially-smoothed", color="#2166AC")
    ax.set_xlabel("pseudotime"); ax.set_ylabel("n spots")
    ax.set_title("Pseudotime distributions")
    fig.legend(fontsize=8, ncol=1, loc="outside right center")
    save_fig(cfg, "03-08-02-unit1-pseudotime-distributions", fig)

    fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(64)))
    ax.scatter(expr_pt, spatial_pt, s=5, alpha=0.4, color="#444444")
    ax.set_xlabel("expression-only pseudotime")
    ax.set_ylabel("spatially-smoothed pseudotime")
    ax.set_title(f"Spearman rho = {rho_two:+.3f}")
    save_fig(cfg, "03-08-02-unit2-pseudotime-scatter", fig)

    # ---- 9. 状态 -----------------------------------------------------------
    improved = I_spatial > I_expr
    status = {
        "dataset_id": cfg["dataset_id"],
        "status": "ok",
        "n_spots": n_obs,
        "smoothing_alpha": alpha,
        "spatial_weights": w_info,
        "spatial_graph": adj_info,
        "root_selection": root_record,
        "morans_I_expression_only": round(I_expr, 4),
        "morans_I_spatially_smoothed": round(I_spatial, 4),
        "morans_I_expected_no_autocorrelation": round(expected, 4),
        "morans_I_gain": round(I_spatial - I_expr, 4),
        "spatial_smoothing_improves_coherence": bool(improved),
        # **§3.6 点名的 StPedf / SpaceFlow / ISORT / Stereopy-TGPI / stLearn
        # 一个都没跑。** 上面跑的是扩散图 + DPT 的内置实现 ——
        # 不记这条，读者会以为空间轨迹用的是文档点名的框架。
        "named_tools": probe_named_tools(
            log=log_warn,
            only=("StPedf", "SpaceFlow", "ISORT", "Stereopy-TGPI", "stLearn")),
        "neighbor_gap_expression_only": round(gap_expr, 4),
        "neighbor_gap_spatially_smoothed": round(gap_spatial, 4),
        "spearman_between_two": round(rho_two, 4),
        "gradient": grad_rows,
        "n_genes_along": len(gene_rows),
        "genes_along_pseudotime": {
            "status": "failed" if genes_along_error else "ok",
            "reason": genes_along_error,
            "n_genes_reported": len(gene_rows),
            "csv": "spatial_trajectory_genes.csv",
        },
        "top_genes": gene_rows[:20],
        "method": ("在空间平滑过的 PCA 嵌入上跑扩散图 + DPT，与直接在原始嵌入上"
                   "跑的结果**定量对比**（Moran's I 与邻居间差值）。"
                   "平滑复用 03_spatial_domains 的 smooth_embeddings，"
                   "空间自相关复用 04_svg 的 morans_i"),
        "named_tools_not_used": {
            "StPedf": ("空间邻近嵌入 + 密度自适应融合。未使用：不在 PyPI 上，"
                       "源码需从论文补充材料取，无法在托管 runner 上稳定安装"),
            "SpaceFlow": ("生成平滑伪空间图。未使用：仅 GitHub 发布，无 PyPI wheel，"
                          "依赖链未固定，装失败会让整条流水线红"),
            "ISORT": ("空间 RNA 速率。未使用：**需要 spliced/unspliced 两套计数**，"
                      "本流水线的输入只有一套"),
            "Stereopy-TGPI": ("多样本时间序列。未使用：需要多个时间点样本，"
                              "本数据集是单个时间点"),
            "stLearn": ("整合组织学图像做轨迹。未使用：需要 H&E 图像像素级处理"
                        "（本流水线有图但未做图像特征提取），且它假定的是"
                        "有形态学梯度的组织"),
        },
        "limitations": [
            "**这不是发育轨迹。** 空间数据没有时间信息，这里报的是沿组织"
            "主导轴的**空间排序**。淋巴结的滤泡↔T 细胞区是空间分区，"
            "不是分化先后",
            f"根是启发式选的（{root_record['method']}）：换根会让整条轴反向，"
            "方向本身不携带生物学结论",
            "空间平滑会**同时抹掉真实的局部结构**：alpha 越大越平滑，"
            "但小的空间异质性也会被平均掉。alpha 是配置项，应做敏感性检查",
            "Moran's I 依赖空间权重矩阵的定义（本流水线用 kNN + 距离阈值）；"
            "换邻居定义数值会变，但两组之间的**相对**高低通常稳定",
            "**空间自相关高不等于轨迹正确**：一条纯空间梯度（如组织边缘"
            "到中心）也能给出很高的 Moran's I，而它与分化无关",
        ],
    }
    write_json(res_dir / "spatial_trajectory_status.json", status)
    log_info(f"空间拟时序完成：Moran's I {I_expr:.4f} → {I_spatial:.4f}"
             f"（{'改善' if improved else '未改善'}）")
    return status


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        run_08_spatial_trajectory(cfg)
        record_step(cfg, "spatial_trajectory", "ok", time.time() - t0)
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "spatial_trajectory", "failed", time.time() - t0, message=str(e))
        raise
