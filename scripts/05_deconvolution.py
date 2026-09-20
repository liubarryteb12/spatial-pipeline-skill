#!/usr/bin/env python3
"""
05_deconvolution.py — spot 组成解卷积（NNLS）

**这是 marker 签名 + NNLS，不是 RCTD/cell2location/DestVI。** 差别必须说清楚：

  RCTD（R）/ cell2location / DestVI（scvi-tools）用**全转录组**的
  scRNA-seq 参考，对每个细胞类型的完整表达谱（所有基因的均值与离散度）
  建模，并用负二项/泊松似然拟合。

  本流水线用**精简 marker 签名 + 非负最小二乘**：只用每种类型的
  一小组标志基因，解一个非负线性混合方程。

后果（全部写进产物的 limitations）：
  - 只能区分**签名差异明显**的类型
  - 比例的**绝对值不可信**，只有**相对的空间趋势**有意义
  - 签名重叠的类型会被合并或互相混淆（如淋巴结的 B_naive 与
    B_germinal_center 共享 MS4A1/CD79A）
  - 不在签名里的类型会被强行分配给最接近的

**为什么用计数而不是 log 后的值。** NNLS 解的是线性混合
`spot 表达 = Σ_k 比例_k × 签名_k`。log 变换破坏了线性关系，
用 log 值解出的"比例"没有意义 —— 而结果看起来仍然像比例。

**重建误差是必报的诊断量。** 误差大的 spot 说明"用这套签名解释不了
这个 spot"，可能是有未列入的类型、或者是技术坏点。
不报误差的话，一个完全解释不了的 spot 也会得到一个看起来正常的组成。
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
import yaml  # noqa: E402
from scipy.optimize import nnls  # noqa: E402

from common import (df_to_records, ensure_dirs, load_config, log_info,  # noqa: E402
                    log_warn, parse_args, probe_named_tools, record_step, save_fig,
                    set_seed,
                    spot_radius_plot_units, write_json, spatial_xy, W_DOUBLE, W_ONE_HALF, W_SINGLE, mm,)


def load_signature(cfg: dict):
    p = Path(__file__).resolve().parent.parent / "assets" / "reference_signatures.yml"
    with open(p, encoding="utf-8") as fh:
        doc = yaml.safe_load(fh)
    key = (cfg.get("deconvolution") or {}).get("signature_set")
    sigs = doc.get("signatures", {})
    if key not in sigs:
        raise KeyError(f"signature_set='{key}' 不在 {p.name} 里；"
                       f"可选: {', '.join(sorted(sigs))}")
    return key, sigs[key]


def marker_score_composition(adata, sig, C: np.ndarray, genes: list) -> tuple:
    """
    **marker 打分法**（不是解卷积）—— 无外部参考时唯一站得住的做法。

    ---
    ## 为什么不用 NNLS（这里踩过一个很典型的坑）

    第一版用 NNLS 解 `y ≈ S^T w`，其中 S 的每一行是一个细胞类型的参考谱。
    没有外部参考时，我把参考谱构造成"只在 marker 基因位置取该类型 marker
    的平均表达，其余位置为 0"。**结果是完全错的**：

      12 种类型的组成全部约等于 0.085 ≈ 1/12，重建误差中位 0.64，
      3836/4025 个 spot 超过误差阈值。

    原因是模型在说"**非 marker 基因应该不表达**"（所有参考谱在那里都是 0），
    而实际上非 marker 基因表达得很正常。残差被这些基因主导，NNLS 只能
    给一个无意义的均匀解 —— **而它看起来仍然像一组比例**。

    更根本的问题是：**没有外部参考时，marker-only NNLS 不是适定问题。**
    方程 `y[markers] = Σ_k w_k · s_k[markers]` 里的 `s_k`（每种类型在
    marker 上的参考表达谱）是未知的。把它从目标数据本身估出来是循环论证。

    ---
    ## 这里做的

    对每种类型算一个**相对富集分数**：该类型 marker 的平均表达，
    减去该 spot 全部基因的平均表达（校正测序深度），再除以该类型 marker
    在全部 spot 上的标准差（把不同类型放到可比的尺度上）。
    然后按 spot 做 softmax 式归一化，让各类型的分数可比。

    **这不给出细胞比例。** 它给出的是"这个 spot 里哪类细胞的标志基因
    更活跃"，可以看空间趋势，不能读成百分比。
    """
    gi = {g: i for i, g in enumerate(genes)}
    types, score_cols, coverage = [], [], {}
    for ct, d in sig["celltypes"].items():
        present = [g for g in d.get("markers", []) if g in gi]
        missing = [g for g in d.get("markers", []) if g not in gi]
        coverage[ct] = {"n_markers": len(d.get("markers", [])),
                        "n_present": len(present), "missing": missing}
        if len(present) < 3:
            continue
        sub = C[:, [gi[g] for g in present]]
        # 该 spot 该类型 marker 的平均表达，减去该 spot 的整体平均（深度校正）
        s = sub.mean(axis=1) - C.mean(axis=1)
        sd = float(s.std())
        score_cols.append(s / sd if sd > 1e-9 else s)
        types.append(ct)

    if not types:
        return None, None, coverage, None

    raw = np.column_stack(score_cols)
    # 按 spot 归一化成"相对权重"：先移到非负，再 L1 归一化
    shifted = raw - raw.min(axis=1, keepdims=True)
    tot = shifted.sum(axis=1, keepdims=True)
    tot[tot <= 0] = 1.0
    props = shifted / tot
    return props, types, coverage, raw


def deconvolve_nnls(counts: np.ndarray, S: np.ndarray) -> tuple:
    """
    对每个 spot 解非负最小二乘。

    **每个 spot 单独解** —— 这是 NNLS 解卷积的定义。
    返回比例矩阵与每个 spot 的重建误差。

    比例做 L1 归一化（和为 1），这样数值可解释为"组成"。
    但注意：归一化后的值**不是细胞比例**，只是相对权重。
    """
    n_spots = counts.shape[0]
    n_types = S.shape[0]
    props = np.zeros((n_spots, n_types))
    errors = np.zeros(n_spots)
    for i in range(n_spots):
        y = counts[i]
        if y.sum() <= 0:
            errors[i] = np.nan
            continue
        # 解 y ≈ S^T @ w, w >= 0
        w, resid = nnls(S.T, y)
        s = w.sum()
        if s > 0:
            props[i] = w / s
        errors[i] = float(np.sqrt(resid ** 2 / max((y ** 2).sum(), 1e-12)))
    return props, errors


def run_05_deconvolution(cfg: dict) -> dict:
    ensure_dirs(cfg)
    set_seed(cfg)
    data_dir = Path(cfg["output"]["data_dir"])
    res_dir = Path(cfg["output"]["results_dir"])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dec = cfg.get("deconvolution") or {}
    mode = dec.get("reference", "builtin")
    status = {"dataset_id": cfg["dataset_id"]}

    if mode == "none":
        status.update({"status": "not_done",
                       "reason": "配置 deconvolution.reference=none —— "
                                 "**没有参考就不能解卷积**，不用代理指标冒充"})
        write_json(res_dir / "deconvolution_status.json", status)
        log_warn(status["reason"])
        return status

    adata = sc.read_h5ad(data_dir / "domains.h5ad")

    # ---- 1. 取计数 ----------------------------------------------------------
    if "counts" not in adata.layers:
        raise RuntimeError(
            "adata.layers['counts'] 不存在 —— 解卷积必须用原始计数。\n"
            "  为什么必须停: NNLS 解的是线性混合 spot = Σ 比例 × 签名；\n"
            "  log 变换破坏线性关系，用 log 值解出的『比例』没有意义 ——\n"
            "  而结果看起来仍然像比例。")
    C = adata.layers["counts"]
    C = C.toarray() if sp.issparse(C) else np.asarray(C)
    C = C.astype(np.float64)
    genes = list(adata.var_names)
    log_info(f"计数矩阵: {C.shape[0]} spot x {C.shape[1]} 基因（HVG）")

    # ---- 2. 参考 ------------------------------------------------------------
    if mode == "builtin":
        sig_key, sig = load_signature(cfg)
        log_info(f"用内置 marker 签名 '{sig_key}'：{len(sig['celltypes'])} 种类型")
        log_info("**marker 打分法，不是解卷积** —— 无外部参考时 NNLS 不适定"
                 "（见 marker_score_composition 的说明）")
        props, types, coverage, raw_scores = marker_score_composition(
            adata, sig, C, genes)
        if props is None:
            status.update({"status": "no_usable_signature",
                           "reason": "没有任何细胞类型的 marker 基因在数据里出现 >=3 个",
                           "coverage": coverage})
            write_json(res_dir / "deconvolution_status.json", status)
            log_warn(status["reason"])
            return status
        errors = None  # 打分法没有重建误差
        ref_desc = {"kind": "builtin_marker_score", "signature_set": sig_key,
                    "signature_note": sig.get("note", ""),
                    "is_deconvolution": False}
        log_info(f"marker 打分完成：{len(types)} 种类型")
        mean_prop = pd.Series(props.mean(axis=0), index=types).sort_values(
            ascending=False)
        log_info("平均相对权重: " + ", ".join(
            f"{k}={v:.3f}" for k, v in mean_prop.items()))

    elif mode == "h5ad":
        ref_path = dec.get("reference_h5ad")
        ct_key = dec.get("celltype_key")
        if not ref_path or not ct_key:
            raise ValueError("deconvolution.reference=h5ad 时必须给 "
                             "reference_h5ad 与 celltype_key")
        ref = sc.read_h5ad(ref_path)
        if ct_key not in ref.obs.columns:
            raise KeyError(f"参考数据的 obs 里没有 '{ct_key}'")
        R = ref.layers["counts"] if "counts" in ref.layers else ref.X
        R = R.toarray() if sp.issparse(R) else np.asarray(R)
        R = R.astype(np.float64)
        common = [g for g in genes if g in set(ref.var_names)]
        if len(common) < 50:
            raise RuntimeError(f"参考与目标只有 {len(common)} 个共同基因，太少")
        ri = {g: i for i, g in enumerate(ref.var_names)}
        gi = {g: i for i, g in enumerate(genes)}
        Rc = R[:, [ri[g] for g in common]]
        C = C[:, [gi[g] for g in common]]
        genes = common
        types = sorted(ref.obs[ct_key].astype(str).unique())
        # **有参考时才是真正的解卷积**：参考谱覆盖全部共同基因，
        # 不是只在 marker 位置非零 —— 那正是第一版错的原因。
        S = np.vstack([Rc[(ref.obs[ct_key].astype(str) == t).values].mean(axis=0)
                       for t in types])
        coverage = {t: {"n_markers": None, "n_present": None, "missing": []}
                    for t in types}
        ref_desc = {"kind": "scrna_reference", "path": str(ref_path),
                    "celltype_key": ct_key, "n_common_genes": len(common),
                    "is_deconvolution": True}
        log_info(f"用 scRNA-seq 参考：{len(types)} 种类型，{len(common)} 个共同基因")

        scale = C.sum(axis=1, keepdims=True)
        scale[scale == 0] = 1.0
        Cn = C / scale * 1e4
        log_info(f"解卷积 {Cn.shape[0]} 个 spot x {len(types)} 种类型（NNLS）")
        props, errors = deconvolve_nnls(Cn, S)
        mean_prop = pd.Series(props.mean(axis=0), index=types).sort_values(
            ascending=False)
        log_info("平均组成: " + ", ".join(
            f"{k}={v:.3f}" for k, v in mean_prop.items()))
    else:
        raise ValueError(f"不支持的 deconvolution.reference: {mode}")

    # ---- 3. 后处理 ----------------------------------------------------------
    min_prop = float(dec.get("min_proportion", 0.01))
    props_f = props.copy()
    props_f[props_f < min_prop] = 0.0
    rs = props_f.sum(axis=1, keepdims=True)
    rs[rs == 0] = 1.0
    props_f = props_f / rs

    max_err = float(dec.get("max_reconstruction_error", 0.5))
    if errors is not None:
        n_unreliable = int(np.nansum(errors > max_err))
        log_info(f"重建误差: 中位 {np.nanmedian(errors):.4f}，"
                 f"最大 {np.nanmax(errors):.4f}；"
                 f"{n_unreliable} 个 spot 误差 >{max_err}")
    else:
        n_unreliable = 0
        log_info("打分法没有重建误差这个诊断量 —— 见 limitations")

    prop_df = pd.DataFrame(props_f, index=adata.obs_names, columns=types)
    prop_df.to_csv(res_dir / "deconvolution_proportions.csv")
    if errors is not None:
        pd.DataFrame({"reconstruction_error": errors,
                      "reliable": errors <= max_err},
                     index=adata.obs_names).to_csv(
            res_dir / "deconvolution_errors.csv")

    # ---- 4. 空间分布 --------------------------------------------------------
    show = mean_prop.head(12).index.tolist()
    ncol = 4
    nrow = int(np.ceil(len(show) / ncol))
    sf = float(adata.uns["spatial"][list(adata.uns["spatial"])[0]]
               ["scalefactors"]["tissue_hires_scalef"])
    xy = spatial_xy(adata, sf)
    fig, axes = plt.subplots(nrow, ncol, figsize=(W_DOUBLE, W_DOUBLE * 3.1 * nrow / (3.0 * ncol)))
    axes = np.atleast_1d(axes).ravel()
    for ax, ct in zip(axes, show):
        s = ax.scatter(xy[:, 0], xy[:, 1], c=prop_df[ct].values, s=4,
                       cmap="magma", vmin=0)
        ax.set_title(f"{ct}\nmean={mean_prop[ct]:.3f}", fontsize=8)
        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(s, ax=ax, shrink=0.75)
    for ax in axes[len(show):]:
        ax.axis("off")
    fig.suptitle("Deconvolved composition — spatial trends only, "
                 "absolute values are not cell fractions")
    save_fig(cfg, "deconvolution_spatial", fig)

    # 组成堆叠（按某个域聚合，看域之间的组成差异）
    if "domain" in adata.obs.columns:
        by_dom = prop_df.groupby(adata.obs["domain"].astype(str).values).mean()
        # 宽度夹在 [单栏半, 双栏]，避免域数多时画出装不进一页的图
        fig, ax = plt.subplots(figsize=(min(W_DOUBLE, max(W_ONE_HALF, 0.55 * len(by_dom) + 2.5)), mm(64)))
        bottom = np.zeros(len(by_dom))
        cmap = plt.get_cmap("tab20")
        for i, ct in enumerate(mean_prop.index):
            ax.bar(range(len(by_dom)), by_dom[ct].values, bottom=bottom,
                   label=ct, color=cmap(i % 20))
            bottom += by_dom[ct].values
        ax.set_xticks(range(len(by_dom)))
        ax.set_xticklabels(by_dom.index, fontsize=8)
        ax.set_xlabel("spatial domain"); ax.set_ylabel("mean relative weight")
        ax.set_title("Composition per spatial domain")
        ax.legend(fontsize=6, ncol=2, loc="center left", bbox_to_anchor=(1.0, 0.5))
        save_fig(cfg, "deconvolution_by_domain", fig)
        by_dom.to_csv(res_dir / "deconvolution_by_domain.csv")

    # ---- 5. 重建误差的空间分布（只有解卷积才有）----------------------------
    if errors is not None:
        fig, ax = plt.subplots(figsize=(W_SINGLE, mm(76)))
        s = ax.scatter(xy[:, 0], xy[:, 1], c=errors, s=5, cmap="magma")
        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"Reconstruction error (median {np.nanmedian(errors):.3f})\n"
                     f"high = signature set cannot explain this spot")
        fig.colorbar(s, ax=ax, shrink=0.8, label="relative error")
        save_fig(cfg, "deconvolution_error_map", fig)

    # ---- 6. 落盘 ------------------------------------------------------------
    is_deconv = bool(ref_desc.get("is_deconvolution"))
    err_block = None
    if errors is not None:
        err_block = {
            "median": round(float(np.nanmedian(errors)), 4),
            "p95": round(float(np.nanpercentile(errors, 95)), 4),
            "max": round(float(np.nanmax(errors)), 4),
            "n_unreliable": n_unreliable,
            "frac_unreliable": round(n_unreliable / len(errors), 4),
        }

    common_limits = [
        "Visium 的 spot 含 1-10 个细胞，所以即使组成准确，"
        "空间分辨率也受限于 spot 尺寸",
        # §3.3 的缺口 —— 与 method 并列，不藏起来
        "**§3.3 点名的 RCTD / cell2location 一个都没跑。** RCTD（`spacexr`）"
        "是 R 包，PyPI 上无同名包；cell2location 有 PyPI 真包，但依赖 "
        "scvi-tools + torch + pyro-ppl + opencv-python，CPU CI 上跑不动。"
        "逐条理由见 `named_tools` 字段。",
    ]
    if is_deconv:
        method_desc = "**NNLS 解卷积**，参考谱来自 scRNA-seq 参考（覆盖全部共同基因）"
        limits = common_limits + [
            "NNLS 没有似然模型，无法给出每个估计的不确定性"
            "（RCTD / cell2location 有）",
            "比例之和被强制为 1，所以一个类型的比例被高估会压低其他类型",
            f"重建误差中位 {err_block['median']:.3f}；误差大的 spot 说明"
            "这套参考解释不了它 —— 看 deconvolution_error_map",
        ]
    else:
        method_desc = ("**marker 打分法（相对富集），不是解卷积。** "
                       "无外部参考时 NNLS 不适定 —— 第一版试过，"
                       "得到 12 种类型全是 0.085≈1/12 的均匀解、"
                       "重建误差中位 0.64，见 marker_score_composition 的说明")
        limits = common_limits + [
            "**不给出细胞比例。** 给出的是『这个 spot 里哪类细胞的标志基因"
            "更活跃』，可以看空间趋势，**不能读成百分比**",
            "各类型分数的可比性依赖 marker 数量和表达水平，"
            "marker 多的类型有系统性优势",
            "打分法**没有重建误差**这个诊断量 —— 无法判断"
            "『这套签名解释不了某个 spot』",
            "签名重叠的类型会互相混淆（淋巴结 B_naive 与 B_germinal_center "
            "共享 MS4A1/CD79A）",
            "要用真正的解卷积：在 config 里设 deconvolution.reference: h5ad "
            "并给一个 scRNA-seq 参考",
        ]

    status = {
        "dataset_id": cfg["dataset_id"],
        "status": "ok",
        "reference": ref_desc,
        "is_deconvolution": is_deconv,
        "n_spots": int(props_f.shape[0]),
        "n_celltypes": len(types),
        "celltypes": types,
        "signature_coverage": coverage,
        "min_proportion": min_prop,
        "max_reconstruction_error": max_err if errors is not None else None,
        "reconstruction_error": err_block,
        "mean_composition": {k: round(float(v), 5) for k, v in mean_prop.items()},
        "method": method_desc,
        # **§3.3 点名的 RCTD / cell2location 一个都没跑。**
        # 上面的 method 是 marker 打分法（is_deconvolution=False）——
        # 不把这条放在 method 旁边，读者会以为"解卷积做过了"。
        "named_tools": probe_named_tools(log=log_warn, only=("RCTD", "cell2location")),
        "limitations": limits,
    }
    write_json(res_dir / "deconvolution_status.json", status)
    return status


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        run_05_deconvolution(cfg)
        record_step(cfg, "deconvolution", "ok", time.time() - t0)
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "deconvolution", "failed", time.time() - t0, message=str(e))
        raise
