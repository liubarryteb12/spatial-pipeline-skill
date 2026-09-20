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
                    log_warn, parse_args, pkg_version, probe_named_tools,
                    record_cross_language, record_decision, record_step,
                    save_fig, set_seed,
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


# ---------------------------------------------------------------------------
# §3.3 点名的 cell2location
# ---------------------------------------------------------------------------
# **它和 NNLS 用同一份参考，所以两者可比** —— 这是选它做交叉验证而不是
# 另找一个工具的理由。差别在似然：NNLS 是最小二乘（高斯噪声），
# cell2location 是负二项（计数噪声 + 每个基因自己的离散度）。
#
# **装得上 ≠ 跑得了。** 上一轮把 cell2location 记成 `deps`（"GPU 导向、
# CPU CI 吃不消"），读它的接口就知道归因错了：真正卡住的是**数据** ——
# `Cell2location(...)` 必须给 `cell_state_df`（每个细胞类型的全转录组
# 后验表达谱），而本流水线默认配的是 `reference: builtin`（marker 签名），
# 根本没有这份参考。所以状态是 `needs_reference`，不是装不上。
#
# **它不会在本仓库的默认 CI 里执行**（默认配置没有 scRNA 参考）。
# 这条路径的意义是：一旦配上 `reference: h5ad` + `celltype_key`，
# 打开 `deconvolution.use_cell2location: true` 就能直接跑，
# 而不是只有一句"没做"。
#
# 代价（写进产物）：CPU 上参考模型 + 映射模型两轮变分推断，
# 实测数量级是**几十分钟**，所以它必须是**可选**的，不能进默认路径。

CELL2LOCATION_DEFAULTS = {
    "max_epochs_reference": 250,
    "max_epochs_mapping": 3000,
    "n_cells_per_location": 8.0,
    "detection_alpha": 20.0,
    "min_cells_per_type": 10,
    "num_samples": 1000,
}


def try_cell2location(C, genes, ref_a, ct_key: str, cfg: dict,
                      obs_names=None, log=log_info):
    """跑 cell2location（文档 §3.3），返回 `(props | None, types | None, info)`。

    参数都是**已经对齐到同一组共同基因**的：
      - `C`     : spot x n_common 的**原始计数**（不是 library-size 归一化后的
                  —— cell2location 的负二项似然要计数）
      - `genes` : `C` 的列名，顺序一致
      - `ref_a` : 参考 AnnData，`X` 是同一组基因的原始计数，
                  `obs[ct_key]` 是细胞类型标签
      - `obs_names` : spot 的 barcode。**必须传** —— 结果按行对应 spot，
                  用整数序号当索引的话，后面和 `adata.obs_names` 对齐时
                  只能靠"顺序没变"这个假设，而那个假设没有任何东西保证。

    只有在 `deconvolution.reference == "h5ad"`（真的有参考）且
    `deconvolution.use_cell2location: true` 时才会被调用。

    **任何异常都吞掉并写进 info** —— 点名工具跑不起来不该让整步失败，
    但必须让人看见（AGENTS 规则 4）。返回的比例按行归一化，
    与 NNLS 那条路的 `props` 口径一致，才可比。
    """
    info = {"attempted": True, "status": None, "reason": "", "section": "§3.3"}
    try:
        from cell2location.models import Cell2location, RegressionModel
    except Exception as exc:  # noqa: BLE001
        info["status"] = "package_missing"
        info["reason"] = (f"cell2location 未安装（{type(exc).__name__}: {exc}）—— "
                          f"PyPI 有真包，但依赖 scvi-tools + torch + pyro-ppl + "
                          f"opencv-python；已回退到同参考的 NNLS")
        log_warn(f"  §3.3 cell2location 未跑：{info['reason'][:110]}")
        return None, None, info

    info["version"] = pkg_version("cell2location")
    c2l = dict(CELL2LOCATION_DEFAULTS)
    c2l.update(((cfg.get("deconvolution") or {}).get("cell2location") or {}))
    info["params"] = {k: c2l[k] for k in CELL2LOCATION_DEFAULTS}

    try:
        # ---- 1. 参考模型：估每个细胞类型的全转录组表达谱 -------------------
        #
        # 与 NNLS 那条路的关键区别：这里用**全部共同基因**，不是只有 marker。
        ref_a = ref_a.copy()
        ref_a.obs["_ct"] = ref_a.obs[ct_key].astype(str)
        ref_a.obs["_batch"] = "ref"
        vc = ref_a.obs["_ct"].value_counts()
        keep = vc[vc >= int(c2l["min_cells_per_type"])].index
        dropped = sorted(set(vc.index) - set(keep))
        if dropped:
            log_warn(f"  cell2location 参考里丢掉细胞数 <{c2l['min_cells_per_type']} "
                     f"的类型：{dropped}（回归会因样本太少而不稳）")
        info["dropped_celltypes"] = dropped
        ref_a = ref_a[ref_a.obs["_ct"].isin(keep)].copy()
        if ref_a.obs["_ct"].nunique() < 2:
            info["status"] = "reference_too_thin"
            info["reason"] = (f"参考里只剩 {ref_a.obs['_ct'].nunique()} 种类型"
                              f"（阈值 min_cells_per_type="
                              f"{c2l['min_cells_per_type']}）—— 无法建回归模型")
            log_warn(f"  §3.3 cell2location 未跑：{info['reason']}")
            return None, None, info

        RegressionModel.setup_anndata(ref_a, labels_key="_ct", batch_key="_batch")
        mod_ref = RegressionModel(ref_a)
        mod_ref.train(max_epochs=int(c2l["max_epochs_reference"]),
                      train_size=1, accelerator="cpu", devices=1)
        ref_a = mod_ref.export_posterior(
            ref_a, sample_kwargs={"num_samples": int(c2l["num_samples"]),
                                  "batch_size": min(2048, ref_a.n_obs)})
        key_varm = "q05_cell_abundance_w_sf"
        if key_varm not in ref_a.varm:
            info["status"] = "failed"
            info["reason"] = (f"参考模型跑完了但 varm 里没有 {key_varm} —— "
                              f"cell2location 的导出接口变了")
            log_warn(f"  §3.3 cell2location 未跑：{info['reason']}")
            return None, None, info
        inf_aver = pd.DataFrame(ref_a.varm[key_varm])
        # 列名：varm 里没有类型名，顺序由 setup_anndata 的 labels 类别决定。
        # **不猜顺序** —— 拿类别名对上，对不上就直接报错。
        cats = [str(c) for c in ref_a.obs["_ct"].astype("category").cat.categories]
        if len(cats) != inf_aver.shape[1]:
            info["status"] = "failed"
            info["reason"] = (f"参考谱有 {inf_aver.shape[1]} 列，但标签有 "
                              f"{len(cats)} 种 —— 列名对不上，不猜")
            log_warn(f"  §3.3 cell2location 未跑：{info['reason']}")
            return None, None, info
        inf_aver.columns = cats
        log_info(f"  参考模型完成：{inf_aver.shape[0]} 基因 x "
                 f"{inf_aver.shape[1]} 种类型")

        # ---- 2. 映射模型：把谱映射到每个 spot ------------------------------
        sp_a = sc.AnnData(X=np.asarray(C, dtype=np.float32))
        sp_a.var_names = list(genes)
        if obs_names is not None:
            # **必须是真 barcode，不能是 0..n-1。** 结果按行对应 spot，
            # 用序号当索引的话，后面和 adata.obs_names 对齐只能靠
            # "顺序没变"这个假设 —— 而那个假设没有任何东西保证。
            sp_a.obs_names = [str(b) for b in obs_names]
        sp_a.obs["_batch"] = "spatial"
        Cell2location.setup_anndata(sp_a, batch_key="_batch")
        mod_map = Cell2location(
            sp_a, cell_state_df=inf_aver,
            N_cells_per_location=float(c2l["n_cells_per_location"]),
            detection_alpha=float(c2l["detection_alpha"]))
        mod_map.train(max_epochs=int(c2l["max_epochs_mapping"]),
                      train_size=1, accelerator="cpu", devices=1)
        sp_a = mod_map.export_posterior(
            sp_a, sample_kwargs={"num_samples": int(c2l["num_samples"]),
                                 "batch_size": sp_a.n_obs})
        if key_varm not in sp_a.obsm:
            info["status"] = "failed"
            info["reason"] = f"映射模型跑完了但 obsm 里没有 {key_varm}"
            log_warn(f"  §3.3 cell2location 未跑：{info['reason']}")
            return None, None, info
        q = sp_a.obsm[key_varm]
        # **行顺序不能猜。** obsm 里放 DataFrame 时，AnnData 允许它的索引
        # 与 obs_names 不同 —— 那意味着行可能已经重排过。这时按位置取
        # `q.values` 会把比例贴到错误的 spot 上，而且**不会报错**。
        if hasattr(q, "index") and not q.index.equals(pd.Index(sp_a.obs_names)):
            info["status"] = "failed"
            info["reason"] = ("obsm 里比例矩阵的索引与 obs_names 不一致 —— "
                              "行顺序无法确定，不猜（按位置取会静默贴错 spot）")
            log_warn(f"  §3.3 cell2location 未跑：{info['reason']}")
            return None, None, info
        types = [str(c).replace(f"{key_varm}_", "") for c in q.columns]
        props = np.clip(np.asarray(q.values, dtype=float), 0, None)
        rs = props.sum(axis=1, keepdims=True)
        rs[rs == 0] = 1.0
        props = props / rs
        info["status"] = "ok"
        info["n_types"] = len(types)
        info["n_genes_in_reference"] = int(inf_aver.shape[0])
        info["mean_props"] = {t: round(float(props[:, i].mean()), 4)
                              for i, t in enumerate(types)}
        log_info(f"  §3.3 cell2location 完成：{props.shape[0]} spot x "
                 f"{len(types)} 种类型")
        return props, types, info
    except Exception as exc:  # noqa: BLE001
        info["status"] = "failed"
        info["reason"] = f"{type(exc).__name__}: {exc}"
        log_warn(f"  §3.3 cell2location 跑失败：{info['reason']}")
        return None, None, info


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

    # §3.3 点名工具 cell2location 的落地状态。**默认值就是 needs_reference** ——
    # 只有 h5ad 分支（真的有带标签的参考）才可能把它改成别的。
    # 放在这里而不是分支里，是为了让"没有参考"这条路径也留下痕迹：
    # 键消失和"状态是 needs_reference"看起来完全不同。
    c2l_info = {
        "attempted": False, "status": "needs_reference", "section": "§3.3",
        "reason": (f"deconvolution.reference = '{mode}' —— **没有带细胞类型标签的 "
                   f"scRNA-seq 参考**。cell2location 的 `Cell2location(...)` 必须"
                   f"给 `cell_state_df`（每个类型的全转录组表达谱），"
                   f"没有参考就跑不了。**这不是装不上** —— PyPI 上有真包"
                   f"（0.1.5），缺的是数据。配上 `reference: h5ad` + "
                   f"`celltype_key` + `use_cell2location: true` 就会走它"),
    }

    if mode == "none":
        # §0.2：`cross_language` 会是空的 —— **必须解释为什么空**，
        # 否则空数组和"忘了记"长得一模一样（两个 Python 仓库都踩过）。
        record_decision(
            cfg, "cross_language",
            "§0.2 跨部分交接：这一轮有没有 Part 2 → Part 3 的单细胞参考？",
            "**没有。** 配置里 `deconvolution.reference = none` —— "
            "本轮不做解卷积，因此也不需要 Part 2 的参考 h5ad",
            evidence="deconvolution_status.json 的 status=not_done",
        )
        status.update({"status": "not_done",
                       "cell2location": c2l_info,
                       "reason": "配置 deconvolution.reference=none —— "
                                 "**没有参考就不能解卷积**，不用代理指标冒充"})
        write_json(res_dir / "deconvolution_status.json", status)
        log_warn(status["reason"])
        return status

    if mode != "h5ad":
        # `builtin` 走的是本仓库自己的 `assets/reference_signatures.yml`
        # （marker 签名），**不是另一个部分的产物**。所以这一轮
        # `cross_language` 为空是设计如此，不是遗漏。
        record_decision(
            cfg, "cross_language",
            "§0.2 跨部分交接：这一轮有没有 Part 2 → Part 3 的单细胞参考？",
            f"**没有。** 配置里 `deconvolution.reference = '{mode}'` —— "
            f"参考来自本仓库的 `assets/reference_signatures.yml`"
            f"（marker 签名），不是 Part 2 导出的 h5ad。"
            f"**所以 `cross_language` 为空是设计如此，不是遗漏**；"
            f"真正的解卷积需要 `reference: h5ad` + `celltype_key`",
            evidence=(f"deconvolution.reference='{mode}'；"
                      f"参考文件 assets/reference_signatures.yml 在本仓库内"),
        )

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

        # ---- §0.2 跨部分交接：这一份参考是 **Part 2 的产物** ----------------
        #
        # `deconvolution.reference: h5ad` 不是"本仓库内的一个文件" ——
        # 它按约定是 Part 2 导出的带细胞类型标签的单细胞 h5ad
        # （Part 2 产物 `clustered.h5ad`，落在它的 data_dir 下，
        # obs 里有细胞类型列）。
        # 所以这一处**必须**走 `record_cross_language`：
        # 不记的话，清单里 `cross_language` 是空的，而空数组和
        # "这一轮没有跨部分交接"长得一模一样。
        #
        # `before` 记 Part 2 那边有什么，`after` 记本步骤真正用上了什么，
        # `lost_fields` 记**没有跟过来**的东西。这三样缺一不可：
        # 只记 "n_genes" 的话，读者不知道参考里的 embedding、
        # 其他 obs 列、以及 HVG 限制都丢了。
        _ref_layers = sorted(ref.layers.keys())
        _ref_obsm = sorted(ref.obsm.keys())
        # **计数层是硬要求，不是偏好。** NNLS 解的是线性混合
        # `spot = Σ 比例 × 签名`，log 变换破坏线性关系（见文件头）。
        # 所以参考谱必须来自原始计数 —— 拿 log 值平均出来的"签名"
        # 解出的比例没有意义，**而它看起来仍然像一组比例**。
        _has_counts = "counts" in ref.layers
        if not _has_counts:
            log_warn("参考 h5ad 里没有 layers['counts']，退回用 .X —— "
                     "如果 .X 是 log 后的值，解出的比例没有意义。"
                     "Part 2 的 clustered.h5ad 带 counts 层，优先用它")
        R = ref.layers["counts"] if _has_counts else ref.X
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
                    "is_deconvolution": True,
                    # 记下**用的是计数层还是退回 .X** —— 不记的话，
                    # "退回 log 值"这个错误前提会一路静默到比例上
                    "used_counts_layer": _has_counts,
                    "matrix_source": "layers['counts']" if _has_counts else ".X",
                    "n_cells_in_reference": int(ref.n_obs),
                    "n_genes_in_reference": int(ref.n_vars)}
        log_info(f"用 scRNA-seq 参考：{len(types)} 种类型，{len(common)} 个共同基因"
                 + ("" if _has_counts else "（**无 counts 层，退回 .X**）"))

        # 交接日志。`record_cross_language` 在 common.py 里定义，
        # 记的是"数据跨过仓库边界时前后各是什么形状、丢了什么"。
        record_cross_language(
            cfg,
            src=f"Part 2 (Python) {Path(ref_path).name}",
            dst="05_deconvolution 参考谱 S（类型 × 共同基因）",
            fmt="h5ad",
            before={
                "n_cells": int(ref.n_obs),
                "n_genes": int(ref.n_vars),
                "obs_columns": sorted(ref.obs.columns.tolist()),
                "layers": _ref_layers,
                "obsm": _ref_obsm,
            },
            after={
                "n_types": len(types),
                "n_common_genes": len(common),
                "n_spots_deconvolved": int(C.shape[0]),
                "celltype_key": ct_key,
                "used_counts_layer": _has_counts,
            },
            lost=(
                # 参考里没跟过来的东西。**逐条列出来**，不要只说"丢了一些字段"
                [c for c in sorted(ref.obs.columns.tolist()) if c != ct_key]
                + [f"layers['{k}']" for k in _ref_layers if k != "counts"]
                + [f"obsm['{k}']" for k in _ref_obsm]
                # 基因交集是**交集**，不是全部：参考若是 HVG 子集，
                # 这里就丢掉了非 HVG 的基因。这一条最容易被忽略，
                # 而它直接决定参考谱覆盖多少基因
                + ([f"非共同基因 {int(ref.n_vars) - len(common)} 个"
                    f"（参考 {int(ref.n_vars)} → 共同 {len(common)}）"]
                   if int(ref.n_vars) != len(common) else [])
                + ([] if _has_counts
                   else ["**参考没有 counts 层，退回 .X** —— "
                         "若 .X 是 log 值，参考谱不是计数谱"])
            ),
            note=("Part 2 → Part 3 的单细胞参考交接。Part 2 侧的产出物是 "
                  "`clustered.h5ad`（落在它的 data_dir 下，obs 里有"
                  "细胞类型列、带 `layers['counts']`）。"
                  "**本仓库只用到计数矩阵与细胞类型标签两样**，"
                  "参考里的 embedding / 其他 obs 列 / 非 counts 的层"
                  "都不参与计算 —— 已逐条列进 lost_fields。"),
        )

        scale = C.sum(axis=1, keepdims=True)
        scale[scale == 0] = 1.0
        Cn = C / scale * 1e4
        log_info(f"解卷积 {Cn.shape[0]} 个 spot x {len(types)} 种类型（NNLS）")
        props, errors = deconvolve_nnls(Cn, S)
        mean_prop = pd.Series(props.mean(axis=0), index=types).sort_values(
            ascending=False)
        log_info("平均组成: " + ", ".join(
            f"{k}={v:.3f}" for k, v in mean_prop.items()))

        # ---- §3.3 点名工具 cell2location（同一份参考，所以可比）-------------
        #
        # **默认不开。** CPU 上两轮变分推断是几十分钟量级，放进默认路径会
        # 让每一轮 CI 都为它付时间。打开后它和 NNLS 用**同一份参考**，
        # 差异只来自似然（负二项 vs 最小二乘），这才是有信息量的对照。
        if dec.get("use_cell2location"):
            ref_sub = sc.AnnData(X=Rc.astype(np.float32))
            ref_sub.obs[ct_key] = ref.obs[ct_key].astype(str).values
            ref_sub.var_names = list(common)
            props_c2l, types_c2l, c2l_info = try_cell2location(
                C, common, ref_sub, ct_key, cfg, obs_names=adata.obs_names)
            if props_c2l is not None:
                pd.DataFrame(props_c2l, index=adata.obs_names,
                             columns=types_c2l).to_csv(
                    res_dir / "deconvolution_proportions_cell2location.csv")
                if set(types_c2l) == set(types):
                    from scipy.stats import spearmanr
                    rho = {}
                    for t in types:
                        r, _ = spearmanr(props[:, types.index(t)],
                                         props_c2l[:, types_c2l.index(t)])
                        rho[t] = None if np.isnan(r) else round(float(r), 4)
                    vals = [v for v in rho.values() if v is not None]
                    c2l_info["vs_nnls"] = {
                        "per_celltype_spearman": rho,
                        "mean_spearman": (round(float(np.mean(vals)), 4)
                                          if vals else None),
                        "note": ("按细胞类型在 spot 上算 Spearman —— "
                                 "**只比空间趋势，不比绝对值**（两条路的"
                                 "比例尺度本来就不同）"),
                    }
                    log_info(f"  cell2location vs NNLS：平均 Spearman "
                             f"{c2l_info['vs_nnls']['mean_spearman']}")
                else:
                    c2l_info["vs_nnls"] = {
                        "compared": False,
                        "reason": (f"细胞类型集合不同（NNLS {len(types)} 种 / "
                                   f"cell2location {len(types_c2l)} 种）"
                                   f"—— 不做逐类型相关"),
                    }
        else:
            c2l_info = {
                "attempted": False, "status": "disabled", "section": "§3.3",
                "reason": ("配置里没有 `deconvolution.use_cell2location: true`。"
                           "cell2location 在 CPU 上要跑两轮变分推断"
                           "（参考模型 + 映射模型，几十分钟量级），"
                           "**默认不放进路径**；已用同一份参考的 NNLS 作主结果"),
            }
            log_info("  §3.3 cell2location：未启用（配置未打开）—— 主结果是同参考的 NNLS")
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
    # **12 个 colorbar 必须共享量程**（评审 3.5：原先各自 vmin=0 + 自动 vmax，
    # 实测量程 0–0.3 / 0–0.6 / 0–0.15 / 0–0.20 四种，面板之间深浅不可比 ——
    # 同一个"偏黄"在两张图里是两倍差异）。统一 vmax = 全部展示类型
    # 联合 p99，并给每个 colorbar 一个说明性 label。
    vmax_shared = float(np.quantile(
        np.concatenate([prop_df[ct].values for ct in show]), 0.99))
    for ax, ct in zip(axes, show):
        s = ax.scatter(xy[:, 0], xy[:, 1], c=prop_df[ct].values, s=4,
                       cmap="magma", vmin=0, vmax=vmax_shared)
        ax.set_title(f"{ct}\nmean={mean_prop[ct]:.3f}", fontsize=8)
        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
        # **colorbar 不加长 label**：12 个面板各挂一条竖排文字会互相重叠
        # （实测第一版把 "Fibroblastic_reticular_cell (shared scale)" 挤进
        # 相邻面板）。共享量程这件事由 suptitle 统一说明，刻度数字足够。
        fig.colorbar(s, ax=ax, shrink=0.75, pad=0.02, fraction=0.046)
    for ax in axes[len(show):]:
        ax.axis("off")
    fig.suptitle("Deconvolved composition — spatial trends only, "
                 "absolute values are not cell fractions\n"
                 f"all panels share one colour scale (0 – {vmax_shared:.2f})")
    save_fig(cfg, "03-05-01-unit1-deconvolution-spatial", fig)

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
        # **x 轴按数值感知排序，不按字符串。** 域标签是 "0".."12" 这类数字串时，
        # sorted() 给 0,1,10,11,12,2,...—— 刻度顺序与读图直觉相反。
        # 能转 int 的按 int 排；混入非数字标签时回退字典序（顺序稳定即可）。
        dom_order = sorted(by_dom.index,
                           key=lambda d: (0, int(d)) if str(d).isdigit() else (1, str(d)))
        by_dom = by_dom.loc[dom_order]
        ax.set_xticks(range(len(by_dom)))
        ax.set_xticklabels(by_dom.index, fontsize=8)
        ax.set_xlabel("spatial domain")
        # **ylabel 必须说清这不是细胞比例。** 本仓库无参考时走 marker 打分法
        # （domain_status.json 里 is_deconvolution: false），数值只是"该域内
        # marker 打分的相对构成"。姊妹图 deconvolution_spatial 有免责句，
        # 这张原先一句都没有 —— 评审认定它才是真正会误导读者的那张。
        ax.set_ylabel("mean marker score per domain (relative;\n"
                      "NOT deconvolved cell fractions)")
        ax.set_title("Marker-score composition per spatial domain\n"
                     "(no reference available: scores, not real deconvolution)")
        ax.legend(fontsize=6, ncol=2, loc="center left", bbox_to_anchor=(1.0, 0.5))
        save_fig(cfg, "03-05-02-unit1-deconvolution-by-domain", fig)
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
        save_fig(cfg, "03-05-03-unit1-deconvolution-error-map", fig)

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
        # §3.3 的落地边界 —— 与 method 并列，不藏起来
        ("**§3.3 点名的 RCTD 没跑**：`spacexr` 是 R 包，PyPI 上无同名包，"
         "本仓库 CI 不装 R。**cell2location 的落地状态见 "
         "`cell2location` 字段** —— 它是 `needs_reference`（装得上，缺的是"
         "带标签的 scRNA 参考），不是装不上。"),
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
        # §3.3 点名工具真的跑了 / 启用了但没跑成 —— 两种都要说出来
        if c2l_info.get("status") == "ok":
            rho = (c2l_info.get("vs_nnls") or {}).get("mean_spearman")
            method_desc += ("；**同时跑了 §3.3 点名的 cell2location**"
                            f"（v{c2l_info.get('version') or '未知'}），"
                            f"与 NNLS 的平均 Spearman rho = {rho}")
            limits.append(
                "cell2location 与 NNLS 的差异只来自似然（负二项 vs 最小二乘），"
                "因为两者用的是**同一份参考** —— 所以 rho 低说明的是"
                "『似然的选择改变了空间趋势』，不是『哪个更对』。"
                "**真值未知，两者都不能当基准**")
        elif dec.get("use_cell2location"):
            limits.append(
                f"**cell2location 已启用但没跑成**（status="
                f"{c2l_info.get('status')}）：{c2l_info.get('reason')}"
                f" —— 上面的结果是同参考的 NNLS，不是 cell2location")
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
        # **§3.3 点名的两个工具各自的落地状态。**
        # RCTD 是 R 包（结构上跑不了）；cell2location 是 needs_reference
        # （装得上，缺数据）—— 两者不能混成一句话。
        "cell2location": c2l_info,
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
