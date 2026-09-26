#!/usr/bin/env python3
"""
07_spatial_communication.py — 空间约束的配体-受体分析

**核心增量：只有空间上接近的细胞才可能通讯。**

普通单细胞流程算配体-受体，隐含假设"样本里所有细胞都可能互相作用" ——
在 5000 个细胞里，A 类细胞和 B 类细胞永远"共现"，所以永远算出信号。
空间数据能加一个真实的约束：**只统计距离在阈值内的 spot 对**。

这里报两类结果：
  1. **空间富集**：把 LR 表达量按"距离内 vs 距离外"分组比较 ——
     同一对 LR，在近邻对上是否比远距离对上更强？
  2. **按域/类型的 LR 强度**：哪些域组合在哪些 LR 上活跃。

**必须报"可用 LR 对数"。** 配体-受体基因可能不在数据里（Visium 的
基因覆盖不全），或者表达太低（零膨胀）。如果 38 对里只有 3 对可用，
那结论建立在 3 对上，必须说清楚 —— 这是 Part 2 实测踩过的坑
（HVG 子集上只有 3/38 可用，全基因集上是 27/38）。

**本工具不做"通讯与否"的二元判定。** 真正的细胞通讯推断需要
CellPhoneDB / CellChat 那样的统计框架（含置换检验与受体复合物建模），
本工具只给"空间约束下的 LR 共表达强度"这个描述性量。
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

from common import (df_to_records, ensure_dirs, finite_round,  # noqa: E402
                    load_config, log_info,
                    log_warn, parse_args, probe_named_tools, record_step, save_fig,
                    set_seed,
                    write_json, spatial_xy, W_ONE_HALF, mm,
                    PAL,)


def _r5(x):
    """四舍五入到 5 位；**非有限值落 `None`**（E-69 同族 Form C）。

    薄壳，逻辑在 `common.finite_round` —— **全仓只有一份收口实现**。
    保留本名字是因为调用点较多，且 `_r5` 是 5 位的语义标记。
    """
    return finite_round(x, 5)


def top_lr_for_plots(lr_df, n: int = 3):
    """从 LR 打分表里选出用于空间图谱的对：**只取 z 有定义的**（E-69）。

    抽成纯函数有两个理由（AGENTS 规则 30.1 / E-63）：

    1. 自检要能调真代码，不能重写一遍判据；
    2. 调用点原先直接 `lr_df[lr_df["z_score"].notna()].head(3)` 然后
       `np.concatenate(...)` —— **全部 z 都算不出来时它返回空表，
       `np.concatenate([])` 抛 `ValueError: need at least one array to
       concatenate`**，「一个量算不出来」的处境变成整步崩溃。

    返回的可能是空表，调用方**必须判 `empty`**。这里不抛错：空是合法状态
    （零模型退化时没有可展示的对），该由调用方决定怎么报。
    """
    return lr_df[lr_df["z_score"].notna()].head(n)


def spatial_z_score(near_mean: float, perm_means: list) -> tuple:
    """由近邻均值与零模型置换样本算 z，**并显式报告「算不出来」**（E-69）。

    返回 `(null_mu, null_sd, z, z_reason)`，四态：

    - 置换样本**全是 nan**（`flat` 为空或 `n_near=0`）：`null_sd=None`、
      `z=None`，原因是「零模型退化」。原实现写
      `null_sd = float(np.nanstd(perm_means)) or 1e-9` —— **`or` 对 nan 无效**：
      `bool(nan)` 是 `True`，所以 `nan or 1e-9` 仍是 `nan`，于是
      `z = (near_mean - null_mu) / nan` 得 nan，四个 `round()` 全部落裸 nan，
      **「零模型抽不出样本」与「z 恰好很小」在产物里长得一样**。
    - 置换标准差**恰好为 0**：`null_sd=0.0`、`z=None`。z 没有定义（不是 0）——
      所有置换给出同一个均值时，分母为 0。
    - `near_mean` 是 nan（阈值内一个 spot 对都没有）：`z=None`。
    - 正常：`z` 是有限浮点。

    **抽成纯函数是为了让自检能调真代码**（AGENTS 规则 30.1 / E-63）。
    """
    null_mu = (float(np.nanmean(perm_means))
               if np.isfinite(perm_means).any() else float("nan"))
    sd_raw = float(np.nanstd(perm_means))
    if not np.isfinite(sd_raw):
        return (null_mu, None, None,
                "零模型的置换样本全部是 nan（`flat` 为空或 `n_near=0`），"
                "`null_sd` 算不出来 —— 这不是「z 很小」，是「这个量在这里"
                "没有定义」，排查方向是看 `n_spot_pairs_within_threshold` "
                "与 `max_distance_um` 是否把阈值设得比最近邻还小")
    if sd_raw == 0.0:
        return (null_mu, 0.0, None,
                "零模型的置换标准差恰好为 0 —— 所有置换给出同一个均值，"
                "z 没有定义（不是 0）")
    if not np.isfinite(near_mean):
        return (null_mu, sd_raw, None,
                "近邻对均值为 nan（`n_near=0`，阈值内一个 spot 对都没有）—— "
                "`z` 算不出来")
    return (null_mu, sd_raw, float((near_mean - null_mu) / sd_raw), None)


def load_lr_pairs(cfg: dict):
    p = Path(__file__).resolve().parent.parent / "assets" / "ligand_receptor.yml"
    with open(p, encoding="utf-8") as fh:
        doc = yaml.safe_load(fh)
    return doc.get("pairs", []), p


def get_expression(adata):
    """
    取表达矩阵。

    **必须用全基因集，不是 HVG。** Part 2 实测：在 2000 个 HVG 上
    38 对 LR 只有 3 对可用；换到全基因集（13714 基因）变成 27 对。
    配体/受体基因大多是低表达基因，几乎不可能进 HVG。
    """
    if adata.raw is not None:
        X = adata.raw.X
        X = X.toarray() if sp.issparse(X) else np.asarray(X)
        return X.astype(np.float64), list(adata.raw.var_names), "adata.raw (全基因集)"
    raise RuntimeError(
        "adata.raw 为空 —— 配体-受体分析需要全基因集。\n"
        "  为什么必须停: 配体/受体大多是低表达基因，几乎不可能进 HVG。\n"
        "  Part 2 实测：2000 个 HVG 上 38 对 LR 只有 3 对可用，\n"
        "  全基因集上是 27 对 —— 用 HVG 会让结论建立在 3 对上。")


def run_07_spatial_communication(cfg: dict) -> dict:
    ensure_dirs(cfg)
    set_seed(cfg)
    data_dir = Path(cfg["output"]["data_dir"])
    res_dir = Path(cfg["output"]["results_dir"])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    c = cfg.get("communication") or {}
    if not c.get("enabled", True):
        status = {"dataset_id": cfg["dataset_id"], "status": "disabled",
                  "reason": "配置 communication.enabled=false"}
        write_json(res_dir / "communication_status.json", status)
        return status

    adata = sc.read_h5ad(data_dir / "domains.h5ad")
    X, genes, xsrc = get_expression(adata)
    log_info(f"表达矩阵: {X.shape[0]} spot x {X.shape[1]} 基因（来源 {xsrc}）")

    pairs, lr_path = load_lr_pairs(cfg)
    log_info(f"配体-受体对: {len(pairs)} 对（{lr_path.name}）")

    gi = {g: i for i, g in enumerate(genes)}

    # ---- 1. 可用性筛查 ------------------------------------------------------
    usable, skipped = [], []
    for pr in pairs:
        lig, rec = pr.get("ligand"), pr.get("receptor")
        if lig not in gi or rec not in gi:
            skipped.append({"pair": f"{lig}-{rec}",
                            "reason": "基因不在数据里",
                            "missing": [g for g in (lig, rec) if g not in gi]})
            continue
        li, ri_ = gi[lig], gi[rec]
        # 至少要在一部分 spot 里有表达，否则全零的列算出来全是 0
        frac_l = float((X[:, li] > 0).mean())
        frac_r = float((X[:, ri_] > 0).mean())
        if frac_l < 0.01 or frac_r < 0.01:
            skipped.append({"pair": f"{lig}-{rec}", "reason": "表达过低",
                            "frac_ligand": round(frac_l, 4),
                            "frac_receptor": round(frac_r, 4)})
            continue
        usable.append({**pr, "ligand_idx": li, "receptor_idx": ri_,
                       "frac_ligand": round(frac_l, 4),
                       "frac_receptor": round(frac_r, 4)})

    log_info(f"可用 LR 对: {len(usable)}/{len(pairs)}"
             f"（跳过 {len(skipped)} 对）")
    if len(usable) / max(len(pairs), 1) < 0.25:
        log_warn(f"**只有 {len(usable)}/{len(pairs)} 对 LR 可用** —— "
                 f"结论建立在这少数对上。常见原因：配体/受体基因"
                 f"不在 Visium 的基因覆盖里，或表达太低（零膨胀）")
    if not usable:
        status = {"dataset_id": cfg["dataset_id"], "status": "no_usable_pairs",
                  "reason": f"{len(pairs)} 对 LR 里没有一对的两个基因都可用",
                  "skipped": skipped[:50]}
        write_json(res_dir / "communication_status.json", status)
        return status

    # ---- 2. 空间距离矩阵（用 kNN，不是全对全）------------------------------
    from scipy.spatial import cKDTree
    xy = spatial_xy(adata)
    tree = cKDTree(xy)
    d, idx = tree.query(xy, k=min(31, len(xy)))
    med = float(np.median(d[:, 1]))
    # Visium spot 中心间距约 100 μm；把 μm 阈值换成像素
    um_per_px = 100.0 / med
    max_um = float(c.get("max_distance_um", 200))
    max_px = max_um / um_per_px
    log_info(f"空间阈值: {max_um} μm = {max_px:.1f} px"
             f"（spot 间距 {med:.1f} px ≈ 100 μm）")

    # 近邻掩码：距离 <= 阈值
    near = d[:, 1:] <= max_px
    n_near = int(near.sum())
    log_info(f"阈值内的 spot 对: {n_near}（平均每 spot {near.sum(1).mean():.1f} 个）")

    # ---- 3. 每对 LR 的空间富集 ---------------------------------------------
    # 统计量：近邻对上的 LR 共表达均值，减去随机邻居对的均值，再标准化。
    # **这个量是描述性的，不是假设检验。**
    #
    # **M6：两个缺陷都在这一段的零模型里。**
    #
    # ① `rng` 原先写在 `for pr in usable:` **循环体内** —— 每对 LR 都用
    #    同一条随机数流。50 个组合、每个组合的零分布**是同一批随机数**，
    #    于是不同 LR 对的 `null_mean`/`null_sd` 之间不是独立的，
    #    `z_score` 的排序里混进了"谁的 prod 分布恰好对上第一串随机数"。
    #    修法：`rng` 提到循环外，用 `rng.spawn()` 给每对一条独立的流
    #    （同时保持"同一个 seed 可复现"）。
    #
    # ② 零分布原先抽的是**单个 spot 的 `prod` 值**（`prod[rng.integers(...)]`），
    #    而 `near_mean` 是**邻居位置上的均值**。这两个量不是一回事：
    #    邻居位置上的 `prod` 之间存在空间自相关（相邻 spot 表达相似），
    #    所以"n_near 个独立抽样"的方差**小于**真实零分布 ——
    #    `null_sd` 偏小 → `z` 系统性偏大。修法：零分布改成从
    #    **全部邻居对**（`nb_vals`，即 kNN 的 30 个邻居，不受距离阈值限制）
    #    里随机抽同样多个，这样零分布保留 `prod` 的空间自相关结构，
    #    检验变成"近邻对比随机邻居"，这才是 `near_mean` 对应的零假设。
    rng = np.random.default_rng(cfg["analysis"]["seed"])
    n_perm = int(c.get("n_permutations", 50))
    lr_scores = []
    for pr in usable:
        L = X[:, pr["ligand_idx"]]
        R = X[:, pr["receptor_idx"]]
        # 归一化表达（每 spot 的总量已由 scanpy 归一化过，这里用 log1p 值）
        prod = L * R
        # 近邻对上的平均（用邻居索引取值）
        nb_vals = prod[idx[:, 1:]]
        near_mean = float(nb_vals[near].mean()) if n_near else np.nan
        # 零模型：**在全部邻居对里随机抽同样多个**（见上面 M6 ② 的理由）。
        # 抽的是邻居对而不是单 spot —— 零分布必须与被检验统计量同构。
        child = rng.spawn(1)[0]
        flat = nb_vals.ravel()
        if flat.size:
            draw = child.integers(0, flat.size, size=(n_perm, int(n_near)))
            perm_means = [float(flat[d].mean()) for d in draw]
        else:
            perm_means = [float("nan")] * n_perm
        null_mu, null_sd, z, z_reason = spatial_z_score(near_mean, perm_means)
        lr_scores.append({
            "ligand": pr["ligand"], "receptor": pr["receptor"],
            "pathway": pr.get("pathway", ""),
            "near_mean": _r5(near_mean),
            "null_mean": _r5(null_mu),
            "null_sd": _r5(null_sd),
            "z_score": None if z is None else round(float(z), 3),
            "z_defined": z is not None,
            "z_note": z_reason,
            "frac_ligand": pr["frac_ligand"], "frac_receptor": pr["frac_receptor"],
        })

    lr_df = pd.DataFrame(lr_scores).sort_values("z_score", ascending=False)
    lr_df.to_csv(res_dir / "communication_lr_scores.csv", index=False)
    # **`z_score` 可能是 `None`，格式化前必须判**（E-69 同族 Form A）：
    # `f"{None:.2f}"` 抛 `TypeError: unsupported format string passed to
    # NoneType.__format__`，会把整步带崩 —— 一个「算不出来」的诊断信息
    # 变成一次崩溃，是最不该发生的连锁。
    _top5_txt = ", ".join(
        f"{r.ligand}-{r.receptor}(z={r.z_score:.2f})" if r.z_score is not None
        else f"{r.ligand}-{r.receptor}(z=undefined)"
        for r in lr_df.head(5).itertuples())
    log_info("空间富集 top5: " + _top5_txt)
    if lr_df["z_score"].isna().any():
        n_undef = int(lr_df["z_score"].isna().sum())
        log_warn(f"{n_undef}/{len(lr_df)} 对 LR 的 z-score **算不出来**"
                 f"（零模型退化，见 communication_status.json 的 z_note）"
                 f"—— 这不是「z 很小」，排序里它们不该被读成「最不富集」")

    # ---- 4. 按域/类型的 LR 强度 ---------------------------------------------
    dom_blocks = []
    if "domain" in adata.obs.columns:
        doms = adata.obs["domain"].astype(str).values
        for pr in usable:
            prod = X[:, pr["ligand_idx"]] * X[:, pr["receptor_idx"]]
            for dm in sorted(set(doms)):
                m = doms == dm
                dom_blocks.append({"domain": dm, "ligand": pr["ligand"],
                                   "receptor": pr["receptor"],
                                   "pathway": pr.get("pathway", ""),
                                   "mean_product": round(float(prod[m].mean()), 5),
                                   "n_spots": int(m.sum())})
    if dom_blocks:
        dom_df = pd.DataFrame(dom_blocks)
        dom_df.to_csv(res_dir / "communication_by_domain.csv", index=False)
        log_info(f"按域的 LR 强度表: {len(dom_df)} 行")

    # ---- 5. 出图 ------------------------------------------------------------
    # **高度公式必须单位一致。** 原写法 max(mm(56), 0.32*len(top)+1.6) 把
    # 英寸当毫米用：mm(56)=2.2in，而 0.32*20+1.6=8.0"=203mm —— 评审量到
    # 203.2mm 高、宽高比 0.67 的根因。现在统一为英寸并夹在 [2.2, 5.4]"。
    top = lr_df.head(20).iloc[::-1]
    fig_h = min(max(2.2, 0.18 * len(top) + 1.2), 5.4)
    fig, ax = plt.subplots(figsize=(W_ONE_HALF, fig_h))
    # **条形颜色有含义就必须有图例。** 红灰蓝三色此前没有任何说明，
    # 蓝色 z=-2 虚线更是画在"通常没有数据"的左侧空白处，读者无从知道它是什么。
    #
    # **颜色一律走 PAL，不写裸字面量**（AGENTS 规则 13）。原写法就地写了
    # `#B2182B` / `#999999` / `#2166AC` —— 与 PAL 的语义色**不一致**
    # （PAL["highlight"] 是 #D55E00 橙，不是红），于是同一张图里"阈值线"
    # 和别处的"阈值线"不同色，读者要重新学一遍配色。用户反馈的
    # "图例布局有问题"正是这套自造配色 + 框内长图例叠加的结果。
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    # 用本仓库 PAL 的语义键（不是 geo 侧的 up/down —— 那套名字这里没有）。
    # 富集=highlight（橙）、耗竭=blue、中性=muted（灰）；三者色相与对比度
    # 已由 check_palette.mjs 校验。
    c_up, c_dn, c_ns = PAL["highlight"], PAL["blue"], PAL["muted"]
    handles = [
        Patch(facecolor=c_up, label="z > +2 (enriched near)"),
        Patch(facecolor=c_ns, label="-2 <= z <= +2"),
        Patch(facecolor=c_dn, label="z < -2 (depleted near)"),
        Line2D([0], [0], color=c_up, ls="--", lw=0.8, label="z = +2"),
        Line2D([0], [0], color=c_dn, ls="--", lw=0.8, label="z = -2"),
    ]
    # **图例放框外右侧、纵向单列**（约定 v2）。原 `frameon=True` 的框内图例
    # 压在条形上（20 条时图例正好盖住中段数据）。
    fig.legend(handles=handles, fontsize=5.5, ncol=1,
               loc="outside right center", frameon=False)
    # **z 可能是 `None`（E-69 同族 Form A）**：`null_sd` 算不出来时 `z_score`
    # 是 `None`，而 `None > 2` 会抛 `TypeError`。排序上 `None` 被 pandas 当作
    # NaN 排在最后 —— 这是想要的（没定义的 z 不该排进 top），但**画图必须
    # 显式跳过**，不能假设它一定是浮点。
    z_vals = top["z_score"]
    n_z_undef = int(z_vals.isna().sum())
    colors = [c_ns if (z is None or not np.isfinite(z))
              else (c_up if z > 2 else (c_dn if z < -2 else c_ns))
              for z in z_vals]
    ax.barh(range(len(top)), z_vals.fillna(0.0), color=colors)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([f"{r.ligand}–{r.receptor}" for r in top.itertuples()],
                       fontsize=7)
    ax.axvline(0, color=PAL["black"], lw=0.6)
    ax.axvline(2, color=c_up, ls="--", lw=0.8)
    ax.axvline(-2, color=c_dn, ls="--", lw=0.8)
    ax.set_xlabel("spatial enrichment z-score (near vs random)")
    # 有 z 没定义时把这件事写在标题里 —— 图上不能只有一根 0 长的灰条。
    _undef_line = (f"\n{n_z_undef} pair(s) z undefined (see z_note)"
                   if n_z_undef else "")
    ax.set_title(f"Ligand–receptor spatial enrichment\n"
                 f"{len(usable)}/{len(pairs)} pairs usable{_undef_line}")
    save_fig(cfg, "03-07-01-unit1-communication-lr-enrichment", fig)

    # 空间表达图：top 3 对
    # **单图原则拆分（D-006）**：top3 面板 -> 3 张独立单图（P7 LR 空间图谱
    # 分解链）。共享量程 0->p99 保留，写进各图 title。
    #
    # **top3 只从 z 有定义的对里选**（E-69 同族 Form A）：`sort_values` 把
    # `z_score=None` 排在最后，正常情况下选不到；但如果**有定义的对不足 3 个**，
    # 未定义的对就会被选进来画一张「z=undefined」的图 —— 那等于把「算不出来」
    # 当成一个发现展示。
    sf = float(adata.uns["spatial"][list(adata.uns["spatial"])[0]]
               ["scalefactors"]["tissue_hires_scalef"])
    xyp = xy * sf
    top3 = top_lr_for_plots(lr_df, 3)
    # **`top3` 可能是空的**（E-69 同族 Form A）：如果**所有** LR 对的 z 都算不出来
    # （零模型退化），过滤后是 0 行，于是 `np.concatenate([])` 抛
    # `ValueError: need at least one array to concatenate` —— 「一个量算不出来」
    # 的处境变成整步崩溃。这不是理论情况：只要 `max_distance_um` 设得比最近邻
    # 距离还小，`n_near=0` 就会让全部 z 无定义。
    if top3.empty:
        log_warn("所有 LR 对的 z-score 都算不出来 —— 跳过 top3 空间图谱"
                 "（没有可展示的对；见 communication_status.json 的 z_note）")
    else:
        _prods = []
        for r in top3.itertuples():
            _prods.append(X[:, gi[r.ligand]] * X[:, gi[r.receptor]])
        vmax_lr = float(np.quantile(np.concatenate(_prods), 0.99))
        DYNAMIC_FIG_BASES = {"02": 3}
        for ui, r in enumerate(top3.itertuples(), start=1):
            li, ri_ = gi[r.ligand], gi[r.receptor]
            prod = X[:, li] * X[:, ri_]
            fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(62)))
            s = ax.scatter(xyp[:, 0], xyp[:, 1], c=prod, s=4, cmap="viridis",
                           vmin=0, vmax=vmax_lr)
            ax.set_title(f"{r.ligand} × {r.receptor}  z={r.z_score:.2f}\n"
                         f"shared scale 0 - {vmax_lr:.2f} (ligand × receptor, log1p)")
        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(s, ax=ax, shrink=0.8, pad=0.02, fraction=0.046,
                     label="ligand × receptor (log1p expression)")
        pair_slug = f"{r.ligand}-{r.receptor}".lower()
        save_fig(cfg, f"03-07-02-unit{ui}-{pair_slug}", fig)

    # ---- 6. 落盘 ------------------------------------------------------------
    status = {
        "dataset_id": cfg["dataset_id"],
        "status": "ok",
        "expression_source": xsrc,
        "n_pairs_total": len(pairs),
        "n_pairs_usable": len(usable),
        "n_pairs_skipped": len(skipped),
        "usable_fraction": round(len(usable) / max(len(pairs), 1), 4),
        "max_distance_um": max_um,
        "max_distance_px": round(max_px, 2),
        "median_spot_spacing_px": round(med, 2),
        "distance_scale": {
            "um_per_px_assumed": round(float(um_per_px), 5),
            "assumption": ("Visium spot 中心间距名义 100 μm（芯片规格）—— "
                           "**这是假定值，不是从本切片估计出来的**；"
                           "`median_spot_spacing_px` 才是观测值"),
            "implied_um_per_px_from_observation": (
                round(float(100.0 / med), 5) if med > 0 else None),
        },
        "null_model": {
            "description": ("零分布 = 从**全部 kNN 邻居对**里随机抽同样多个"
                            "（不是从单个 spot 抽）—— 保留 prod 的空间"
                            "自相关结构"),
            "rng": "每对 LR 独立 spawn（原先每对重用同一条随机流）",
            "n_permutations": n_perm,
        },
        "n_spot_pairs_within_threshold": n_near,
        # **「算不出来」必须与「算出来很小」分开报**（E-69 同族 Form A）。
        # 原来 `z_score` 一律是浮点，读者无法区分 `z=-0.3`（真测出来）与
        # `z=nan`（零模型退化）—— 后者在排序里落到最末，被读成「最不富集」。
        "n_z_defined": int(lr_df["z_score"].notna().sum()),
        "n_z_undefined": int(lr_df["z_score"].isna().sum()),
        "z_undefined_note": (
            "z 未定义的对：零模型 50 次置换全部为 nan（`n_near=0` 或 `flat` 为空）"
            "或置换标准差恰为 0。逐对的取值与原因见 `top_enriched[].z_note`；"
            "**这不是「z 很小」，这些对不能参与富集排序**"
            if int(lr_df["z_score"].isna().sum()) else None),
        "top_enriched": df_to_records(lr_df.head(15)),
        "skipped_examples": skipped[:20],
        "method": ("空间约束的配体-受体共表达强度："
                   "统计阈值内 spot 对上的配体×受体表达均值，"
                   "与随机 spot 对比较得 z-score"),
        "not_a_call": ("**这不是『通讯与否』的判定。** 真正的细胞通讯推断需要 "
                       "CellPhoneDB / CellChat 那样的统计框架"
                       "（含置换检验与受体复合物建模）。本工具只给"
                       "空间约束下的描述性强度"),
        # §3.5 点名的 CellChat 没跑 —— 和 not_a_call 并列，别让读者以为
        # "既然有 not_a_call 说明替代方案已经上了"
        "named_tools": probe_named_tools(log=log_warn, only=("CellChat",)),
        "limitations": [
            f"**只有 {len(usable)}/{len(pairs)} 对 LR 可用** —— "
            "结论建立在这些对上；不可用的原因是基因不在 Visium 覆盖里"
            "或表达太低",
            "配体×受体表达量是**代理量**，不代表蛋白水平的信号传递",
            "**共表达 ≠ 通讯。** 两个基因在同一个 spot 里高，可能是"
            "同一个细胞表达了两者（自分泌），也可能是两个细胞紧邻",
            "Visium 的 spot 含 1-10 个细胞，所以『阈值内』是 spot 层面"
            "的接近，不是细胞接触",
            "z-score 用 50 次随机采样估计零分布，精度有限；"
            "且没有做多重检验校正（组合数是 O(配体×受体)，z 的排序不能"
            "直接读成显著性）",
            "**零模型与统计量同构但仍是近似**：零分布从全部 kNN 邻居对里抽，"
            "保住了 prod 的空间自相关，但『邻居对』本身的图结构"
            "（六边形、边界 spot 度数少）没有完全复制",
            "受体复合物（如 IL2 受体的 α/β/γ 三聚体）被简化成单个受体基因",
        ],
    }
    write_json(res_dir / "communication_status.json", status)
    return status


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        run_07_spatial_communication(cfg)
        record_step(cfg, "spatial_communication", "ok", time.time() - t0)
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "spatial_communication", "failed", time.time() - t0,
                    message=str(e))
        raise
