#!/usr/bin/env python3
"""
03_spatial_domains.py — 空间域识别

**空间域 ≠ 表达簇。** 这是本步骤存在的全部理由。

普通的 Leiden 聚类只看表达相似度，所以同一个空间连续区域可能被切成
几块、或者两个不相邻的区域被合成一个簇 —— 因为聚类不知道 spot 在哪。
空间域要求"域内空间连续"，做法是把每个 spot 的表达与其**空间邻居**
平均后再聚类。

这里**两种都跑，并把差别量化**：
  - `expression_only`：不平滑，等价于普通表达聚类
  - `spatial`：平滑后聚类

量化指标：
  - **空间连贯性**：每个域在邻居图上的连通分量数（越少越连贯）
  - **邻居同域率**：相邻 spot 属于同一域的比例（越高越平滑）
  - **域数**：平滑通常会减少域数（把碎块合并）

**平滑不是无代价的。** 平滑强度过高会把真实的微小结构（如一个小的
生发中心）抹平。所以扫描多个平滑强度，让选择有依据。

出图用 squidpy 把域叠到 H&E 上 —— 这是唯一能判断"域划分是否对应
真实组织学结构"的方法。
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
from scipy.sparse.csgraph import connected_components  # noqa: E402

from common import (df_to_records, ensure_dirs, load_config, log_info,  # noqa: E402
                    log_warn, parse_args, record_step, save_fig, set_seed,
                    spot_radius_plot_units, write_json, spatial_xy, W_DOUBLE, mm, PAL,)


def spatial_neighbor_graph(adata, n_neighbors: int = 6):
    """
    建空间邻居图。

    用 **kNN + 距离阈值**，不是纯 kNN。理由：Visium 网格在组织边缘会有
    "半个六边形"的情况，纯 kNN 会把边缘 spot 和距离很远的 spot 连起来，
    而距离阈值能挡掉。

    阈值取"最近邻距离中位数的 1.6 倍" —— Visium 六边形网格里，
    6 个直接邻居的距离是 d，次近邻是 √3·d ≈ 1.73d，所以 1.6 倍
    正好把直接邻居包进来、把次近邻挡在外面。
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
    dists = d[:, 1:].ravel()
    mask = dists <= thresh
    adj = sp.coo_matrix((np.ones(int(mask.sum())), (rows[mask], cols[mask])),
                        shape=(len(xy), len(xy))).tocsr()
    # 对称化
    adj = adj.maximum(adj.T)
    return adj, {"n_neighbors_k": int(k - 1), "median_nn_dist": round(med, 2),
                 "distance_threshold": round(thresh, 2),
                 "mean_degree": round(float(adj.getnnz(axis=1).mean()), 2)}


def smooth_embeddings(adata, adj, alpha: float, key: str = "X_pca",
                      out_key: str = "X_pca_smooth") -> np.ndarray:
    """
    把每个 spot 的嵌入与其空间邻居平均。

    `alpha` 是平滑强度：0 = 完全用自身，1 = 完全用邻居均值。
    用行归一化的邻接矩阵做加权平均。

    **注意邻接矩阵要行归一化**，否则度数高的 spot 会被放大。
    """
    X = np.asarray(adata.obsm[key], dtype=float)
    if alpha <= 0:
        return X.copy()
    deg = np.asarray(adj.sum(axis=1)).ravel()
    deg[deg == 0] = 1.0
    norm_adj = sp.diags(1.0 / deg) @ adj
    nb = norm_adj @ X
    return (1.0 - alpha) * X + alpha * nb


def domain_metrics(adata, labels: np.ndarray, adj) -> dict:
    """
    量化一个域划分的"空间性"。

    - `neighbor_same_frac`：相邻 spot 属于同一域的比例。
      **这是核心指标** —— 必须与 `random_baseline_same_frac`
      （Σ 各域占比²）比较才有意义。绝对值的解读依赖于域数和域大小分布。
    - `n_components_total`：所有域在邻居图上的连通分量总数。
    - `fragmented_domains`：被切成多块的域的数量。

    **`fragmented_domains` 不是缺陷指标，解读要小心。** 它假设
    "一个域 = 一块连续区域"，而这个假设对**有重复结构的组织不成立**：
      - 淋巴结的滤泡、生发中心本身就是散布的多个斑块
      - 肿瘤的癌巢、三级淋巴结构同理
    实测淋巴结数据：不平滑时 9/9 个域都是多块，平滑后 11/13 ——
    **那不是"聚类没做好"，那是滤泡本来就有很多个。**
    只有对"应当连续"的组织（如脑的层状结构、上皮分层）才能把它当缺陷看。
    """
    labels = np.asarray(labels).astype(str)
    uniq = sorted(set(labels))
    coo = adj.tocoo()
    same = float((labels[coo.row] == labels[coo.col]).mean())

    # 随机基线：按域占比算期望同域率
    p = np.array([(labels == u).mean() for u in uniq])
    baseline = float((p ** 2).sum())

    n_comp_total, n_frag = 0, 0
    for u in uniq:
        m = labels == u
        sub = adj[m][:, m]
        nc, _ = connected_components(sub, directed=False)
        n_comp_total += nc
        if nc > 1:
            n_frag += 1
    return {
        "n_domains": len(uniq),
        "neighbor_same_frac": round(same, 4),
        "random_baseline_same_frac": round(baseline, 4),
        "spatial_enrichment": round(same - baseline, 4),
        "n_components_total": int(n_comp_total),
        "fragmented_domains": int(n_frag),
        "domain_sizes": {u: int((labels == u).sum()) for u in uniq},
    }


def run_03_spatial_domains(cfg: dict) -> dict:
    ensure_dirs(cfg)
    set_seed(cfg)
    data_dir = Path(cfg["output"]["data_dir"])
    res_dir = Path(cfg["output"]["results_dir"])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    adata = sc.read_h5ad(data_dir / "normalized.h5ad")
    log_info(f"读入 {adata.n_obs} spot x {adata.n_vars} HVG")

    dom = cfg["domains"]
    n_spat = int(dom.get("n_spatial_neighbors", 6))
    adj, graph_info = spatial_neighbor_graph(adata, n_spat)
    log_info(f"空间邻居图: 平均度数 {graph_info['mean_degree']}，"
             f"最近邻距离中位 {graph_info['median_nn_dist']}，"
             f"阈值 {graph_info['distance_threshold']}")

    sc.pp.neighbors(adata, n_neighbors=int(cfg["reduce"]["n_neighbors"]),
                    n_pcs=int(cfg["reduce"]["n_pcs"]),
                    random_state=cfg["analysis"]["seed"])
    sc.tl.umap(adata, random_state=cfg["analysis"]["seed"])

    res_used = float(dom.get("resolution", 0.8))
    alpha_used = float(dom.get("smoothing", 0.5))

    # ---- 1. 对照：不平滑 vs 平滑 --------------------------------------------
    sc.tl.leiden(adata, resolution=res_used, key_added="domain_expr_only",
                 flavor="igraph", n_iterations=2, directed=False,
                 random_state=cfg["analysis"]["seed"])
    m_expr = domain_metrics(adata, adata.obs["domain_expr_only"].values, adj)
    log_info(f"不平滑: {m_expr['n_domains']} 域，邻居同域率 "
             f"{m_expr['neighbor_same_frac']:.3f}（随机基线 "
             f"{m_expr['random_baseline_same_frac']:.3f}），"
             f"{m_expr['fragmented_domains']} 个域被切碎")

    adata.obsm["X_pca_smooth"] = smooth_embeddings(adata, adj, alpha_used)
    sc.pp.neighbors(adata, n_neighbors=int(cfg["reduce"]["n_neighbors"]),
                    use_rep="X_pca_smooth", key_added="spatial_nn",
                    random_state=cfg["analysis"]["seed"])
    sc.tl.leiden(adata, resolution=res_used, key_added="domain",
                 neighbors_key="spatial_nn",
                 flavor="igraph", n_iterations=2, directed=False,
                 random_state=cfg["analysis"]["seed"])
    m_spat = domain_metrics(adata, adata.obs["domain"].values, adj)
    log_info(f"平滑({alpha_used}): {m_spat['n_domains']} 域，邻居同域率 "
             f"{m_spat['neighbor_same_frac']:.3f}，"
             f"{m_spat['fragmented_domains']} 个域被切碎")

    # ---- 2. 平滑强度扫描 ----------------------------------------------------
    scan = []
    for a in [0.0, 0.25, 0.5, 0.75, 0.9]:
        emb = smooth_embeddings(adata, adj, a)
        adata.obsm["_scan"] = emb
        sc.pp.neighbors(adata, n_neighbors=int(cfg["reduce"]["n_neighbors"]),
                        use_rep="_scan", key_added="_scan_nn",
                        random_state=cfg["analysis"]["seed"])
        sc.tl.leiden(adata, resolution=res_used, key_added="_scan_lab",
                     neighbors_key="_scan_nn",
                     flavor="igraph", n_iterations=2, directed=False,
                     random_state=cfg["analysis"]["seed"])
        m = domain_metrics(adata, adata.obs["_scan_lab"].values, adj)
        scan.append({"smoothing": a, **{k: v for k, v in m.items()
                                        if k != "domain_sizes"}})
    scan_df = pd.DataFrame(scan)
    scan_df.to_csv(res_dir / "smoothing_scan.csv", index=False)
    log_info("平滑扫描: " + ", ".join(
        f"{r.smoothing}->{r.n_domains}域/同域率{r.neighbor_same_frac:.3f}"
        for r in scan_df.itertuples()))

    fig, axes = plt.subplots(1, 2, figsize=(W_DOUBLE, mm(64)))
    axes[0].plot(scan_df["smoothing"], scan_df["neighbor_same_frac"], "o-",
                 color=PAL["primary"], label="observed")
    axes[0].plot(scan_df["smoothing"], scan_df["random_baseline_same_frac"], "s--",
                 color=PAL["highlight"], label="random baseline")
    axes[0].set_xlabel("smoothing strength α")
    axes[0].set_ylabel("neighbor same-domain fraction")
    axes[0].set_title("Spatial coherence vs smoothing")
    axes[0].legend(fontsize=7)
    axes[1].plot(scan_df["smoothing"], scan_df["n_domains"], "o-", color=PAL["primary"])
    axes[1].set_xlabel("smoothing strength α")
    axes[1].set_ylabel("number of domains")
    axes[1].set_title("Domain count vs smoothing")
    save_fig(cfg, "smoothing_scan", fig)

    # ---- 3. 分辨率扫描 ------------------------------------------------------
    rscan = []
    for r in dom.get("resolution_scan", [0.2, 0.4, 0.6, 0.8, 1.0, 1.5]):
        sc.tl.leiden(adata, resolution=float(r), key_added=f"_r{r}",
                     neighbors_key="spatial_nn",
                     flavor="igraph", n_iterations=2, directed=False,
                     random_state=cfg["analysis"]["seed"])
        m = domain_metrics(adata, adata.obs[f"_r{r}"].values, adj)
        rscan.append({"resolution": float(r), "n_domains": m["n_domains"],
                      "neighbor_same_frac": m["neighbor_same_frac"],
                      "fragmented_domains": m["fragmented_domains"]})
    pd.DataFrame(rscan).to_csv(res_dir / "domain_resolution_scan.csv", index=False)

    # ---- 4. Marker 与注释 ---------------------------------------------------
    sc.tl.rank_genes_groups(adata, "domain", method="wilcoxon", use_raw=True,
                            random_state=cfg["analysis"]["seed"])
    top_n = int((cfg.get("analysis") or {}).get("top_markers", 25))
    frames = []
    for g in adata.obs["domain"].cat.categories:
        df = sc.get.rank_genes_groups_df(adata, group=g).head(top_n)
        df.insert(0, "domain", str(g))
        frames.append(df)
    markers = pd.concat(frames, ignore_index=True)
    markers.to_csv(res_dir / "domain_markers.csv", index=False)
    log_info(f"域 marker 表: {len(markers)} 行")

    top3 = (markers.sort_values(["domain", "scores"], ascending=[True, False])
            .groupby("domain", observed=True).head(3)["names"].unique().tolist())
    # **上限按图宽算，不是随手取 40。** 双栏 183 mm 下每个基因约 7 mm，
    # 再多标签就挤成一片。原来取 40 会画出 370 mm 宽的图 —— 装不进任何
    # 期刊的一页。
    top3 = [g for g in top3 if g in adata.raw.var_names][:24]
    if top3:
        sc.pl.dotplot(adata, top3, groupby="domain", use_raw=True, show=False,
                      standard_scale="var")
        fig = plt.gcf()
        # scanpy 自己按基因数定尺寸，这里拉回标准双栏宽。
        # 高度 96 mm 是实测值：80 mm 时域标签顶出画布 +4.2%，88 mm 时 +2.6%
        fig.set_size_inches(W_DOUBLE, mm(96))
        fig.suptitle("Top markers per spatial domain")
        save_fig(cfg, "domain_markers_dotplot", fig)

    # ---- 4b. 域的组织学标签（用 marker 签名打分）----------------------------
    # **这不是"域 = 某一种细胞"。** Visium 的 spot 含 1-10 个细胞，
    # 一个域是"细胞组成相近的一片区域"。这里做的是给域一个**可读的标签**。
    #
    # ---
    # ## 为什么要按类型做 z-score（第一版错在这里）
    #
    # 第一版直接比"该域里各类型 marker 的平均表达"，结果 13 个域里
    # **12 个都被标成 Plasma_cell** —— 因为浆细胞的 marker
    # （MZB1 / JCHAIN / IGHG1 / IGKC …）在淋巴结里表达量本来就极高，
    # 绝对表达一比就压倒所有其他类型。
    #
    # **这和之前 NNLS 那次是同一类错误：一个看起来像答案、但不含信息的标签。**
    #
    # 正确的问法是"**哪个域对某个类型相对最富集**"，而不是"哪个类型
    # 表达最高"。所以先对每个类型跨全部域做 z-score —— 这样消掉了
    # "某些类型的 marker 天生高表达"这个系统性偏差，
    # 每个类型都站在自己的尺度上比。
    annot = {}
    try:
        import yaml
        sig_path = (Path(__file__).resolve().parent.parent / "assets"
                    / "reference_signatures.yml")
        sig_key = ((cfg.get("deconvolution") or {}).get("signature_set")
                   or (cfg.get("analysis") or {}).get("celltype_markers"))
        with open(sig_path, encoding="utf-8") as fh:
            sigs = yaml.safe_load(fh).get("signatures", {})
        sig = sigs.get(sig_key)
        if sig:
            rawX = adata.raw.X
            rawX = rawX.toarray() if sp.issparse(rawX) else np.asarray(rawX)
            rgenes = list(adata.raw.var_names)
            rgi = {g: i for i, g in enumerate(rgenes)}
            doms = adata.obs["domain"].astype(str).values
            dom_list = sorted(set(doms), key=lambda x: int(x) if x.isdigit() else x)

            # 矩阵 (域 × 类型) 的原始平均 marker 表达
            used_ct, M = [], []
            for ct, d in sig["celltypes"].items():
                present = [g for g in d.get("markers", []) if g in rgi]
                if len(present) < 3:
                    continue
                col = [float(rawX[doms == dm][:, [rgi[g] for g in present]].mean())
                       for dm in dom_list]
                used_ct.append(ct)
                M.append(col)
            M = np.asarray(M, dtype=float)          # (n_ct, n_dom)

            # **按类型（行）做 z-score** —— 消掉"某些类型 marker 天生高表达"
            mu = M.mean(axis=1, keepdims=True)
            sd = M.std(axis=1, keepdims=True)
            sd[sd < 1e-12] = 1.0
            Z = (M - mu) / sd                       # (n_ct, n_dom)

            for j, dm in enumerate(dom_list):
                col = Z[:, j]
                order = np.argsort(-col)
                top_i, second_i = int(order[0]), int(order[1]) if len(order) > 1 else int(order[0])
                z_top, z_second = float(col[top_i]), float(col[second_i])
                annot[dm] = {
                    "label": used_ct[top_i],
                    "z_score": round(z_top, 4),
                    "runner_up": used_ct[second_i],
                    "runner_up_z": round(z_second, 4),
                    # margin 小时标签不该被当结论（沿用 Part 2 的做法）
                    "z_margin": round(z_top - z_second, 4),
                    "assignment_confident": bool(z_top - z_second > 0.5),
                    "raw_mean_expression": round(float(M[top_i, j]), 4),
                    "n_spots": int((doms == dm).sum()),
                    "top5": [{"celltype": used_ct[int(i)],
                              "z_score": round(float(col[int(i)]), 4)}
                             for i in order[:5]],
                }
            log_info("域标签: " + ", ".join(
                f"{k}->{v['label']}(z {v['z_score']:.2f}, Δ{v['z_margin']:.2f})"
                for k, v in annot.items()))
            n_unc = sum(1 for v in annot.values() if not v["assignment_confident"])
            if n_unc:
                log_warn(f"{n_unc}/{len(annot)} 个域的标签 z_margin <= 0.5 —— "
                         f"这些标签不该被当结论")
            write_json(res_dir / "domain_annotation.json",
                       {"signature_set": sig_key, "domains": annot,
                        "scoring": ("按类型跨域做 z-score 后取最高 —— "
                                    "问的是『哪个域对这个类型相对最富集』，"
                                    "不是『哪个类型表达最高』"),
                        "why_zscore": ("直接比绝对平均表达时，13 个域里 12 个都被"
                                       "标成 Plasma_cell —— 浆细胞的 marker "
                                       "（MZB1/JCHAIN/IGHG1/IGKC）在淋巴结里"
                                       "表达量本来就极高，绝对表达一比就压倒"
                                       "所有其他类型。z-score 消掉了这个偏差"),
                        "caveat": ("**标签是提示，不是结论。** Visium 的 spot 含 "
                                   "1-10 个细胞，域是『组成相近的一片区域』，"
                                   "不是某一种细胞。z_margin 小时标签不可信")})
    except Exception as e:  # noqa: BLE001
        log_warn(f"域标签打分跳过: {type(e).__name__}: {e}")

    # ---- 5. 叠到 H&E 上（唯一能判断域是否对应组织学的方法）-------------------
    lib = list(adata.uns["spatial"].keys())[0]
    entry = adata.uns["spatial"][lib]
    img = entry["images"]["hires"]
    sf = float(entry["scalefactors"]["tissue_hires_scalef"])
    r_plot = spot_radius_plot_units(entry["scalefactors"], "hires")
    xy = spatial_xy(adata, sf)

    fig, axes = plt.subplots(1, 3, figsize=(W_DOUBLE, mm(56)))
    for ax, key, title in zip(
            axes,
            ["domain_expr_only", "domain", None],
            [f"Expression-only clustering ({m_expr['n_domains']} clusters)",
             f"Spatial domains, α={alpha_used} ({m_spat['n_domains']} domains)",
             "H&E reference"]):
        ax.imshow(img, alpha=1.0 if key is None else 0.55)
        if key is not None:
            cats = adata.obs[key].astype(str).values
            uniq = sorted(set(cats), key=lambda x: int(x) if x.isdigit() else x)
            cmap = plt.get_cmap("tab20")
            for i, u in enumerate(uniq):
                m = cats == u
                ax.scatter(xy[m, 0], xy[m, 1], s=8, color=cmap(i % 20),
                           label=u, linewidths=0)
            ax.legend(fontsize=6, markerscale=2.2, loc="upper right",
                      ncol=2, framealpha=0.8)
        ax.set_title(title)
        ax.set_xticks([]); ax.set_yticks([])
    # **suptitle 要折行。** constrained layout 不会给文字换行 —— 一行长标题
    # 会把整张图撑得比 183 mm 宽，而 savefig.bbox="standard" 下多出来的
    # 部分直接裁掉。_content_overflow() 就是靠这条抓到的。
    fig.suptitle("Spatial domains overlaid on H&E\n"
                 "(the only way to judge whether domains match real histology)")
    save_fig(cfg, "domains_on_he", fig)

    # ---- 6. 落盘 ------------------------------------------------------------
    adata.obs[["domain", "domain_expr_only"]].to_csv(res_dir / "spatial_domains.csv")
    out = data_dir / "domains.h5ad"
    adata.write_h5ad(out)
    log_info(f"已写出 {out}")

    status = {
        "dataset_id": cfg["dataset_id"],
        "n_spots": int(adata.n_obs),
        "resolution": res_used,
        "smoothing_alpha": alpha_used,
        "spatial_graph": graph_info,
        "expression_only": m_expr,
        "spatial": m_spat,
        "smoothing_scan": df_to_records(scan_df),
        "resolution_scan": rscan,
        "improvement": {
            "neighbor_same_frac_delta": round(
                m_spat["neighbor_same_frac"] - m_expr["neighbor_same_frac"], 4),
            "fragmented_domains_delta": int(m_spat["fragmented_domains"] -
                                            m_expr["fragmented_domains"]),
            "n_domains_delta": int(m_spat["n_domains"] - m_expr["n_domains"]),
        },
        "note": ("**空间域 ≠ 表达簇。** 不平滑的聚类不知道 spot 在哪，"
                 "会把连续区域切碎或把不相邻区域合并。"
                 "`neighbor_same_frac` 是判断空间性的核心指标，"
                 "它要与 `random_baseline_same_frac` 比较才有意义"),
        "limitations": [
            "平滑会抹平真实的微小结构（如小的生发中心）；α 过高时域数减少但可能丢细节",
            "Visium 的 spot 直径 55 μm，含 1-10 个细胞 —— 域的空间分辨率受此限制，"
            "**不能说『某个域是某一种细胞』**，只能说这个区域的细胞组成不同",
            "域边界的位置有 ±1 个 spot 的不确定性",
            "**`fragmented_domains` 高不等于聚类失败。** 淋巴结的滤泡、"
            "肿瘤的癌巢本身就是散布的多个斑块 —— 同一个域出现在多个不相邻"
            "位置是生物学事实。只有对『应当连续』的组织（脑的层状结构、"
            "上皮分层）才能把这个指标当缺陷看",
        ],
        "status": "ok",
    }
    write_json(res_dir / "domain_status.json", status)
    return status


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        run_03_spatial_domains(cfg)
        record_step(cfg, "spatial_domains", "ok", time.time() - t0)
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "spatial_domains", "failed", time.time() - t0, message=str(e))
        raise
