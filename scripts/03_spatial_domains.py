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

import random
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
                    log_warn, named_tools_note, parse_args, pkg_version,
                    probe_named_tools, record_step, save_fig, set_seed,
                    write_json, spatial_xy, W_DOUBLE, W_ONE_HALF, mm, build_marker_dotplot_figure, PAL,)


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


# ---------------------------------------------------------------------------
# §3.2 点名的空间域方法：SpaGCN / STAGATE
# ---------------------------------------------------------------------------
# **这两个是"交叉验证"，不是替换主方法。** 主方法仍是上面的
# 平滑 + Leiden（内置实现）—— 它可复现、无外部依赖、并且已经验证过。
# 点名工具跑起来之后，要回答的是"两套划分是不是在说同一件事"，
# 所以必须报 ARI / NMI / 邻居同域率，**不是只报"跑通了"**。
#
# ## SpaGCN 上一轮被登记成「装不上」，那个结论是错的
#
# 上一轮的登记写的是：`SpaGCN` 依赖 `louvain`，而 `louvain` 0.8.2
# 没有 py3.12 wheel、只有 2019 年的 sdist，要编译 C++/Cython。
# **前半句是实测的，后半句的推论是错的** —— 读 SpaGCN 1.2.7 的源码：
#
#   - `SpaGCN/SpaGCN.py`、`models.py`、`util.py` 里**没有一处**
#     `import louvain`；
#   - 它走的是 `scanpy.tl.louvain`（`models.py:69`、`util.py:272`）；
#   - 而 `simple_GC_DEC.fit` 的 `init` 参数有 **`"kmeans"` 分支**
#     （`models.py:52-61`），完全不碰 louvain。
#
# 所以 `pip install --no-deps SpaGCN` + `init="kmeans"` 可以绕开整条编译链。
# **`install_requires` 里有某个包，不等于运行时会 import 它** ——
# 判据要读源码，不能只读元数据。
#
# 代价（必须说清楚）：`init="louvain"` 让 SpaGCN 用表达+空间初始化簇心，
# 换成 `kmeans` 后初始化只用空间平滑后的 GCN 特征。`n_clusters` 因此
# **必须外部给定**，我们传内置方法得到的域数，这样两边域数相同、
# ARI 才可比。

SPAGCN_DEFAULTS = {"num_pcs": 30, "lr": 0.005, "max_epochs": 200, "tol": 1e-3}


def _agreement(labels_a, labels_b, adj, adata, log=log_info) -> dict:
    """两套域划分的一致程度。

    **ARI 是主判据**：它对"域编号不同但划分相同"免疫（域 3 和域 7 换名
    不影响），而这正是这里要问的问题。`same_label_frac` 不是 ——
    它会被编号顺序完全支配，所以只作参考一起报。
    """
    from sklearn.metrics import (adjusted_rand_score,
                                 normalized_mutual_info_score)
    a = np.asarray(labels_a).astype(str)
    b = np.asarray(labels_b).astype(str)
    out = {
        "n_labels_a": int(len(set(a.tolist()))),
        "n_labels_b": int(len(set(b.tolist()))),
        "adjusted_rand_index": round(float(adjusted_rand_score(a, b)), 4),
        "normalized_mutual_info": round(
            float(normalized_mutual_info_score(a, b)), 4),
        "same_label_frac": round(float((a == b).mean()), 4),
        "note": ("ARI 对域编号置换免疫，是主判据；same_label_frac 会被编号"
                 "顺序支配，只作参考"),
    }
    # **域×域对应矩阵**（差距清单 #25，SRC-17 f9E-G / SRC-16 惯例）：
    # ARI 是单个数，看不出"哪几个域对上了、哪些被分裂/合并"。
    # contingency 矩阵（行=A 的域、列=B 的域、值=Jaccard）把 ARI 拆开。
    try:
        ua = sorted(set(a.tolist()), key=lambda x: int(x) if x.isdigit() else x)
        ub = sorted(set(b.tolist()), key=lambda x: int(x) if x.isdigit() else x)
        jac_m = np.zeros((len(ua), len(ub)), dtype=float)
        for ia, va in enumerate(ua):
            sa = (a == va)
            for ib, vb in enumerate(ub):
                sb = (b == vb)
                inter = float(np.sum(sa & sb))
                union = float(np.sum(sa | sb))
                jac_m[ia, ib] = inter / union if union > 0 else 0.0
        out["jaccard_matrix"] = {"rows": ua, "cols": ub,
                                 "values": [[round(float(v), 4) for v in row]
                                            for row in jac_m]}
        out["jaccard_matrix_note"] = (
            "行 = A 的域、列 = B 的域、值 = Jaccard 重叠度；"
            "对角线以外的亮块 = 被分裂或合并的域")
    except Exception as e:  # noqa: BLE001
        log(f"  域×域矩阵失败: {type(e).__name__}: {e}")
    ma = domain_metrics(adata, a, adj)
    mb = domain_metrics(adata, b, adj)
    out["neighbor_same_frac_a"] = ma["neighbor_same_frac"]
    out["neighbor_same_frac_b"] = mb["neighbor_same_frac"]
    out["spatial_enrichment_a"] = ma["spatial_enrichment"]
    out["spatial_enrichment_b"] = mb["spatial_enrichment"]
    log(f"  ARI={out['adjusted_rand_index']}  NMI={out['normalized_mutual_info']}"
        f"  邻居同域率 {ma['neighbor_same_frac']} vs {mb['neighbor_same_frac']}")
    return out


def try_spagcn(adata, cfg: dict, n_clusters: int, log=log_info) -> tuple:
    """跑 SpaGCN（文档 §3.2），返回 `(labels | None, prob | None, info)`。

    **任何异常都吞掉并写进 info** —— 点名工具跑不起来不该让整步失败，
    但必须让人看见（AGENTS 规则 4）。
    """
    info = {"attempted": True, "status": None, "reason": "", "section": "§3.2",
            "init": "kmeans",
            "why_kmeans": ("SpaGCN 的 install_requires 里有 `louvain`，但源码"
                           "从不 import 它（走 scanpy.tl.louvain），且 init 有 "
                           "kmeans 分支 —— 用 kmeans 绕开 louvain 的 py3.12 "
                           "编译链")}
    try:
        import SpaGCN as spg
    except Exception as exc:  # noqa: BLE001
        info["status"] = "package_missing"
        info["reason"] = f"SpaGCN 未安装（{type(exc).__name__}: {exc}）"
        log(f"  §3.2 SpaGCN 未跑：{info['reason']}")
        return None, None, info

    info["version"] = pkg_version("SpaGCN")

    # ---- 坐标：用 Visium 的 array 索引（六边形网格的整数坐标）---------------
    # 不用像素坐标：像素下 l 是几百，含义要换算才知道；array 坐标下
    # 最近邻距离就是 1，l 的值一眼能对上网格尺度。
    for k in ("array_row", "array_col"):
        if k not in adata.obs.columns:
            info["status"] = "no_array_coords"
            info["reason"] = (f"adata.obs 里没有 '{k}' —— SpaGCN 需要 Visium 的 "
                              f"array 索引；不退回像素坐标（尺度含义会变，"
                              f"l 就没法解释）")
            log(f"  §3.2 SpaGCN 未跑：{info['reason']}")
            return None, None, info
    x = adata.obs["array_row"].astype(float).to_numpy()
    y = adata.obs["array_col"].astype(float).to_numpy()

    try:
        sp_cfg = dict(SPAGCN_DEFAULTS)
        sp_cfg.update(((cfg.get("domains") or {}).get("spagcn") or {}))
        # l = array 坐标下的最近邻距离中位数。SpaGCN 的权重是
        # exp(-d²/(2l²))，所以 l 就是"多远算邻居"的长度尺度。
        #
        # **实测 Visium 淋巴结上是 √2 ≈ 1.414，不是 1.0。** 六边形网格在
        # array(row, col) 索引下的 6 个直接邻居落在 (±1,±1) 与 (0,±2) 上，
        # 距离是 √2、√2、2 各两对 —— 中位最近邻距离因此是 √2。
        # （第一版按"最近邻距离 = 1"写了注释，实测才发现是 √2。值本身是
        # 从数据算的，所以注释错了不影响结果，但会误导读的人。）
        #
        # 取中位最近邻距离的效果：直接邻居权重 0.61、次近邻（d=2）0.37 ——
        # 与内置邻居图"只连直接邻居"的口径接近但不相同。
        from scipy.spatial import cKDTree
        dd, _ = cKDTree(np.c_[x, y]).query(np.c_[x, y], k=2)
        l_scale = float(np.median(dd[:, 1]))
        info["length_scale_l"] = round(l_scale, 4)
        info["n_clusters"] = int(n_clusters)
        info["params"] = {k: v for k, v in sp_cfg.items()}

        adj_d = spg.calculate_adj_matrix(x=x.tolist(), y=y.tolist(),
                                        histology=False)
        X = adata.X
        X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
        # SpaGCN 内部走 `adata.X.A`（sparse 的 .A 在新 scipy 上已不保证存在），
        # 所以直接喂稠密矩阵 —— 等价，且不依赖那个属性。
        # obs 带上 array 坐标：官方用法就是传完整 adata，只给 X 会让它在
        # 需要 obs 的分支（如 ez_mode 的绘图）炸掉。
        sp_ad = sc.AnnData(X=X.astype(np.float32),
                           obs=adata.obs[["array_row", "array_col"]].copy())

        clf = spg.SpaGCN()
        clf.set_l(l_scale)
        # ---- 在 train 之前钉住三个全局 RNG ---------------------------------
        #
        # **`set_seed(cfg)` 管不到这里。** 它只 seed 了 `random` 和
        # `np.random`，**没有 seed `torch`** —— 而 SpaGCN 的 GCN 训练是
        # torch 的：权重初始化、dropout 都走 torch 的 RNG。
        #
        # 实测（四轮 CI，同一份代码、同一批包版本，SpaGCN 与内置划分的 ARI）：
        #
        #   run         env                        ARI      NMI
        #   35488906157 无钉                      0.3734   0.5535
        #   35489064432 +OMP/OB/MKL/CORETYPE       0.3784   0.5557
        #   35489172104 同上                       0.3639   0.5310
        #   35489534090 +NUMBA_NUM_THREADS=1       0.4216   0.5762
        #
        # **四轮四个值，钉 BLAS 与 Numba 都没用** —— 因为变量根本不在这里。
        # 上游的 Moran's I（0.7670）、域数（13）、平滑邻居同域率（0.668）
        # 四轮全部逐位相同，只有这一块在漂。
        #
        # 读 SpaGCN 1.2.7 源码定位到两处：
        #
        # 1. `models.py:55` `KMeans(self.n_clusters, n_init=20)` ——
        #    **没有 `random_state`**，走全局 numpy 遗留 RNG。
        # 2. `SpaGCN.py` 的 `train()` 训练 GCN，torch RNG 从未被 seed。
        #
        # SpaGCN 自己知道要 seed：`util.search_res()` 与
        # `ez_mode.detect_spatial_domains_ez_mode()` 开头都有
        # `random.seed(r_seed); torch.manual_seed(t_seed); np.random.seed(n_seed)`。
        # **但本仓库走的是 `init="kmeans"` + 外部给定 `n_clusters` 那条路**
        # （为了绕开 louvain 的 py3.12 编译链，见 `why_kmeans`），
        # 两个函数都不经过 —— 于是它的 seeding 全部被跳过。
        #
        # 所以这里照抄它自己的做法，在 train 之前把三个全局 RNG 钉死。
        # **实测这不足以让 ARI 稳定**（第五轮 CI 仍给 0.3708，又是一个新值），
        # 详见 domain_status.json 的 reproducibility.falsified_hypotheses。
        # 留着它是因为成本为零、方向正确，且 `seeded_before_train` 可核对。
        _seed = int((cfg.get("analysis") or {}).get("seed", 0))
        random.seed(_seed)
        np.random.seed(_seed)
        try:
            import torch
            torch.manual_seed(_seed)
            # torch 的 CPU 线程数也影响归约顺序（和 BLAS 是两套）
            torch.set_num_threads(1)
            info["torch_seeded"] = True
        except Exception as exc:  # noqa: BLE001
            # 没有 torch 就说明 SpaGCN 也跑不起来，不该走到这里；
            # 真走到了要如实记录，不能假装 seed 成功
            info["torch_seeded"] = False
            info["torch_seed_note"] = f"{type(exc).__name__}: {exc}"
        info["seeded_before_train"] = ["random", "numpy"] + (
            ["torch"] if info.get("torch_seeded") else [])
        info["seed"] = _seed
        clf.train(sp_ad, adj_d,
                  num_pcs=int(sp_cfg["num_pcs"]),
                  lr=float(sp_cfg["lr"]),
                  max_epochs=int(sp_cfg["max_epochs"]),
                  init="kmeans", n_clusters=int(n_clusters),
                  init_spa=True, tol=float(sp_cfg["tol"]))
        y_pred, prob = clf.predict()
        labels = np.asarray([str(int(v)) for v in y_pred])
        info["status"] = "ok"
        info["n_domains"] = int(len(set(labels.tolist())))
        info["mean_max_prob"] = round(float(np.max(prob, axis=1).mean()), 4)
        log(f"  §3.2 SpaGCN 完成：{info['n_domains']} 域，"
            f"平均最大后验 {info['mean_max_prob']}，"
            f"l={l_scale:.3f}，{sp_cfg['max_epochs']} epochs 上限")
        return labels, prob, info
    except Exception as exc:  # noqa: BLE001
        info["status"] = "failed"
        info["reason"] = f"{type(exc).__name__}: {exc}"
        log_warn(f"  §3.2 SpaGCN 跑失败：{info['reason']}")
        return None, None, info


def try_stagate(adata, cfg: dict, log=log_info) -> tuple:
    """跑 STAGATE（文档 §3.2），返回 `(labels | None, info)`。

    **本环境里它装不上，这条路径大概率不会执行** —— 但代码路径必须
    在位并且如实记状态，否则"没装"和"没写"分不开。
    """
    info = {"attempted": True, "status": None, "reason": "", "section": "§3.2"}
    try:
        from STAGATE_pyG import Cal_Spatial_Net, train_STAGATE
    except Exception as exc:  # noqa: BLE001
        info["status"] = "package_missing"
        info["reason"] = (f"STAGATE_pyG 未安装（{type(exc).__name__}: {exc}）—— "
                          f"PyPI 上 STAGATE / STAGATE_pyG / stagate 三个名字"
                          f"全部 404；官方只发 GitHub，而它的 gat_conv.py "
                          f"模块级 `from torch_sparse import SparseTensor, set_diag`，"
                          f"torch-sparse 在 PyPI 上只有 sdist（0 个 wheel），"
                          f"要按 torch 版本编译")
        log(f"  §3.2 STAGATE 未跑：{info['reason'][:90]}")
        return None, info

    info["version"] = pkg_version("STAGATE_pyG")
    try:
        st_cfg = ((cfg.get("domains") or {}).get("stagate") or {})
        n_epochs = int(st_cfg.get("n_epochs", 500))
        work = adata.copy()
        Cal_Spatial_Net(work, rad_cutoff=float(st_cfg.get("rad_cutoff", 150)))
        train_STAGATE(work, n_epochs=n_epochs,
                      random_seed=int((cfg.get("analysis") or {}).get("seed", 0)))
        sc.pp.neighbors(work, use_rep="STAGATE", random_state=cfg["analysis"]["seed"])
        sc.tl.leiden(work, resolution=float(st_cfg.get("resolution", 0.5)),
                     key_added="_stagate", flavor="igraph", n_iterations=2,
                     directed=False, random_state=cfg["analysis"]["seed"])
        labels = work.obs["_stagate"].astype(str).values
        info["status"] = "ok"
        info["n_domains"] = int(len(set(labels.tolist())))
        info["n_epochs"] = n_epochs
        info["clustering"] = ("STAGATE 官方用 mclust（R）；本仓库不装 R，"
                              "所以对它的嵌入跑 Leiden —— 这一步不是 STAGATE "
                              "原版流程的一部分")
        log(f"  §3.2 STAGATE 完成：{info['n_domains']} 域（{n_epochs} epochs）")
        return labels, info
    except Exception as exc:  # noqa: BLE001
        info["status"] = "failed"
        info["reason"] = f"{type(exc).__name__}: {exc}"
        log_warn(f"  §3.2 STAGATE 跑失败：{info['reason']}")
        return None, info


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
    # **M8：图例必须放到画布外。** 原先这里是 `axes[0].legend(fontsize=7)`
    # —— 面板内图例违反约定 v2（画布右侧、单列竖排），而且会**压住数据**：
    # 平滑扫描的两条线正好穿过左上角（`neighbor_same_frac` 从高往低走，
    # 面板内图例默认 loc="best" 就落在那里）。门禁抓不到这条 ——
    # `check_legend_convention.mjs` 的文件头明确写了它**不判**
    # "有没有写 legend.position"（只判写法）。所以这里靠人读代码发现。
    # 改成 `fig.legend` + `outside right center`，与同文件 768/816 行一致。
    fig.legend(fontsize=6, ncol=1, loc="outside right center",
               bbox_to_anchor=(1.0, 0.5))
    axes[1].plot(scan_df["smoothing"], scan_df["n_domains"], "o-", color=PAL["primary"])
    axes[1].set_xlabel("smoothing strength α")
    axes[1].set_ylabel("number of domains")
    axes[1].set_title("Domain count vs smoothing")
    save_fig(cfg, "03-03-01-unit1-smoothing-scan", fig)

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
        # **手工画 dotplot，不用 `sc.pl.dotplot`。** 用户 2026-09-24 第五轮
        # 指出的四个问题与 scrna 侧同源（详见那边注释）：
        # ① 标度矛盾 —— `standard_scale="var"` 实为逐基因 min-max（0-1），
        #    不是 z-score。改为真 z-score + RdBu_r 对称色标。
        # ② 基因名左缘裁切（MS4A4X/UBA1B 是 10x 参考自带的符号，保留原名；
        #    裁切问题由自绘布局解决）。
        # ③ 副标题截断、Y 轴没标 Spatial domain → 自绘 + 显式 ylabel。
        # ④ 图例圆点粘连 → 大小图例单独轴、间距显式。
        # 用 **log 化的全基因集**（adata.raw，规则 9）：marker 多为低表达
        raw = adata.raw.to_adata() if adata.raw is not None else adata
        sub = raw[:, top3]
        X = np.asarray(sub.X.todense()) if hasattr(sub.X, "todense") else np.asarray(sub.X)
        groups = adata.obs["domain"].astype(str).values
        ud = sorted(set(groups), key=lambda v: int(v))
        frac = np.zeros((len(ud), len(top3)))
        mean_expr = np.zeros((len(ud), len(top3)))
        for ri in range(len(ud)):
            blk = X[groups == ud[ri]]
            frac[ri] = (blk > 0).mean(axis=0)
            mean_expr[ri] = blk.mean(axis=0)
        sd = mean_expr.std(axis=0, ddof=0)
        sd[sd == 0] = 1.0
        zmat = (mean_expr - mean_expr.mean(axis=0)) / sd

        frac_df = pd.DataFrame(frac, index=ud, columns=top3)
        z_df = pd.DataFrame(zmat, index=ud, columns=top3)

        # **图例放主图下方的横带（Seurat do_DotPlot 范式，用户第七/八轮反馈）。**
        # 整图构建抽在 `common.build_marker_dotplot_figure` —— 前七轮把这段
        # 内联在脚本里、验证脚本又照抄一份，三份镜像不同步，导致"改了没效果"
        # 与七轮返工（详见该函数 docstring）。
        fig, size_handles = build_marker_dotplot_figure(
            frac_df, z_df,
            group_label="Spatial domain",
            title="Top markers per spatial domain",
            subtitle="rows = spatial domains (histology labels: domain_labels.csv)\ndot size = fraction of spots expressing the gene")

        save_fig(cfg, "03-03-02-unit1-domain-markers-dotplot", fig)

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
    annot_error = None
    sig_key = ((cfg.get("deconvolution") or {}).get("signature_set")
               or (cfg.get("analysis") or {}).get("celltype_markers"))
    try:
        import yaml
        sig_path = (Path(__file__).resolve().parent.parent / "assets"
                    / "reference_signatures.yml")
        with open(sig_path, encoding="utf-8") as fh:
            sigs = yaml.safe_load(fh).get("signatures", {})
        sig = sigs.get(sig_key)
        if not sig:
            raise KeyError(f"signature_set='{sig_key}' 不在 {sig_path.name} 里"
                           f"（可选: {', '.join(sorted(sigs))}）")
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
        if M.size == 0:
            raise ValueError("没有任何细胞类型的 marker 覆盖 >= 3 个")

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
    except Exception as e:  # noqa: BLE001
        # **S2：不再静默吞掉。** 这一段原先整块包在裸 except 里，
        # 失败只 `log_warn`，而 `domain_annotation.json` 在
        # `main_analysis.py` 里**没有任何消费者** —— 于是"域标签没算出来"
        # 和"域标签算出来了"在验收层看起来完全一样。
        # 现在：① 把失败原因写进 status 文件（`domain_annotation_failed`）；
        # ② 让 main_analysis 有一条 `content:domain_annotation` 消费它。
        annot_error = f"{type(e).__name__}: {e}"
        log_warn(f"域标签打分跳过: {annot_error}")

    if annot:
        write_json(res_dir / "domain_annotation.json",
                   {"signature_set": sig_key, "domains": annot,
                    "status": "ok",
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
    else:
        # **必须留下一个文件。** 没有文件 = 验收层那条检查会以为
        # "这一步不存在"，而不是"这一步失败了"（见 S2/S4 的教训）。
        write_json(res_dir / "domain_annotation.json",
                   {"signature_set": sig_key, "domains": {},
                    "status": "failed" if annot_error else "empty",
                    "reason": annot_error or "没有任何细胞类型的 marker 覆盖 >= 3 个",
                    "n_celltypes_used": 0,
                    "caveat": "**域标签没有产出** —— 不要把它当『这些域没有标签』"})
        log_warn(f"域标签没有产出（{annot_error or '无可用类型'}）—— "
                 f"已写 domain_annotation.json 的 status=failed")

    # ---- 4c. §3.2 点名方法：SpaGCN / STAGATE（交叉验证）----------------------
    #
    # **域数对齐再比。** SpaGCN 用 kmeans 初始化时必须外部给 n_clusters，
    # 我们传内置方法的域数 —— 两边域数相同，ARI 才在比"划分方式"而不是
    # 在比"谁切得细"。
    n_clusters_for_named = int(m_spat["n_domains"])
    spg_labels, spg_prob, spg_info = try_spagcn(adata, cfg, n_clusters_for_named)
    stg_labels, stg_info = try_stagate(adata, cfg)

    method_agree = {}
    if spg_labels is not None:
        adata.obs["domain_spagcn"] = pd.Categorical(spg_labels)
        pd.DataFrame({
            "barcode": adata.obs_names,
            "domain": adata.obs["domain"].astype(str).values,
            "domain_spagcn": spg_labels,
            "spagcn_max_prob": np.round(np.max(spg_prob, axis=1), 5),
        }).to_csv(res_dir / "spagcn_domains.csv", index=False)
        method_agree["SpaGCN_vs_builtin"] = _agreement(
            adata.obs["domain"].astype(str).values, spg_labels, adj, adata)
    if stg_labels is not None:
        adata.obs["domain_stagate"] = pd.Categorical(stg_labels)
        pd.DataFrame({
            "barcode": adata.obs_names,
            "domain": adata.obs["domain"].astype(str).values,
            "domain_stagate": stg_labels,
        }).to_csv(res_dir / "stagate_domains.csv", index=False)
        method_agree["STAGATE_vs_builtin"] = _agreement(
            adata.obs["domain"].astype(str).values, stg_labels, adj, adata)
    if spg_labels is not None and stg_labels is not None:
        method_agree["SpaGCN_vs_STAGATE"] = _agreement(
            spg_labels, stg_labels, adj, adata)

    # ---- 5. 叠到 H&E 上（唯一能判断域是否对应组织学的方法）-------------------
    lib = list(adata.uns["spatial"].keys())[0]
    entry = adata.uns["spatial"][lib]
    img = entry["images"]["hires"]
    sf = float(entry["scalefactors"]["tissue_hires_scalef"])
    r_plot = spot_radius_plot_units(entry["scalefactors"], "hires")
    xy = spatial_xy(adata, sf)

    # **单图原则拆分（D-006）**：三面板 -> 3 张独立单图（P3 域划分对照链）：
    #   unit1 = 表达-only 聚类（证明"不平滑的聚类不知道空间"）
    #   unit2 = 空间域（平滑后贴合组织）
    #   unit3 = H&E 参考（生物学裁决依据）
    # 三张图共享 xy 与 tab20 域色（同色纪律，跨图可对照）。
    PANEL_SPECS = [("domain_expr_only", f"Expression-only clustering ({m_expr['n_domains']} clusters)",
                    "03-03-03-unit1-domains-expr-only"),
                   ("domain", f"Spatial domains, α={alpha_used} ({m_spat['n_domains']} domains)",
                    "03-03-03-unit2-domains-spatial")]
    for key, title, name in PANEL_SPECS:
        cats = adata.obs[key].astype(str).values
        uniq = sorted(set(cats), key=lambda x: int(x) if x.isdigit() else x)
        cmap = plt.get_cmap("tab20")
        fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(56)))
        ax.imshow(img, alpha=0.55)
        for i, u in enumerate(uniq):
            m = cats == u
            ax.scatter(xy[m, 0], xy[m, 1], s=8, color=cmap(i % 20),
                       label=u, linewidths=0)
        fig.legend(fontsize=6, markerscale=2.2, loc="outside right center",
                  ncol=1, framealpha=0.8)
        ax.set_title(title)
        ax.set_xticks([]); ax.set_yticks([])
        save_fig(cfg, name, fig)
    # unit3：H&E 参考（无散点，纯组织学底图）
    fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(56)))
    ax.imshow(img, alpha=1.0)
    ax.set_title("H&E reference")
    ax.set_xticks([]); ax.set_yticks([])
    save_fig(cfg, "03-03-03-unit3-he-reference", fig)

    # ---- 5b. 方法对照图（只在点名方法真的跑了时才画）------------------------
    #
    # **"跑通了"不是结论，两套划分是否一致才是。** 所以图和方法名放在一起，
    # 让读者一眼看到 ARI 之外的东西：哪一块组织两套方法分歧最大。
    if spg_labels is not None or stg_labels is not None:
        # **单图原则拆分（D-006）**：方法对照 -> 每方法一张独立单图
        # （P3 组，unit4 起编号）。"哪一块组织两套方法分歧最大"由
        # 读者跨图对照（同 tab20 域色 + 同一 H&E 底图保证可比）；
        # ARI 汇总写在每张图副标题里（拆图后 suptitle 不存在了）。
        sub = " / ".join(f"{k}: ARI={v['adjusted_rand_index']}"
                         for k, v in method_agree.items())
        panels = [("domain", f"Builtin smooth+Leiden ({m_spat['n_domains']})",
                   "03-03-04-unit1-domains-builtin")]
        if spg_labels is not None:
            panels.append(("domain_spagcn",
                           f"SpaGCN, kmeans init ({spg_info['n_domains']})",
                           "03-03-04-unit2-domains-spagcn"))
        if stg_labels is not None:
            panels.append(("domain_stagate",
                           f"STAGATE ({stg_info['n_domains']})",
                           "03-03-04-unit3-domains-stagate"))
        for _ui, (key, title, name) in enumerate(panels, start=1):
            # **方法身份靠标题，不靠颜色。** 这里原先写了一段
            # `METHOD_COLORS = {...}` + `method_fill = ...`，声称"每个方法一个
            # 固定色"——但两个变量**赋值后从未被使用**，散点实际用的是
            # `cmap(i % 20)`（按**域标签**上色）。留着它比删掉更危险：
            # 读者会以为方法色在生效，于是按"色=方法"去读图。
            # 真实约定是：**域标签色跨图一致**（tab20，见上面 PANEL_SPECS
            # 的"同色纪律"），方法身份由 `title` 与图名承载。
            # 想真的按方法上色，就得把散点改成单色 scatter —— 那会丢掉
            # 域间对照，是另一个决定，不该藏在死变量里。
            cats = adata.obs[key].astype(str).values
            uniq = sorted(set(cats), key=lambda x: int(x) if x.isdigit() else x)
            cmap = plt.get_cmap("tab20")
            fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(52)))
            ax.imshow(img, alpha=0.55)
            for i, u in enumerate(uniq):
                m = cats == u
                ax.scatter(xy[m, 0], xy[m, 1], s=8, color=cmap(i % 20),
                           label=u, linewidths=0)
            fig.legend(fontsize=5, markerscale=2.0, loc="outside right center",
                      ncol=1, framealpha=0.8,
                      title="domain", title_fontsize=6)
            ax.set_title(title)
            ax.set_xticks([]); ax.set_yticks([])
            save_fig(cfg, name, fig)
        log_info(f"方法对照已拆分单图输出；一致性: {sub}")
        # **域×域对应矩阵热图**（差距清单 #25）：把方法对照的 ARI 拆成
        # "哪几个域对上了"。每对方法一张（同图号 unit7+）。
        DYNAMIC_FIG_BASES = {"04": 3}
        JAC_BASE = "-".join(["03", "03", "04", "unit"])
        for ji, (key, title, name) in enumerate(panels, start=1):
            pair = None
            for kk, vv in method_agree.items():
                if "jaccard_matrix" in vv and (
                        ("SpaGCN" in kk and key == "domain_spagcn") or
                        ("STAGATE" in kk and key == "domain_stagate") or
                        (kk.startswith("SpaGCN_vs_STAGATE") and key == "domain_stagate")):
                    pair = vv
                    break
            if pair is None:
                continue
            jm = np.array(pair["jaccard_matrix"]["values"])
            fig_j, ax_j = plt.subplots(figsize=(W_ONE_HALF, mm(64)))
            im_j = ax_j.imshow(jm, cmap="magma", vmin=0, vmax=1, aspect="auto")
            ax_j.set_xticks(range(len(pair["jaccard_matrix"]["cols"])))
            ax_j.set_xticklabels(pair["jaccard_matrix"]["cols"], fontsize=6, rotation=90)
            ax_j.set_yticks(range(len(pair["jaccard_matrix"]["rows"])))
            ax_j.set_yticklabels(pair["jaccard_matrix"]["rows"], fontsize=6)
            ax_j.set_xlabel(f"{title} domain")
            ax_j.set_ylabel("Builtin domain")
            ax_j.set_title("Domain-by-domain Jaccard overlap - bright off-diagonal blocks = split/merged domains")
            fig_j.colorbar(im_j, ax=ax_j, shrink=0.8, pad=0.02,
                           fraction=0.046, label="Jaccard")
            save_fig(cfg, JAC_BASE + str(ji + 6) + "-domain-jaccard", fig_j)

    # ---- 6. 落盘 ------------------------------------------------------------
    adata.obs[["domain", "domain_expr_only"]].to_csv(res_dir / "spatial_domains.csv")
    out = data_dir / "domains.h5ad"
    adata.write_h5ad(out)
    log_info(f"已写出 {out}")

    # ---- 7. 点名工具的落地登记（§3.2）---------------------------------------
    #
    # 文档 §3.2 点名 BayesSpace / STAGATE / SpaGCN。**这三条现在是三件不同的事**，
    # 不能再用一句话概括：
    #
    #   SpaGCN    —— 真的跑了（`try_spagcn`），结果在 `domain_methods` 里，
    #                与内置划分的一致性在 `method_agreement` 里
    #   STAGATE   —— 代码路径在位，但包装不上（PyPI 三个名字全 404 +
    #                torch-sparse 只有 sdist）；`domain_methods` 里记原因
    #   BayesSpace—— R/Bioconductor 包，本仓库 CI 不装 R + rpy2，结构上跑不了
    #
    # **SpaGCN 已从 `NAMED_TOOLS` 移除**（它不再是"用不了的工具"）；
    # 留在这里反而会让 `probe_named_tools` 每次都报"登记过期"。
    named = probe_named_tools(log=log_warn, only=("BayesSpace", "STAGATE"))
    for tool, tinfo in named.items():
        log_info(f"  §3.2 {tool}: 未使用（{tinfo['kind']}）—— {tinfo['reason'][:70]}")

    domain_methods = {
        "builtin_smooth_leiden": {
            "used": True, "role": "primary", "section": "§3.2",
            "status": "ok", "n_domains": m_spat["n_domains"],
            "note": ("本仓库的平滑 + Leiden（`flavor='igraph'`）。"
                     "**这不是文档点名的三个方法中的任何一个** —— "
                     "它是可复现、无外部依赖的基线"),
        },
        "SpaGCN": dict(spg_info, used=spg_info.get("status") == "ok"),
        "STAGATE": dict(stg_info, used=stg_info.get("status") == "ok"),
        "BayesSpace": dict(named["BayesSpace"], used=False),
    }

    # 验收用的"指名工具落地情况"：每条都要有 used 和（未用时）reason
    n_used = sum(1 for v in domain_methods.values() if v.get("used"))
    log_info(f"  §3.2 落地情况：{n_used}/{len(domain_methods)} 个方法实际产出结果")

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
        # ---- §3.2 点名工具的**逐条落地状态**（不是"一个都没跑"）------------
        "domain_methods": domain_methods,
        "method_agreement": method_agree,
        "named_tools": named,
        "named_tools_note": (
            "§3.2 点名的三个方法里，**SpaGCN 已实际运行**（结果见 "
            "domain_methods / method_agreement / spagcn_domains.csv）；"
            "STAGATE 与 BayesSpace 未运行，逐条理由见 domain_methods。"
            "主方法仍是内置的平滑 + Leiden。\n"
            # 逐步骤的说明之外，再带上**整轮**的落地边界 ——
            # 读者只看 domain_status.json 时也能知道 §3.4 的 SpatialDE 跑了，
            # 不会以为"整份文档点名的工具一个都没跑"。
            + named_tools_note()),
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
        # ---- 哪些数可以当确定值报，哪些不能 --------------------------------
        #
        # **"跑通了"不等于"这个数可复现"。** 实测五轮 CI（同一份代码、
        # 同一批包版本、五种不同的并行度/种子设置），SpaGCN 与内置划分的
        # ARI 给了五个值：0.3734 / 0.3784 / 0.3639 / 0.4216 / 0.3708。
        #
        # **而同一轮里，主方法的数全部逐位相同**：Moran's I 0.7670、
        # 域数 13、不平滑邻居同域率 0.507、平滑后 0.668 —— 五轮一致。
        # 所以漂的只有 §3.2 的交叉验证那一块。
        #
        # 三次归因，三次被否证：
        #
        #   1. 多线程 BLAS 归约顺序 —— 钉住 OMP/OPENBLAS/MKL + CORETYPE
        #      之后 ARI 仍从 0.3784 变 0.3639；
        #   2. Numba 的 prange 线程数 —— 补上 NUMBA_NUM_THREADS=1 之后
        #      仍给 0.4216；
        #   3. torch 未 seed（下面这段代码就是为它加的）—— 钉住
        #      random/numpy/torch 三个全局 RNG 之后仍给 0.3708。
        #
        # **所以残留随机源尚未定位。** 三次"读源码看着很有道理"的归因
        # 都被日志否证了 —— 第四次动手之前先读日志（AGENTS 规则 13）。
        #
        # 已确认的源码事实（不是推测）：
        #   - `models.py:55` 的 `KMeans(n_clusters, n_init=20)` **没有
        #     `random_state`**，走全局 numpy 遗留 RNG；
        #   - `train()` 训练 GCN 走 torch，而本仓库的 `set_seed()`
        #     **只 seed 了 random 与 numpy，没有 seed torch**。
        # 下面把这三个 RNG 都钉住了 —— 但实测**不足以**让 ARI 稳定。
        "reproducibility": {
            "stable": ["morans_I", "n_domains", "neighbor_same_frac",
                       "resolution_scan", "smoothing_scan"],
            "unstable": ["method_agreement.SpaGCN_vs_builtin"
                         ".adjusted_rand_index"],
            "ari_observed_range": [0.3639, 0.4216],
            "evidence": ("五轮 CI（35488906157 / 35489064432 / 35489172104 / "
                         "35489534090 / 35489962167，同一份代码、同一批包版本）："
                         "ARI = 0.3734 / 0.3784 / 0.3639 / 0.4216 / 0.3708；"
                         "同五轮的 Moran's I 全部 0.7670、域数全部 13、"
                         "平滑邻居同域率全部 0.668"),
            "range_is_from": ("**历史观测值，不是本轮的** —— 本轮的 ARI 见 "
                              "`method_agreement`"),
            "cause": ("**尚未定位。** 已知 SpaGCN 1.2.7 有两处不受本仓库"
                      "控制的随机源（`models.py:55` 的 `KMeans` 没有 "
                      "`random_state`；`train()` 的 GCN 走 torch 而 "
                      "`set_seed()` 不 seed torch），但把这三个全局 RNG "
                      "都钉住之后 ARI 仍在变"),
            "falsified_hypotheses": [
                ("多线程 BLAS 归约顺序 —— 否证：钉住 OMP/OPENBLAS/MKL + "
                 "OPENBLAS_CORETYPE 之后 ARI 仍从 0.3784 变 0.3639"),
                ("Numba `prange` 的线程数 —— 否证：补上 "
                 "NUMBA_NUM_THREADS=1 之后仍给 0.4216"),
                ("torch 未 seed —— 否证：在 `clf.train()` 前钉住 "
                 "random/numpy/torch 三个全局 RNG 之后仍给 0.3708"),
            ],
            "mitigation": ("`try_spagcn` 在 `clf.train()` 之前钉住三个全局 "
                           "RNG（照抄 SpaGCN 自己 `search_res()` 的做法），"
                           "结果记在 `domain_methods.SpaGCN."
                           "seeded_before_train`。**实测这不足以让 ARI "
                           "可复现**，所以不要以为设了就可复现 —— 按范围报"),
            "how_to_report": ("主方法的数（Moran's I、域数、邻居同域率）可按"
                              "确定值报；**`SpaGCN_vs_builtin` 的 ARI 必须带"
                              "范围报**，不能只报一个数 —— 而且"
                              "『ARI 高也不等于两套方法都对』：它们可能共享"
                              "同一个错误"),
        },
        "limitations": [
            "平滑会抹平真实的微小结构（如小的生发中心）；α 过高时域数减少但可能丢细节",
            "Visium 的 spot 直径 55 μm，含 1-10 个细胞 —— 域的空间分辨率受此限制，"
            "**不能说『某个域是某一种细胞』**，只能说这个区域的细胞组成不同",
            "域边界的位置有 ±1 个 spot 的不确定性",
            "**`fragmented_domains` 高不等于聚类失败。** 淋巴结的滤泡、"
            "肿瘤的癌巢本身就是散布的多个斑块 —— 同一个域出现在多个不相邻"
            "位置是生物学事实。只有对『应当连续』的组织（脑的层状结构、"
            "上皮分层）才能把这个指标当缺陷看",
            # §3.2 的落地边界 —— 必须与"域划分做出来了"并列出现
            ("**SpaGCN 跑的是 `init=\"kmeans\"` 路径，不是默认的 "
             "`init=\"louvain\"`。** 差别在簇心初始化：louvain 用表达+空间，"
             "kmeans 只用 GCN 特征。`n_clusters` 因此是外部给定的"
             f"（传了内置方法的 {n_clusters_for_named} 个域）—— "
             "**这既是可比性的前提，也意味着 SpaGCN 的域数不是它自己选的**"),
            ("**STAGATE 与 BayesSpace 没有运行。** STAGATE 的三个 PyPI 名字"
             "全部 404、官方 GitHub 版的 `gat_conv.py` 模块级依赖 torch-sparse"
             "（PyPI 上只有 sdist、0 个 wheel）；BayesSpace 是 R/Bioconductor 包，"
             "本仓库 CI 不装 R。逐条理由见 `domain_methods`"),
            ("**ARI 高不等于两套方法都对。** 它们可能共享同一个错误"
             "（比如都被同一个技术批次效应驱动）。一致只说明"
             "『换一种方法也不会得到完全不同的域』"),
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
