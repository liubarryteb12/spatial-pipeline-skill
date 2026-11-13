"""spatial_stats：Moran's I / Geary's C 的自实现分块算子（Line 3 公共核心库）。

为什么不用 squidpy：
  1. squidpy 1.8.3 已移除 morans_i / gearys_c（只剩 neighbors 系列），装了也调不到；
  2. 即便有，整矩阵稠密化在 12,560 spot x 16,765 基因上是 1.68 GB float64，
     实测直接触发 ArrayMemoryError 并把手机打崩。

数学定义（与 squidpy / PySAL 一致）：
  Moran: I = (n / W) * (x' W x) / (x' x)
  Geary: C = ((n - 1) / (2W)) * sum_ij w_ij (x_i - x_j)^2 / sum_i (x_i - mean)^2

  x 为已中心化的表达向量，W 为空间权重矩阵，W = sum_ij w_ij，n 为 spot 数。

Geary 分子的分块技巧（关键实现细节，单测会验证与定义式双循环严格相等）：
  sum_ij w_ij (x_i - x_j)^2 = sum_i d_i x_i^2 - 2 x' W x + sum_j e_j x_j^2
  其中 d_i 为 W 的第 i 行和、e_j 为第 j 列和。
  这样循环体内只需要 O(n_spot x 块宽) 的稠密阵，不必构造「边 x 基因」的大矩阵
  （后者在 12,560 spot x 6 邻域上是 7.5 万边 x 1.6 万基因 = 1.57 GB）。
"""

import numpy as np
import scipy.sparse as sps
from scipy.spatial import cKDTree

GENE_CHUNK = 1000


def build_knn_graph(coords, k=6):
    """按坐标连最近 k 个点（对称化）。返回 CSR 权重矩阵，权重为 1。

    对称化用 maximum 而不是加法：两个点只要互为近邻就保留一条边，
    避免某些点因为「被更多人选中」而权重虚高。
    """
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("coords 必须是 (n, 2)，实际 %r" % (coords.shape,))
    n = coords.shape[0]
    if n < 2:
        raise ValueError("点数 %d < 2，无法建图" % n)
    k_eff = min(int(k) + 1, n)
    dist, idx = cKDTree(coords).query(coords, k=k_eff)
    dist = np.atleast_2d(dist)
    idx = np.atleast_2d(idx)
    rows = np.repeat(np.arange(n), k_eff - 1)
    cols = idx[:, 1:].ravel()
    vals = dist[:, 1:].ravel()
    conn = sps.csr_matrix((np.ones_like(vals), (rows, cols)), shape=(n, n))
    return conn.maximum(conn.T)


def morans_i_dense(W, X):
    """参考实现：把 X 当稠密矩阵直接按定义算，用于小规模单测对拍。"""
    W = sps.csr_matrix(W, dtype=np.float64)
    X = np.asarray(X, dtype=np.float64)
    if X.ndim == 1:
        X = X[:, None]
    n = W.shape[0]
    w_sum = float(W.sum())
    Xc = X - X.mean(axis=0, keepdims=True)
    denom = (Xc * Xc).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        num = (W @ Xc * Xc).sum(axis=0)
        return (n / w_sum) * num / denom


def gearys_c_dense(W, X):
    """参考实现：按定义式显式遍历每条边，用于小规模单测对拍。"""
    W = sps.coo_matrix(W, dtype=np.float64)
    X = np.asarray(X, dtype=np.float64)
    if X.ndim == 1:
        X = X[:, None]
    n = W.shape[0]
    w_sum = float(W.sum())
    Xc = X - X.mean(axis=0, keepdims=True)
    denom = (Xc * Xc).sum(axis=0)
    num = np.zeros(X.shape[1], dtype=np.float64)
    for i, j, w in zip(W.row, W.col, W.data):
        diff = Xc[i, :] - Xc[j, :]
        num += w * diff * diff
    with np.errstate(invalid="ignore", divide="ignore"):
        return ((n - 1) / (2.0 * w_sum)) * num / denom


def morans_i(W, X, gene_chunk=GENE_CHUNK):
    """分块版 Moran's I，返回形状 (n_genes,) 的数组。"""
    return _chunked(W, X, gene_chunk)[0]


def gearys_c(W, X, gene_chunk=GENE_CHUNK):
    """分块版 Geary's C，返回形状 (n_genes,) 的数组。"""
    return _chunked(W, X, gene_chunk)[1]


def _chunked(W, X, gene_chunk=GENE_CHUNK):
    W = sps.csr_matrix(W, dtype=np.float64)
    n = W.shape[0]
    w_sum = float(W.sum())
    if w_sum == 0:
        raise ValueError("权重矩阵全零，统计量无定义")
    d_row = np.asarray(W.sum(axis=1)).ravel()
    d_col = np.asarray(W.sum(axis=0)).ravel()
    if sps.issparse(X):
        Xc_all = X.tocsc()
    else:
        Xc_all = sps.csc_matrix(np.asarray(X, dtype=np.float64))
    n_genes = Xc_all.shape[1]
    means = np.asarray(Xc_all.mean(axis=0)).ravel()
    I = np.full(n_genes, np.nan, dtype=np.float64)
    C = np.full(n_genes, np.nan, dtype=np.float64)
    for lo in range(0, n_genes, gene_chunk):
        hi = min(lo + gene_chunk, n_genes)
        Xc = np.asarray(Xc_all[:, lo:hi].todense(), dtype=np.float64)
        Xc -= means[lo:hi]
        denom = (Xc * Xc).sum(axis=0)
        WX = W @ Xc
        num_moran = (WX * Xc).sum(axis=0)
        sqx = Xc * Xc
        with np.errstate(invalid="ignore", divide="ignore"):
            I[lo:hi] = (n / w_sum) * num_moran / denom
            sq_term = (d_row[:, None] * sqx).sum(axis=0) - 2.0 * num_moran \
                      + (d_col[:, None] * sqx).sum(axis=0)
            C[lo:hi] = ((n - 1) / (2.0 * w_sum)) * sq_term / denom
    return I, C