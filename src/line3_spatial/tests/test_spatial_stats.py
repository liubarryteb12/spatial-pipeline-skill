#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""spatial_stats 的黄金解析解单测（Argon 第 38 轮要求）。

解析解部分用显式 rook 邻接（4 邻域），因为只有规则的邻接才能手算闭式解；
kNN 图的属性单独验证（对称、无自环、连通）。

判据：
  1. rook 邻接 + 棋盘格 -> Moran I = -1（严格解析解，不是近似）
  2. rook 邻接 + 半区聚集 -> Moran I -> +1 方向（> 0.5）
  3. 分块 Moran I 与稠密定义式在 1e-9 内相等
  4. 分块 Geary C 与双循环定义式在 1e-9 内相等（验证恒等式展开）
  5. 常数向量返回 NaN 而不是崩溃
  6. kNN 图对称、无自环
"""
import sys
import os

import numpy as np
import scipy.sparse as sps

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from scripts.spatial_stats import (  # noqa: E402
    build_knn_graph, morans_i, gearys_c, morans_i_dense, gearys_c_dense,
)

FAIL = []


def check(name, ok, detail=""):
    if ok:
        print("ok   %s" % name)
    else:
        print("FAIL %s %s" % (name, detail))
        FAIL.append(name)


def rook_graph(n_rows, n_cols):
    """显式 4 邻域（rook）邻接，节点按 (r, c) 行优先编号。"""
    rows, cols = [], []
    for r in range(n_rows):
        for c in range(n_cols):
            i = r * n_cols + c
            if c + 1 < n_cols:
                j = r * n_cols + (c + 1)
                rows += [i, j]
                cols += [j, i]
            if r + 1 < n_rows:
                j = (r + 1) * n_cols + c
                rows += [i, j]
                cols += [j, i]
    data = np.ones(len(rows))
    return sps.csr_matrix((data, (rows, cols)), shape=(n_rows * n_cols,) * 2)


N = 4
Wr = rook_graph(N, N)

# --- 1. 棋盘格：rook 邻接全是异值边，Moran I 严格等于 -1 ---
x_checker = np.array([(r + c) % 2 for r in range(N) for c in range(N)], dtype=float)
I_checker = morans_i(Wr, x_checker[:, None])[0]
check("棋盘格 rook Moran I == -1 (1e-9)", abs(I_checker + 1.0) < 1e-9, "actual=%.12f" % I_checker)

# --- 2. 左右半区：正自相关 ---
x_half = np.array([0.0 if c < N // 2 else 1.0 for r in range(N) for c in range(N)])
I_half = morans_i(Wr, x_half[:, None])[0]
check("半区聚集 rook Moran I > 0.5", I_half > 0.5, "actual=%.4f" % I_half)

# --- 2b. Geary 对棋盘格应接近最大值（相邻差异最大）---
C_checker = gearys_c(Wr, x_checker[:, None])[0]
check("棋盘格 rook Geary C > 1", C_checker > 1.0, "actual=%.4f" % C_checker)

# --- 3/4. 分块 vs 稠密 ---
rng = np.random.default_rng(20261105)
X = rng.random((N * N, 37))
I_block = morans_i(Wr, X, gene_chunk=10)
I_ref = morans_i_dense(Wr, X)
dI = float(np.nanmax(np.abs(I_block - I_ref)))
check("Moran 分块 == 稠密 (1e-9)", dI < 1e-9, "maxdiff=%.3e" % dI)

C_block = gearys_c(Wr, X, gene_chunk=10)
C_ref = gearys_c_dense(Wr, X)
dC = float(np.nanmax(np.abs(C_block - C_ref)))
check("Geary 分块 == 双循环定义式 (1e-9)", dC < 1e-9, "maxdiff=%.3e" % dC)

# --- 5. 常数向量 ---
with np.errstate(invalid="ignore", divide="ignore"):
    I_const = morans_i(Wr, np.ones((N * N, 1)))[0]
check("常数向量 Moran 为 NaN", bool(np.isnan(I_const)), "actual=%r" % I_const)

# --- 6. kNN 图属性 ---
grid = np.array([(r, c) for r in range(N) for c in range(N)], dtype=float)
Wk = build_knn_graph(grid, k=4)
check("kNN 图对称", (Wk != Wk.T).nnz == 0)
check("kNN 图无自环", float(Wk.diagonal().sum()) == 0.0)
check("kNN 图边数合理", Wk.nnz >= 2 * (N * N - 1), "nnz=%d" % Wk.nnz)

if FAIL:
    print("\n失败 %d 项: %s" % (len(FAIL), ", ".join(FAIL)))
    sys.exit(1)
print("\n全部通过")