---
name: spatial-transcriptomics-pipeline
description: Run an end-to-end Visium spatial transcriptomics pipeline — fetch matrix plus spatial data, spot QC with tissue-aware thresholds, normalization, spatially-smoothed domain detection overlaid on H&E, Moran's I spatial variable genes, marker-based spot composition, neighborhood enrichment with permutation testing, and spatially-constrained ligand-receptor analysis — and ship it as a GitHub Actions workflow that uploads results as an artifact. Use when the user asks for spatial transcriptomics, Visium analysis, spatial domains, spatial variable genes, Moran's I, spot deconvolution, neighborhood enrichment, niche analysis, spatial cell-cell communication, or a 10x Genomics spatial workflow.
---

# 空间转录组（Visium）分析流水线

## 这个 skill 做什么

从 10x Visium 的**两半数据**（表达矩阵 + 空间信息）出发，跑到
空间域、空间高变基因、spot 组成、空间邻域与空间通讯，
全程在 GitHub Actions 上跑，产物作为 artifact 下载。

参考数据集：**正常人淋巴结**（10x 官方 Visium，4035 spot）。

## 为什么需要单独的流水线（不是把单细胞流程改改）

**空间数据的价值全在"位置"上，而位置有三类特有的坑：**

1. **数据是两半的。** 只读表达矩阵，你有一个普通单细胞数据集，
   但**以为自己有空间数据**。流水线强制要求两半都给。

2. **坐标错了不报错。** `obsm['spatial']` 列序写反（转置）时，
   散点图仍然画出"看起来像组织形状"的点阵。只有叠到 H&E 上、
   或做定量网格几何检查才发现。

3. **"表达簇"不等于"空间域"。** 不看位置的聚类会把连续区域切碎、
   把不相邻区域合并。必须把表达与空间邻居平均后再聚类。

## 流程

| 步骤 | 脚本 | 做什么 |
|---|---|---|
| 0 | `00_fetch.py` | 下载矩阵 + 空间信息，建 `obsm['spatial']` 与 `uns['spatial']` |
| 1 | `01_qc.py` | spot 质控（组织特异阈值）、过滤、**组织连通性检查** |
| 2 | `02_normalize.py` | HVG、标准化、PCA、**PC 的空间投影诊断** |
| 3 | `03_spatial_domains.py` | **空间域**（平滑 vs 不平滑对照）+ 域标签 |
| 4 | `04_svg.py` | **Moran's I / Geary's C** + 置换检验 |
| 5 | `05_deconvolution.py` | spot 组成（marker 打分 / 真解卷积） |
| 6 | `06_niche.py` | **邻域富集**（置换检验 z-score）+ 共现曲线 |
| 7 | `07_spatial_communication.py` | **空间约束**的配体-受体分析 |
| — | `main_analysis.py` | 串起来 + 26 项验收检查 |

## 怎么用

### 跑官方淋巴结数据

```bash
pip install -r requirements.txt
python scripts/main_analysis.py --config assets/config.lymph_node.yml
```

### 跑自己的数据

1. 在 `assets/datasets.yml` 里加一项（需要 `matrix_url` + `spatial_url`）
2. 复制 `assets/config.lymph_node.yml`，改 `dataset_id` 与 `source.dataset`
3. **改 `qc.max_pct_mt`** —— 阈值是组织特异的：
   PBMC 5-10 / 淋巴结 20 / 心肌 30-50 / 肝 30+
4. 跑

### 云端

推送到 GitHub，`.github/workflows/spatial_analysis.yml` 自动跑，
产物在 Actions 页面的 artifact 里。

## 产物

**17 张图**（`results/<dataset_id>/figures/`），关键的几张：

- `domains_on_he.png` —— 空间域叠在 H&E 上。**这是判断域划分是否
  对应真实组织学结构的唯一方法。**
- `pca_on_tissue.png` —— PC 的空间投影。有空间结构才说明主成分
  抓到了组织学差异而不是技术噪声。
- `smoothing_scan.png` —— 平滑强度 vs 空间连贯性。让"选 α=0.5"有依据。
- `svg_top_genes.png` —— top 空间高变基因的空间分布。
- `niche_enrichment_celltypes.png` —— 邻域富集 z-score 热图。
- `deconvolution_spatial.png` —— 各细胞类型的相对权重空间分布。

**24 个结果文件**（JSON/CSV），每个都带 `method` 与 `limitations`。

## 实测结果（淋巴结，本地与 CI 一致）

- 4035 spot × 36601 基因 → QC 后 **4025 spot**，组织连通性 1 块
- 空间域：不平滑 9 域（邻居同域率 0.536，随机基线 0.133）；
  平滑 α=0.5 后 **13 域**（同域率 **0.672**）
- 域标签 10 种，`B_germinal_center`(z=+3.20)、`Plasma_cell`(+2.50)、
  `T_cell`(+2.78)、`Smooth_muscle`(+2.57)
- SVG：4000 个基因里 **2890 个** BH 校正后 p<0.05；top 是
  `IGKC`、`CCL21`、`FDCSP` —— **正是淋巴结里有已知空间结构的基因**
- 邻域富集：`FRC ↔ T_cell` **z=13.8**（T 细胞区支架）、
  `B_naive ↔ FDC` **z=9.4**（滤泡）、`B_naive ↔ FRC` **z=−28.5**
  （**B 细胞滤泡与 T 细胞区互斥** —— 淋巴结的分区结构）
- 空间通讯：**62/63 对 LR 可用**；top 富集是
  `CCL21-CCR7`(z=3.31)、`CCL19-CCR7`(z=2.81) —— **定义淋巴结
  T 细胞区招募的两条趋化因子轴**
- 全流程 **约 2 分钟**（本地）

## 这些结论是怎么被独立验证的

流水线没有注入任何淋巴结的先验知识，但独立复现了它的经典结构：

1. **SVG top 基因**是 `CCL21`（FRC 分泌的 T 区趋化因子）、
   `FDCSP`（生发中心滤泡树突细胞）、免疫球蛋白基因（浆细胞）
2. **邻域富集**发现 B 细胞滤泡与 T 细胞区**互相排斥**
   （z=−28.5 / −24.6），而 FRC 与 T 细胞**紧密共定位**（z=13.8）
3. **空间通讯** top 是 CCL19/CCL21–CCR7 —— 就是驱动上述分区的机制

三条独立分析指向同一个已知生物学，这是实现正确性的强证据。

## 必须知道的限制

**这些限制写进了产物，不是免责声明，是结论的适用范围。**

| 分析 | 限制 |
|---|---|
| 空间域 | spot 直径 55 μm 含 1-10 个细胞 —— **不能说"某个域是某一种细胞"**，只能说这个区域的细胞组成不同。域边界有 ±1 spot 不确定性 |
| 域标签 | `z_margin <= 0.5` 时标签不可信（实测 7/13 个域如此） |
| SVG | Moran's I 依赖空间权重矩阵；置换次数限制了最小 p 值 |
| 组成 | `builtin` 模式是 **marker 打分法，不是解卷积**（`is_deconvolution: false`）；**不给出细胞比例**，只有相对空间趋势 |
| 邻域 | 用 argmax 硬分配细胞类型，丢失 spot 内混合信息；z-score 没做多重检验校正 |
| 通讯 | **共表达 ≠ 通讯**；不做"通讯与否"的判定 |

## 目录结构

```
scripts/         00_fetch … 07_spatial_communication, main_analysis
scripts/lib/     common.py（spatial_xy 等）
assets/          datasets.yml, config.*.yml, reference_signatures.yml,
                 ligand_receptor.yml
tools/           check_py_syntax.mjs, check_figures.mjs,
                 check_artifact_paths.mjs, verify_spatial_alignment.py
references/      methods.md, troubleshooting.md
```

## 工程规则

见 `AGENTS.md`。14 条规则全部来自实测踩过的坑，
每条都写了"为什么"和"错误长什么样"。
