# 线3 / 单细胞空间转录组 + 虚拟敲除线

> 术语见 `GLOSSARY.md`，隔离约束见 `docs/isolation.md`。
> 本目录只放线3 代码。虚拟敲除是本线**末端模块**，不是独立一条线。

## 输入

下载得到的单细胞空间转录组数据（10x Visium 等）。

## 前置步骤

1. 下载
2. QC
3. 标准化
4. 降维
5. 空间域识别 / 细胞分割注释

## 下游步骤（九项）

1. 空间域识别
2. 空间域 marker / 空间差异表达
3. 空间自相关 / 空间表达模式：Moran's I、Geary's C、SpatialFeaturePlot
4. 空间邻域 / 共现分析
5. 区域组成差异分析
6. 空间通讯分析
7. 解卷积 / 细胞类型注释
8. 目标基因空间定位
9. 虚拟敲除（本线末端模块）

## 工具分工（线内混用）

| 工具 | 语言 | 用途 |
|---|---|---|
| Squidpy | Python | 空间自相关 / 邻域 / 共现 |
| Seurat / Semla | R | 空间域与可视化 |
| BayesSpace | R | 空间域识别 |
| RCTD / cell2location | R / Python | 解卷积 |
| NICHES | Python | 空间通讯 |

虚拟敲除工具同线2 候选（scTenifoldKnk / Geneformer / CellOracle），在本线内部跑。

## 解卷积参考来源（硬约束）

空间线做解卷积时，参考**不得取自线2**。
只允许两个来源：

- `public_references/` 下的公共参考
- 线3 自有单细胞分辨率数据

## 虚拟敲除约束

- 优先单细胞分辨率空间平台。
- 若是 spot 级空间数据，直接虚拟敲除需谨慎：可按空间域聚合或按子集运行。
- 输出 DR genes（**DR 基因不是 DEG**）。
- 敲除后只做本线内：DR 基因富集、映射回空间坐标看分布、与线内区域 DEG / 空间域 marker 交集。

## 产出去向

算力节点 `<BIO_DATA_ROOT>/projects/<项目>/line3_spatial/{raw_data,processed,results,figures}/`
（脚本 `bin/line3_spatial/` 与环境声明 `envs/line3_spatial/` 全局一份，所有项目共用）。
台账 `line_id` = `line3_spatial`。

## 本目录结构（本地代码本体，与算力节点侧数据目录一一对应）

本地只放代码；数据、产物、环境全部在算力节点 `<BIO_DATA_ROOT>/projects/<项目>/line3_spatial/`（代码与环境在顶层全局一份）。

```
src/line3_spatial/
├── README.md            # 本文件
├── scripts/             # 线内分析脚本（.py / .R 可混用，但只在本线内）
└── checkpoint/          # 桥上检查点脚本，命名 checkpoint_*.py / checkpoint_*.R
```

桥接检查点（Python 出 `.h5ad` 给 R 读之类）必须校验四项：基因名、矩阵方向、细胞名、embedding。
约束与校验规则见 `docs/isolation.md`。

## 状态

**步骤 1 `ingest_spatial` 已完成**（2026-10-13，界面验收 exit 0，台账 completed/passed）。

数据源 GSE192741（Guilliams et al. 2022 Cell，肝巨噬细胞空间蛋白基因组研究）：15 样本 =
10 小鼠（GPL24247 / mm10）+ 5 人（GPL24676 / hg19），组覆盖 StSt/SD/WD/Healthy/Steatotic
五类（grouping.txt 5 个声明组合全覆盖；其"共 10 个样本"与实际 15 个不符，如实记 warning
不阻断，按实际 15 个全量读入）。

生产脚本 `src/line3_spatial/scripts/ingest_spatial`（bash 包装 + python3 heredoc）：
输入按文件名特征四分流（RAW tar / GPL24247 soft / GPL24676 soft / grouping），
CSC→CSR 转置（索引压 int32），var_names 用 Ensembl ID（symbol 有重复：小鼠 35 / 人 91），
obs_names 加样本前缀防 barcode 跨样本碰撞。自校验内建：tar 成员集完整无意外、
soft 元数据齐全、平台-物种-基因组三方一致、矩阵 barcode 全落 positions in_tissue=1、
矩阵 spot 数 == soft "number of spots"、基因 ID 唯一且同物种跨样本同序、计数非负、
组映射完整、ADT barcode 可对齐（BH001 93 抗体，引号 + 无 -1 后缀，重合 1698/1708）。
产物 4 个：`spatial_counts_mouse.h5ad`（17290 spots × 31053 基因，10 样本）、
`spatial_counts_human.h5ad`（6889 spots × 32738 基因，5 样本）、
`samples_meta.tsv`（15 行 29 列）、`ingest_report.txt`。

坐标口径（下游步骤必须遵守）：`obsm["spatial"]` 为 pxl 像素坐标，**仅限样本内组织图叠加**；
跨样本共同坐标框是 `obs["array_row"]/["array_col"]`（scalefactors 两档扫描设置不同：
spot_diameter 21.39 vs 147.16，像素坐标跨样本不可比）。`uns["spatial"][sample_id]`
只存 lowres 图 + 4 个 scalefactors；hires 未入 h5ad（12MB × 15，需要时从 tar 重取）。
BH001 的 ADT 只记录台账，不并入基因表达矩阵（模态不同）。

**步骤 2 `qc_spots` 已完成**（2026-10-13，界面验收 exit 0，台账 completed/passed）。

生产脚本 `src/line3_spatial/scripts/qc_spots`：输入契约 17 项（obs 16 列 + var 3 列 +
uns 11 键 + obsm/uns spatial 双在位 + CSR 格式 + 物种唯一），指标重算
（ncount/nfeature/pct_mito）为判定唯一口径并与 ingest 落盘列对账（口径漂移即失败）。

阈值口径（2026-10-13 算力节点三轮探针实测定稿，固定不随结果调整）：

- `ncount >= 100`（spot 总计数，低于此视为空点/残点）
- `nfeature >= 50`（spot 检出基因数）
- `pct_mito <= 50`（线粒体计数占比，仅滤病态 spot）
- 基因表达 spot 数 >= 3（噪声基因）

定稿依据（探针实测，GSE192741 17290 小鼠 + 6889 人源 spot）：

- 沉积数据已由 Space Ranger 预过滤（仅组织下 spot 进矩阵），低指标 spot 极少：
  ncount<100 或 nfeature<50 小鼠 7 个 / 人源 10 个。
- **mito 阈值定 50 而非 20 的关键证据**：mito>20% 的 spot 小鼠 42.9% / 人源仅 6.6%，
  且人源集中在 JBO014 单样本（掉 19.63%）、不聚在捕获区边缘（边缘占比与全体持平），
  是样本级质量特征而非边缘损伤；滤了会在组织上打洞破坏邻域分析，故阈值定 50
  只滤病态（小鼠 1 个 / 人源 2 个），样本级偏高记 warning 不过滤（BH001 已正确触发）。
- 基因 <3 spot：小鼠 2 个 / 人源 0 个；被滤基因 UMI 占比仅 0.001-0.004%。
  34% 基因（小鼠 10500/31053）在任何 spot 零表达，属参考基因组死基因，正常。

自校验内建：每样本保留比例 >= 75%（阈值设计失误早炸而非静默毁数据）、保留矩阵
无负值、obs 列不丢、obs_names 唯一、uns['samples'] 与 obs 一致、uns['spatial']
样本集与 obs 一致、保留矩阵指标与阈值一致（过滤生效证明）。任一不过非零退出。

产物 4 个：`spatial_qc_mouse.h5ad`（17283 spots × 18628 基因）、
`spatial_qc_human.h5ad`（6877 spots × 19613 基因）、`spot_qc_metrics.tsv`
（全量 24179 spot，含被滤带去留原因，可审计）、`qc_report.txt`。obs 追加
`total_counts`/`n_genes_by_counts`/`pct_mito` 三列，`uns['qc']` 记口径键；
`uns['spatial']` 图与 scalefactors 随矩阵保留（只丢 spot 不丢样本，图不受影响）。

环境：Python 3.12.3 / scanpy 1.12.4 / anndata 0.13.4 / pandas 2.3.3（manifest 已由误抄的
R 修正，坑 #106 同模式）。写 h5ad 显式 `ad.settings.allow_write_nullable_strings = False`（坑 #107）。

**步骤 3 `normalize_hvg` 已完成**（2026-10-13，界面验收 exit 0，台账 completed/passed）。

生产脚本 `src/line3_spatial/scripts/normalize_hvg`：输入契约 32 项（obs 20 列 + var 3 列 +
uns 12 键 + obsm/uns spatial 双在位 + CSR 格式 + 物种唯一），**逐物种处理**（小鼠 /
人源基因宇宙不同，禁止跨物种选取 HVG）。

归一化口径：`sc.pp.normalize_total(target_sum=1e4) + log1p`（Scanpy 标准口径，与 Seurat
LogNormalize(type='count') 数值等价），原始整数计数存 `layers['counts']`；mito 基因为
**内源基因不剔除**（不同于 line2 剔除 ERCC spike-in），入选 HVG 数如实报告
（小鼠 10 / 人源 5）。

HVG 定稿（算力节点探针 + 生产双跑一致）：seurat flavor（dispersion 口径），
n_top_genes=2000，按物种全数据集选取；小鼠 top20 为免疫球蛋白为主（Igkc/Igha/Ighm 等，
与肝巨噬细胞研究主题一致），人源 top20 为 C7/FBLN1/PTGDS 等肝组织高表达基因，生物学合理。

自校验内建（任一不过非零退出）：输入契约 32 项、归一化后每 spot 库大小精确 1e4
（rtol 1e-4 / atol 1e-2）、log1p 后 X 全有限非负、归一化不改变稀疏结构（nnz 前后一致）、
HVG 数精确 2000、**layers['counts'] spot 合计与 obs['total_counts'] 逐 spot 对账**、
obs 20 列不丢。回读验证 147 项全过（含 tsv 与 h5ad HVG 集合按 index 对账、
disp_norm 降序、per-sample spots 数）。

产物 5 个：`spatial_norm_mouse.h5ad`（17283 spots × 18628 基因，layers['counts'] 保留）、
`spatial_norm_human.h5ad`（6877 spots × 19613 基因）、`hvg_genes_mouse.tsv`、
`hvg_genes_human.tsv`（各 2000 行，means/dispersions/dispersions_norm）、`norm_report.txt`。
uns 追加 `normalization`/`hvg`/`seed`/`source_h5ad`。

**步骤 4 `dim_reduce` 已完成**（2026-10-13，界面验收 exit 0，台账 completed/passed）。

生产脚本 `src/line3_spatial/scripts/dim_reduce`：逐物种在 HVG 2000 子集上
scale(max_value=10) → PCA(50) → neighbors(15,30PC) → UMAP，X_pca/X_umap 与
neighbors 图贴回全基因对象（X 保持归一化 log 值不动，layers['counts'] 保留）。

口径决策（算力节点探针 + 生产双跑一致）：

- **不做跨样本批次校正**：下游空间域识别逐样本做（像素坐标跨样本不可比，
  见坐标口径），嵌入仅作质控与无监督聚类参考；样本/组分离度如实入报告。
- 实测分离度（PCA10 silhouette，采样 5000 spots）：mouse sample_id 0.1105 /
  group 0.3894，human sample_id 0.2560 / group 0.2360。样本不主导方差
  （mouse 0.11 低），跨样本嵌入可比性可接受。
- 方差解释：mouse PC1 6.421% / cum30 16.822%，human PC1 2.376% / cum30 10.302%
  （空间数据表达异质性高、无单一主导轴，属正常）。

自校验内建（任一不过非零退出）：输入契约（步骤 3 产物 32 项 + 未做过降维）、
HVG 数精确 2000、X_pca 形状 (n,50) 与全有限、X_umap 全有限、方差解释率非空、
layers['counts'] nnz 不变、obsm spatial / uns spatial 不丢、obs 20 列不丢。
回读验证 54 项全过（含 connectivities 每 spot 有邻居、counts 层与
obs['total_counts'] 对账、per-sample spots 与步骤 3 一致）。

产物 3 个：`spatial_dimred_mouse.h5ad`（17283 spots × 18628 基因，
obsm X_pca/X_umap + obsp 邻居图）、`spatial_dimred_human.h5ad`
（6877 spots × 19613 基因）、`dimred_report.txt`（方差解释 + silhouette
诊断 + 计时 + 口径说明）。uns 追加 `dimred` 口径键。

**步骤 5 `spatial_domains` 已完成**（2026-10-13，界面验收 exit 0，台账 completed/passed）。

生产脚本 `src/line3_spatial/scripts/spatial_domains`：**逐样本**在自身 X_pca 上重算
neighbors(15,30PC) → Leiden(resolution=0.5, igraph, n_iterations=2)，域编号按样本
命名空间隔离（`obs['domain'] = "<sample_id>_d<n>"`）。

工具路线决策（dry-run 证据）：README 首选 BayesSpace（R）运行环境未装；squidpy
dry-run 会连带**降级 numpy 2.5.3→2.4.6** 并拖入 dask/distributed/geopandas/pyarrow
等 60+ 包（坑 #107 场景），故走 scanpy 原生零新依赖方案，决策证据落报告
`tool_route` 行。

resolution 定稿依据（算力节点探针双跑）：0.5 vs 1.0 对比，逐样本空间连贯性
（spot 的样本内像素最近邻与自身同域比例）mouse mean 0.5163 vs 0.3627、
human 0.6951 vs 0.4647；碎片域（<30 spots）mouse 4 个 vs 11 个。取 0.5：
域更贴合空间结构、碎片更少。

关键口径（坑 #114）：**多样本像素坐标不是全局坐标系**（各样本范围重叠，
全数据集 KD 树 91.5% 最近邻落在他样本），空间连贯性必须逐样本算；
后续 Moran's I、邻域共现等像素空间指标同此口径。

结果：mouse 57 域（10 样本，每样本 5-7 域，连贯性 mean 0.5163 / min 0.4143）；
human 18 域（5 样本，每样本 3-5 域，连贯性 mean 0.6951 / min 0.5349）。

自校验内建：输入契约（步骤 4 产物 32 项 + 未做过域识别）、每 spot 均有域标签、
域命名空间隔离（前缀须等于本 spot 的 sample_id）、每样本域数 > 1、连贯性在 0-1、
X_pca/X_umap 不被改动、counts 层 nnz 不变、spatial 双不丢、obs 20 列不丢。
回读验证 43 项全过（含 X_pca 与步骤 4 逐元素一致、per-sample spots 不变、
uns 与 obs 对账）。

产物 3 个：`spatial_domains_mouse.h5ad`、`spatial_domains_human.h5ad`
（obs 加 `domain` 列，uns 加 `domains` 口径键含逐样本域数与连贯性）、
`domains_report.txt`（域数/大小/连贯性/resolution 敏感性/工具路线）。

**步骤 6 `spatial_markers` 已完成**（2026-10-13，界面验收 exit 0，台账 completed/passed）。

生产脚本 `src/line3_spatial/scripts/spatial_markers`：逐物种对 domain 分组跑
wilcoxon 秩和检验（scanpy 原生，pts=True 同时取 pct_nz_group/pct_nz_reference
两列检出率），n_genes=200，key_added=domain_markers。

关键口径（算力节点三轮探针实测定稿）：

- **marker names 是 Ensembl ID**（ingest 口径 var_names=ID、symbol 在 var 列），
  肝巨噬 panel 用 symbol 直配命中为空即因此；tsv 导出前经
  `var['gene_symbol'].to_dict()` 映射出 symbol 列，0 缺失（symbol 有重复属已知口径，
  按 ID 逐行取）。
- 生物学合理性（mouse JBO009_d2 巨噬域）：Cd5l rank7 / C1qb rank17 / C1qa rank35 /
  Clec4f rank185 / Tmsb4x rank1 / Cd74 rank5 / Lyz2 rank34 同域；human MARCO
  （JBO022_d2）/ C1QB / C1QA 弱命中方向一致（人肝 KC 标记 spot 级混合细胞信号弱，
  属正常）。每域 top1 呈现肝 zonation：mouse Mup3/Elovl3/Cyp2c69/Isg15/Hamp/Apom，
  human CYP2A6/CES1/ALB/APOA1，与 Guilliams 2022 肝分区发现一致。
- 计时：wilcoxon pts=True mouse 29.9s / human 8.2s（pts=False 24.9s/8.1s，
  成本可接受保留 pts）。

自校验内建（任一不过非零退出）：输入契约（步骤 5 产物键集合 obs 21 + var 7 +
uns 21 + domain 列 + 未做过 marker）、marker 表行数 = 域数 x 200（mouse 11400 /
human 3600）、top20 行数 = 域数 x 20（1140 / 360）、symbol 映射 0 缺失、泛免疫基因
sanity（mouse Tmsb4x/Cd5l/C1qb 与 human CD74 必须在 top200，探针实测组内排名
1/7/17/20 远离边界）、X nnz 不变、counts 层 nnz 不变、X_pca 不被改动、domain 列
不被改动、spatial 双不丢、obs 列不丢。回读验证 43 项全过（含 tsv 与
uns['domain_markers'] 按组 top20 对账、每组 scores 降序、uns['markers'] 参数键）。

产物 5 个：`spatial_marker_mouse.h5ad`（17283 spots × 18628 基因，uns 加
`domain_markers` top200 全量排名 + `markers` 口径键）、`spatial_marker_human.h5ad`
（6877 spots × 19613 基因）、`domain_markers_mouse.tsv`（1140 行）、
`domain_markers_human.tsv`（360 行，列：group, rank, ensembl_id, gene_symbol,
scores, logfoldchanges, pvals, pvals_adj, pct_nz_group, pct_nz_reference）、
`marker_report.txt`（域数/marker 行数/panel 命中明细/每域 top1/计时）。

**步骤 7 `spatial_autocorr` 已完成**（2026-10-13，界面验收 exit 0，台账 f5d8bea8
completed/passed）。

生产脚本 `src/line3_spatial/scripts/spatial_autocorr`：逐物种、逐样本在像素坐标上建
KNN(k=6) 空间权重（cKDTree → 对称化 w+w.T → 行标准化），样本内逐基因算 Moran's I
与 Geary's C（批量稀疏乘全基因一次算完），基因级取跨样本 mean/min/max 与
n_valid_samples。矩阵取 X（归一化 log 值，squidpy/spatial_stats 惯例口径）。

口径定稿依据（算力节点两轮探针实测）：

- **逐样本**（坑 #114）：跨样本像素坐标范围重叠，禁在全数据集建树。
- 计时：批量路线 mouse 全基因 11.5s / human 4.9s（无需子集）。
- 对称化前后 top10 排序一致，定稿对称化对齐 squidpy/scanpy symmetrize 惯例。
- Moran's I 分布合理：mouse mean 0.0164 / max 0.5949；human max 0.7374。
- 生物学 sanity：mouse Mup17(0.595)/Cyp2e1(0.590)/Alb(0.590)/Mup20/Fabp1 等肝
  zonation 基因 I 0.56-0.60 z 57-62；human MT-ND4(0.734)/MT-CO1/ALB(0.698) 等
  I 0.64-0.74 z 58-65。known genes：Apoa2 0.521 / Mup3 0.557 / Hamp 0.449
  （zonation 高）vs Cd5l 0.157 / C1qb 0.120 / MARCO 0.101 / C1QB 0.092
  （免疫弥散低），方向与生物学一致；Geary's C 同向互证（zonation 0.43-0.55
  vs 免疫 0.84-0.91，均 <1 即全基因组整体聚集、zonation 聚集更强）。

推断口径：**Moran's I 用 Cliff-Ord 随机化假设的正态近似解析 z**（esda norm 模式
同式，E[I]=-1/(n-1)，Var[I] 含 S0/S1/S2 矩，零置换成本）；z 逐样本算（E[I]/sd
依赖样本 n）再跨样本平均。**Geary's C 只报统计量不做解析推断**（其随机化方差依赖
基因四阶矩，跨基因解析式与 esda permutation 模式不等价），方向解读以 Moran's I
为准，报告已注明。nan 基因（全 spot 零表达 denom=0）实测 0 个（步骤 2 QC 已滤
min_cells>=3）；部分样本零表达基因 mouse 4236/18628 / human 3451/19613 如实记
n_valid_samples 与 partial_zero_genes，不填 0 不删行。

自校验内建（任一不过非零退出）：输入契约（步骤 6 产物键集合 obs 21 + var 7 +
uns 21 含 markers + 未做过自相关）、每样本权重行和=1、I 理论范围 [-1,1]、
C 理论范围 [0,2]、zonation 基因 I 全部高于免疫弥散基因（两物种分别硬断言）、
X nnz 不变、counts 层 nnz 不变、X_pca 不被改动、domain 列不被改动、obs/var
列不丢。回读验证 38 项全过（含朴素双重循环教科书式核对单样本单基因 I 落在
min..max 内、tsv 按 i_mean 降序、不变量与上游产物逐项一致）。

产物 5 个：`spatial_autocorr_mouse.h5ad`（17283 spots × 18628 基因，var 加
morans_i_mean/min/max、morans_z_mean、gearys_c_mean/min/max、
autocorr_n_valid_samples 八列，uns 加 `autocorr` 口径键）、
`spatial_autocorr_human.h5ad`（6877 spots × 19613 基因）、
`gene_autocorr_mouse.tsv` / `gene_autocorr_human.tsv`（全基因按 morans_i_mean
降序，列：var_index, ensembl_id, gene_symbol, 八列统计 + n_valid）、
`autocorr_report.txt`（口径/sanity 明细/top10/partial_zero 统计/计时）。

**步骤 8 `spatial_cooccurrence` 已完成**（2026-10-13，界面验收 exit 0，台账 94eadd75
completed/passed）。

生产脚本 `src/line3_spatial/scripts/spatial_cooccurrence`：逐物种、逐样本在像素坐标上建
KNN(k=6) 空间权重（**与步骤 7 同一权重定义**：有向邻接 A → 对称化 A+A.T（互向近邻
权重为 2）→ 行标准化，跨步骤可比），对 domain 标签算共现矩阵
C[di,dj] = onehot_di^T W onehot_dj（one-hot 稀疏乘一次算全类型对），
100 次标签置换（seed 注入，保持几何权重只打乱标签）得每 domain 同域富集 z
（squidpy gr.nhood_enrichment 同思路，零新依赖自实现）。

口径定稿依据（算力节点探针实测）：

- 计时：one-hot 稀疏乘 + 100 次置换 mouse 全 10 样本 0.6s / human 0.3s。
- 共现矩阵对角占优：所有样本所有 domain 同域邻域富集 z 全为正（min 2.89 /
  mean 14-48），与步骤 5 Leiden 空间域的空间连贯性自洽。
- same_frac 0.09-0.90：JBO006 有碎片域 0.089（与步骤 5 碎片域记录一致），
  JBO014 等连贯样本 0.64-0.90。

自校验内建（任一不过非零退出）：输入契约（步骤 7 产物键集合 obs 22 + var 15 +
uns 24 含 autocorr + 未做过共现）、每样本类型数 = 样本内 domain 数、共现行和 > 0、
same_frac 在 [0,1]、同域富集 z 全部为正（空间域聚集性，探针实测远离 0 稳定可断）、
X nnz 不变、counts 层 nnz 不变、X_pca 不被改动、domain 列不被改动、obs/var 列不丢。
回读验证 47 项全过（含朴素双重循环教科书式核对单样本共现矩阵对角，与生产稀疏实现
diff 4.74e-05 纯 4 位小数舍入；spot 级两列与 tsv 逐行对账 75 行全过；
top_partner 与 uns 归一化矩阵逐项一致）。

产物 5 个：`spatial_cooccurrence_mouse.h5ad`（17283 spots × 18628 基因，obs 加
`domain_same_frac`/`domain_same_z` 两列 spot 级贴回，uns 加 `cooccurrence` 口径键
含逐样本归一化共现矩阵与 z 矩阵）、`spatial_cooccurrence_human.h5ad`
（6877 spots × 19613 基因）、`domain_cooccurrence_mouse.tsv`（57 行）/
`domain_cooccurrence_human.tsv`（18 行，列：sample_id, domain, n_spots,
same_frac, same_z, top_partner_domain, top_partner_frac）、
`cooccurrence_report.txt`（口径/两物种汇总/top partner 示例/计时）。

**步骤 9 `spatial_composition` 已完成**（2026-10-13，界面验收 exit 0，台账 9272b61a
completed/passed；界面首跑台账 22e09761 failed 是 fill 工具未把值写进 React state
（坑 #118 复发，原生 setter 重填后通过））。

生产脚本 `src/line3_spatial/scripts/spatial_composition`：逐物种逐样本算**域组成的
rank 对齐规模谱**（每样本内各域 spot 占比降序序列，ties 按域名排序保证确定性）
+ 样本级多样性（Shannon H / 有效域数 exp(H) / top1 / top3）+ 组间检验
（condition_class 与 condition 两级）。

口径定稿依据（算力节点探针实测，`cn_debug/s4_composition_probe.py`）：

- **domain 身份样本私有不可跨样本直比**：每样本的 domain 编号来自该样本自己的
  Leiden 聚类（`<sample_id>_d<n>`），逐域 spot size 在全样本上全 0（域间大小无
  可比性），Bray-Curtis 在域身份空间全 1.0（任意两样本"完全不同"）。按域身份比
  组成无意义，改为 rank 对齐的域规模谱 + 样本级多样性，这是本步骤最重要的口径
  选择，下游解读不得绕过。
- 小样本统计口径（与 line2 DEG 教训同款）：显著性只作参考，重点报效应量与置换 p。
  精确置换枚举：2 组且 C(n,n1) <= 20000 时全枚举（mouse C(10,4)=210 /
  human C(5,2)=10），无蒙特卡洛误差、结果与 seed 无关；非 2 组回退 seed 注入
  10000 次随机置换（add-one 校正）。human 2v3 仅 C(5,2)=10 种分配，p 最小可能值
  0.1，永远不可能 p<0.05（数据设计属性，已写入报告）。
- 教科书 KW 独立实现：H = 12/(n(n+1))·Σ(Rj²/nj) − 3(n+1) + tie 校正，
  p = chi2.sf(H, k−1)，不调 scipy.stats.kruskal（防同源互证，坑 #119 教训）。

实测结果：mouse condition_class（nafld 4 vs stst 6）逐 rank 0/7 显著、
置换 p=0.1095（exact，210 分配）；mouse condition 级 SD 单样本如实跳过；
human 两级（steatotic 3 vs non_steatotic 2 / Steatotic 3 vs Healthy 2）
逐 rank 0/5、置换 p=0.2000（exact，10 分配）。两组均无显著组成差异，
如实报告（小样本设计属性）。

自校验内建（任一不过非零退出）：输入契约（步骤 8 键集合：obs 22 列含
domain_same_frac/domain_same_z + var 15 列 + uns 24 键含 cooccurrence +
未做过 composition）、每样本未舍入组成谱行和=1（1e-9 紧容差）、Shannon
在 [0, log K]、与步骤 8 uns 样本记录对账（domain 数/n_spots）、rank 1..K
完整、落盘舍入行和 1e-3、domain_prop 在 [0,1]、spot 级 rank 与所属域一致、
X nnz / counts 层 nnz / X_pca / domain 列 / obs-var 列不变量。测试跑连修三坑：
categorical `value_counts()` 把未观测类别以 count 0 列出（坑 #120，修为
`astype(str)` 后统计，自校验正确拦截无错误产物落盘）；舍入行和容差（步骤 8
同课复发，紧校验放未舍入值上）；uns 禁 list-of-dicts（per_rank 改 dict 键控
"rank1".."rankK"，None 数值改 nan / 跳过原因空串）。回读验证 308 项全过
（`cn_debug/s4_composition_readback.py`）：所有参考值从步骤 8 产物独立重算
（朴素 Counter 组成谱 / 教科书 KW 自带 tie 校正不调 scipy / 朴素双重循环 L1 +
itertools 精确枚举置换），A 键集合 / B tsv 结构 / C 朴素 vs tsv 逐样本 /
D spot 级两列 / E uns samples / F 组间检验 / G 不变量 / H 报告 token。

产物 5 个：`spatial_composition_mouse.h5ad`（17283 spots × 18628 基因，obs 加
`domain_rank`/`domain_prop` 两列 spot 级贴回，uns 加 `composition` 口径键含
per_rank / sample_metrics / permutation / samples）、`spatial_composition_human.h5ad`
（6877 spots × 19613 基因）、`sample_composition_mouse.tsv`（10 行）/
`sample_composition_human.tsv`（5 行）、`composition_report.txt`（口径/两物种
汇总/组间检验/小样本设计属性说明）。

**步骤 10 `spatial_communication` 已完成**（2026-10-13，界面验收 exit 0，台账 c0005d86
completed/passed）。

生产脚本 `src/line3_spatial/scripts/spatial_communication`：LR 资源取 **liana consensus**
（本地已装 1.10.0，离线零新依赖；资源盘点 4620 对，cellphonedb 1223 / cellchatdb 1912 /
icellnet 738 可作备选；squidpy / CellPhoneDB / NICHES / COMMOT / SpatialDM 均未装）。
配体/受体基因用**大写归一匹配**：liana 资源是人类命名（大写），mouse symbol 是 Title
case（Tgfb1/Csf1），直接匹配 mouse 只命中 2/4620，大写归一后 2213/4620；human 2418/4620。
重复 symbol 取首个 var 索引（与步骤 6 marker 映射口径一致）。

空间打分口径（探针 1-5 定稿）：逐物种、逐样本像素 KNN(k=6) 邻域权重（与步骤 7/8/9
同一 build_knn_weights），对每对 LR 算 domain 对级强度
score(p; di, dj) = sum_i sum_j lv_ip * w_ij * rv_jp * 1[dom_i=di] * 1[dom_j=dj]，
再除以 domain 对邻域质量 mass(di,dj)，得**单位质量 LR 强度**，消除域规模差异
（与步骤 9 rank 对齐精神一致）。向量化实现（探针定稿 human 1.5s / mouse 6-8s；
逐 pair 循环 253s/588s 不可接受，提速 100 倍）：
rec_side = W @ (rv * mask_dj)，score[:, di, dj] += einsum('ip,ip->p', lv*mask_di, rec_side)。

统计口径同步骤 9：样本级指标（lr_top_score / lr_mean_score / lr_n_active）KW + eta2
+ L1 置换检验（condition_class 与 condition 两级；2 组且 C(n,n1)<=20000 精确枚举）。
spot 级贴回 lr_score_recv（每 spot 作为受体端收到的最大单位质量 LR 强度，
max over pairs × 源域；**先收集 max 再一次写入**，双重循环直接赋值会被后写覆盖，
回读对账抓到此点并已修）。

实测结果：mouse condition_class（nafld 4 vs stst 6）置换 p=0.0524（exact 210 分配）；
mouse condition 级 SD 单样本如实跳过；human 两级置换 p=0.6000（exact 10 分配）。
两组均无显著通讯差异，如实报告。生物学合理性：mouse top 跨域对 APOA2-LRP1
（载脂蛋白代谢，肝主题）、human top HP-ASGR1/HP-ASGR2（肝细胞去唾液酸糖蛋白受体，
Guilliams 2022 肝巨噬细胞主题正中）。

自校验内建：输入契约（步骤 9 键集合 obs 24 含 domain_rank/domain_prop + var 15 +
uns 25 含 composition + 未做过通讯）+ 权重行和=1 + score 非负 + 样本内观测域
对角质量>0（doms 是全物种并集，坑 #120 同源）+ 与步骤 9 uns 记录对账（n_spots/
样本内 domain 数）+ spot 级非负有限 + X/counts/X_pca/domain/obs/var 不变量 +
参考值对账（单样本 score 最大组合朴素双重循环，两物种六位小数精确一致）。
tsv 增落未舍入 raw 列（lr_top_score_raw / lr_mean_score_raw）：tsv 的 4 位小数
舍入会引入额外 ties，使教科书 tie 校正后的 KW H 与生产不一致，属舍入伪差异非口径差异。

产物 5 个：`spatial_communication_mouse.h5ad`（17283 spots × 18628 基因，obs 加
lr_score_recv，uns 加 `communication` 口径键含 method/n_lr_pairs_*/domains/
domain_pair_mass/top_pairs/group_tests）、`spatial_communication_human.h5ad`
（6877 × 19613）、`sample_communication_mouse.tsv`（10 行，12 列）/
`sample_communication_human.tsv`（5 行）、`communication_report.txt`。

**步骤 12 `spatial_localization`（目标基因空间定位）已完成**（2026-10-13，界面验收
exit 0，台账 38b675f2 completed/passed）。

生产脚本 `src/line3_spatial/scripts/spatial_localization`：目标基因候选池 = 步骤 6
rank_genes_groups 的**各域 marker 前 TOPK_MARKER(20) 名去重**（探针定性：wilcoxon
top1 在 spot 噪声数据上含大量泛表达基因 frac>=0.99 无定位信息量，取 top20 让打分
排序自行筛选）。逐基因逐样本算域富集（domain 身份样本私有，跨样本累加无意义，
同步骤 9 教训）：enrich(sample) = 域内表达 spot 占比 / 该基因全 spot 占比；
self_consistent(sample) = 该基因 marker 域（样本内）恰为 enrich 最高域。
目标基因排序 = self_consistent 样本数降序 → best_enrich 降序 → frac_all 升序
（泛表达基因 frac 高自然沉底，不用硬阈值拍脑袋）。

实测结果：mouse 候选 363 基因，top 目标 Hba-a1（自洽 7/10，enrich 2.30）/
Hbb-bt（6/10，4.04）/Gsn/Spp1（5/10；Spp1 与步骤 6 巨噬域 marker 一致）/
Cyp2c69/Cyp2c40（5/10，肝 zonation 门脉周基因）；human 候选 190 基因，top
MTRNR2L1（3/5）/MT-ND4L/MT-ND4（线粒体）/PTGDS（2/5，步骤 3 已知人源 HVG top）/
DGAT2（NAFLD 脂代谢主题）。spot 级贴回 target_gene_expr（top 目标基因表达值）。

自校验内建：输入契约（步骤 10 键集合 obs 25 含 lr_score_recv + var 15 + uns 26
含 communication + 未做过定位）+ 候选池非空 + 候选基因全映射 var_names + enrich
非负 + 自洽样本数 <= 总样本数 + spot 级非负有限 + X/counts/X_pca/domain/obs/var
不变量。回读验证 62 项全过（`cn_debug/s4_loc_readback.py`）：候选池朴素重算、
前 20 基因 enrich/自洽/计数逐点一致、tsv 排序键独立验证、targets 与 tsv 前 N
一致、spot 级与 top 基因列逐元素一致、不变量。回读首跑 4 FAIL 是回读自身把 uns
的 target 键按字符串字典序排（target1/target10/target2...，target10 插到第二位），
修为按数值后缀排序后全过（产物本身正确）。

产物 5 个：`spatial_localization_mouse.h5ad`（17283 spots × 18628 基因，obs 加
target_gene_expr，uns 加 `localization` 口径键含 method/topk_marker/n_candidates/
n_targets/top_target_gene/targets）、`spatial_localization_human.h5ad`
（6877 × 19613）、`gene_localization_mouse.tsv`（363 行）/ gene_localization_human.tsv
（190 行，列：ensembl_id/gene_symbol/n_expr/frac_all/n_samples_top/n_samples_total/
best_enrich/best_example/n_marker_domains）、`localization_report.txt`。

**步骤 13 `spatial_knockout`（虚拟敲除）已完成**（2026-10-13，界面验收 exit 0，
台账 1414fbdf completed/passed）。

生产脚本 `src/line3_spatial/scripts/spatial_knockout`：敲除对象 = **逐样本 spot 级
子集运行**（aims.txt：spot 级数据可按空间域聚合或子集运行）。域聚合路线
（样本 x 域 pseudo-bulk，mouse 57 / human 18 个）实测废弃：human 聚合后
scTenifold QC 只剩 1 个 pseudo-cell，`pc_net` 抛 `n_comp must be < number of
non-constant genes (0)`；mouse 侧虽能跑但两物种路线不对称。基因宇宙 = HVG 2000
+ 敲除基因（与 line2 步骤 9 同口径；全基因 18628 x ~1700 spot 建共表达网络
实测不可行）。敲除基因 = Spp1/SPP1（线内证据：步骤 12 target5 + 步骤 6 巨噬
域 marker 同域；不取 line2 TNF，三线不交叉）。**scTenifold mito QC 失真
（坑 #122）**：其 mito 比率按喂入矩阵库计算，HVG 子矩阵分母缩水使 mouse
比率放大约 5 倍，修为全 counts 口径预剔 mito 比 >= 0.1 的 spot（阈值与匹配
规则同 scTenifold）+ `qc_kws={"max_mito_ratio": 1.0}` 关库内过滤。

物种级显著 DR = 跨样本并集（adj p<0.05 且非敲除基因；DR 不是 DEG，计算模拟
非实验验证）。排序 = 显著样本数降序再最小 adj_p 升序。KO 表达占比 <5% 或
mito 预剔后 <50 spot 的样本跳过并如实记录。spot 级贴回 `obs['ko_dr_expr']`
（显著 DR 基因 counts 均值），供映射回空间坐标看分布。

实测结果：mouse 9/10 样本跑通 72.8s，54 个物种级显著 DR（top：Mup11 6 样本
adj_p=9.8e-44 / Rnase4 6 样本 / mt-Co1 5 样本 / Mup1/Rgn；Mup/Rnase4/Rgn 是
肝 Kupffer 细胞表型）；human 4/5 样本跑通 25.2s，14 个（top：MT-ATP6 3 样本 /
SERPINA1 2 样本 / FGL1/CYP2A6/NNMT/SDS 各 1 样本，肝细胞功能基因）。
跳过样本：mouse CAP002（KO%=0.0157）/ human JBO015（KO%=0.0351）。
注意：mito 基因大量进 top 是 scTenifold 保留 mito 于网络的默认行为（不同于
scRNA 流程剔除 ERCC），已写报告说明。

自校验内建：输入契约（步骤 12 键集合 obs 26 + var 15 + uns 27 + 未做过敲除）+
KO 基因 counts 有表达 + 逐样本行数/DR 表结构七列 + 敲除基因不进显著 DR +
ko_dr_expr 有限 + X/counts/X_pca/X_umap/spatial/domain/obs/var 不变量。
回读验证 119 项全过（`cn_debug/s4_ko_readback.py`）：top_dr 独立重算、
ko_dr_expr 按 counts 独立重算、逐样本行数对账 uns、显著并集数、敲除基因
互斥、不变量全套。

产物 5 个：`spatial_knockout_mouse.h5ad`（17283 spots × 18628 基因，obs 加
ko_dr_expr，uns 加 `knockout` 口径键含 top_dr dict 键控 + sample_dr_stats_json）、
`spatial_knockout_human.h5ad`（6877 × 19613）、`ko_dr_genes_mouse.tsv`
（10341 行 = 9 样本 DR 表，列含 species/ko_gene/sample_id/Gene/Distance/Z/FC/
p-value/adjusted p-value/is_knockout_gene/significant）、`ko_dr_genes_human.tsv`
（3885 行）、`knockout_report.txt`。

**下一步：步骤 11（解卷积/细胞类型注释）参考来源需用户决策**（线3 raw_data 无
单细胞分辨率参考，硬约束不得取自 line2，仅允许 public_references/ 或线内自有
数据）。除步骤 11 外，line3_spatial 九项下游清单（空间域识别/marker/自相关/
邻域共现/区域组成/空间通讯/解卷积/目标基因定位/虚拟敲除）已落地 8 项。
像素空间指标一律逐样本（坑 #114）。
