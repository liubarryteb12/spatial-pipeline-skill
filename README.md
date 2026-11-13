# 空间转录组（Visium）流水线

<p align="center">
  <a href="https://github.com/liubarryteb12/spatial-pipeline-skill/actions/workflows/spatial_analysis.yml"><img alt="Spatial pipeline" src="https://github.com/liubarryteb12/spatial-pipeline-skill/actions/workflows/spatial_analysis.yml/badge.svg?branch=main"></a>
  <a href="https://github.com/liubarryteb12/spatial-pipeline-skill/actions/workflows/docs_check.yml"><img alt="Docs check" src="https://github.com/liubarryteb12/spatial-pipeline-skill/actions/workflows/docs_check.yml/badge.svg?branch=main"></a>
  <img alt="Python" src="https://img.shields.io/badge/Python-3.10%2B-3776AB">
  <img alt="Visium" src="https://img.shields.io/badge/Platform-10x%20Visium-1F9E89">
</p>

> **L3 分析线：空间转录组。** 从 Visium 的「表达矩阵 + 空间信息」两半数据出发，产出空间域、空间高变基因、spot 组成、邻域富集、空间通讯与空间拟时序，并在 GitHub Actions 上自动运行、打包为 artifact。

## 三线体系

| 线 | 仓库 | 数据 | 核心问题 |
|---|---|---|---|
| L1 | [`geo-normal-pipeline-skill`](https://github.com/liubarryteb12/geo-normal-pipeline-skill) | bulk 芯片 | 差异表达、通路富集、共表达、预后签名 |
| L2 | [`scrna-pipeline-skill`](https://github.com/liubarryteb12/scrna-pipeline-skill) | 单细胞 | 细胞类型、轨迹、通讯、调控 |
| **L3** | **本仓库：`spatial-pipeline-skill`** | Visium 空间转录组 | 空间域、空间高变基因、邻域、空间通讯 |

三条线共用同一套工程约定：配置驱动、状态落盘、验收判据、静态门禁。本仓库负责把空间特有的几何与邻域关系处理正确。

## 这个仓库是什么

- **框架，不是焊死的流水线。** 步骤、参数、产物与验收项都随分析需求变化。加一步、换一种方法、关掉某个可选步骤，都是预期用法。
- **空间信息是一等公民。** 坐标、缩放系数、H&E 图像与表达矩阵一起进入流程，不做“只读矩阵然后当单细胞处理”的降级。
- **两种划分都跑。** 空间平滑域与普通表达聚类并列输出，并用 ARI、NMI、邻居同域率量化差异。
- **云端仓库是唯一真源。** 本 skill 不需要安装，也不依赖某台机器上的目录，按需从 GitHub 拉取即可。

## 它解决什么问题

空间转录组里，“位置”本身就是信息，但位置也最容易悄悄出错。本流水线把三类空间特有错误挡在前面。

### 1. 只读了矩阵，没有空间信息

Visium 的数据是两半：表达矩阵 + `spatial.tar.gz`（坐标、缩放、H&E 图像）。只读矩阵，你有一个普通的单细胞数据集，却以为自己有空间数据，而且不会报错。`00_fetch.py` 强制要求两半同时给。

### 2. 坐标转置了，但看不出来

`obsm['spatial']` 的列序是 `(x, y) = (pxl_col_in_fullres, pxl_row_in_fullres)`。写反了散点图仍然画出一个“看起来像组织形状”的点阵，只是转置了。`tools/verify_spatial_alignment.py` 用网格几何判据定量检查：正确映射下对角线相关系数为 2.0000、交叉项为 0.1915；转置假设下几何分从 +1.81 翻成 -1.81。

### 3. 把表达簇当成了空间域

普通 Leiden 聚类不知道 spot 在哪，会把空间上连续的区域切成几块，或把不相邻的区域合成一个簇。流水线两种都跑并量化差别：不平滑 11 域（邻居同域率 0.507，随机基线 0.111），平滑后 13 域（同域率 0.668）。

此外还跑 SpaGCN 做交叉验证（`init="kmeans"`，域数对齐到 13），并报告它与内置划分的 ARI、NMI、邻居同域率。“跑通了”不是结论，两套划分是否在说同一件事才是。

## 分析流程

| 步骤 | 脚本 | 主要内容 |
|---|---|---|
| 00 取数 | `00_fetch.py` | 下载矩阵与空间信息，写入 `obsm['spatial']` 与 `uns['spatial']` |
| 01 质控 | `01_qc.py` | spot 质控、组织特异阈值、组织连通性检查 |
| 02 标准化 | `02_normalize.py` | HVG、标准化、PCA、PC 的空间投影诊断 |
| 03 空间域 | `03_spatial_domains.py` | 平滑 vs 不平滑对照、SpaGCN 交叉验证、域标签 |
| 04 SVG | `04_svg.py` | Moran's I、Geary's C、置换检验、SpatialDE 交叉验证 |
| 05 解卷积 | `05_deconvolution.py` | marker 打分、NNLS 解卷积、cell2location 边界 |
| 06 邻域 | `06_niche.py` | 邻域富集置换检验、共现曲线 |
| 07 空间通讯 | `07_spatial_communication.py` | 空间约束的配体-受体分析 |
| 08 空间拟时序 | `08_spatial_trajectory.py` | 空间感知与朴素排序的定量对比 |

每一步都是独立可执行脚本，可以单独重跑。`main_analysis.py` 负责编排与验收。

## 快速开始

```bash
pip install -r requirements.txt
python scripts/main_analysis.py --config assets/config.lymph_node.yml
```

参考配置约 2 分钟跑完，产出 20 张图与 24 个结果文件。推送到 GitHub 后，`.github/workflows/spatial_analysis.yml` 会自动运行并上传 artifact。

### 作为 Agent Skill 使用

```bash
./use.sh                          # 拉取或更新到 ~/.cache/dsh-skills/，打印路径
./use.sh --register               # 需要本机 agent 直接发现时加
./use.sh --ref v1.0               # 钉住某一版
SKILL=$(./use.sh --print-path)    # 只取路径，便于脚本使用
./use.sh --clean                  # 清理缓存副本并撤销注册
```

拉取后会校验 `SKILL.md` 存在。远端改名或换结构时会明确报错，不会安静地给一个空目录。

## 参考数据集与实测结果

参考数据是 10x Genomics 官方的正常人淋巴结 Visium 数据。流水线没有注入任何淋巴结先验，但独立复现了它的经典结构。

| 分析 | 结果 | 对应生物学 |
|---|---|---|
| SVG top | `IGKC`、`CCL21`、`FDCSP` | 浆细胞、T 细胞区趋化因子、生发中心 FDC |
| 邻域富集 | `FRC` 与 `T_cell`，z=13.8 | 成纤维网状细胞是 T 细胞区的支架 |
| 邻域富集 | `B_naive` 与 `FDC`，z=9.4 | 滤泡树突细胞在 B 细胞滤泡内 |
| 邻域富集 | `B_naive` 与 `FRC`，z=-28.5 | B 细胞滤泡与 T 细胞区互斥 |
| 空间通讯 | `CCL21-CCR7` z=3.31，`CCL19-CCR7` z=2.81 | 驱动 T 细胞区招募的趋化因子轴 |

其他基线：4035 spot × 36601 基因，QC 后 4025 spot，组织连通性 1 块；4000 个基因中 2890 个 BH 校正后 p<0.05；空间通讯 62/63 对配体-受体可用；验收 26 项，required、content、honesty 三类均为 0 失败。

> ARI 是一个范围，不是一个数。五轮 CI（同一份代码、同一批包版本）给过 0.3639 到 0.4216。主方法的数（Moran's I 0.7670、域数 13、邻居同域率）五轮逐位相同，可以当定值用。残留随机源的定位记录见 [`references/module0.md`](references/module0.md) §9.1 与 `domain_status.json` 的 `reproducibility` 字段。

## 点名工具的落地边界

文档点名了一批具名工具。本仓库不是“全都没做”，也不是“全都做了”，逐条状态写在 `domain_status.json`、`deconvolution_status.json`、`svg_status.json` 里。

| 工具 | 状态 |
|---|---|
| **SpaGCN** | 已跑，`init="kmeans"` 绕开 louvain 的 Python 3.12 编译链 |
| **SpatialDE** | 已跑，带 scipy `misc.derivative` 垫片 |
| **cell2location** | 可安装，但缺带标签的 scRNA 参考，状态为 `needs_reference` |
| STAGATE | 装不上：PyPI 三个名字全 404，官方版模块级依赖 torch-sparse（无 wheel） |
| BayesSpace、RCTD、SPARK-X、CellChat | R/Bioconductor 包，CI 不装 R |
| StPedf、SpaceFlow、ISORT、Stereopy-TGPI、stLearn | 不在 PyPI、依赖链跑不动或名字被顶 |

`needs_reference` 与“装不上”是两件事。把“缺数据”写成“装不上”，下一个人就会去折腾安装，方向完全错。理由与判据见 [`references/module0.md`](references/module0.md) §5。

## 产物

```text
results/<dataset_id>/
  figures/                     20 张 PNG，图名格式见 AGENTS 规则 23
  acceptance_report.json       26 项验收检查
  *_status.json                每步的 method、limitations 与回退原因
  figure_overflow.json         内容溢出记录，验收层会读并判红
data/<dataset_id>/
  domains.h5ad                 最终态，含 raw 全基因集、counts 层与域
  spatial_alignment_check.json 坐标对齐验证
```

只上传 `domains.h5ad`。每个 h5ad 带完整 36601 基因矩阵与 `layers['counts']`，单个 200 到 320 MB；`domains.h5ad` 是最终态，包含前面所有步骤的信息。中间态的排除在 `tools/check_artifact_paths.mjs` 里显式声明，有意排除与“忘了加清单”必须长得不一样。

## 质量门禁

```bash
node tools/check_py_syntax.mjs
python tools/check_py_names.py .
node tools/check_artifact_paths.mjs
node tools/check_legend_convention.mjs
node tools/check_figures.mjs results/lymph_node/figures
node tools/check_fig_names.mjs
node tools/check_fig_sizes.mjs results/lymph_node/figures
node tools/check_fig_dpi.mjs results/lymph_node/figures
node tools/check_palette.mjs
node tools/check_doc_refs.mjs
```

CI 还会跑 `check_legend_convention.mjs --selftest`、`check_figures.mjs --selftest`、`selftest_summary_fields.py` 与 `check_doc_refs.mjs --selftest`。门禁自带反向标定用例，避免“检查了但什么都没查”的假阴性。

## 适用范围与限制

- **Visium 的 spot 直径 55 μm，含 1 到 10 个细胞。** 不能说“某个域是某一种细胞”，只能说这个区域的细胞组成不同。
- **`builtin` 组成模式是 marker 打分法，不是解卷积。** 产物里的 `is_deconvolution: false`。它不给出细胞比例，只有相对空间趋势。要做真解卷积，配 `reference: h5ad` 与 scRNA-seq 参考。
- **共表达不等于通讯。** 空间通讯分析不做“通讯与否”的判定，只给空间约束下的共表达强度。
- **坐标对齐检查检测不到镜像。** 镜像保持相关系数绝对值不变。若怀疑镜像，需要人工核对 `aligned_fiducials.jpg` 里的基准框方位。
- **`fragmented_domains` 高不等于聚类失败。** 淋巴结的滤泡本身就是散布的多个斑块。

这些不是免责声明，是结论的适用范围。

## 目录结构

```text
scripts/         00-08 步骤脚本 + main_analysis.py
scripts/lib/     common.py（spatial_xy、日志、配置、绘图工具）
assets/          datasets.yml、config.lymph_node.yml、参考签名、配体受体库
tools/           静态门禁、坐标对齐验证、产物清单交叉检查
references/      methods.md、module0.md、troubleshooting.md、agents-detail.md
docs/            figure_groups.md
AGENTS.md        工程规则，全部来自实测踩过的坑
SKILL.md         Agent Skill 元数据
```

## 文档导航

| 文档 | 内容 |
|---|---|
| [`AGENTS.md`](AGENTS.md) | 工程规则、验收判据、踩坑记录 |
| [`references/methods.md`](references/methods.md) | 方法学依据与参数 |
| [`references/module0.md`](references/module0.md) | 模块零运行时规范与工具边界 |
| [`references/troubleshooting.md`](references/troubleshooting.md) | 故障排查 |
| [`docs/figure_groups.md`](docs/figure_groups.md) | 图组清单与命名 |

## 赞赏支持

如果这个仓库帮你节省了时间，欢迎请作者喝杯咖啡。赞赏完全自愿，不影响代码使用与后续更新。

<p align="center">
  <img src="assets/sponsor-qr.png" alt="伊霖的赞赏码" width="260">
</p>

> 感谢支持，我会尽快更新代码，维护代码。

## 授权

本仓库暂未附 LICENSE 文件。如需商用或二次分发，请先联系作者确认。
