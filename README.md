# L3 空间转录组与虚拟敲除流水线

> 本仓库是 **L3 生产分析线**的源码镜像。完整技术说明、输入契约、步骤细节与实测基线见
> [`src/line3_spatial/README.md`](src/line3_spatial/README.md)。

## 三线体系

| 线 | 仓库 | 数据 | 核心问题 |
|---|---|---|---|
| L1 | [`geo-normal-pipeline-skill`](https://github.com/liubarryteb12/geo-normal-pipeline-skill) | bulk 表达矩阵 | 差异表达、共表达、富集、PPI、机器学习 |
| L2 | [`scrna-pipeline-skill`](https://github.com/liubarryteb12/scrna-pipeline-skill) | 单细胞计数 | 细胞类型、轨迹、通讯、调控、虚拟敲除 |
| **L3** | **本仓库：`spatial-pipeline-skill`** | Visium 空间转录组 | 空间域、空间统计、通讯、定位、虚拟敲除 |

三条线默认互不交叉。线内允许按工具分工混用 R 与 Python，线间不互相取参考、不合并数据集、不跨线验证。约束全文见 [`src/line3_spatial/docs/isolation.md`](src/line3_spatial/docs/isolation.md)。

## 这套框架是什么

- **模块注册表**：根目录 [`workflows.json`](workflows.json) 是本线 stage、mod、依赖产物、成本估算与约束的唯一事实源。
- **步骤母版**：`src/line3_spatial/scripts/` 下每个文件就是一个可执行步骤，文件名即 step 名。
- **跨语言桥**：`src/bridge_guard.sh` 在 Python 与 R 之间传 `.h5ad` / `.rds` 前校验基因名、矩阵方向、细胞名、embedding 四项断言。
- **环境声明**：`src/line3_spatial/envs/manifest.txt` 记录主语言、版本、环境名与关键依赖，分析中途不更换。
- **检查点**：`src/line3_spatial/checkpoint/` 放跨语言检查点脚本。
- **共享规范**：`src/line3_spatial/GLOSSARY.md`、`src/line3_spatial/docs/` 保存术语、隔离、目录与编码规范。

## 分析阶段

| 阶段 | 模块 | 脚本 |
|---|---|---|
| 1 数据接入 | 空间数据装配 | `ingest_spatial` |
| 2 质控与归一化 | Spot QC、归一化与高变基因 | `qc_spots`、`normalize_hvg` |
| 3 降维与空间域 | 降维、空间域识别、域 Marker | `dim_reduce`、`spatial_domains`、`spatial_markers` |
| 4 空间统计 | 空间自相关、域共现、样本组成、空间通讯 | `spatial_autocorr`、`spatial_cooccurrence`、`spatial_composition`、`spatial_communication` |
| 5 可选支线 | 细胞类型注释、目标基因定位、富集、富集图件 | `celltype_annotation`、`spatial_localization`、`enrichment`、`enrichment_figures` |
| 6 虚拟敲除 | 空间虚拟敲除、DR 后分析 | `spatial_knockout`、`dr_post_analysis` |
| 7 图件汇总 | 线 3 图件包 | `spatial_figures` |

虚拟敲除是本线末端模块，产出 **DR genes**，不是 DEG。DR 基因与差异表达基因两个概念不混写。

## 输入

- Visium 表达矩阵：`mtx` 目录或 h5。
- `tissue_positions_list.csv`、`scalefactors_json.json`、组织图像。
- 双平台 soft 注释与 `grouping.txt` 分组表。
- ADT 数据可选，单独记录，不并入基因表达矩阵。
- 项目数据目录：`<BIO_DATA_ROOT>/projects/<project>/line3_spatial/raw_data/`。

同一数据集里可能有多个物种与平台。脚本按平台与基因组分别校验，再统一到样本级。

## 坐标口径

- `obsm["spatial"]` 是像素坐标，只用于样本内组织图叠加。
- 跨样本共同坐标框是 `obs["array_row"]` 与 `obs["array_col"]`。
- 不同样本的 `scalefactors` 可能不同，像素坐标跨样本不可比。

## 环境

`src/line3_spatial/envs/manifest.txt` 是环境正本：

- Python 3.12.3
- scanpy 1.12.4、anndata 0.13.4
- pandas 2.3.3、numpy 2.5.4、scipy 1.18.1
- liana 1.10.0、scTenifoldpy 0.5.1、gseapy 1.3.1
- matplotlib 3.11.2、h5py 3.16.0、pillow 13.0.0

改版本号视为新分析项目，不做热切换。

## 产出

产物落在 `<BIO_DATA_ROOT>/projects/<project>/line3_spatial/`：

```text
raw_data/    原始数据
processed/   中间 h5ad、坐标与域标签
results/     统计表、报告、DR 基因表
figures/     单图输出
```

每一步输出 `*_report.txt`，记录口径、阈值、诊断量、计时与产物清单。自校验不过不落盘。

## 运行方式

生产环境由任务引擎执行，不依赖 GitHub Actions：

1. 运维层把 `src/line3_spatial/scripts/` 部署到算力节点的 `bin/line3_spatial/`，并注入 PATH。
2. 任务引擎把工作目录切到 `<line>/processed/`，传入 `--input`、`--output`、`--seed`。
3. 引擎注入 `BIOANALYS_SEED`；缺失即退出码 3，不设默认值。
4. 解卷积参考只取 `public_references/` 或本线自有单细胞分辨率数据，不从线 2 读取。

退出码约定：

| 码 | 含义 |
|---|---|
| 0 | 成功，产物与台账齐 |
| 1 | 自校验失败（计算或产物不合口径） |
| 2 | 参数错误或输入契约不满足 |
| 3 | `BIOANALYS_SEED` 未注入 |
| 4 | 同步骤锁互斥 |
| 5 | 启动前内存余量不足 |

本仓库只保留源码与规范，不包含数据、产物、凭证或设备身份信息。路径中的 `<BIO_DATA_ROOT>` 由部署环境提供。

## 目录结构

```text
workflows.json              L3 模块注册表
src/
  bridge_guard.sh           跨语言桥接四项断言
  line3_spatial/
    README.md               生产说明与实测基线
    GLOSSARY.md             领域术语表
    docs/                   隔离、目录、编码规范
    envs/manifest.txt       环境与版本声明
    checkpoint/             桥接检查点
    scripts/                步骤母版
assets/sponsor-qr.png       赞赏码
```

## 迁移说明

2026-11-13 起，本仓库只保留线 3 生产母版。旧版 Scanpy / Squidpy / SpaGCN 全流程脚本、旧 CI 与旧图表模板已移除。生产行为以 `src/line3_spatial/` 与 `workflows.json` 为准。

## 赞赏支持

如果这个仓库帮你节省了时间，欢迎请作者喝杯咖啡。赞赏完全自愿，不影响代码使用与后续更新。

<p align="center">
  <img src="assets/sponsor-qr.png" alt="伊霖的赞赏码" width="260">
</p>

> 感谢支持，我会尽快更新代码，维护代码。

## 授权

本仓库暂未附 LICENSE 文件。如需商用或二次分发，请先联系作者确认。
