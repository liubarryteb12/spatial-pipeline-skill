# spatial-pipeline-skill

## 这是一个**框架**，不是一条焊死的流水线

本仓库提供的是**生信分析的骨架与判据**：数据门禁、方法学约定、验收项、
产物清单。**具体跑什么由输入数据和配置文件决定**，步骤本身可增删 ——
加一步、换一种方法、关掉某个可选步骤，都是预期用法，不是"改坏了"。

所以「这个仓库能做什么」的答案在 `assets/config.*.yml` 和验收项里，
**不在目录结构里**。`scripts/` 中没有任何一处硬编码某个疾病或某个平台。

后续会有一个独立的「流水线编排模块」，让使用者挑选分析模块并串起来，
再与本仓库对接。**那部分不在本仓库职责范围内** —— 本仓库只负责把每一步做对。

## 怎么拿到它：云端仓库是唯一真源

本 skill **不需要"安装"**，也不依赖任何一台机器上的目录。真源是 GitHub 仓库：

    https://github.com/liubarryteb12/spatial-pipeline-skill

要用的时候从云端拉下来：

```bash
./use.sh                          # 拉取/更新到 ~/.cache/dsh-skills/，打印路径
./use.sh --register               # 需要本机 agent 直接发现它时才加
./use.sh --ref v1.0               # 钉住某一版
SKILL=$(./use.sh --print-path)    # 只取路径，便于脚本里用
./use.sh --clean                  # 清掉缓存副本并撤销注册
```

**拉取后会校验 `SKILL.md` 存在。** 远端改名或换结构时会明确报错，
而不是安静地给一个空目录 —— 实测过：`git clone` 失败时后面的步骤照样会跑，
最后就是靠这道校验拦住的。

> `use.sh` 里每个可能失败的步骤都显式 `|| die`，**不依赖 `set -e`**。
> 实测（bash 5.3）在 `resolved="$(pull)"` 这种「函数在命令替换里」的结构下，
> 函数内部的失败不一定会中止外层脚本。出错的路径必须自己说出来。

端到端的 Visium 空间转录组流水线，跑在 GitHub Actions 上，产物作为 artifact 下载。

参考数据集：**正常人淋巴结**（10x Genomics 官方 Visium，4035 spot）。

---

## 快速开始

```bash
pip install -r requirements.txt
python scripts/main_analysis.py --config assets/config.lymph_node.yml
```

约 2 分钟跑完，产出 17 张图 + 24 个结果文件。

**云端**：推送到 GitHub，`.github/workflows/spatial_analysis.yml` 自动跑。

---

## 它解决什么问题

空间转录组的数据里，"位置"本身就是信息 —— 但位置也是最容易悄悄出错的地方。

这个流水线把**三类空间特有的错误**显式地挡在前面：

### 1. 只读了矩阵，没有空间信息

Visium 的数据是**两半**：表达矩阵 + `spatial.tar.gz`（坐标/缩放/H&E 图像）。
只读矩阵，你有一个普通的单细胞数据集，但**以为自己有空间数据** ——
而且不会报错，因为矩阵本身完全正常。

`00_fetch.py` 强制要求两半都给。

### 2. 坐标转置了，但看不出来

`obsm['spatial']` 的列序是 `(x, y) = (pxl_col_in_fullres, pxl_row_in_fullres)`。
写反了散点图**仍然画出一个"看起来像组织形状"的点阵**，只是转置了。

`tools/verify_spatial_alignment.py` 用**网格几何**判据定量检查：
正确映射下 array(row,col) 与像素轴的相关系数矩阵对角线是 2.0000、
交叉项 0.1915；转置假设下几何分从 **+1.81 翻成 −1.81**。

### 3. 把"表达簇"当成了"空间域"

普通的 Leiden 聚类不知道 spot 在哪，会把空间上连续的区域切成几块、
把不相邻的区域合成一个簇。

流水线**两种都跑并量化差别**：不平滑 9 域（邻居同域率 0.536，
随机基线 0.133）→ 平滑后 13 域（同域率 **0.672**）。

---

## 流程

```
00_fetch      下载矩阵 + 空间信息 → obsm['spatial'], uns['spatial']
01_qc         spot 质控（组织特异阈值）+ 组织连通性检查
02_normalize  HVG + 标准化 + PCA + PC 的空间投影诊断
03_domains    空间域（平滑 vs 不平滑对照）+ 域标签
04_svg        Moran's I / Geary's C + 置换检验
05_deconvo    spot 组成（marker 打分 / 真解卷积）
06_niche      邻域富集（置换检验）+ 共现曲线
07_comm       空间约束的配体-受体分析
08_spatial_traj 空间拟时序（空间感知 vs 朴素的定量对比）
```

每一步都是独立的可执行脚本，可以单独重跑。

---

## 实测结果（淋巴结）

流水线**没有注入任何淋巴结的先验知识**，但独立复现了它的经典结构：

| 分析 | 结果 | 对应的已知生物学 |
|---|---|---|
| SVG top | `IGKC`、`CCL21`、`FDCSP` | 浆细胞、T 细胞区趋化因子、生发中心 FDC |
| 邻域富集 | `FRC ↔ T_cell` **z=13.8** | 成纤维网状细胞是 T 细胞区的支架 |
| 邻域富集 | `B_naive ↔ FDC` **z=9.4** | 滤泡树突细胞在 B 细胞滤泡内 |
| 邻域富集 | `B_naive ↔ FRC` **z=−28.5** | **B 细胞滤泡与 T 细胞区互斥** |
| 空间通讯 | `CCL21-CCR7` z=3.31、`CCL19-CCR7` z=2.81 | 驱动 T 细胞区招募的趋化因子轴 |

三条独立分析指向同一个已知生物学 —— 这是实现正确性的强证据。

其他数字：

- 4035 spot × 36601 基因 → QC 后 **4025 spot**，组织连通性 1 块
- 域标签 10 种：`B_germinal_center`(z=+3.20)、`T_cell`(+2.78)、
  `Smooth_muscle`(+2.57)、`Plasma_cell`(+2.50)…；**7/13 个域
  `z_margin <= 0.5` 被标为不可信**
- SVG：4000 个基因里 **2890 个** BH 校正后 p<0.05
- 空间通讯：**62/63 对 LR 可用**（用全基因集；Part 2 在 HVG 上只有 3/38）
- 验收：**26 项检查，required / content / honesty 三类全部 0 失败**

---

## 产物

```
results/lymph_node/
  figures/                     17 张 PNG
    domains_on_he.png          空间域叠在 H&E 上（判断域是否对应组织学）
    pca_on_tissue.png          PC 的空间投影
    smoothing_scan.png         平滑强度 vs 空间连贯性
    svg_top_genes.png          top 空间高变基因
    niche_enrichment_*.png     邻域富集 z-score 热图
    deconvolution_spatial.png  各类型相对权重的空间分布
  acceptance_report.json       26 项验收检查
  *_status.json                每步的 method + limitations
data/lymph_node/
  domains.h5ad                 最终态（含 raw 全基因集 + counts 层 + 域）
  spatial_alignment_check.json 坐标对齐验证
```

**只上传 `domains.h5ad`。** 每个 h5ad 带完整 36601 基因矩阵 + `layers['counts']`，
单个 200-320 MB；`domains.h5ad` 是最终态，包含前面所有步骤的信息。
四个全传约 1 GB，其中约 700 MB 冗余。中间态的排除在
`tools/check_artifact_paths.mjs` 里**显式声明**（`INTERMEDIATE` 白名单）——
有意排除和"忘了加清单"必须长得不一样。

---

## 限制

**这些不是免责声明，是结论的适用范围。**

- **Visium 的 spot 直径 55 μm，含 1-10 个细胞。** 所以不能说"某个域是
  某一种细胞"，只能说这个区域的细胞组成不同。
- **`builtin` 组成模式是 marker 打分法，不是解卷积**
  （产物里 `is_deconvolution: false`）。不给出细胞比例，只有相对空间趋势。
  要用真解卷积：配 `reference: h5ad` + scRNA-seq 参考。
- **共表达 ≠ 通讯。** 空间通讯分析不做"通讯与否"的判定，
  只给空间约束下的共表达强度。
- **坐标对齐检查检测不到镜像** —— 镜像保持 |相关系数| 不变。
- **`fragmented_domains` 高不等于聚类失败** —— 淋巴结的滤泡本身就是
  散布的多个斑块。

---

## 目录

```
scripts/         8 个步骤 + main_analysis
scripts/lib/     common.py（spatial_xy、日志、配置）
assets/          datasets.yml, config.lymph_node.yml,
                 reference_signatures.yml, ligand_receptor.yml
tools/           check_py_syntax.mjs     语法 + 4 条仓库规则
                 check_figures.mjs       像素级检查空白图
                 check_artifact_paths.mjs 产物与 artifact 清单一致性
                 verify_spatial_alignment.py  坐标对齐定量验证
references/      methods.md, troubleshooting.md
AGENTS.md        14 条工程规则（全部来自实测踩过的坑）
```

---

## 工程原则

1. **静态检查在前，装依赖在后。** 语法错误 1 秒能发现，不该等 3 分钟。
2. **"没做"必须和"做了没问题"长得不一样。** 做不了就写
   `status: not_done` + 原因，不用代理指标冒充。
3. **判据必须有区分力。** 一个所有假设都通过的检查等于没检查。
4. **优化后必须验证结果一致。** 只看"跑通了"不够。
5. **每轮 CI 只修日志里明确显示的问题。**

详细展开见 `AGENTS.md`。
