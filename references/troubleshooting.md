# 排错手册

按"症状 → 原因 → 修法"组织。所有条目都来自实际踩过的坑。

---

## 数据获取

### `HTTP Error 403: Forbidden` 从 10x CDN 下载

**症状**：同样的 URL 用浏览器打开是 200，用 Python 是 403。
错误信息只有 `HTTP Error 403: Forbidden`。

**原因**：10x 的 CDN（`cf.10xgenomics.com`）对 `Python-urllib/3.x`
这个 User-Agent 直接拒绝。

**修法**：带一个浏览器 UA。`00_fetch.py` 的 `download()` 已经带了
Chrome UA 与 `Accept: */*`。

**为什么难查**：报错信息完全看不出是 UA 的问题，
**很容易误判成"链接失效了"**。

---

### `ValueError: 配置的 source 段必须同时给 matrix_url 与 spatial_url`

**原因**：只给了矩阵，没给空间信息。或者写了 `source.dataset`
但那个键不在 `assets/datasets.yml` 里。

**修法**：确认 `datasets.yml` 里的条目有 `matrix_url` **和** `spatial_url`。

**不要试图绕过这个检查。** 只读矩阵 = 你有一个普通单细胞数据集，
但以为自己有空间数据。

---

### 坐标读出来少一个 spot / 列名全错

**原因**：Space Ranger v2 的 `tissue_positions.csv` **有表头**，
v1 的 `tissue_positions_list.csv` **没有**。把 v2 当 v1 读，
表头被当成数据行。

**修法**：`common.read_tissue_positions()` 两种都处理。
它按文件内容判断，不靠文件名。

---

## 坐标与图像

### spot 叠到 H&E 上明显错位 / 转置

**症状**：散点图单独看"像组织形状"，但叠到 H&E 上对不上。

**原因**（四种，都不会报错）：
1. `obsm['spatial']` 的 x/y 顺序反了（转置）
2. 轴向反了（镜像）
3. 忘了乘 `tissue_hires_scalef`（坐标是全分辨率尺度）
4. 用了 fullres 坐标配 hires 图

**修法**：跑 `python tools/verify_spatial_alignment.py --config <cfg>`。

它用**网格几何**判据：正确映射下 array(row,col) 与像素轴的相关系数
矩阵对角线是 2.0000、交叉项 0.1915；转置假设下几何分从 **+1.81
翻成 −1.81**。

**注意这个工具检测不到镜像** —— 镜像保持 |相关系数| 不变。
若几何分正常但仍怀疑镜像，必须人工核对 `aligned_fiducials.jpg`
里的基准框方位。

---

### 对齐检查"通过"了，但其实什么都没验证

**症状**：`verify_spatial_alignment.py` 的第一版对四个假设
（含转置、镜像）都给出相同的分数，然后报"通过"。

**原因**：判据是"spot 中心处是不是组织像素"，而实测整张 hires 图的
组织像素占比是 **1.0000**（组织铺满整帧，没有白色背景）——
任何点都"落在组织上"，判据失去区分力。

**教训**：**判据必须有区分力，否则它给的是虚假的安心。**
现在这个判据仍然报出来，但会标注 `discriminative: false`。

---

### spot 半径画出来不对

**原因**：`spot_diameter_fullres` 是全分辨率尺度下的直径，
画在 hires 图上要乘 `tissue_hires_scalef`，画在 lowres 图上要乘
`tissue_lowres_scalef`。

**修法**：用 `common.spot_radius_plot_units(scalefactors, "hires")`。

---

## 质控

### 过滤后 spot 少得离谱

**原因**：`max_pct_mt` 用的是别的组织的阈值。

阈值是**组织特异**的：PBMC 5-10 / 淋巴结 20 / 心肌 30-50 / 肝 30+。

**修法**：看 `qc_status.json` 的 `filter_chain`，找出是哪一条滤掉了大部分。
`01_qc.py` 在过滤后 spot 少于 `min_spots` 时会直接抛异常并打印过滤链。

---

### `pct_counts_mt` 全是 0

**原因**：物种判断错了。人是 `MT-`，小鼠是 `mt-`。

**为什么危险**：过滤形同虚设，但**日志里看不出任何异常**。

**修法**：`01_qc.py` 在匹配到 0 个线粒体基因时直接抛异常，
错误信息里会打印它用的前缀和物种判断。

---

### 过滤后组织被割裂

**症状**：`qc_status.json` 的 `connectivity.fragmented` 是 `true`。

**原因**：滤掉的 spot 集中在组织某一区域，把组织切成了几块。

**后果**：空间邻域/niche/通讯**只在块内有效**。

**修法**：检查阈值是不是过严；或者组织本身就有分离的区域
（那种情况下要在报告里说明）。

---

## 空间域

### 所有域都被标成同一种细胞类型

**症状**：13 个域里 12 个都被标成 `Plasma_cell`。

**原因**：直接比"该域里各类型 marker 的平均表达"。浆细胞的 marker
（MZB1/JCHAIN/IGHG1/IGKC）在淋巴结里表达量本来就极高，
绝对表达一比就压倒所有其他类型。

**修法**：按类型跨域做 z-score，再取最高。
问的是"**哪个域对这个类型相对最富集**"，不是"哪个类型表达最高"。

**这类错误的共同特征**：**一个看起来像答案、但不含信息的标签。**

---

### 域数随平滑强度增加（反直觉）

**症状**：α=0 时 9 个域，α=0.5 时 13 个域。

**原因**：平滑让嵌入更平滑，同一分辨率下 Leiden 能分出更多
互不重叠的区域。

**这不是 bug。** 看 `smoothing_scan.csv` 里的 `neighbor_same_frac`
（0.537 → 0.672）才是有意义的指标。

---

### `fragmented_domains` 很高

**症状**：9 个域里 9 个都是多块的。

**先别改代码。** 淋巴结的滤泡、肿瘤的癌巢**本身就是散布的多个斑块** ——
同一个域出现在多个不相邻位置是生物学事实。

只有对"应当连续"的组织（脑的层状结构、上皮分层）才能当缺陷看。

---

## 空间高变基因

### 置换检验太慢

**原因**：在全基因集（36601 基因）上做置换。

**修法**：按方差粗筛（`svg.max_genes`，默认 4000）后做置换检验。
实测 4000 个基因 × 100 次置换 = 48 秒。

**注意**：粗筛会漏掉低方差但空间结构强的基因。
这个取舍写在 `svg_status.json` 的 `gene_set_note` 里。

---

### `ModuleNotFoundError: No module named 'skmisc'`

**原因**：`seurat_v3` 口味的 HVG 需要 `scikit-misc`。

**修法**：`pip install scikit-misc`。

**如果不装**：流水线会降级到 `seurat` 口味，并在
`normalize_status.json` 里记录 `hvg_flavor_requested` vs
`hvg_flavor_used` 与 `hvg_fallback`。

**为什么必须记录**：两种口味选出的基因集合不同，
下游 PCA/空间域/SVG 全部跟着变。

---

## 组成 / 解卷积

### 所有细胞类型的组成几乎相同（≈ 1/n）

**症状**：12 种类型的组成全部约等于 0.085 ≈ 1/12，
重建误差中位 0.64。

**原因**：用 NNLS 解卷积，但参考谱被构造成"只在 marker 基因位置
非零，其余为 0"。模型在说"**非 marker 基因应该不表达**"，
而实际上它们表达得很正常。残差被这些基因主导，
NNLS 只能给一个无意义的均匀解。

**根本原因**：**没有外部参考时，marker-only NNLS 不是适定问题。**
方程 `y[markers] = Σ_k w_k · s_k[markers]` 里的 `s_k` 是未知的，
从目标数据本身估出来是循环论证。

**修法**：`builtin` 模式用 marker 打分法（产物里
`is_deconvolution: false`）；要真解卷积就配 `reference: h5ad`
+ 一个 scRNA-seq 参考。

---

### 重建误差普遍很高

**症状**：`deconvolution_status.json` 里
`reconstruction_error.frac_unreliable` 很高。

**原因**（真解卷积模式下）：
- 参考的细胞类型覆盖不全（目标里有参考里没有的类型）
- 参考与目标批次差异大
- 共同基因太少

**修法**：看 `deconvolution_error_map.png` —— 误差高的 spot 是否
集中在某个区域？那通常说明那个区域有未列入的类型。

---

## 空间通讯

### 可用的配体-受体对很少

**症状**：`communication_status.json` 里 `usable_fraction` 很低。

**原因**：用了 HVG 子集。配体/受体大多是低表达基因，几乎不进 HVG。
实测 Part 2 在 2000 个 HVG 上 38 对里只有 3 对可用；
全基因集上是 27 对。

**修法**：确保 `adata.raw` 存在（全基因集）。
`get_expression()` 在 `.raw` 为空时直接抛异常，不会悄悄退回 HVG。

---

### YAML 解析报 `expected alphabetic or numeric character, but found '*'`

**原因**：YAML 里以 `*` 开头的**裸标量**被当成别名（alias）引用。

例如：
```yaml
note: **方向反了**（CTLA4 在 T 细胞上）
```

**修法**：加引号：
```yaml
note: "**方向反了**（CTLA4 在 T 细胞上）"
```

---

## 绘图

### 图上标题显示成一个个方框（豆腐块）

**原因**：图标签里写了中文，而 matplotlib 用 DejaVu Sans，
**没有 CJK 字形**。

**为什么难查**：代码不报错、CI 是绿的、`check_figures.mjs` 也只看到
"有墨迹"。**只有打开图才发现标题不可读。**

**修法**：图上标签用英文；中文解释放注释和 JSON 产物里。
`tools/check_py_syntax.mjs` 会检查这一条。

---

### `TypeError: ... got an unexpected keyword argument 'return_fig'`

**原因**：scanpy 1.12 移除了 `sc.pl.highly_variable_genes(return_fig=True)`。

**修法**：用 `plt.gcf()` 拿当前 figure。

---

### `TypeError: Axes.boxplot() got an unexpected keyword argument 'labels'`

**原因**：matplotlib 3.9 把 `labels` 改名 `tick_labels`，3.11 直接报错。

**修法**：`requirements.txt` 钉上限；版本敏感的调用加 try/except 回退。

---

## CI

### artifact 里缺文件，但 job 是绿的

**原因**：验收检查的是 **runner 工作目录**里的文件，
artifact 的 `path:` 是一个**完全独立的列表**。两者不一致时：
验收通过、job 绿，但用户下载 artifact 后核心产物不在里面。

**修法**：`node tools/check_artifact_paths.mjs`。

**注意**：显式声明的中间态（`INTERMEDIATE` 白名单）不算漏 ——
有意排除和"忘了加清单"必须长得不一样。

---

### 装了依赖之后才发现语法错误

**原因**：静态检查没跑在装依赖之前。

**修法**：workflow 里 `check_py_syntax.mjs` 必须在 `pip install` 之前。
语法错误 1 秒能发现，不该等 3 分钟。

---

### `harmonypy` 编译失败

**原因**：Python 3.13/3.14 上还没有 wheel，要从源码编译，
常因缺 numpy 头文件失败。

**修法**：CI 用 **Python 3.12**（本流水线不做多样本整合，
但同类包都有这个问题）。

---

### 本地和 CI 的聚类标签不一样

**症状**：本地 clusters `['6','2','7']`，CI `['3','1','0']`。

**原因**：**聚类的索引标签在不同平台上不可移植**，即使固定了随机种子。
线程调度（`OMP_NUM_THREADS` 等）与 OpenBLAS 的运行时 SIMD 内核分派
都会造成浮点末位差异，累积到聚类就变成标签顺序不同。

**关键**：**簇的数量与细胞类型构成是一致的。** 不一致的只是标签编号。

**修法**：不要跨平台比对标签编号。要比就比簇数、
比每个簇的 marker、比邻域富集的结构。

**不要试图"修"这个。** 它来自浮点非结合性，不是 bug。
