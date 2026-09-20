# 空间转录组流水线的工程规则

> **治理层（2026-09-20 起）**：本仓库是 `scientific_agent_skill` 工作区三仓库之一，
> 受工作区治理层约束：一切产物只写工作区内；任务先登记在
> `governance/02_TASKLIST.md`；推送前跑 `node governance/hooks/pre-push.mjs`；
> checkpoint 台账见 `governance/04_CHECKPOINT_PLAN.md`；行为准则
> `governance/01_SPEC_v1.0.md`。冲突按规范 §7.4 报告裁决。

姊妹项目 `scrna-pipeline-skill/AGENTS.md` 与
`geo-normal-pipeline-skill/AGENTS.md` 的规则在这里同样适用。
本文只写**空间数据特有的**部分。

---

## 1. Visium 的数据是两半，只给矩阵就退化成单细胞

表达矩阵（`filtered_feature_bc_matrix.h5`）**不含**空间信息。
空间信息在 `spatial.tar.gz` 里：spot 坐标、缩放系数、H&E 图像。

只读矩阵 = 你有一个普通的单细胞数据集，但以为自己有空间数据。
`00_fetch.py` 强制要求 `matrix_url` 与 `spatial_url` **同时给**，
缺一个就报错，不给"降级成单细胞"的选项。

## 2. `obsm['spatial']` 的列序是 `(x, y) = (pxl_col_in_fullres, pxl_row_in_fullres)`

**写反了不会报错。** 散点图仍然画出一个"看起来像组织形状"的点阵，
只是转置了。只有把点叠到 H&E 图上、或者做定量的网格几何检查，
才能发现不对。

所以坐标一律通过 `common.spatial_xy()` 取，不直接访问
`adata.obsm["spatial"]`。`tools/check_py_syntax.mjs` 会拦住直接访问。

注意 `imshow(img)` 用的是 `(row, col)` 顺序，而 scatter 用 `(x, y)` ——
叠图时先 `imshow` 再 `scatter(x, y)`，不要反过来。

## 3. 坐标验证要用有区分力的判据

`tools/verify_spatial_alignment.py` 的第一版用"spot 中心处是不是组织像素"
作为判据。**对淋巴结数据完全无效** —— 实测整张 hires 图的组织像素占比是
**1.0000**（组织铺满整帧，没有白色背景），所以四个假设（含转置、镜像）
都得 1.000 分，检查报"通过"而实际上什么都没验证。

**判据必须有区分力，否则它给的是虚假的安心。**

现在用两条真正有区分力的判据：
- **网格几何**（主判据）：array(row,col) 与像素轴的相关系数矩阵。
  实测正确映射下对角线 2.0000、交叉项 0.1915；转置假设下翻成
  −1.8085。**转置立刻暴露。**
- **网格规整度**：最近邻距离变异系数（实测 0.031）。尺度用错会暴露。

**这个工具检测不到镜像** —— 镜像保持 |相关系数| 不变。要判镜像需要看
`aligned_fiducials.jpg` 里的基准框方位，本工具没做。若怀疑镜像，
必须人工核对。

## 4. QC 阈值是组织特异的

`max_pct_mt`：PBMC 5-10 / 淋巴结 20 / 心肌 30-50 / 肝 30+。
用 PBMC 的阈值套到心肌上会滤掉大部分 spot。

而且 **spot 是 1-10 个细胞的混合物**，所以 `pct_counts_mt` 的分布比
单细胞**窄得多**（多个细胞的线粒体信号被平均掉了）。

过滤后必须检查**组织是否被割裂**（`check_connectivity`）——
单细胞里滤掉细胞没有几何后果，空间数据里会留洞。

## 5. 小鼠的线粒体前缀是 `mt-`，人是 `MT-`

用错前缀会让 `pct_counts_mt` 全 0，**过滤形同虚设，而日志里看不出任何异常**。
`01_qc.py` 在匹配到 0 个线粒体基因时直接抛异常。

## 6. 空间域 ≠ 表达簇

不平滑的 Leiden 聚类不知道 spot 在哪，会把连续区域切碎、
把不相邻区域合并。`03_spatial_domains.py` 两种都跑并量化差别：

- `neighbor_same_frac`：相邻 spot 同域比例。**必须与
  `random_baseline_same_frac`（Σ 域占比²）比较才有意义。**
  实测：不平滑 11 域 / `0.507` vs 基线 `0.111`；平滑后 13 域 /
  `0.668` vs `0.096`。

  > **这几个数改过。** 早期版本是"不平滑 9 域 / 0.536 vs 0.133；
  > 平滑后 0.672"，平滑后基线还写过 `0.099`（当前产物是 **0.0958**）。
  > **文档里写死的数值一定会过时**（见规则 21.3）——
  > 引用前先看 `domain_status.json`，不要从文档抄。

**`fragmented_domains` 高不等于聚类失败。** 淋巴结的滤泡、肿瘤的癌巢
本身就是散布的多个斑块 —— 同一个域出现在多个不相邻位置是生物学事实。
只有对"应当连续"的组织（脑的层状结构、上皮分层）才能当缺陷看。

## 7. 域标签要按类型做 z-score，不能比绝对表达

第一版直接比"该域里各类型 marker 的平均表达"，结果 13 个域里
**12 个都被标成 `Plasma_cell`** —— 浆细胞的 marker（MZB1/JCHAIN/IGHG1/IGKC）
在淋巴结里表达量本来就极高，绝对表达一比就压倒所有其他类型。

**这和"没有参考时用 NNLS 解卷积"是同一类错误：一个看起来像答案、
但不含信息的标签。**

正确的问法是"**哪个域对某个类型相对最富集**"：按类型跨域做 z-score
（消掉"某些类型 marker 天生高表达"的偏差），再取最高。
改后 13 个域得到 11 种不同标签，且 `z_margin <= 0.5` 的 8 个被标为不可信。

> **这两个数改过：** 早先写的是"10 种不同标签、7 个不可信"，
> 而近三轮 CI 都是 **11 种 / 8 个**。域标签依赖平滑后的 Leiden 划分，
> 划分一动它就动 —— **引用前看 `domain_status.json`，不要从文档抄**（规则 21.3）。

## 8. 没有外部参考时，marker-only NNLS 不是适定问题

`05_deconvolution.py` 的第一版用 NNLS 解 `y ≈ S^T w`，其中 S 的每一行是
一个细胞类型的参考谱；没有外部参考时把参考谱构造成"只在 marker 基因
位置取平均表达，其余为 0"。**结果完全错**：

12 种类型的组成全部约等于 0.085 ≈ 1/12，重建误差中位 0.64，
3836/4025 个 spot 超过误差阈值。

原因是模型在说"**非 marker 基因应该不表达**"（所有参考谱在那里都是 0），
而实际上它们表达得很正常。残差被这些基因主导，NNLS 只能给一个
无意义的均匀解 —— **而它看起来仍然像一组比例**。

更根本的问题是：方程 `y[markers] = Σ_k w_k · s_k[markers]` 里的
`s_k`（每种类型在 marker 上的参考表达谱）是未知的，从目标数据本身
估出来是循环论证。

所以 `builtin` 模式改成 **marker 打分法**，并在产物里写明
`is_deconvolution: false`。要用真正的解卷积，配 `reference: h5ad` +
一个 scRNA-seq 参考（那时参考谱覆盖全部共同基因，问题才适定）。

## 9. 需要特定基因的分析必须用 `adata.raw`（全基因集）

配体、受体大多是低表达的，几乎不进 HVG。Part 2 实测：2000 个 HVG 上
38 对 LR 只有 3 对可用；全基因集上 27 对。

空间这边实测：全基因集上 **62/63 对可用**。

**假阴性看起来和"真的没有"一模一样**，所以宁可报错也不要悄悄退回 HVG。
`07_spatial_communication.py` 的 `get_expression()` 在 `adata.raw` 为空时
直接抛异常。

## 10. 共表达 ≠ 通讯

`07_spatial_communication.py` 算的是"空间约束下的配体×受体共表达强度"，
**不是"通讯与否"的判定**。真正的细胞通讯推断需要 CellPhoneDB / CellChat
那样的统计框架（含置换检验与受体复合物建模）。

两个基因在同一个 spot 里高，可能是同一个细胞表达了两者（自分泌），
也可能是两个细胞紧邻 —— 而 Visium 的 spot 含 1-10 个细胞。

产物里必须有 `not_a_call` 字段说明这一点。

## 11. 图的标签里不能用中文

matplotlib 用 DejaVu Sans，**没有 CJK 字形** —— 图上显示成一个个方框
（豆腐块），而代码不报错、CI 是绿的、`check_figures.mjs` 也只看到"有墨迹"。
只有打开图才发现标题不可读。

中文解释一律放注释和 JSON 产物里；图上只用英文。
`tools/check_py_syntax.mjs` 会检查。

## 12. 中间态不是交付物，但排除必须显式

每个 h5ad 都带着完整的 36601 基因矩阵 + `layers['counts']`，
单个 200-320 MB。`domains.h5ad` 是最终态，**包含**前面所有步骤的信息。

四个全上传约 1 GB，其中约 700 MB 冗余。所以只上传最终态 ——
但排除在 `tools/check_artifact_paths.mjs` 里以 `INTERMEDIATE` 白名单
**显式声明**。有意排除和"忘了加进清单"必须长得不一样。

## 13. 每轮 CI 只修日志里明确显示的问题

不要凭猜测改代码。先读日志，再改，一次只改日志支持的那一处。

## 14. 优化后必须验证结果一致

`spatial_xy()` 重构涉及 7 个文件、12 处调用。改完立刻对比重构前后的
`communication_lr_scores.csv` 与 `svg_results.csv` —— **0 行差异**。
只看"跑通了"是不够的。

## 15. 图幅按毫米，宽度夹在标准栏宽内

参考规范：K-Dense `scientific-visualization` skill（样式文件已 vendored 到
`assets/publication.mplstyle`，与 `scrna-pipeline-skill` 逐字节一致）。

**期刊栏宽是按毫米规定的**，英寸是排版软件内部单位。写英寸时"这图多宽"
要靠换算才知道，写毫米时一眼能对上投稿要求。

| 常量 | 值 | 用途 |
|---|---|---|
| `W_SINGLE` | 89 mm | 单栏 |
| `W_ONE_HALF` | 136 mm | 一栏半 |
| `W_DOUBLE` | 183 mm | 双栏（通栏）|

- **宽度随类别数增长的图必须夹住**：`min(W_DOUBLE, max(W_ONE_HALF, ...))`。
  无上限增长会画出装不进任何期刊一页的图 —— 实测修之前最宽的
  `domain_markers_dotplot` 是 370 mm。
- `save_fig()` 默认**不再用 tight bbox**。`bbox_inches="tight"` 会**改变物理
  输出尺寸**，让上面的毫米约定失效。溢出改由 `_content_overflow()` 检测并告警。
- **`set_seed()` 末尾会 `apply_style()`** —— rcParams 在**图创建时**就被读取，
  在 `save_fig()` 里设样式已经太晚。
- 图上文字一律英文（见规则 11）。

**两处与参考规范的偏差**（已写在 `.mplstyle` 文件头）：

1. `figure.constrained_layout.use: True` 全局开启。参考规范要求逐图 opt-in，
   但本仓库有 19 张图、9 个脚本，逐处改容易漏。
2. `font.sans-serif: DejaVu Sans, Arial, Helvetica`。参考规范首选 Arial，
   但 Ubuntu CI 上没有 Arial 而 Windows 上有 —— 会导致 CI 与本地渲染出
   **不同的字形**，破坏可复现性。DejaVu Sans 随 matplotlib 分发，处处一致。

**constrained layout 不会自动折行长标题。** 实测 `domains_on_he` 的单行
suptitle 超出 183 mm 宽 2.9%，被静默裁掉（`savefig.bbox: standard` 下文件照样
生成、`check_figures.mjs` 也照样报"有墨"）。长标题必须自己换行。

## 16. 步骤崩溃不能靠旧文件冒充成功

`main_analysis.py` 在每步开跑前**先删掉该步的状态文件**，并在验收里
逐步骤检查 `status == "ok"`。

只检查"产物文件在不在"是不够的：步骤崩溃时旧文件还在，产物检查会通过，
而本轮实际上什么都没产出。姊妹项目 Part 2 实测踩过 —— `07_grn` 因漏 import
崩溃，而它的状态文件是上一轮的，验收照样全绿。

## 17. 静态检查要挡住"未定义名字"

`py_compile` **只做编译，看不出未定义名字**。漏 import 一个 `W_SINGLE`
时它照样报"语法通过"，要等运行时才炸 —— 实测因此白跑一整轮流水线。

`tools/check_py_syntax.mjs` 现在会检查：用到的 `W_SINGLE` / `W_ONE_HALF` /
`W_DOUBLE` / `mm` / `PAL` / `PAL_CYCLE` / `apply_style` 是否都 import 了。

名单是**写死的**，不是"common 导出的所有名字"。后者会把函数参数名当成用法
（`alignment.py` 里 `def verify_alignment(adata, log_info=None)` 的 `log_info`），
要正确处理得做作用域分析 —— 那是重写一个 linter。同时会剥掉注释和字符串再扫，
避免"名字只出现在注释里"的误报。

## 18. 每轮运行必须留下可追溯的运行清单（模块零）

参考规范：三大部分整合文档的「模块零」（§0.2–§0.4）。姊妹项目
`scrna-pipeline-skill/AGENTS.md` 规则 16 有完整说明，这里只写空间特有的。

`common.py` 的清单层产出 `results/<dataset_id>/run_manifest.json`。
空间这边有两点不同：

1. **`KEY_PACKAGES` 包含大量"本仓库没装"的工具** —— STAGATE / SpaGCN /
   BayesSpace / cell2location / Bering / BOMS / SpaceFlow / ISORT 等。
   它们会记成 `null`，这是**有意为之**：空间方法的可选项比单细胞多得多，
   "哪些没装"本身就是结论适用范围的一部分。
   全部省略键会让清单看起来"该有的都有"。
2. **`record_input` 在 `run_acceptance` 开头跑，不在 `run_steps` 开头。**
   可选步骤（去卷积、空间轨迹）这轮有没有产物，要等步骤跑完才知道。

**人工复核节点默认 `pending`，不算失败**（`cell_segmentation` /
`domain_number` / `deconv_reference` / `spatial_traj_direction`）——
和规则 16 同一条理由：判成 FAIL 会让每个 job 都红，反而没人看。

**`init_manifest` 必须清掉上一轮**，否则上轮的清单冒充本轮，
比没有清单更糟（同规则 16）。

## 19. 老包要兼容垫片，但垫片必须可验证、可复算、如实记录

文档 §3.4 点名 SpatialDE。**SpatialDE 1.1.3 是 2019 年的包，
在当代依赖上有两处独立的不兼容**（都实测确认，不是推测）：

| # | 位置 | 症状 | 处理 |
|---|---|---|---|
| 1 | `base.py:12` `from scipy.misc import derivative` | scipy 1.12 移除了 `scipy.misc.derivative` → `ImportError`，**装得上但导不进来** | 垫片补回该名字（3 点中心差分，Vandermonde 解权重）|
| 2 | `base.py:432` → `util.py:19` `pv = pv.ravel()` | 传入的是 pandas Series，**Series 没有 `ravel`** → `AttributeError` | 绕开 `run`，直接用 `base.dyn_de` + `base.get_mll_results`，多重检验校正改用本仓库的 BH |

**两条硬要求：**

1. **垫片必须用已知解析导数的函数对拍**，不能只看"import 成功了"。
   `_central_diff_weights(3,1)` 必须是 `[-1/2, 0, 1/2]`、`(3,2)` 必须是
   `[1,-2,1]`；`sin`/`exp`/`x^3`/`x^4` 的相对误差要 <1e-4。
   **垫片错了只会让结果悄悄偏，而不会报错。**
2. **必须记录本轮是否真的打了垫片**（`scipy_misc_derivative_shimmed`）。
   已经存在时返回 `False`，不无条件声称"我修了"。

**绕开 `run` 的代价要说清楚**：qval 从 Storey q-value 变成 BH。
这反而更好 —— Moran's I 那条路也是 BH，**两条路的校正口径一致才可比**。

**只用 `SpatialDE.base.X`，不能用 `SpatialDE.X`。**
`SpatialDE/__init__.py` 只导出 `dyn_de` / `run` / `model_search` /
`fit_patterns` / `spatial_patterns` 五个名字，`get_l_limits` 与
`get_mll_results` 都在 `base` 里没被提上来。

**代价与取舍**：SpatialDE 每基因要拟合一个高斯过程，实测 4025 个 spot 上
150 个基因约 2.6 分钟。所以只跑「Moran's I 前 100 + 按种子随机抽 50 作背景」。
**子集抽样会引入选择偏差**（前 100 个是 Moran's I 挑的），
所以一致性必须**分 top 组和背景组各报一次** ——
合并成一个数会把"两边都认为强"和"两边都认为弱"平均掉。
**背景组的 rho 才是有信息量的那个数**（实测 `0.7983` ~ `0.7984`，
top 组 `0.8713`）。

> **背景组那个数改过：** 早先写的是 `0.7954`，近三轮 CI 是
> `0.7984` / `0.7983` / `0.7983`。**top 组三轮逐位相同（0.8713），
> 背景组在第 4 位小数上漂** —— 因为背景组是"按种子随机抽 50 个"，
> 抽到的集合对 Moran's I 排序的边缘很敏感。
> **引用前看 `svg_results.csv`，不要从文档抄**（规则 21.3）。


## 20. `install_requires` 里有某个包，不等于运行时会 import 它

`SpaGCN` 上一轮被登记成"装不上"，理由是：

> PyPI 有真包（1.2.7），但依赖 `louvain` —— 该包最新版 0.8.2 没有
> py3.12 wheel，只有 sdist，要从 2019 年的 C++/Cython 源码编译。

**前半句是实测的，后半句的推论是错的。** 读 SpaGCN 1.2.7 的源码：

- `SpaGCN/SpaGCN.py`、以及 SpaGCN 包内的 `models.py`、`util.py`
  里**没有一处** `import louvain`；
- 它走的是 `scanpy.tl.louvain`（包内 `models.py:69`、`util.py:272`）；
- `simple_GC_DEC.fit` 的 `init` 参数有 **`"kmeans"` 分支**
  （包内 `models.py:52-61`），完全不碰 louvain。

所以 `pip install --no-deps SpaGCN==1.2.7` + `init="kmeans"` 绕开了整条
编译链，SpaGCN 现在**真的在跑**。代价（已写进产物的 `limitations`）：
簇心初始化从"表达+空间"变成"只用 GCN 特征"，`n_clusters` 因此必须
外部给定 —— 传内置方法的域数，两边域数相同 ARI 才可比。

**规则：判断一个包能不能用，判据是源码里的 import，不是元数据。**
元数据只说明"pip 会不会去装它"。这一条和规则 19（SpatialDE 垫片）
是同一件事的两面：**装得上但导不进来**（SpatialDE）与
**元数据说装不上但根本不需要**（SpaGCN），都不能靠读 PyPI 页面判断。

反过来也成立：**上游 `setup.py` 里没写的依赖不代表不需要。**
`STAGATE_pyG` 的 `install_requires = ["requests"]`，而它的
`gat_conv.py:10` 是模块级 `from torch_sparse import SparseTensor, set_diag`。
**两个方向都会骗人。**

### 落地登记要能表达"缺数据"这一种

`NAMED_TOOLS` 的 `kind` 原来有四类（`r_package` / `not_on_pypi` /
`deps` / `name_taken`），这一轮加了第五类 **`needs_reference`**：

| kind | 含义 | 例子 |
|---|---|---|
| `needs_reference` | **包装得上，缺的是数据** | `cell2location`（§3.3） |

`cell2location` 0.1.5 的依赖（scvi-tools + torch + pyro-ppl + opencv-python）
都装得上，卡住的是 `Cell2location(...)` 需要 `cell_state_df` ——
每个细胞类型的全转录组后验表达谱，而默认配置是 marker 签名。

**把"缺数据"写成"装不上"会让下一个人去折腾安装，方向完全错。**
所以 `probe_named_tools` 对 `needs_reference` 这类**不报"登记过期"** ——
它可 import 是符合预期的。

### 点名工具跑了，就必须量化它与主方法的一致性

**"跑通了"不是结论。** 一个点名工具跑出 13 个域，如果不说它和内置划分
的 ARI，读者只能看到一个孤立的划分。所以 `main_analysis.py` 里有一条
honesty 检查：`domain_methods` 里任何 `used: True` 的工具，
`method_agreement` 里必须有对应的 `<tool>_vs_builtin`。

ARI 是主判据（对域编号置换免疫）；`same_label_frac` 只作参考 ——
它会被编号顺序完全支配。**ARI 高也不等于两套方法都对**：
它们可能共享同一个错误（比如被同一个技术批次效应驱动）。

### 登记位置不唯一，验收要把两张表合并起来看

§3.2 的 SpaGCN 跑起来之后从 `named_tools` 搬到了 `domain_methods`
（它不再是"用不了的工具"）。所以 `named_tools:*` 这条验收检查必须
**把 `named_tools` 与 `domain_methods` 合并**再比对点名的工具名单 ——
只认 `named_tools` 会把"已经跑了的工具"报成"缺登记"，
**正好把好事判成坏事**。

### 按步检查查不出"整段没登记"——要再做一次全局清点

上面那条是**按步**查的：每个状态文件声明自己负责哪几个工具，检查它们在不在。
它查不出"**没有任何一步声明过某个工具**"。

实测踩到：`Bering` 与 `BOMS` 的 `section` 是 **§3.1**，而 §3.1 对应的步骤
（`01_qc` / `02_normalize`）**没有 `named_tools` 字段** —— 于是这两个工具
只活在 `common.NAMED_TOOLS` 这个 Python 常量里，**任何产物里都看不到它们**。
按步检查全部通过。

**这是"缺口没登记"最隐蔽的形态：不是登记错了，是根本没登记，
而验收看起来是绿的。**

所以 `main_analysis.py` 另有三条：

| 检查 | 判据 |
|---|---|
| `honesty:named_tools_registry` | `NAMED_TOOLS` 每条都有 `kind` 与 `reason` |
| `honesty:named_tools_never_probed` | 清点哪些工具**没有任何一步探测过**（可见，不判失败）|
| （落盘） | 全表探测结果写进清单 `params.named_tools_probe` |

第三条不能省：前两条只是 CI 日志里的两行，**日志会滚掉，清单不会**。
不加 `only=` 过滤 —— 全局清点过滤掉谁都会漏。

---

## 21. `set_seed()` 不 seed torch，而第三方工具的随机性可能在 torch 里

姊妹项目 `scrna-pipeline-skill/AGENTS.md` 规则 20 记了同一类问题的另一面
（那边是 `pynndescent` 的 Numba 并行）。空间这边实测踩到的是 **torch**。

**证据（五轮 CI，同一份代码、同一批包版本，只改并行度/种子设置）：**

| run | 并行度 / 种子设置 | SpaGCN vs 内置 ARI | NMI |
|---|---|---|---|
| 35488906157 | 无钉 | 0.3734 | 0.5535 |
| 35489064432 | +OMP/OPENBLAS/MKL + CORETYPE | 0.3784 | 0.5557 |
| 35489172104 | 同上 | 0.3639 | 0.5310 |
| 35489534090 | +`NUMBA_NUM_THREADS=1` | 0.4216 | 0.5762 |
| 35489962167 | +`random`/`numpy`/`torch` 三个全局 RNG | **0.3708** | 0.5451 |

**五轮五个值。三次归因、三次被否证 —— 变量不在这里。**

而同一五轮里，**主方法的数全部逐位相同**：Moran's I（空间平滑后）
`0.7670`、域数 `13`、不平滑邻居同域率 `0.507`、平滑后 `0.668`。

**"某个量在变"不等于"所有量都在变"** —— 先看哪些**没**变，
范围一下就缩小了。这一步比"设了种子"有用得多。

### 21.1 三次归因都被否证了

| # | 归因 | 依据 | 否证 |
|---|---|---|---|
| 1 | 多线程 BLAS 归约顺序 | geo 规则 12 的经验 | 钉住 OMP/OPENBLAS/MKL + CORETYPE 后仍从 0.3784 变 0.3639 |
| 2 | Numba `prange` 线程数 | SpaGCN 走 numpy，怀疑并行归约 | 补 `NUMBA_NUM_THREADS=1` 后仍给 0.4216 |
| 3 | **torch 未 seed** | `set_seed()` 确实没 seed torch（读源码确认） | 钉住三个全局 RNG 后仍给 **0.3708** |

**第 3 条的依据是真的**（`set_seed()` 确实只 seed 了 `random` 与
`numpy`），**但"依据为真"不等于"它就是变量"。** 三次都是"读源码
看着很有道理"就动手，三次都被日志否证。

> **残留随机源尚未定位。** 第四次动手之前先读日志（规则 13）。
> 上面那五轮 ARI 摆在一起，说明问题不在"哪个 RNG 没设"这个层面 ——
> 继续加环境变量只是碰运气。

### 21.2 已确认的源码事实（不是推测）

```python
def set_seed(cfg):
    random.seed(seed)
    np.random.seed(seed)      # <-- 没有 torch.manual_seed(seed)
    apply_style(cfg)
```

SpaGCN 1.2.7 的两处随机源：

1. `models.py:55` `KMeans(self.n_clusters, n_init=20)` —— **没有
   `random_state`**，走全局 numpy 遗留 RNG；
2. `train()` 训练 GCN 走 **torch**（权重初始化 + dropout）。

**SpaGCN 自己知道要 seed。** `util.search_res()` 与
`ez_mode.detect_spatial_domains_ez_mode()` 开头都有

```python
random.seed(r_seed); torch.manual_seed(t_seed); np.random.seed(n_seed)
```

**但本仓库走的是 `init="kmeans"` + 外部给定 `n_clusters` 那条路**
（规则 20 为了绕开 louvain 的 py3.12 编译链选的），两个函数都不经过 ——
**它的 seeding 全部被跳过。**

### 21.3 处理

在 `clf.train()` 之前照抄 SpaGCN 自己的做法把三个全局 RNG 钉死：

```python
random.seed(seed); np.random.seed(seed)
torch.manual_seed(seed); torch.set_num_threads(1)
```

并把结果记进 `domain_methods.SpaGCN.seeded_before_train` ——
**没有 torch 时要如实记 `torch_seeded: False`，不能假装 seed 成功。**

**留着它**：成本为零、方向正确、`seeded_before_train` 可核对。
**但不要以为设了就可复现** —— 实测不足以让 ARI 稳定。

### 21.4 报数

**ARI 是范围，不是定值。** `domain_status.json` 的 `reproducibility`
字段写明哪些量能按定值报（Moran's I、域数、邻居同域率）、
哪些必须带范围报（`SpaGCN_vs_builtin` 的 ARI，实测
`0.3639` ~ `0.4216`），以及**被否证的三条假设**。

**别把数值写死在 README 里。** 实测 README 里"不平滑 9 域、同域率 0.536、
基线 0.133"早就对不上了（现在是 11 域 / 0.507 / 0.111）——
**文档里写死的数一定会过时**，而读者不会知道它过时了。

### 21.5 诊断字段必须扛得住日志截断

排查这个 ARI 时卡了很久，一部分原因是**日志里根本看不到需要的字段**：
`汇总产物` 那一步原来打的是 `json.dumps(v)[:500]`，而 `domain_methods`
一长就从中间被切断 —— `seeded_before_train` / `seed` 恰好落在第 500 个
字符之后。**只能靠猜，而我已经猜错三次了。**

现在改成**按字段名挑**：短标记与数值逐条打全，长散文（`reason` /
`note` / `why_kmeans`）才截断，并把 `reproducibility` 一起打出来。

> **规则：日志的截断位置不该由"字典有多长"决定，该由"哪些字段能定位
> 问题"决定。** 一个 `[:500]` 会让下一轮排查同样卡住 ——
> 而 `results/` 的 artifact 有 14 天保留期，日志滚得更快。

---

## 22. 文档改动走独立的 `docs_check.yml`

`spatial_analysis.yml` 有 `paths:` 过滤（只跑 `scripts/` `tools/` `assets/`
`requirements.txt`），**`*.md` 的改动不触发它** —— 所以文档里的死链接
以前**没有任何门禁能挡住**，CI 每次都是绿的。

实测就踩到了：AGENTS 与 `references/module0.md` 让读者去看 `models.py` /
`util.py` / `setup.py`，**但没说那是哪个包的** —— 它们是 SpaGCN 与
STAGATE_pyG 包内的文件，不在本仓库里。

现在 `tools/check_doc_refs.mjs` + `.github/workflows/docs_check.yml` 兜住这一类，
**约 20 秒**、不装依赖、不跑分析。两个逃生舱写在工具文件头：

| 逃生舱 | 标记 | 为什么必须写 |
|---|---|---|
| 指向**已删**的文件 | 已删 / 已移除 / 不再存在 / 曾经 / 当时的 / deleted | 读者要能区分"历史"和"笔误" |
| 指向**第三方包**源码 | 包内 / 上游 / 源码 / site-packages / 该包 | 读者要知道去哪个包里翻 |

**两个逃生舱都故意做成"要写一句话"的** —— 静默豁免会让检查退化成没有检查。

> **为什么文档不并进主 workflow：** 那样改一个错别字要跑完整的空间流水线
> （约 25 分钟），而且会被主流水线的偶发失败牵连（规则 21 的 ARI 五轮五个值
> 就是实测记录）—— 文档改动因无关原因判红，反而让"每次推送 CI 必须绿"失效。

## 23. 出图名要带代码坐标

参考规范：geo 仓库 AGENTS 规则 30（同一套约定，阶段号换成 `03`）。

**命名格式五个字段，用 `-` 连起来：**

```
<阶段>-<模块>-<图>-unit<单元>-<名称>
 03      03     03   unit1     domains-on-he
```

| 字段 | 取值 | 来源 |
|---|---|---|
| 阶段 | `03` | 空间转录组 = 第三部分（geo `01` / scrna `02`）|
| 模块 | 两位数字 | **脚本文件名的前两位**（`04_svg.py` → `04`）|
| 图 | 两位数字，从 `01` 起 | 该脚本内第几张图 |
| 单元 | `unit1` 起 | 同一张图里的功能单元 |
| 名称 | 小写连字符 slug | 图的内容 |

**为什么写成门禁而不是靠人记：** 图名和脚本序号是**两处**，而"图名里的
模块号写错了"没有任何东西能发现 —— 图照样生成、CI 照样绿、验收照样过，
只是读者按图名去 `scripts/` 里找代码时会**找错文件**。
`tools/check_fig_names.mjs` 查四条：格式合规、**模块号与脚本文件名一致**、
图号/单元号连续、全仓库无重名。它**只认字符串字面量**；经辅助函数传名的
调用会逐条列出（可见）但不判失败。

**图名改了，引用它的地方也要改。** 本仓库 `run_acceptance` 的图检查是
glob 整个 `figures/` 数张数、不引用具体文件名（规则 12 的中间态排除
同一条思路），所以验收不会因为改名而静默变 false —— 但 README / SKILL
里的产物清单会。**`check_doc_refs.mjs` 管不到这一类**：它查的是
源码/配置文件名，图名是**运行产物**、本地不存在。改图名时手工 grep
一遍 `\.png` / `\.pdf`。

> 三个仓库的 `check_fig_names.mjs` 是**同一份文件**（geo 那份只多识别
> R 侧 `save_pdf` 的调用写法）。阶段号表写在文件头的 `PART_BY_REPO`，
> 仓库目录名认不出时它直接报错退出，不会静默放行。

## 24. 定量面板必须有可读的色标；颜色有含义就必须有图例

两条出图层约定，都是出图质量评审抓出来的（评审报告存档在工作区的
治理层 reviews/ 目录，不在本仓库内）：

**24.1 多面板共享尺度的定量图，色标必须存在且共享。**
实测 `svg_top_genes`：30 个面板逐基因各自 min-max 归一化、零 colorbar ——
每个面板自己的"最暗"都是 0、"最亮"都是自己最大值，**面板间完全不可比、
绝对值读不出**。修法：`vmin=0` 固定（log1p 后 0 = 不表达），`vmax` 用
top 基因的**全局 p99**（绘图分位数、不是统计阈值），逐面板显式传
`norm=Normalize(vmin, vmax)`，整图一个共享 colorbar（ScalarMappable +
`ax=list(axes)`，constrained layout 自动让位），suptitle 写明色标范围。

**24.2 条形/散点颜色编码了阈值或类别，就必须有图例。**
实测 `communication_lr_enrichment`：红（z>2）/灰/蓝（z<−2）三色条无任何说明，
z=−2 蓝虚线画在没有数据的左侧空白处，读者无从知道它是什么。修法：
Patch/Line2D 显式图例（含两条阈值线），title 里不再用 "`|z|>2 dashed`"
这种没有图例时的凑合话。

**顺带两条同源教训（同一次修复实测）：**

- **figsize 公式的单位必须一致。** `max(mm(56), 0.32*len+1.6)` 把英寸当毫米比，
  20 条目给 8.0"=203mm —— 图幅超限不是设计出来的，是单位混用算出来的。
  `mm()` 的结果只能和毫米比较，或换算成英寸后再夹。
- **colorbar 会挤压面板标题。** constrained layout 给 colorbar 分的空间
  会把多面板图的标题压劈/截断（实测 `spatial_pseudotime_maps` 三面板
  逐个受损）。colorbar 显式给 `fraction≈0.046, pad≈0.02, shrink≈0.8`，
  标题压成短行、统计量移到第二行。

## 25. `PAL` 的键名不能靠记；多面板别给每个 colorbar 挂长标签

**25.1 `PAL["xxx"]` 的键必须真存在 —— 本地门禁已挡住这一类。**
实测（2026-09-20）：给 `svg_stat_distribution` 加参考线时写了
`PAL["up"]` —— 那是 **geo（R 侧）**的语义键；本仓库 PAL 只有
`highlight` / `primary` / `muted` / 定性色。本地静态检查全绿（`PAL`
这个**名字**确实 import 了），CI 跑到 svg 步骤才 `KeyError: 'up'`、整步判红。
`PAL` 键是**本仓库自己的字典，本地完全查得到**，所以
`check_py_syntax.mjs` 现在会从 `common.py` 解析 PAL 的字面量键集合再扫
所有脚本，键不存在就本地判红（scrna 侧同款检查器，两份逐字节相同）。

> **实现上踩到的坑（值得记）：** 第一版在 `stripComments()` 的结果上扫，
> 而它把字符串字面量换成 `""` —— `PAL["up"]` 变成 `PAL[""]`，
> **检查器把要检查的东西本身擦掉了**，负向验证时照样报"通过"。
> 现在改为在**原始源码**上正则。与 `check_r_syntax` 当年
> "stripLiterals 擦掉隐式拼接"同一个坑：**新写检查器时，必须用
> "故意塞一个错"验证它真的会响**（本条已做双向验证：正向绿、错键红）。

**25.2 多面板共享量程时，别给每条 colorbar 挂长 label。** 实测
`deconvolution_spatial` 给 12 个面板各挂 `"<类型名> (shared scale)"`，
长类型名（Fibroblastic_reticular_cell）的竖排文字挤进相邻面板 ——
修共享量程却引入新重叠。共享量程由 suptitle 统一说明，刻度数字足够。

