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

> **`z_margin` 空着不等于「不确定」（E-68，2026-09-26）。**
> `domain_annotation.json` 的每个域现在带一个 `margin_state` 四态：
> `ok`（差距 > 0.5，可信）/ `low_margin`（有第二名但差距 <= 0.5，
> 两条路真的分不开，去比对两条注释路）/ `single_celltype`
> （**候选细胞类型只有一个，没有第二名可比**，`z_margin` 是 `null`）/
> `margin_undefined`（数值算不出来）。后两者合并计数为
> `n_margin_undefined`，与 `n_low_margin` **分开报** ——
> 旧写法在只有一个候选类型时把第二名取成第一名自己，`z_margin` 恒为
> `0.0` ⇒ `0.0 > 0.5` 为假 ⇒ 记成 `assignment_confident: false`，
> 于是日志说「N/M 个域的标签 z_margin <= 0.5 —— 这些标签不该被当结论」，
> **而其中可能一个域的 margin 都没算出来**。两者排查方向相反：
> 前者去比对两条路，后者去**补签名基因**。
> `assignment_confident` 因此是三态（`true` / `false` / `null`），
> 不要再当布尔量用。验收层 `content:domain_annotation` 会把这个分布
> 打进 `acceptance.json` —— 状态写出来不算数，**有人读**才算（同 L6 第二半）。

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
  输出尺寸**，让上面的毫米约定失效。溢出改由 `_content_overflow()` 检测并
  **落盘成 `figure_overflow.json`**（供验收判红，见规则 28.3）。
  门禁层的另一半是 `check_figures.mjs` 的内容贴边检查（规则 27）——
  **两者不能互相替代**：贴边是**事后**发现，落盘 + 验收判红是**当轮**发现。
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

### 17.1 写死的名单只覆盖了它自己 —— 现在是真正的逐作用域分析（E-61）

上面那份名单**只有 10 个名字**，而 `common.py` 实际导出 **68 个** —— 挡住的全是
"恰好被列进去的"，没列进去的照旧漏到 CI。实测 `03_spatial_domains.py` 漏 import
`spot_radius_plot_units`，本地报"全部通过"，CI 跑 25 分钟到 H&E 叠图段才
`NameError`。**根因不是"名单短了一点"，而是判据的输入域与它要防的缺陷不匹配。**

现在 `tools/check_py_names.py` 用标准库 `symtable` 做逐作用域分析（判四件事，
见上表）。**三个实现坑**（不安全构造必须用 AST 判、行号必须按作用域定位、
提示要扣掉 common 自己 import 的名字）、geo 判"**不适用**不是失败"、
`PY_GATES` 接线与 13 项标定：`references/agents-detail.md#原规则-171`。

### 17.2 手工删 import 会连带删掉同一行的活名字

E-61 的直接触发动作是**手工清理死 import**：`unused_imports.py` 正确报出
死名字 `plot_marker_dotplot`，但我在执行删除时把同一 import 行里相邻的
`spot_radius_plot_units`（它**有**调用点，不在死名单里）一起删了。

**扫描器没错，是执行删除这一步错了。** 所以：删 import 时逐名核对
"这个名字在文件里还有引用吗"，不要按行删。死名字扫描器的输出是**名单**，
不是**待删行号**。

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
| 1 | SpatialDE **包内** `base.py:12` `from scipy.misc import derivative` | scipy 1.12 移除了 `scipy.misc.derivative` → `ImportError`，**装得上但导不进来** | 垫片补回该名字（3 点中心差分，Vandermonde 解权重）|
| 2 | SpatialDE **包内** `base.py:432` → `util.py:19` `pv = pv.ravel()` | 传入的是 pandas Series，**Series 没有 `ravel`** → `AttributeError` | 绕开 `run`，直接用 `base.dyn_de` + `base.get_mll_results`，多重检验校正改用本仓库的 BH |

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

**规则：** 判一个包"装不上"之前，**读它的源码看运行路径上到底 import 谁** —— 依赖表里有的包可能一次都没用到。`SpaGCN` 曾被登记成"装不上"（PyPI 有真包 1.2.7，但依赖 `louvain`，该包没有 py3.12 wheel）；**前半句是实测的，后半句的推论是错的** —— 读 SpaGCN 1.2.7 源码：包内 `SpaGCN.py` / `models.py` / `util.py` **没有一处** `import louvain`（走的是 `scanpy.tl.louvain`），且 `simple_GC_DEC.fit` 的 `init` 参数有 `"kmeans"` 分支可绕开。**"实测的前半句 + 未验证的推论"合成一条读起来完全合理的结论，是这类记录最危险的形态。**

完整原文（含源码证据与排查顺序）：`references/agents-detail.md#原规则-20`

## 21. `set_seed()` 不 seed torch，而第三方工具的随机性可能在 torch 里

**规则：** **"我设了种子"不等于"结果可复现"** —— `set_seed()` 覆盖不到第三方库内部的迭代求解器、并行归约与 GPU 内核。实测：同一份代码、同一批包版本，**只改并行度/种子设置**，SpaGCN 与内置方法的 ARI/NMI 就逐轮漂（五轮 CI 对比表见原文）。空间这边踩到的是 **torch**（姊妹项目 `scrna-pipeline-skill/AGENTS.md` 规则 20 记的是同一类问题的另一面：`pynndescent` 的 Numba 并行）。并行度有三个独立旋钮（BLAS 线程数 / `OPENBLAS_CORETYPE` / `NUMBA_NUM_THREADS`），**设它们是对的，但设了不保证可复现**。

完整原文（含五轮 CI 对比表与三个旋钮）：`references/agents-detail.md#原规则-21`

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


## 26. 图例一律图框外右侧、纵向排列（用户约定 v2，2026-09-23）

**图例不能画在图框（panel）里面，也不能放在顶部** —— 框内会压住数据点，
顶部会把主图压扁变形（实测校准图 / UMAP / 去卷积条形图被压得很扁）。

- **matplotlib**：`fig.legend(loc="outside right center", ncol=1)`。
  `loc="outside ..."` **只对 `fig.legend()` 有效**，传给 `ax.legend()` 会报
  `ValueError: 'outside' option ... only works for figure legends`（实测
  run 35749568552 因此崩了整个 job）。
  所有 axes 级图例必须改成 `fig.legend(...)`。
- **`ncol=1` 强制纵向单列** —— 多图例时默认可能横排，必须显式指定。
- constrained layout 会自动为框外图例让出空间；**改完要亲读**确认没被裁掉
  （`savefig.bbox: standard` 下溢出是静默裁，见规则 13 同类问题）。

**门禁**：`node tools/check_legend_convention.mjs`
（三仓同一份，静态扫源码；CI 在"静态检查（不装依赖）"那一步跑）。
它抓两类**确定错**：`fig.legend` 缺 `loc="outside..."` 或缺 `ncol=1`、
`ax.legend` 用了 `loc="outside..."`。

> **为什么需要门禁而不是靠记：** `check_py_syntax.mjs` 只查未定义名字，
> 看不见图例位置。一张图例压在数据点上的图，在
> "文件存在 / 有墨迹 / 图名合规 / 配色合规 / 图幅合规"眼里**全都是合格的** ——
> 这正是工作区治理层错误台账（`governance/15_ERROR_LEDGER.md`，**不在本仓库内**）E-06「门禁本身有盲区」的又一例。


## 27. 门禁要查"内容贴边"，不能只查"有没有墨"（Q-26，2026-09-25）

`tools/check_figures.mjs` 原来只查两件事：**有没有墨**（`frac >= MIN_INK`）、
**是不是糊死**（颜色数）。

**被裁掉的图照样有墨** —— 这正是它当年漏掉姊妹项目那张超宽标题图的
原因（台账 E-49：`_content_overflow()` 检测到 `width_overflow_frac=0.3523`
并打了 WARN，但没人读；`savefig.bbox: standard` 下标题超出画布的部分
被静默切掉，文件照常生成、门禁照样报"有墨"）。

裁切的物理后果是**墨迹延伸到画布边缘**，所以现在量**非背景像素外接框
到左右边的距离**：

```js
const EDGE_MIN_PX = 3;
const edgeBad = !darkBg && margins &&
                (margins.left < EDGE_MIN_PX || margins.right < EDGE_MIN_PX);
```

命中时打 `[贴边]` 并计入失败，OK 行尾追加实际边距（`边距 L/R 13/10`）
—— 边距是**可读的量**，不是"通过/不通过"一个比特。

**阈值 3 px 是实测标定的**（scrna 34 张 + spatial 76 张）：

- scrna 边距分布 `{0:1, 10:9, 11:2, 12:6, 13:13, 14:2, 41:1}` ——
  `L<3 或 R<3`、`L==0 或 R==0`、`L<6 或 R<6` **都只命中那一张**，零误伤。
- **本仓库 76 张全部通过，零误伤。**

两条刻意的取舍：

1. **只判左右，不判上下。** 实测姊妹项目
   `02-05-05-unit1-pseudotime-by-cluster` 上边距 = 0 —— 那是布局取舍，
   不是裁切。单看某一边会误伤。
2. **暗底图跳过。** 整幅都是"墨"，外接框必满幅，判了全是假阳性。

> **本文件与 `scrna-pipeline-skill/tools/check_figures.mjs` 必须逐字相同**
> （两仓同一份，见该文件头注释）；改一侧必须同步另一侧并比对哈希。
> 当前两侧 SHA256 已比对一致。
>
> **另一半是验收层的落盘 + 判红**（规则 28.3）。贴边门禁是**事后**发现，
> 落盘 + 验收判红是**当轮**发现，两者不能互相替代 —— 一张在 CI 上被
> 静默裁掉的图，要等 artifact 下载下来跑一遍门禁才被发现，而那时
> 这一轮早就判绿了。


## 28. 验收层四条补强（Q-27，2026-09-25）

姊妹项目 `scrna-pipeline-skill/AGENTS.md` 规则 25 是同一批补强的另一半。
**本仓库的口径与它不同**（规则 23：图检查是 glob 数张数、不引用具体文件名），
所以**没有照抄**，而是按空间侧结构重新设计：

- **28.1 步骤函数自己返回的 `status` 必须有人读。** `run_steps` 把 `fn(cfg)` 的
  返回值记进 `results[name]["result_status"]`，**而全仓没有一个消费者** ——
  于是 `08_spatial_trajectory.py` 返回 `bad_root` / `not_applicable` /
  `missing_pca` 时验收照样记 `status="ok"`。**"步骤崩了"与"步骤跑完了但结果是
  没做成"是两件事**，后者最容易读成成功。现在 `step_result:<sid>` 消费它，
  只把 `failed`/`error`/`fail` 判红。
  > **坑：`chk(cid, kind, ok, detail, severity="required")` 的第二个位置参数是
  > `kind`（归类），不是 severity。** 传错**不报错**，只是把"只可见"记成
  > "必须通过"。**同名同型同位置数相邻，传错不报错。**
- **28.2 嵌套 `status` 的 `failed` 必须判红。** `domain_methods.STAGATE.status` /
  `cell2location.status` / `spatialde.status` 都是嵌套的，旧验收层只看顶层 ——
  **里面崩了、顶层还是 `ok`**。实现 `_iter_nested_status(obj, path=())` 递归扫
  **全部** `*status.json`，只判 `failed`/`error`/`fail`（其余取值判红会让每个 job
  都红，反而没人看）。本仓实测嵌套分布 `{ok: 11, needs_reference: 1,
  package_missing: 1}` —— 零命中、零误伤。
- **28.3 溢出检测必须落盘，不能只打 WARN。** E-49 的形态：`_content_overflow()`
  **正确检测到** `width_overflow_frac=0.3523`、**正确打了 WARN**，然后**没有任何
  人读**。**"检测到了"不等于"有人会知道"。** 现在 `_record_figure_overflow()`
  落盘到 `results/<dataset_id>/figure_overflow.json`，验收层读它并**判红**。
  三条约束：**跑前必须删旧文件**（否则图修好了旧记录还在）、**只累积不覆盖**、
  记录可被修复。
- **28.4 图验收从"数张数"换成"声明过的每张在不在"。** 旧判据 `len(figs) >= 8`
  是**纯计数**，于是"该有的图没有"这一整类没有任何检查看得见。声明**从源码扫
  出来**（`declared_figures()`）不手抄表（手抄的表会漂，漂移方向恰好是"新加的
  图不在表里"），扫时要**排除 `main_analysis.py` 自己**。三条判据分工：
  `figures:declared` 逐名核存在性 / `figures:dynamic` 声明了槽位的动态图名每组
  **至少 1 张** / `figures:count >= 8` 保留为**下限兜底**。
  > **"判据因为输入域重叠而永远为真"的坑，正向标定看不出来：** `03-03-04` 前缀
  > **同时**有静态图（unit1/2/3）与动态图（unit8 Jaccard）。按前缀计数而不扣掉
  > 静态声明图，动态循环整段没跑时计数仍是 2（全是静态的），**判据永远不响**。
- **28.5 双向标定：正向零误报 + 反向逐类注入，缺一不可。** 正向（真实 artifact
  干跑）只说明"没误报"；**一个永远返回 True 的判据在干净 artifact 上也是零命中**。
  逐类注入实测：删 1 张静态声明图 → 只 `declared` 红；删整组动态图 → 只
  `dynamic` 红；写 `figure_overflow.json` → 只 `overflow` 红；嵌套
  `STAGATE.status=failed` → 只 `status:nested` 红；删全部图 → 三条都红。
  正向：76 张图干跑，新增 4 条检查全 PASS，`n_checks` 59 → 63，**没有一条既有
  检查消失**。（干跑把所有写盘函数 patch 成 no-op —— 否则会重写 artifact 里的
  `acceptance_report.json`，变成"自己改自己的证据"。）

完整原文（含注入表与全部实测证据）：`references/agents-detail.md#原规则-28`

## 29. 叠 H&E 的图必须 `set_aspect("equal")`，且**不能**跟着别的脚本 `invert_yaxis()`

**这条是补上的规则缺口**：`03_spatial_domains.py` 的三处 H&E 叠图从仓库建立起
就没有 `set_aspect("equal")`（`git log -S` 零命中），而
`assets/publication.mplstyle:108` 是 `image.aspect: auto`。后果是
1921x2000 的 hires 图被铺进 136mm x 56mm 的 axes —— **横向拉伸 2.585 倍**，
组织被拉扁、**fiducial 点阵的圆点被画成横椭圆**（实测渲染产物里 w/h 中位 3.000，
源图里是圆）。

### 29.1 为什么不能照抄同仓其它脚本

全仓 `set_aspect("equal"); ax.invert_yaxis()` 出现在 `01_qc.py` / `02_normalize.py` /
`04_svg.py` / `05_deconvolution.py` / `07_spatial_communication.py` ——
**这些脚本都没有 `imshow`**（纯 scatter，坐标是 y 向下的原始像素，所以要 invert）。
`ax.imshow(img)` **默认就已经把源图朝向显示对了**：实测 `ax.yaxis_inverted() == True`，
按网格点与源图灰度做相关，同点 **0.9894**，上下翻转 0.7251、左右翻转 0.7239；
**显式 `invert_yaxis()` 之后同向相关反而掉到 0.7199**。

> **一句话**：`imshow` 的 y 轴本来就是反的，再 invert 一次就把图上下翻过来了。
> 叠图段只加 `set_aspect("equal")`，**不要**加 `invert_yaxis()`。

### 29.2 图幅要跟着宽高比走，不能写死

等比例之后 `136mm x 56mm` 只剩 31% 利用率。正确做法是让画布高由源图宽高比推出 ——
`lib/common.py:fit_fig_to_aspect(fig, ax, width_mm, aspect)` 用
`fig.get_tightbbox()` 解方程（迭代 3 次、容差 0.01mm，第 2 轮即收敛）。

**这里有一个量纲坑**：`get_tightbbox()` 返回**英寸**，`ax.get_window_extent()` 返回**像素**，
直接相减得 `deco_w = -120.3 mm` —— **负数不报错，只给出荒谬的画布尺寸**。
必须先乘 `fig.dpi` 再减。

**还有一个"在方画布上量边距"的坑**：不能用 `fig_h - bb_h` 估装饰边距 ——
等比例没占满的那部分长宽比余量会被当成装饰，结果把画布锁死在方形
（合成横长/竖长用例全败，fig 恒 136x136）。必须用 tightbbox 解方程。

实测（修后）：真实图 136mm → fig 136x133.2、axes w/h **0.960500 逐位等于源图**、
利用率 **85.0%**；无图例 94.0%；合成横长/方形/竖长/极宽/极高 6 种用例
83.1%~95.9%，全部无溢出。

### 29.3 spot 大小要由尺度因子推，不要写魔法数

硬编码 `s=8` 与源图无关。用 `spot_radius_plot_units(scalefactors, "hires")`
拿到**数据单位**半径，再用 `lib/common.py:marker_area_pt2(fig, ax, radius_data)`
换算成 `scatter(s=)` 的面积（内部 `ax.transData.transform` 量 2r 个数据单位的
显示长度，乘 `72/fig.dpi` 得点数）。

**换算必须在 `fig.canvas.draw()` 之后做**（constrained layout 定稿后
transform 才是最终的）；`ax.get_window_extent()` 在 `savefig` 前后不变（实测 True），
所以可以在画图前一次算好。

> `references/troubleshooting.md:98` 早就点名了这个修法，但 03 一直没照做 ——
> **文档写了正确做法、代码做的是相反的事**，是 E-58 那一类缺陷的同族。

### 29.4 门禁覆盖

`check_legend_convention.mjs` 只能看图例；等比例**没有自动门禁**（像素级形变
需要渲染后比对，成本高）。所以本规则靠**人读图 + 像素级探针**守：
`D:\tmp\_r04\verify_he_fix.py` 会把渲染结果里 fiducial 暗斑的 w/h 量出来，
**圆的就该是 1.000，旧产物是 3.000**。

## 30. 门禁要带内建自检，且自检本身必须被反向标定（E-62 / E-63，2026-09-26）

**没有自检的门禁只能证明"它没报错"，不能证明"它检查了"。** 两条实测：

| 台账 | 门禁 | 缺陷形态 |
|---|---|---|
| E-62 | `check_legend_convention.mjs` | Python 侧没抹注释 → 注释里一个 `fig.legend(` 让括号配平**一路吞到文件尾**，其后所有真调用一个都没查，门禁照样打绿 |
| E-63 | `check_figures.mjs` | WARN 落盘分支从落地起**一次都没执行过**，里面有两个必崩的错（`INK_FAIL_MIN` 未定义、报告路径用了循环变量）|

两条的共同点：**假阴性**。门禁的失败方式不是"报错"，而是"什么都不报" ——
而"什么都没发现"与"检查通过了"在输出上完全一样。**假阴性比假阳性危险得多**：
假阳性会被人骂着修掉，假阴性会被当成绿。

### 30.1 自检要调真代码，不能自己重写一遍逻辑

E-63 的自检第一版有 9 个用例，**用例 8 自己另写了一遍路径拼接**
（`join(dirname(join(dir, "figures")), "warn_report.json")`），没调真代码。
反向标定把缺陷注回去（`outPath` 改回 `join(dir, ...)`）→ **自检仍然通过**。

**抽函数**才解决：`checkDirs()` / `buildWarnReport()` / `writeWarnReport()`
三个纯函数，`main()` 与自检**都调它们**。抽完再注一次缺陷 →
`ReferenceError: dir is not defined`、exit=1。

> 同 E-62 的「检查器要检查的东西，与检查器描述自己要检查什么，在纯文本上
> 无法区分」是同一个坑的两种形态：**自检里重实现一遍被测逻辑，等于没测。**

### 30.2 反向标定：逐个把原缺陷注回去，确认自检真的会红

**正向通过证明不了任何事** —— 一个永远返回 True 的用例在干净产物上也是绿的。

| 注入 | 期望 | 实测 |
|---|---|---|
| `inkFailBlank: MIN_INK` → `INK_FAIL_MIN` | 回归 1/2 红 | ✅ `ReferenceError` |
| `outPath` 改回 `join(dir, ...)` | 回归 2 红 | ✅ 第一版**不红**（假自检）→ 抽函数后红 |
| 尾斜杠处理删掉 | 回归 3 红 | ✅ |

### 30.3 自检必须接进 CI 与 pre-push —— 没人跑的自检是同一类缺陷

写了 `--selftest` 却只在本地手敲，等于又造了一个"从未执行过的分支"。
现在三处都接：

| 位置 | 内容 |
|---|---|
| `.github/workflows/spatial_analysis.yml`（本仓）/ `scrna-pipeline-skill/.github/workflows/scrna_analysis.yml` | 「静态检查（不装依赖）」那一步跑 `check_legend_convention.mjs --selftest` + `check_figures.mjs --selftest` |
| `geo-normal-pipeline-skill/.github/workflows/geo_analysis.yml` | 同一步跑 `check_legend_convention.mjs --selftest`（geo 无 `check_figures.mjs`）|
| `governance/hooks/pre-push.mjs` | 新增 `SELFTESTS` 表，在**所有**静态门禁之后跑，日志标签是 `（自检）` |

> **跨仓引用必须写仓名前缀。** 三个仓的目录结构很像（都有
> `.github/workflows/`），只写文件名（不带仓名）读者不知道去哪个仓找。
> `check_doc_refs.mjs` 的判据 C 就是为此加的 —— 它以前**两条判据都不进**，
> 于是这类引用**完全不被检查**（E-62 / E-63 同一类：检查器"没报错"
> 不等于"检查了"）。

**日志标签必须区分"带镜像目录"与"带 `--selftest`"** —— 两者都走 `extraArgs`，
但一个是拿真实产物判、一个是拿合成用例判，长得一样就没法排查。

**`SELFTESTS` 的接线也做了反向标定**：把 `check_figures.mjs` 自检里
"全白判红"用例的条件改成 `false` → pre-push 输出
`✗ spatial-pipeline-skill tools/check_figures.mjs 未通过`、
`PRE-PUSH 未通过：1 项判红 —— 禁止 push。` —— **说明接线真的会拦，
而不是只在日志里多打一行 `✓`。**

> **接了线但从不失败的检查，与没接线是一样的。** 每加一条自检，都要问
> "我怎样让它红一次"——答不上来就说明它现在是个装饰。

### 30.4 判据的"通过数"必须能看见 0 —— 否则死代码与"没有这类输入"无法区分（E-64）

补判据 C（跨仓引用）时，主体写在 `if (!isA && !isB) continue` **之后** ——
而跨仓 token 正是"两条判据都不进"的那一类，于是**判据 C 的分支永远走不到**。
更糟的是报告那行是 `if (nCheckedC) console.log(...)` 守卫的：计数恒 0 时
**连打印都不打印**，输出看起来与本仓没有跨仓引用**完全一样**。
三仓复跑全绿、exit=0、毫无异常迹象 —— **假阴性顺手把自己藏了起来**。

**两条规则：**

1. **新判据要放在早退分支之前** —— 分流顺序错了分支就是死代码，
   而**死代码不报错、只是永远不执行**。
2. **计数行不能加 `if (n)` 守卫** —— **0 也是信息**。"检查了 0 条"与
   "根本没检查"必须在输出上长得不一样。

**修后计数变化本身就是证据**：geo 102→107 条（C=5）、scrna 113→139 条（C=24）、
spatial 103→117 条（C=14）—— geo/scrna 的跨仓引用**以前一条都没被检查过**。
`--selftest` 11 个用例在 `mkdtempSync` 造的隔离 workspace 里自建两个真姊妹仓
+ 一个**伪装成仓库的非 git 目录**（`releases/`）；其中回归 2 就是"判据 C 计数为 0
时那行照样打印"。

> **自检里"看起来像仓库但不是仓库"的干扰项是必需的** —— 首版把 `releases/`
> 当成姊妹仓，提示给出 `releases/scrna_analysis.yml`，**连路径都是错的**。
> 修法是只认有 `.git` 的目录。**指错地方比不指地方更糟**（同 E-61 防复发④）。

### 30.5 门禁的"检查范围"要和"它守护的动作"对齐（E-65）

`governance/hooks/pre-push.mjs` 的 `changesOf(repo)` 原来只读
`git status --porcelain`（**只含未提交改动**）。而 pre-push 是 `git push`
的钩子，**它唯一被调用的时刻就是"已经提交、还没推送"** —— 那时 porcelain
为空，三仓全走 `无改动，跳过`，**[4/6] 静态门禁段一条都没跑**，
而打印的是 `PRE-PUSH 通过（0 条提醒）。可以 push。`

**这不是边角，是主路径**：正常情况下它每次都在空转，给出虚假的安心。
触发它的是"门禁说通过、而门禁自己（`check_doc_refs.mjs`）说 exit=1"
这个**直接矛盾**。修法是取并集 —— 除未提交改动外，再加
`git log --name-only --pretty=format: @{upstream}..HEAD`（**已提交未推送**）；
没有 upstream 时那一段抛错被吞（"不适用"不是失败）。

> **一个门禁段被整段跳过时，不能打印"通过"。** `无改动，跳过` 用的是
> `ok()`（绿勾），它和"查过了没问题"在输出上一样。写门禁前先跑一次
> "什么都没改"的路径 —— 如果它空转时也说通过，那它有改动时说的通过
> 也不可信。**反向标定的输出是 `exit=0` 而三仓全被跳过** ——
> 它的失败方式不是失败，是**沉默**。
>
> **正确的提交顺序是：改完 → 跑 pre-push → 提交 → push。**

## 31. 「写出来了」不等于「有人读」：三种形态与三个守卫（E-69，2026-09-26）

**规则：一个字段/一条判据的价值不在于它被算出来，而在于有人消费它；
而"算不出来"和"算出来很小"必须在产物里长得不一样。**

E-68 是同族第一次自查（规则 30 那批门禁的连带产物），E-69 是拿它的判据
**回头扫全仓**抓到的第二次 —— 第三次跨仓抓到东西。三种形态：

### 31.1 Form A：`nan` 参与比较会静默变成 `False`

`nan > x` / `nan < x` **都是 `False`，且不报错**。于是"这个量算不出来"
被读成"这个量很小 / 没改善"。三处现场：

| 现场 | 旧写法 | 后果 |
|---|---|---|
| `05_deconvolution.py` 重建误差 | `np.nanmean` / `np.nanpercentile` / `np.nanmax` 在**全 nan** 时返回 nan 并继续算 | `frac_unreliable` 的分母用了 `len(errors)` 而不是有限值个数 |
| `07_spatial_communication.py` z 分数 | `null_sd = float(np.nanstd(perm_means)) or 1e-9` | 全 nan 时 `or` **不触发** → `z = nan` 写进 top5；`denom` 全零 → `100.0/1e-9 = 1e11` 假放大 |
| `08_spatial_trajectory.py` Moran's I | `improved = I_spatial > I_expr` | `nan > nan` 为 `False` → 状态声称"平滑损害了一致性"，而**两个量都没算出来** |

**修法：抽纯函数，让"算不出来"有名字。**

- `summarize_morans_pair(i_expr, i_spatial, ndigits=4) -> {"defined","improved","gain","expression_only","spatially_smoothed","note"}`
- `spatial_z_score(near_mean, perm_means) -> (null_mu, null_sd, z, z_reason)`
- `summarize_reconstruction_error(errors, max_err) -> dict`（含 `error_defined` / `n_unreliable` / `frac_unreliable`）
- `common.finite_round(x, n)` **三态**：nan/inf/None → `None`；有限 → 四舍五入；**`0.0` 要保留，判空用 `is not None`**

**为什么必须抽函数**：抽出来才能被标定脚本**直接调**（139 项纯函数断言）。
内联写法只能靠"跑整条流水线看产物"，而那是 25 分钟一轮。

### 31.2 Form B：只写不读的状态字段

全仓 **15 个字段没有任何消费者** —— 坏值和"根本不存在"在验收层长得一样。
现在 12 条探针注册表（`main_analysis.py` 的
`dict(cid, base, f, path, kind, severity, opt, fn, good, bad)`）：
`status:input_is_counts` / `coords_dropped` / `hvg_flavor` / `full_gene_counts` /
`svg_gene_subset` / `svg_gene_selection` / `svg_genes_dropped` /
`proportions_truncated` / `deconv_matrix_source` / `niche_spot_alignment` /
`smoothing_improves` / `morans_I_defined`。

- `fn` 返回 `None` ⇒ **不适用**（PASS），不是失败。
- `opt=False` = **无条件写**；缺失 ⇒ 判红「产生端不再写了」。
  `opt=True` = 条件写；缺失 ⇒ PASS。

**三个守卫（都是踩出来的）：**

1. **`chk(cid, kind, ok, detail, severity="required")` 的第二个位置参数是
   `kind`，不是 severity。** 传错**不报错**，只是把"只可见"记成"必须通过"。
2. **父状态白名单**：`_PROBE_PARENT_OK = (None, "ok")` + `_PARENT_NOT_EXECUTED`
   把 12 个非 ok 状态字面量映射成中文说明；**未知状态跳过但打印原始值**
   —— 静默放行会让"新增了一种没见过的状态"退化成没有检查。
3. **`_dig_present(obj, path) -> (found, value)`** —— `_dig()` 对"键不存在"
   与"键存在但值为 `None`"返回同一个 `None`，而 `hvg_fallback=None`
   的意思是**没有回退**（好事）。两者必须能区分。

### 31.3 Form C：恒真判据

`manifest:human_review` 的 `ok` 曾写死 `True`、`required: False` ——
**"一个都没登记"被写成"全部已确认"**。修法：`_n_hr > 0` 才可能为真；
`human_review_confirmed` 只列 `status in ("confirmed","overridden","not_needed")`
的节点；**默认 `pending` 不算失败**（判红会让每个 job 都红，反而没人看）。
**两仓同形**（本仓 `main_analysis.py` + `scrna-pipeline-skill/scripts/main_analysis.py`）。

### 31.4 标定：正向 139 项 + 反向 10/10，缺一不可

- 正向 `D:\tmp\_q28\calib_e69_probes.py` **139 项 / 0 失败**；
  `calib_e69.py` 100 项（纯函数）；`calib_e69_consumer.py` 36 项（验收消费端）。
- 反向 `D:\tmp\_q28\neg_e69_probes.py`：11 类注入 → **符合预期 10 类 /
  不符合 0 类**，`exit=0`，末尾 `源码已还原: True（main=88962, traj=25549）`。
- **判红必须伴随非空 FAIL 摘要。** 只有 `exit != 0` 不算证据 ——
  `subprocess.run` 少 `cwd` / `PYTHONIOENCODING` 时，子进程按 cp936 读
  UTF-8 源码抛 `UnicodeDecodeError`，**崩溃的非零退出被误读成"判红"**。

### 31.5 三个反复踩到的实现坑

1. **消费端标定看不见产生端。** 把产生端的 `morans_I_defined` 键改名后
   12 条探针**全绿** —— 因为消费端读的是标定脚本构造的 JSON 夹具，
   **产生端写什么它根本不知道**。补的判据必须是**源码级**：
   `08_spatial_trajectory.py:404` 必须是 `"morans_I_defined": _mi_pair_defined`。
2. **`_code_only`（剥 COMMENT+STRING）不能用来查字典键名** ——
   键名本身就是 STRING token，被一起剥掉，判据**永远找不到**。
   查键名要用只剥 COMMENT 的 `_no_comment`。
   **同一个文件里两种剥离策略各服务一条判据。**
3. **tokenize 的 token 是无空格拼接** —— `"key":value` 匹配得上，
   `"key": value` **永远匹配不上**。

### 31.6 收口

`common.write_json` 加 `allow_nan=False` + `_scrub_nonfinite`（递归
dict/list/tuple/numpy；**dict 的键必须是 `str`，否则 `default` 回调崩**）
—— 以前会写出裸 `NaN`（非法 JSON 字面量）；`lib/alignment.py` 是唯一
绕过 `write_json` 的写盘点，已收口。`common.read_json` 改成**损坏时抛异常**
（原来 `except Exception: return None` 把"文件不存在"和"文件坏了"混成一种），
另新增 `read_json_or_none`。

> **`np.float64` 是 `float` 的子类，`np.float32` 不是** ——
> 后者会一路走到 `json.dump` 的 `default` 回调。

台账：`governance/15_ERROR_LEDGER.md` E-69
（同族前两次：E-68 自查同族两处、E-61 静态检查名单只覆盖它自己）

## 32. 循环里建的图必须在循环体内保存（E-70，2026-09-26）

**规则：`for` 里 `plt.subplots()` 建的图，保存调用必须在**同一个循环体**里 ——
缩进错一格不会报错，只会少出图。**

现场：`scripts/07_spatial_communication.py` 画 top3 配体受体对。E-69 把整段包进
`if top3.empty: ... else:` 时给外层几行加了 4 个空格，**而循环体的 6 行
（含 `save_fig`）留在了原来的 8 空格** —— 于是它们跑到 `for` 外面去了。
后果（实测 CI artifact）：循环建了 3 张图、**一张都没保存**；循环结束后 `ui=3`、
`r` 是最后一行，所以只出 `03-07-02-unit3-icam1-itgal`，
**E-68 基线 76 张 → 本轮 74 张**。**而数据侧完全正确**：
`communication_status.json` 的 `n_z_defined=62 / n_z_undefined=0`，
`top_enriched` 前三名正是 CXCL12/CXCR4 z=1.939、CCL21/CCR7 z=1.901、
ICAM1/ITGAL z=1.569 —— **要画的三对就是它们，图却只出了一张。**

**四套门禁为什么全绿：** `check_figures.mjs` 只看"有没有墨 / 贴不贴边"；
`check_fig_names.mjs` 的账目判据是 `if (deficit > declared)`，**把 `declared`
（`DYNAMIC_FIG_BASES` 槽位和）当上限用** —— 本例 `deficit=1`、`declared=3`
⇒ **`1 > 3` 为假**，结构上不可能看见这一类；验收层 `figures:dynamic` 只要求
"每组 ≥1 张"；`py_compile` / `check_py_names.py` 语法与名字全对。
**只有"跨 artifact 比对图名集合 + 亲读像素"才看得见。**

**门禁：** `tools/check_py_names.py` 新增第三条判据（同文件、不新建检查器）。
对每个 `For` / `AsyncFor` / `While`，在**同作用域**子树（`walk_same_scope`，
**不下潜**嵌套 `def` / `lambda` / `class`）里数建图调用
（`MAKE_FIG_TAILS = {"subplots", "figure"}` 或名字以 `sc.pl.` 开头）与保存调用
（`SAVE_FIG_TAILS = {"save_fig", "savefig", "close"}` **或**
`save_helper_names(tree)` —— 本文件里"函数体内含保存调用"的函数名）。
`makes and not saves` ⇒ 判红。
**`save_helper_names` 是标定抓出来的必要修正**：没有它时"循环里建图后调同文件
`_draw(i)`、`_draw` 内部 `save_fig`"会被**误报**，而**误报会让门禁被关掉**
（同 E-64 首版把 `releases/` 当姊妹仓、E-61 行号指向另一段正确代码）。
**已知局限：跨模块的保存帮助函数仍会判红。**

**双向标定 20 项全过**（`D:\tmp\_e70\calib_gate.py`）：11 条合成用例 + 反向注入
真仓库（把 `pair_slug` + `save_fig` 退回 8 空格 ⇒ 判红、带非空摘要、报
`07_spatial_communication.py:391`、提示含 `L388`、**只有 1 处命中**；还原后逐字节
一致）。**判红必须断言摘要非空** —— 少 `cwd` / `PYTHONIOENCODING` 时子进程按
cp936 读 UTF-8 源码抛 `UnicodeDecodeError`，**崩溃的非零退出会被误读成"判红"**。

**验收层同步把"少出了几张"摆出来（判据不变）：** `figures:dynamic` 详情追加
`（产出/槽位：03-03-04 1/3, 03-04-01 30/30, ...）` —— 用**有缺陷的那一轮**
artifact 干跑，`03-07-02 1/3` 直接印在验收详情里，而 `n_checks` 与基线逐位相同
79、**消失/新增的检查 id 均为空**（只加详情、不加判据）。

台账：`governance/15_ERROR_LEDGER.md` E-70（任务行 `governance/02_TASKLIST.md` R-04h）
