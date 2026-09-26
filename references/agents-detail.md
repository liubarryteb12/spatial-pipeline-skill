# 长文规则原文归档（AGENTS.md 精简版对应的完整原文）

> 本文件是原本写在 `AGENTS.md` 里的长文规则叙事的**未删节归档**：`AGENTS.md` 对每条规则只保留一份精简版，
> 因为 harness 在 **65536 字节**处截断注入的指令文件，长文会让排在后面的规则**整条进不去**。
> **这里没有任何内容被取代** —— 需要完整的事故经过与实测证据时读本文件。
> 持久的事故记录仍然是 `governance/15_ERROR_LEDGER.md`（该仓库在工作区，位于 `..\governance\`）。

## 原规则 20. `install_requires` 里有某个包，不等于运行时会 import 它

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

## 原规则 21. `set_seed()` 不 seed torch，而第三方工具的随机性可能在 torch 里

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

1. SpaGCN **包内** `models.py:55` `KMeans(self.n_clusters, n_init=20)` —— **没有
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

## 原规则 17.1. 写死的名单只覆盖了它自己 —— 现在是真正的逐作用域分析（E-61）

上面那份名单**只有 10 个名字**，而 `common.py` 实际导出 **68 个**。
于是它挡住的全是"恰好被列进去的"，没列进去的照旧漏到 CI。
实测踩到（2026-09-26）：`03_spatial_domains.py` 漏 import
`spot_radius_plot_units`，本地 `check_py_syntax.mjs` 报"全部通过"，
CI 跑 25 分钟到 §3.2 的 H&E 叠图段才
`NameError: name 'spot_radius_plot_units' is not defined`，
04–08 步全没跑。

**根因不是"名单短了一点"，而是判据的输入域与它要防的缺陷不匹配** ——
漏 import 的可以是任意一个导出。把名单补成 68 个只是把同一种漂移往后推一次
（`common.py` 下次加函数又会漂）。

规则 17 当年放弃做作用域分析的理由是"那是重写一个 linter"。
**那个顾虑是对的，但作用域分析不必自己写**：标准库 `symtable` 就是
CPython 编译器的符号表，按作用域给出每个名字是 local / parameter /
imported / free（闭包）/ global，正是这里需要的。

现在 `tools/check_py_names.py` 做**四件事**：

| 检查 | 判据 |
|---|---|
| 未定义名字 | 某名字在某作用域**被引用**，却既非该作用域局部绑定、也非闭包自由变量、模块层也没有、也不是内置 → 运行时必然 `NameError` |
| 幽灵 import | `from common import X` 而 common 模块级没有 `X` |
| 不安全构造 | 出现 `import *` / `globals()` / `exec` / `eval` / `vars` / `locals` → **判红退出**，而不是假装通过 |
| 循环建图未保存（E-70） | 某个 `For` / `AsyncFor` / `While` 的**同作用域**子树里建了图（`subplots` / `figure` / `sc.pl.*`），却**没有**保存调用（`save_fig` / `savefig` / `close` / 本文件内"函数体里含保存调用"的帮助函数）→ 判红。见规则 32 |

**两仓该文件逐字节相同**（SHA256 `712AEA6995962ED70801EF9B39FCA3CF5749B026DB7E0EFF392AA1310A45D8EE`，17939 字节），改一侧必须同步并比对哈希。

**三个实现上的坑（都实测踩过）：**

1. **不安全构造必须用 AST 判，不能扫子串。** 第一版扫裸子串，于是本文件
   自己的模式元组 `("import *", "globals()", "exec(")` 命中了它自己 ——
   **检查器把自己的源码判红**。"检查器要检查的东西"与"检查器描述自己要检查
   什么"在文本上无法区分，只有语法结构能区分（同规则 25.1：`stripComments`
   把字符串换成 `""` 后把要检查的东西本身擦掉了）。确实需要时在同一行写
   `# py-names: unsafe-ok —— <理由>` 豁免（**故意做成要写一句话的**，
   静默豁免会让检查退化成没有检查）。
2. **行号必须按作用域定位，不能全文件找首次出现。** 第一版在整棵 AST 上找
   第一个同名 `Name`，于是 `plt` 被报成 `common.py:859` —— 那是**另一个
   函数**（`plot_marker_dotplot`）里的 `plt`，那个函数自己 `import
   matplotlib.pyplot as plt`，行号指向一处**没问题的代码**。真正有问题的在
   `fix_dotplot_legends` 的 L1086。**指向错的行号比不指行号更糟**：下一个人
   会去读一段正确的代码然后困惑。
3. **"common 提供这个名字"的提示要扣掉 common 自己 import 进来的。**
   `import numpy as np` 被删时报出「common 导出过 np —— 加进
   `from common import (...)`」，而 `np` 出现在 common 命名空间里只是因为
   **common 自己也 import 了 numpy**。照着改会写出 `from common import np`
   —— 把 common 的内部依赖当接口用。两个集合分工：含 import 的用于**幽灵
   import 判据**（`from common import np` 运行期确实能成功，不算幽灵），
   只有 `def`/赋值的才用于**提示**。

**标定（`D:\tmp\_r04\calib_names.py`，13 项全过）：** 正向两仓零命中；
geo（纯 R）判为**不适用且不判红**；传错目录判红；反向逐类注入 —— 删 common
import（报 `log_info`）、删第三方 import（报 `np`，证明不限于 common 导出）、
幽灵 import、`import *`、`globals()`、`eval()` 各自只让它自己那条响；逃生舱
写了理由后放行；a7bd1bd 夹具精确报出
`scripts/03_spatial_domains.py:754 spot_radius_plot_units`。

> **"没有 Python 文件"是"不适用"，不是"失败"。** geo 是纯 R 仓库
> （`scripts/` 下 0 个 `.py`），跑到这里必须放行 —— 判红会让 pre-push 在 geo
> 上永远红，而那条告警与 geo 的改动毫无关系。但"仓库根目录传错了"必须
> 与它长得不一样，所以分开检查（`scripts/` 和 `tools/` 都不存在才判红）。

**门禁接线：** CI 在两个 workflow 的「静态检查（不装依赖）」那一步跑；
`governance/hooks/pre-push.mjs` 的 `PY_GATES` 表（Python 门禁用解释器跑，
不是 `process.execPath`；本机没有 `python` 时记为**提醒**而不是通过 ——
"没跑成"和"跑过了没问题"必须长得不一样）。

## 原规则 28. 验收层四条补强（Q-27，2026-09-25）

姊妹项目 `scrna-pipeline-skill/AGENTS.md` 规则 25 是同一批补强的另一半。
**本仓库的口径与它不同**（规则 23：本仓库的图检查是 glob 数张数、
不引用具体文件名），所以**没有照抄**，而是按空间侧的结构重新设计。

### 28.1 步骤函数自己返回的 `status` 必须有人读（验收层）

`run_steps` 把 `fn(cfg)` 的返回值记进了 `results[name]["result_status"]`，
**而全仓没有一个消费者** —— 于是 `08_spatial_trajectory.py` 返回
`bad_root` / `not_applicable` / `missing_pca` 时，验收照样记 `status="ok"`。

**这不是"步骤崩了"，而是"步骤跑完了、但结果是『没做成』"** ——
恰恰是最容易读成成功的一种。现在 `step_result:<sid>` 消费它，
判据仍只把 `failed`/`error`/`fail` 判红；`bad_root` 之类的取值是
**设计如此地没做成**，只可见（`honesty`）、不阻断。

> **实现踩到的坑：** `chk(cid, kind, ok, detail, severity="required")` 的
> **第二个位置参数是 `kind`（归类），不是 severity**。第一版把
> `"required" if bad else "honesty"` 传在第二个位置，于是
> `bad_root` 这些**本该只可见**的取值全被记成 required —— 注入测试里
> 7 类 `result_status` 只有 `failed` 该红，实测却红了一片。
> **`kind` 与 `severity` 同名同型同位置数相邻，传错不报错。**

### 28.2 嵌套 `status` 的 `failed` 必须判红（验收层）

`status` 不只在顶层。`domain_status.json` 的
`domain_methods.STAGATE.status`、`deconvolution_status.json` 的
`cell2location.status`、`svg_status.json` 的 `spatialde.status` 都是
**嵌套**的，而旧验收层只看顶层 —— **里面崩了、顶层还是 `ok`**。

姊妹仓库实测（E-48）：顶层分布 `{ok: 7, not_configured: 1}`，
嵌套分布里躺着唯一一条 `failed`（`figsize` 三元素元组让整段抛异常），
而它所属文件的顶层 `status` 是 `ok`。**顶层全绿、里面已经崩了。**

实现是 `_iter_nested_status(obj, path=())` 递归产出
`(路径, 值, 同级 reason)`，扫**全部** `*status.json`（不只可选步骤那几个）。
判据只把 `failed`/`error`/`fail` 判红 —— 把 `needs_reference` /
`package_missing` / `not_done` / `not_applicable` / `disabled` /
`missing_domains` / `bad_root` 判红会让每个 job 都红，反而没人看。

**本仓库实测嵌套分布 `{ok: 11, needs_reference: 1, package_missing: 1}` ——
零命中、零误伤。**

### 28.3 溢出检测必须落盘，不能只打 WARN（产物级）

E-49 的形态：`_content_overflow()` **正确检测到了**
`width_overflow_frac=0.3523`、**正确打了 WARN**，然后**没有任何人读**。

**"检测到了"不等于"有人会知道"。** 现在 `save_fig()` 里的
`_record_figure_overflow(cfg, name, bad)` 把溢出**落盘**到
`results/<dataset_id>/figure_overflow.json`，验收层读它并**判红（required）**。

三条配套约束：

1. **与状态文件同理，跑前必须删旧文件** —— 否则图修好了旧记录还在，
   验收把已修好的图判红（规则 16 的同一条道理）。
2. **只累积、不覆盖** —— 同一次运行多张图溢出要全部记下。
3. **记录可被修复** —— 标题改短后不应再新增记录。

### 28.4 图验收从"数张数"换成"声明过的每张在不在"（验收层）

旧判据是 `len(figs) >= 8` —— **纯计数**，于是"该有的图没有"这一整类
问题没有任何检查看得见：姊妹仓库实测 5 张图因 `figsize` 三元素元组
从未产出过，而验收 70 项全绿（E-48）。

声明**从源码扫出来**（`declared_figures()`），不手抄表 —— 手抄的表会漂移，
漂移方向恰好是"新加的图不在表里"，等于把盲区原样再造一遍。
扫的时候要**排除 `main_analysis.py` 自己**（它的检查项里也有图名字符串，
不排除会把"检查项的名字"当成"声明的图"）。

三条判据分工：

| 检查 | 判据 |
|---|---|
| `figures:declared` | 逐名核存在性；条件图不适用时进合法豁免 |
| `figures:dynamic` | 声明了槽位的**动态图名**每组**至少产出 1 张** |
| `figures:count` | `>= 8`，**保留为下限兜底**（声明扫描本身失效时还能拦住） |

**动态图名为什么要单独一条**：名字运行期才拼得出来（`DYNAMIC_FIG_BASES`
声明槽位数），不能逐名检查；但**整组 0 张**说明那段循环整段没跑 ——
槽位是上限不是精确值，少出合法。

> **这里有一个"判据因为输入域重叠而永远为真"的坑，正向标定看不出来：**
> `03-03-04` 这个前缀**同时**有静态图（unit1/unit2/unit3）和动态图
> （unit8 的 Jaccard 矩阵）。按前缀计数而不扣掉静态声明图的话，
> 动态循环整段没跑时计数仍是 2（全是静态的），**判据永远不会响**。
> 只有注入"整组动态图消失"才会暴露 —— 见规则 29。

### 28.5 双向标定：正向零误报 + 反向逐类注入

**这两半缺一不可。** 正向（拿真实 artifact 干跑）只能说明"没误报"；
**一个永远返回 True 的判据在干净 artifact 上也是零命中**。
所以逐类注入，确认每一类缺陷都真的会让**它自己那条**判据变红：

| 注入 | 应红的判据 | 实测 |
|---|---|---|
| 删 1 张静态声明图 | `figures:declared` | ✅ 只它红 |
| 删整组动态图（30 张） | `figures:dynamic` | ✅ 只它红 |
| 写 `figure_overflow.json` | `figures:overflow` | ✅ 只它红 |
| 嵌套 `STAGATE.status = failed` | `status:nested` | ✅ 只它红 |
| 删全部图 | `count` + `declared` + `dynamic` | ✅ 三条都红 |

正向：真实 artifact（`lymph_node`，76 张图）干跑 —— 新增 4 条检查
**全部 PASS**，`n_checks` 59 → 63，**没有一条既有检查消失**。
（干跑只调 `run_acceptance()` 并把所有写盘函数 patch 成 no-op ——
否则会重写 artifact 里的 `acceptance_report.json`，那就变成
"自己改自己的证据"了。）

## 原规则 32. 循环里建的图必须在循环体内保存（E-70，2026-09-26）

**规则：`for` 里 `plt.subplots()` 建的图，保存调用必须在**同一个循环体**里 ——
缩进错一格不会报错，只会少出图。**

现场（本仓，2026-09-26）：`scripts/07_spatial_communication.py` 画 top3 配体受体对。
E-69 把整段包进 `if top3.empty: ... else:` 时给外层几行加了 4 个空格
（`_prods` / `for r in top3.itertuples()` / `vmax_lr` / `DYNAMIC_FIG_BASES` /
`for ui, r in enumerate(top3.itertuples(), start=1):`），**而循环体的 6 行
（`ax.set_aspect` / `invert_yaxis` / `set_xticks` / `fig.colorbar` /
`pair_slug` / `save_fig`）留在了原来的 8 空格** —— 于是它们跑到 `for` 外面去了。

后果（实测 CI artifact）：循环建了 3 张图、**一张都没保存**；循环结束后
`ui=3`、`r` 是最后一行，所以只有 `03-07-02-unit3-icam1-itgal` 被画出来，
**`cxcl12-cxcr4` 与 `ccl21-ccr7` 两张凭空消失，另外两张图泄漏未关闭**。
E-68 基线 76 张 → 本轮 74 张。

**而数据侧完全正确**：`communication_status.json` 的
`n_z_defined=62 / n_z_undefined=0`，`top_enriched` 前三名正是
CXCL12/CXCR4 z=1.939、CCL21/CCR7 z=1.901、ICAM1/ITGAL z=1.569 ——
**要画的三对就是它们，图却只出了一张。**

### 32.1 四套门禁为什么全绿

| 门禁 | 为什么看不见 |
|---|---|
| `check_figures.mjs` | 只看"有没有墨 / 贴不贴边" —— 少一张图与多一张图它都看不见 |
| `check_fig_names.mjs` | 账目判据是 `if (deficit > declared)`（`deficit = 调用数 - 字面量数`、`declared = DYNAMIC_FIG_BASES` 槽位和）。本例 `deficit=1`、`declared=3` ⇒ **`1 > 3` 为假** —— 它把 `declared` 当**上限**用，结构上不可能看见这一类 |
| 验收层 `figures:dynamic` | 只要求"每组 ≥1 张"（规则 28.4：槽位是上限不是精确值，少出合法）|
| `py_compile` / `check_py_names.py` | 语法与名字全对 |

**只有"跨 artifact 比对图名集合 + 亲读像素"才看得见**（E-69 推送后亲读
spatial 产物 74 vs 76）。**这正是"改完自己去看渲染像素"这条纪律的价值所在。**

### 32.2 门禁：`tools/check_py_names.py` 新增第三条判据

同文件、不新建检查器（避免两仓同步 + 接线 + 改文档引用的额外面积）。
对每个 `For` / `AsyncFor` / `While`，在**同作用域**子树（`walk_same_scope`，
**不下潜**嵌套 `FunctionDef` / `Lambda` / `ClassDef`）里数建图调用
（`MAKE_FIG_TAILS = {"subplots", "figure"}`，或名字以 `sc.pl.` 开头）与
保存调用（`SAVE_FIG_TAILS = {"save_fig", "savefig", "close"}`，**或**
`save_helper_names(tree)` —— 本文件里"函数体内含保存调用"的函数名）。
`makes and not saves` ⇒ 判红，打印
`{rel}:{建图行}  循环体里建了图，但同一个循环体内没有保存（循环在 L{循环行}）`。

**`save_helper_names` 是标定抓出来的必要修正**：没有它时
"循环里 `plt.subplots()` 之后调同文件 `_draw(i)`、`_draw` 内部 `save_fig`"
会被**误报**（`walk_same_scope` 不下潜嵌套 `def`）—— **误报会让门禁被关掉**
（同 E-64 首版把 `releases/` 当姊妹仓、E-61 行号指向另一段正确代码）。

**已知局限（写下来，不假装没有）**：**跨模块**的保存帮助函数仍会判红。

### 32.3 双向标定（`D:\tmp\_e70\calib_gate.py`，20 项全过）

11 条合成用例（`tempfile.TemporaryDirectory()` 造仓 + `shutil.copy2` 门禁本体）：
循环建图+循环外保存 **红** / 循环体内保存 绿 / 图建在循环外循环内保存 绿 /
保存在循环里调的同文件 `def` 里 绿 / `plt.close(fig)` 算保存 绿 /
`sc.pl.*` 不保存 **红** / `while` 建图不保存 **红** / `fig.savefig` 绿 /
save 在 `if` 里（仍属循环体）绿 / 删掉 `from common import save_fig` 仍 **红**
（证明 E-61 判据没被破坏）/ 跨模块帮助函数 **红**（已知局限）。
反向注入真仓库：把 `pair_slug` + `save_fig` 两行退回 8 空格 ⇒ 判红、带非空摘要、
报 `07_spatial_communication.py:391`、提示含 `L388`、**只有 1 处命中**（未误伤
其它脚本）；`finally` 还原后回到绿且**逐字节一致**。

> **判红必须断言摘要非空。** `subprocess.run` 少 `cwd` / `PYTHONIOENCODING`
> 时子进程按 cp936 读 UTF-8 源码抛 `UnicodeDecodeError`，
> **崩溃的非零退出会被误读成"判红"**（同规则 31.4）。

### 32.4 验收层同步把"少出了几张"摆出来（判据不变）

`figures:dynamic` 的判据**不改**（少出确实合法），但详情里追加
`（产出/槽位：03-03-04 1/3, 03-04-01 30/30, ...）` —— 用**有缺陷的那一轮
artifact 干跑**，`03-07-02 1/3` 直接印在验收详情里，而 `n_checks` 与基线
逐位相同 79、**消失/新增的检查 id 均为空**（只加详情、不加判据）。
**"检测到了"不等于"有人会知道"**（同规则 28.3）。

台账：`governance/15_ERROR_LEDGER.md` E-70（任务行 `governance/02_TASKLIST.md` R-04h）
