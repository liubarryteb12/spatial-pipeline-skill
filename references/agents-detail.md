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

