# 模块零：运行清单

参考规范：三大部分整合文档「模块零：语言与运行时规范」。
姊妹项目 `geo-normal-pipeline-skill` / `scrna-pipeline-skill` 的
`references/module0.md` 是同一份文档，接口**同名同义**，
三部分的清单可以并排读。

---

## 0. 为什么要有这一层

**没有运行清单的分析结果不是结果。**

半年后拿到一份 `svg_results.csv`，如果不知道当时装的是哪个版本的
SpatialDE、输入的 h5ad 是哪个哈希、随机种子是多少、Moran's I 的
置换跑了几次，那份 CSV 就**无法被复现，也无法被质疑** ——
而不可质疑的结论没有价值。

这一层不产出任何生物学结论，它只回答一个问题：
**"这份结果是在什么条件下跑出来的？"**

产物：`results/<dataset_id>/run_manifest.json`

**为什么不塞进 `state.json`：** `state.json` 记的是"这一步跑没跑成"，
每步重写；manifest 记的是"本轮是在什么条件下跑出来的"，是证据，
写入后不该再变。混在一起会让后者被前者覆盖。

---

## 1. 规范条款与接口的对应

| 条款 | 要求 | 接口 | 落到 manifest 的哪个字段 |
|---|---|---|---|
| §0.3 | `pip freeze` **全量**输出 | `capture_versions()` | `versions` |
| §0.3 | 关键工具**逐个**记版本 | `capture_versions()` 的 `KEY_PACKAGES` | `key_versions` |
| §0.3 | 所有随机过程固定种子并记录 | `record_params()`（`params` 里含 `seed`）+ `m["seed"]` | `params`、`seed` |
| §0.4 | 输入数据哈希 | `record_input()` | `inputs` |
| §0.4 | 全部关键参数 | `record_params()` | `params` |
| §0.4 | Agent 决策链 | `record_decision()` | `decisions` |
| §0.4 | 人工干预记录 | `record_human_review()` | `human_review` |
| §0.2 | 跨语言转换前后维度、丢失字段 | `record_cross_language()` | `cross_language` |

---

## 2. 接口清单

| 函数 | 作用 | 备注 |
|---|---|---|
| `MANIFEST_NAME` | `"run_manifest.json"` | 三个仓库一致 |
| `KEY_PACKAGES` | 文档点名的关键工具 | 见 §4 |
| `NAMED_TOOLS` | 点名但用不了的工具 + 理由 | 见 §5 |
| `manifest_path(cfg)` | 清单路径 | |
| `read_manifest(cfg)` | 读清单；文件不存在返回 `{}` | |
| `init_manifest(cfg, language="python")` | **开新一轮，清掉上一轮** | 见 §3.1 |
| `capture_versions(cfg, ...)` | 全量已装包 + `KEY_PACKAGES` 逐个 | |
| `record_input(cfg, path, ...)` | 输入文件的 sha256 与字节数 | |
| `record_params(cfg, params)` | 参数（`update` 合并，可多次调用） | |
| `record_decision(cfg, node, q, a, evidence)` | 决策链 | |
| `record_human_review(cfg, node, required, status, note)` | 人工复核节点 | |
| `record_cross_language(cfg, src, dst, fmt, before, after, lost, ...)` | 跨语言转换 | |
| `manifest_summary(cfg)` | 供验收用的摘要 | 见 §3.3 |
| `probe_named_tools(log, only)` | 把 `NAMED_TOOLS` 整理成可写进状态 JSON 的表 | |
| `named_tools_note()` | 一句话说明为什么多数点名工具没用上 | |

---

## 3. 四条不能省的约定

### 3.1 `init_manifest` 必须清掉上一轮

上轮的清单留在那里冒充本轮，**比没有清单更糟** —— 它看起来是证据，
实际是过期证据。这和"每步开跑前先删自己的状态文件"是同一个道理
（`STEP_STATUS_FILES` 那段注释）。

### 3.2 输入哈希在**步骤跑完之后**才登记

可选步骤这轮有没有产物，跑完才知道；在开头登记会把"上轮残留"记成本轮输入。

### 3.3 `manifest_summary` 同时报两个缺失数

```
inputs_missing           所有缺失的输入（可见，供人看）
inputs_missing_required  只有 required=True 的缺失（供验收判 FAIL）
```

**第一版只有一个数，把可选输入缺失也算成失败**，结果每轮都红。
可选输入缺失是**正常状态**（`niche_labels.csv` 没配、外部参考不存在），
判成 FAIL 会让红灯失去意义。但它必须**可见** —— 所以两个数都要有。

### 3.4 未装的工具要记成 `None`，不能省略键

```json
"key_versions": { "SpatialDE": "1.1.3", "BayesSpace": null }
```

`"BayesSpace": null` 和"没有 BayesSpace 这个键"是两件事：前者是
"**查过了，装不上**"，后者是"**没查**"。省略会让读者分不清。

### 3.5 `pip freeze` 不用 `subprocess` 抓

**沙箱下管道捕获会 EPERM**（`child_process` 默认 `stdio: 'pipe'`）。
用 `importlib.metadata` 枚举已装包。
包名归一化走 **PEP 503**（`_norm_pkg`）：`spatialde` / `SpatialDE` /
`spatial-de` 是同一个包，不归一会让 `key_versions` 里出现查不到的键。

---

## 4. `KEY_PACKAGES`：三部分的清单

文档点名的工具**逐个列出**，装了的记版本、没装的记 `null`。

| 部分 | 数量 | 内容 |
|---|---|---|
| Part 1 | 15 | `GEOquery` `limma` `WGCNA` `clusterProfiler` `GSVA` `glmnet` `survival` `survminer` `timeROC` `rms` `STRINGdb` + `TRRUST` `ChEA3` + `scTenifoldKnk` `PerturbNet` `RegVelo` |
| Part 2 | 25 | `scanpy` `anndata` `scvi-tools` `cellbender` `harmonypy` `scvelo` `celltypist` `pyscenic` `liana` `doubletdetection` `scrublet` + 拟时序 `palantir` `scfates` `cytotrace` + R 包 `monocle3` `slingshot` `cellchat` `soupx` `scdblfinder` + `scTenifoldKnk` `PerturbNet` `RegVelo` |
| **Part 3（本仓库）** | 16 | `SpatialDE` `SpatialDE2` `spacexr` `BayesSpace` `SPARK-X` `cell2location` `STAGATE` `SpaGCN` `SpaceFlow` `stLearn` `ISORT` `Bering` `BOMS` + LIANA 等 |

**Part 3 的清单里为什么有 R 包（`BayesSpace` / `RCTD` / `CellChat`）：
** 规范说 Part 3 是"Python (+R via rpy2)"，但**本仓库 CI 里没有 R**，
所以这条路径实际不存在。这些键**本来就该是 `null`** —— 但必须留着，
否则读者不知道"是没查还是不该有"。确切原因写在 `NAMED_TOOLS` 里（见 §5）。

---

## 5. `NAMED_TOOLS`：点名工具"为什么没用上"的登记

**这是本仓库最容易被误读的地方。** 域划分、解卷积、通讯、轨迹四步
都有产出、都有图 —— 看起来"§3 做完了"，而文档点名的方法
**一个都没用上**，跑的全是内置实现。

所以把"为什么没用上"逐条落盘。**判据是"理由写了没有"，
不是"工具跑了没有"。** 将来某个工具能装了，验收应该依然 PASS
（理由变成"已装"），而不是因为 `available=False` 就变红 ——
那会把"如实记录"惩罚成失败。

四类 `kind`：

| kind | 含义 | 本仓库的例子 |
|---|---|---|
| `r_package` | R/Bioconductor 包，CI 无 rpy2 | `BayesSpace`（§3.2）`RCTD`（§3.3）`CellChat`（§3.5） |
| `not_on_pypi` | 真包不在 PyPI | `STAGATE`（§3.2）`StPedf`（§3.6） |
| `deps` | PyPI 有真包，依赖链跑不动 | `SpaGCN` `cell2location` `SpaceFlow` `Stereopy-TGPI` `stLearn` `Bering` `BOMS` |
| `name_taken` | **PyPI 上那个名字是另一个不相干的包** | `ISORT`（§3.6） |

### `name_taken` 是最危险的一类

因为 `pip install` 会**成功**，装进来的是完全无关的东西。实测元数据：

```
pip install ISORT -> PyCQA/isort，Python 的 import 排序工具
pip install sparkx -> 高能物理的碰撞相对论运动学
pip install edgeR  -> "Redirect Microsoft Edge to your preferred browser"
pip install slingshot -> "Index Migration for ElasticSearch"
```

`ISORT` 是 §3.6 点名的空间轨迹工具。**这个名字绝对不能进
`requirements.txt`。** 装不上会立刻报错，**装错了要到 import 或
跑出结果才发现** —— 而那时"结果"可能已经在图上看着挺像回事了。

**判断依据是 summary / author / project_urls，不是"名字存不存在"。**
反例：`SingleR` 在 PyPI 上就是真的 BiocPy/singler。

### 几条实测的依赖链细节

- **`SpaGCN`** 1.2.7 是真包，但依赖 `louvain` —— **该包最新版 0.8.2
  没有 py3.12 wheel，只有 sdist**，要从 2019 年的 C++/Cython 源码编译。
- **`STAGATE`** / `STAGATE_pyG` / `stagate` 三个名字在 PyPI **全部 404**。
- **`Stereopy`** 1.6.2 的 `requires_python` 是 `'<3.9,>=3.8'` ——
  **上限卡死**，与 CI 的 3.12 不兼容。
- **`stLearn`** 1.4.1 依赖 `numpy>=2.4.0` / `scipy>=1.17.0` /
  `scanpy>=1.12.0` / `zarr>=3.1` / spatialdata 全家桶 + torch +
  torchvision + geopandas + dask（20+ 个），与本仓库钉的区间冲突。
- **`BOMS`** 依赖 `mkl` + `mkl-service`（conda 时代的 Intel MKL 绑定）。

---

## 6. 人工复核节点

| 节点 | required | 说明 |
|---|---|---|
| `cell_segmentation` | FALSE | 细胞分割参数调整（Visium 上不适用） |
| `domain_number` | TRUE | 空间域数量的确定 |
| `deconv_reference` | TRUE | 去卷积参考数据的选择 |
| `spatial_traj_direction` | TRUE | 空间拟时序轨迹方向的验证 |

**默认 `pending`，不算失败。** 自动化流水线不能替人签字 ——
把未确认的节点记成已确认，等于把复核节点变成摆设。
但必须**可见**：验收里作为"可见但不阻断"的项列出。

---

## 7. 清单必须和结果同时可及

`run_manifest.json` 落在 `results/<dataset_id>/` 下，随 artifact 一起上传
（artifact 包含整个数据集目录）—— **它必须和结果同时可及，
否则追溯链是断的。** 只把清单写在 CI 日志里不行：日志会滚掉，文件不会。

---

## 8. 怎么读一份清单

```python
import json, pathlib
m = json.loads(pathlib.Path("results/lymph_node/run_manifest.json").read_text())
m["language"]        # "python"
m["seed"]            # 随机种子
m["inputs"]          # 每个输入文件的 sha256 与字节数
m["key_versions"]    # 点名工具逐个的版本（null = 查过了装不上）
m["params"]["named_tools"]   # 点名工具的缺口登记
m["decisions"]       # Agent 决策链：问题 / 结论 / 证据
m["human_review"]    # 人工复核节点及状态
m["cross_language"]  # 跨语言转换记录
```

**验收读的是 `manifest_summary(cfg)["inputs_missing_required"]`，不是
`inputs_missing`** —— 后者包含可选输入的缺失，那是正常状态。

---

## 9. 本地与 CI 的数值差异

清单会记下 `versions`，所以**"同一个 commit 在两台机器上给出不同数字"
是可以解释的**，不是玄学。实测：

| 量 | 本地（Python 3.14 / scipy 1.18） | CI（Python 3.12） |
|---|---|---|
| Moran's I（空间平滑后） | `0.7671` | `0.7670` |
| 邻居同类比例（空间平滑后） | `0.5355` | `0.5355` |
| 邻居差（空间平滑后） | `0.0358` | `0.0358` |

末位差异是**版本差异，不是回归**。所以比对数值前先比对 `key_versions` ——
清单存在的意义之一就是让这类差异**可归因**，而不是每次都要重新怀疑代码。

**反过来也要注意：** 本地与 CI 一致的量（域数 9 → 13、邻居差 0.0517 →
0.0358、通讯 62/63 对可用）才是可以当结论用的。报数时区分这两类。
