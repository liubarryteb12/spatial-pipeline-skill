"""
lib/common.py — 空间转录组流水线公共库

与 `scrna-pipeline-skill` / `geo-normal-pipeline-skill` 同一套约定。
本文只写空间特有的部分，其余见那两个仓库的 AGENTS.md。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import sys
import time
from pathlib import Path

import numpy as np
import yaml

# ============================================================================
# 确定性 —— **必须在 import scanpy/numba 之前设**
# ============================================================================
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("OPENBLAS_CORETYPE", "Haswell")
os.environ.setdefault("PYTHONHASHSEED", "0")


# ============================================================================
# 日志
# ============================================================================
_ORCHESTRATED = False


def set_orchestrated(flag: bool = True) -> None:
    global _ORCHESTRATED
    _ORCHESTRATED = flag


def _log(level: str, msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {level:<7} {msg}", flush=True)


def log_info(msg: str) -> None:
    _log("INFO", msg)


def log_warn(msg: str) -> None:
    _log("WARN", msg)


def log_error(msg: str) -> None:
    _log("ERROR", msg)


# ============================================================================
# 参数与配置
# ============================================================================
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="空间转录组流水线",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True, help="配置文件路径")
    p.add_argument("--steps", default=None, help="只跑指定步骤，逗号分隔")
    return p.parse_args(argv)


def load_config(path: str) -> dict:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"配置文件不存在: {path}")
    with open(p, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    if not isinstance(cfg, dict):
        raise ValueError(f"配置文件顶层必须是映射: {path}")
    if not cfg.get("dataset_id"):
        raise ValueError("配置缺少 dataset_id")

    did = str(cfg["dataset_id"])
    out = cfg.setdefault("output", {})
    out.setdefault("results_dir", f"results/{did}")
    out.setdefault("data_dir", f"data/{did}")
    out.setdefault("figures_dir", f"results/{did}/figures")

    ana = cfg.setdefault("analysis", {})
    ana.setdefault("seed", 20260919)
    # **300 dpi 是投稿图的底线。** 原来默认 150：183 mm 宽的图只有 1080 px，
    # 放大或印刷后字形和细线都发虚 —— 而"发虚"从图注上完全看不出来，
    # 文件大小也正常。
    # **注意这只是兜底**：`assets/config.lymph_node.yml` 里显式写了
    # `figure_dpi`，而 `setdefault` 只在键缺失时才生效 —— 所以那份配置
    # 不会用到这个默认值。**两处都要对。**（姊妹项目 scrna 实测踩过：
    # 只改了这里的默认值，artifact 里的 PNG 仍然是 150 dpi。）
    ana.setdefault("figure_dpi", 300)
    return cfg


def ensure_dirs(cfg: dict) -> None:
    for k in ("results_dir", "data_dir", "figures_dir"):
        Path(cfg["output"][k]).mkdir(parents=True, exist_ok=True)


def set_seed(cfg: dict) -> int:
    seed = int(cfg["analysis"]["seed"])
    random.seed(seed)
    np.random.seed(seed)
    # **样式必须在这里应用，不能等到 save_fig。** rcParams 只在 figure
    # 创建时被读取 —— 等图建好了再 plt.style.use，那张图仍然是旧样式。
    # 每个步骤脚本第一句都是 set_seed(cfg)，所以这里是唯一的正确位置。
    apply_style(cfg)
    return seed


def write_json(path, obj) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)

    def default(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.bool_,)):
            return bool(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, Path):
            return str(o)
        return str(o)

    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, ensure_ascii=False, indent=2, default=default)


def read_json(path):
    p = Path(path)
    if not p.exists():
        return None
    try:
        with open(p, encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:  # noqa: BLE001
        return None


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ============================================================================
# 步骤状态
# ============================================================================
def state_path(cfg: dict) -> Path:
    return Path(cfg["output"]["results_dir"]) / "state.json"


def read_state(cfg: dict) -> dict:
    return read_json(state_path(cfg)) or {"steps": []}


def record_step(cfg: dict, step_id: str, status: str, seconds: float = None,
                message: str = "", required: bool = True) -> None:
    st = read_state(cfg)
    st.setdefault("steps", [])
    st["steps"] = [s for s in st["steps"] if s.get("id") != step_id]
    entry = {"id": step_id, "status": status, "required": bool(required),
             "message": message}
    if seconds is not None:
        entry["seconds"] = round(float(seconds), 1)
    st["steps"].append(entry)
    st["dataset_id"] = cfg["dataset_id"]
    st["updated_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    write_json(state_path(cfg), st)


# ============================================================================
# 运行清单（模块零规范层）
# ============================================================================
# 规范来源：用户整合文档「模块零：语言与运行时规范」。
#   §0.2 跨语言接口 —— 只走 CSV；每次转换记录维度/metadata/丢失字段
#   §0.3 版本记录   —— pip freeze 全量 + 文档点名的关键工具单独记版本
#   §0.3 随机种子   —— 所有随机过程固定种子并记录
#   §0.4 运行日志   —— 输入数据哈希、软件版本、关键参数、决策链、
#                      人工干预记录、跨语言转换记录
#
# 产物：results/<dataset_id>/run_manifest.json
#
# **为什么不塞进 state.json：** state.json 记的是"这一步跑没跑成"，每步重写；
# manifest 记的是"本轮是在什么条件下跑出来的"，是证据，写入后不该再变。
# 混在一起会让后者被前者覆盖。
#
# **诚实性要求（AGENTS.md 规则 4）：** 没做的分析、没装的工具、没确认的
# 复核节点，都要在 manifest 里留下痕迹，不能因为"不影响结论"就不写。
# `decisions` 记的是**实际的选择**，不是"应该怎么做"。

MANIFEST_NAME = "run_manifest.json"

# 文档 §3 点名的工具。**没装的记 None，不省略键** —— 键消失和"版本是 None"
# 看起来完全不同，后者才说明"这个工具本该有但没装"。
KEY_PACKAGES = [
    # §3 核心 Python 包
    "scanpy", "anndata", "squidpy", "cell2location", "STAGATE", "SpaGCN",
    "SpaceFlow", "stLearn", "ISORT", "Bering", "BOMS",
    # §3.4 空间可变基因
    "SpatialDE", "SpatialDE2",
    # §3 点名的 R 包（本仓库无 rpy2 路径，正常就是 None）
    "spacexr", "BayesSpace", "SPARK-X",
]


# ---- 文档点名、但本仓库用不了的工具 -----------------------------------------
#
# 规范 §3.2 点名 BayesSpace / STAGATE / SpaGCN，§3.3 点名 RCTD / cell2location，
# §3.4 点名 SPARK-X，§3.5 点名 CellChat，§3.6 点名 StPedf / SpaceFlow / ISORT /
# Stereopy-TGPI / stLearn。
#
# **SpaGCN 已经从这张表里移除了 —— 它现在真的在跑**（`try_spagcn`）。
# 剩下的每一条都把"为什么没用上"记下来，判据全部是实测的。
#
# ## 五种不同的"用不了"，不能混为一谈
#
#   r_package       R/Bioconductor 包，本仓库的 Python CI 里没有 rpy2 → 结构上装不了
#   not_on_pypi     真包不在 PyPI，只能从 GitHub 装（且通常还依赖 torch-geometric）
#   deps            PyPI 上有真包，但依赖链在当前 CI 上跑不动
#   needs_reference 包**装得上**，缺的是数据（如解卷积需要带标签的 scRNA 参考）
#   name_taken      **PyPI 上那个名字是另一个不相干的包** ← 最危险的一类
#
# `needs_reference` 是这一轮新加的：`cell2location` 上一轮被记成 `deps`
# （"GPU 导向，CPU CI 吃不消"），但真读它的接口就会发现**它不是装不上的问题** ——
# 它缺的是一份带细胞类型标签的参考，而本流水线配的是 marker 签名。
# 把"缺数据"写成"装不上"会让下一个人去折腾安装，方向完全错。
#
# ## `name_taken` 为什么单列一类
#
# `pip install <名字>` 会**成功**，装进来的是完全无关的东西。实测：
#
#   pip install edgeR     -> "Redirect Microsoft Edge to your preferred browser"
#   pip install ISORT     -> PyCQA/isort，Python 的 import 排序工具
#   pip install slingshot -> ElasticSearch 索引迁移
#   pip install sparkx    -> 高能物理的碰撞相对论运动学
#
# 这比"装不上"糟得多：装不上会立刻报错，装错了要到 import 或跑出结果才发现，
# 而那时"结果"可能已经在图上看着挺像回事了。所以这些名字**不能出现在
# requirements.txt 里**，哪怕文档点名了。
#
# 反例：`SingleR` 在 PyPI 上是 **BiocPy/singler**，是 R 那个算法的官方
# Python 绑定（作者 Aaron Lun）—— 名字对得上，不是顶名的。所以判断依据是
# **summary / author / project_urls**，不是"名字存在与否"。
NAMED_TOOLS = {
    # ---- §3.2 空间域 --------------------------------------------------------
    # **SpaGCN 已经不在这张表里了。** 它上一轮被登记为 `deps`（理由是
    # "依赖 louvain，而 louvain 没有 py3.12 wheel"）—— 那个理由是
    # **读元数据推出来的，读源码就能推翻**：SpaGCN 1.2.7 的源码里没有一处
    # `import louvain`，它走 `scanpy.tl.louvain`，而 `init="kmeans"` 分支
    # 根本不碰它。现在 `pip install --no-deps SpaGCN` + `init="kmeans"`
    # 真的在跑（见 `03_spatial_domains.py` 的 `try_spagcn`）。
    # **留在表里会让 `probe_named_tools` 每轮都报"登记过期"。**
    "BayesSpace": dict(
        kind="r_package", section="§3.2",
        reason=("Bioconductor R 包（`BayesSpace`），PyPI 上无同名包"
                "（实测 `BayesSpace` 404）；本仓库 CI 不装 R + rpy2，"
                "结构上跑不了"),
    ),
    "STAGATE": dict(
        kind="not_on_pypi", section="§3.2",
        reason=("PyPI 上 `STAGATE` / `STAGATE_pyG` / `stagate` 三个名字"
                "**全部 404**（实测 pypi.org/pypi/<name>/json）。官方只发 GitHub"
                "（`RucDongLab/STAGATE_pyG`），而它的 `setup.py` 里 "
                "`install_requires = [\"requests\"]` —— **元数据完全没写真实依赖**："
                "`STAGATE_pyG/gat_conv.py:10` 是模块级 "
                "`from torch_sparse import SparseTensor, set_diag`。"
                "`torch-sparse` 0.6.18 在 PyPI 上**只有 sdist、0 个 wheel**"
                "（要按 torch 版本编译 C++ 扩展），而且 `gat_conv.py` 还用 "
                "`from torch_geometric.typing import OptPairTensor, Adj, Size, NoneType` "
                "—— 这些名字在 PyG>=2.4 已移除，等于同时钉死 torch-sparse 与 "
                "PyG<2.4。**两条都要满足才跑得起来，代价远大于它作为交叉验证的价值**"),
    ),
    # ---- §3.3 解卷积 --------------------------------------------------------
    "RCTD": dict(
        kind="r_package", section="§3.3",
        reason="`spacexr` 是 R 包（Bioconductor/GitHub dmcable/RCTD），PyPI 上无同名包",
    ),
    "cell2location": dict(
        kind="needs_reference", section="§3.3",
        reason=("PyPI 有真包（0.1.5，BayraktarLab），依赖 scvi-tools + torch + "
                "pyro-ppl + opencv-python —— **都装得上**（torch 本来就要为 "
                "SpaGCN 装）。**但装得上不等于跑得了**：cell2location 的 "
                "`RegressionModel` 需要一份**带细胞类型标签的 scRNA-seq 参考**"
                "（`cell_state_df`），而本流水线的配置是 "
                "`deconvolution.reference: builtin`（marker 签名）—— "
                "没有匹配的参考，硬跑只会得到一个不含信息的均匀组成"
                "（这正是 `05_deconvolution.py` 里 NNLS 那一版踩过的坑）。"
                "所以状态是 **needs_reference，不是装不上**；"
                "配上 `reference: h5ad` + `celltype_key` 就会走它"),
    ),
    # ---- §3.4 空间可变基因 --------------------------------------------------
    "SpatialDE2": dict(
        kind="not_on_pypi", section="§3.4",
        reason=("PyPI 上 **404**（实测 `https://pypi.org/pypi/SpatialDE2/json`）—— "
                "它是 SpatialDE 的继任实现，没有独立发行版，"
                "代码在 GitHub 仓库里（`Teichlab/SpatialDE` 的 `SpatialDE2` 分支/"
                "包内子模块），只能从源码装。本仓库用的是 **SpatialDE 1.1.3**"
                "（PyPI 有 wheel），所以 §3.4 的落地工具是前者不是后者"),
    ),
    "SPARK-X": dict(
        kind="r_package", section="§3.4",
        reason=("SPARK-X 是 R 包（xzhoulab/SPARK），PyPI 上无同名包。"
                "**注意 `pip install sparkx` 会成功但装错东西** —— PyPI 上的 "
                "`sparkx` 2.2.0 是高能物理的碰撞相对论运动学包，"
                "与空间可变基因毫无关系。所以这个名字不能进 requirements.txt"),
    ),
    # ---- §3.5 空间通讯 ------------------------------------------------------
    "CellChat": dict(
        kind="r_package", section="§3.5",
        reason="R 包（GitHub JinmiaoChenLab/CellChat），PyPI 上无同名包",
    ),
    # ---- §3.6 空间轨迹 ------------------------------------------------------
    "StPedf": dict(
        kind="not_on_pypi", section="§3.6",
        reason="PyPI 上无此包（`StPedf` / `stpedf` 都 404）",
    ),
    "SpaceFlow": dict(
        kind="deps", section="§3.6",
        reason=("PyPI 有真包（1.0.4，wheel 是 py3-none-any），但依赖 "
                "torch-geometric + torch-sparse + torch-scatter —— 后两者在 "
                "PyPI 上只有 sdist（0 个 wheel），要按 torch 版本编译，"
                "与 STAGATE 撞在同一道墙上"),
    ),
    "ISORT": dict(
        kind="name_taken", section="§3.6",
        reason=("**PyPI 上的 `ISORT` 是 PyCQA/isort —— Python 的 import 排序工具，"
                "和空间轨迹的 ISORT 毫无关系。** 真 ISORT 无 PyPI 发行版"),
    ),
    "Stereopy-TGPI": dict(
        kind="deps", section="§3.6",
        reason=("PyPI 有 `Stereopy` 1.6.2，但 `requires_python = '<3.9,>=3.8'` —— "
                "**与 CI 的 3.12 不兼容**（是上限卡死，不是下限）"),
    ),
    "stLearn": dict(
        kind="deps", section="§3.5/§3.6",
        reason=("PyPI 有真包（1.4.1，`requires_python='>=3.12'`），但依赖 "
                "`numpy>=2.4.0` / `scipy>=1.17.0` / `scanpy>=1.12.0` / "
                "`zarr>=3.1` / spatialdata 全家桶 + torch + torchvision + "
                "geopandas + dask（20+ 个），与本仓库钉的 scanpy/numpy 区间冲突"),
    ),
    # ---- 细胞分割（§3.1 的平台分支）------------------------------------------
    "Bering": dict(
        kind="deps", section="§3.1",
        reason="PyPI 有真包（0.1.2，KANG-BIOINFO/Bering），但依赖 torch + torch-geometric",
    ),
    "BOMS": dict(
        kind="deps", section="§3.1",
        reason=("PyPI 有真包（1.1.0），但 wheel 只到 cp310（无 cp312），"
                "依赖 `mkl` + `mkl-service` —— conda 时代的 Intel MKL 绑定，"
                "pip 环境下不可靠"),
    ),
}


def pkg_version(name: str):
    """单个发行版的版本号；查不到返回 None（**不编造**）。

    与 `capture_versions` 的区别：那个是全量枚举 + 写清单，这个只回答
    "某个点名工具现在装的是哪个版本"，给状态 JSON 里的
    `domain_methods` / `named_tools` 用。

    名字按 PEP 503 归一化后再查 —— `STAGATE_pyG` 与 `stagate-pyg` 是同一个
    发行版，直接拿原名查会漏。
    """
    from importlib import metadata as _md
    for cand in (name, _norm_pkg(name), name.replace("_", "-").lower()):
        try:
            return _md.version(cand)
        except Exception:  # noqa: BLE001
            continue
    return None


def probe_named_tools(log=None, only=None) -> dict:
    """把 NAMED_TOOLS 整理成可写进状态 JSON 的登记表。

    **不尝试 import** —— 这一节的结论是"没装/装不了"，去 import 只会
    反复报同一个 ImportError。真正的可用性判断在 `importlib.util.find_spec`，
    这里只做一次，用来区分"登记说装不了，但实际上环境里有"（那说明
    登记过期了，值得报出来）。

    返回 {tool: {available, kind, section, reason}}。
    """
    import importlib.util

    # 名字 -> 真正会被 import 的模块名（不总是一样）
    import_name = {
        "SpaceFlow": "spaceflow", "stLearn": "stlearn",
        "cell2location": "cell2location", "Bering": "Bering", "BOMS": "boms",
        "ISORT": "isort", "STAGATE": "STAGATE", "BayesSpace": None,
        "RCTD": None, "CellChat": None, "StPedf": None, "SPARK-X": None,
        "SpatialDE2": None,
        "Stereopy-TGPI": "stereopy",
    }
    out = {}
    for tool, meta in NAMED_TOOLS.items():
        if only is not None and tool not in only:
            continue
        mod = import_name.get(tool)
        avail = False
        if mod:
            try:
                avail = importlib.util.find_spec(mod) is not None
            except (ImportError, ValueError):
                avail = False
        if avail and log and meta["kind"] != "needs_reference":
            # 登记说过不了、环境里却有 —— 登记过期了，必须报出来。
            # **`needs_reference` 不算"说过不了"**：那一类的意思正是
            # "包装得上，缺的是数据"，所以它可 import 是符合预期的，
            # 报"登记过期"是误报。
            log_warn(f"工具 {tool} 登记为不可用，但环境里能 import —— 登记需要更新")
        out[tool] = {
            "available": bool(avail),
            "kind": meta["kind"],
            "section": meta["section"],
            "reason": meta["reason"],
        }
    return out


def named_tools_note() -> str:
    """一句话说明文档点名的工具里哪些没用上、为什么。"""
    return ("文档 §3.2–§3.6 点名的工具**大部分没有用上**："
            "要么是 R/Bioconductor 包（CI 无 rpy2），要么不在 PyPI，"
            "要么依赖链（torch-sparse / torch-scatter / 旧版 mkl）"
            "在当前 CPU CI 上跑不动，要么缺的是数据而不是包"
            "（cell2location 需要带标签的 scRNA 参考），"
            "要么 **PyPI 上那个名字是另一个不相干的包**。"
            "**两个例外，它们真的跑了：**"
            "SpaGCN（§3.2，`init=\"kmeans\"` 绕开了 louvain 的 py3.12 编译链，"
            "状态在 `domain_status.json` 的 `domain_methods` 里）和 "
            "SpatialDE（§3.4，带 scipy `misc.derivative` 垫片，"
            "状态在 `svg_status.json` 的 `spatialde` 里）—— "
            "两者都不在缺口登记里，因为缺口登记记的是**用不了**的工具。"
            "逐条理由见各步状态 JSON 的 `named_tools` / `domain_methods` 字段。")



def manifest_path(cfg: dict) -> Path:
    return Path(cfg["output"]["results_dir"]) / MANIFEST_NAME


def read_manifest(cfg: dict) -> dict:
    return read_json(manifest_path(cfg)) or {}


def init_manifest(cfg: dict, language: str = "python") -> dict:
    """建立本轮清单骨架。**会清掉上一轮的内容** —— 清单描述的是本轮。"""
    m = {
        "dataset_id": cfg.get("dataset_id"),
        "language": language,
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "seed": (cfg.get("analysis") or {}).get("seed"),
        "versions": {},
        "key_versions": {},
        "inputs": [],
        "params": {},
        "decisions": [],
        "human_review": [],
        "cross_language": [],
    }
    write_json(manifest_path(cfg), m)
    return m


def _manifest_append(cfg: dict, key: str, entry) -> None:
    m = read_manifest(cfg)
    m.setdefault(key, [])
    m[key].append(entry)
    m["dataset_id"] = cfg.get("dataset_id")
    m["updated_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    write_json(manifest_path(cfg), m)


def _norm_pkg(name: str) -> str:
    """PEP 503 归一化：包名大小写与 -/_/. 不敏感。"""
    return str(name).strip().lower().replace("_", "-").replace(".", "-")


def capture_versions(cfg: dict, key_packages=None, extra: dict = None) -> dict:
    """§0.3 版本记录。

    全量走 `importlib.metadata` 枚举已安装发行版，**不起子进程** ——
    管道捕获输出在受限沙箱里会 EPERM，而这里拿到的信息与 `pip freeze` 等价。

    文档点名的关键工具单独放进 `key_versions`：全量 freeze 有几百行，
    关键工具淹没在里面。
    """
    from importlib import metadata as _md

    full = {}
    try:
        for d in _md.distributions():
            try:
                n = d.metadata["Name"]
            except Exception:
                continue
            if n:
                full[_norm_pkg(n)] = d.version
    except Exception as exc:
        log_warn(f"枚举已安装包失败（{exc}）—— versions 会不完整")

    key = {}
    for p in (key_packages if key_packages is not None else KEY_PACKAGES):
        key[p] = full.get(_norm_pkg(p))
    for k, v in (extra or {}).items():
        key[k] = v

    m = read_manifest(cfg)
    m["versions"] = dict(sorted(full.items()))
    m["key_versions"] = key
    m["n_packages"] = len(full)
    m["python"] = sys.version.split()[0]
    try:
        import platform
        m["platform"] = platform.platform()
    except Exception:
        pass
    m["dataset_id"] = cfg.get("dataset_id")
    m["updated_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    write_json(manifest_path(cfg), m)

    missing = sorted(k for k, v in key.items() if v is None)
    if missing:
        log_warn(f"关键工具未安装（{len(missing)}/{len(key)}）：{', '.join(missing)}")
    else:
        log_info(f"关键工具全部就位（{len(key)} 个），共记录 {len(full)} 个已安装包")
    return key


def record_input(cfg: dict, path, label: str = "", required: bool = True) -> dict:
    """§0.4 输入数据哈希。

    文件不存在时**记 missing 而不是抛异常** —— 调用点未必知道某个输入
    这轮会不会产生（可选步骤的产物就是）。`required=True` 时 missing
    会在验收里被看见。
    """
    p = Path(path)
    entry = {"label": label or p.name, "path": str(p), "required": bool(required)}
    if p.exists() and p.is_file():
        entry["sha256"] = sha256_file(p)
        entry["bytes"] = p.stat().st_size
        entry["status"] = "present"
    else:
        entry["status"] = "missing"
    _manifest_append(cfg, "inputs", entry)
    return entry


def record_params(cfg: dict, params: dict) -> None:
    """§0.4 关键参数完整记录（含随机种子）。"""
    m = read_manifest(cfg)
    m.setdefault("params", {})
    m["params"].update(params or {})
    m["seed"] = (cfg.get("analysis") or {}).get("seed")
    m["dataset_id"] = cfg.get("dataset_id")
    m["updated_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    write_json(manifest_path(cfg), m)


def record_decision(cfg: dict, node: str, question: str, answer, evidence: str = "") -> None:
    """§0.4 Agent 决策链：从原始问题到最终结论的每一步推理。

    `evidence` 要写**支持这个选择的实际数字**，不是"因为这是通行做法"。
    没有量化依据的决策也要记，但 evidence 就写"没有量化依据"。
    """
    _manifest_append(cfg, "decisions", {
        "node": node, "question": question, "answer": answer,
        "evidence": evidence, "at": time.strftime("%Y-%m-%dT%H:%M:%S"),
    })


def record_human_review(cfg: dict, node: str, required: bool = True,
                        status: str = "pending", note: str = "") -> None:
    """§0.4 人工干预记录。

    status：pending（需确认，未确认）/ confirmed / overridden（人推翻了
    自动结果，note 写改成什么）/ not_needed（本数据集不涉及）。

    **默认 pending 而不是 confirmed。** 自动化流水线不能替人签字 ——
    把未确认的节点默认记成已确认，等于把复核节点变成摆设。
    """
    _manifest_append(cfg, "human_review", {
        "node": node, "required": bool(required), "status": status,
        "note": note, "at": time.strftime("%Y-%m-%dT%H:%M:%S"),
    })


def record_cross_language(cfg: dict, src: str, dst: str, fmt: str,
                          before: dict = None, after: dict = None,
                          lost=None, tool: str = "", note: str = "") -> None:
    """§0.2 跨语言转换记录。

    文档要求记录转换前后维度、metadata 字段数、丢失字段清单。桥接工具限定
    zellkonverter / anndata2ri，**禁止 sceasy**（维护状态差、metadata 丢失
    风险高）。

    本仓库与姊妹仓库之间只走 CSV，所以正常路径下 `before`/`after` 是
    行列数与列名集合；真正发生对象级转换时才填 `tool`。
    """
    _manifest_append(cfg, "cross_language", {
        "src": src, "dst": dst, "format": fmt, "tool": tool,
        "before": before or {}, "after": after or {},
        "lost_fields": sorted(lost or []),
        "note": note, "at": time.strftime("%Y-%m-%dT%H:%M:%S"),
    })


def manifest_summary(cfg: dict) -> dict:
    """给验收用的一行摘要。

    **必需项缺失和可选项缺失要分开报。** 可选项本来就允许不存在，把它算进
    "缺失"会让每个没有该文件的数据集都判失败 —— 那是把"设计如此"当成
    "出错了"。两者都必须**可见**，但只有必需项判失败。
    """
    m = read_manifest(cfg)
    if not m:
        return {"present": False}
    ins = m.get("inputs") or []
    miss = [i for i in ins if i.get("status") == "missing"]
    return {
        "present": True,
        "n_versions": len(m.get("versions") or {}),
        "n_inputs": len(ins),
        "inputs_missing": sorted(i.get("label") or "?" for i in miss),
        "inputs_missing_required": sorted(
            i.get("label") or "?" for i in miss if i.get("required")
        ),
        "n_decisions": len(m.get("decisions") or []),
        "human_review_pending": sorted(
            h["node"] for h in (m.get("human_review") or [])
            if h.get("status") == "pending"
        ),
        "n_cross_language": len(m.get("cross_language") or []),
    }


# ============================================================================
# 出图样式与调色板
# ============================================================================
# 规范来源：scientific-agent-skills/skills/scientific-visualization
#   assets/publication.mplstyle 与 assets/color_palettes.py（K-Dense，MIT）。
# 本仓库把样式文件放在 assets/publication.mplstyle，两处偏离写在该文件头部。
#
# **颜色不能是唯一线索。** Okabe-Ito 只是把颜色本身做成色盲友好；
# 分类图上仍要加 marker / 线型 / 直接标注，否则灰度打印就全糊了。
PAL = {
    # Okabe-Ito（Wong, Nature Methods 8:441, 2011）里**能做实心标记的**几个。
    #
    # **原来的八个标准色里有两个用不了**，理由是实测的，不是偏好：
    #   * `yellow` #F0E442 —— 白底对比度只有 **1.32:1**，实心圆点在白底上
    #     几乎看不见（`tools/check_palette.mjs` 会拦下来）；
    #   * `sky_blue` #56B4E9 —— 与 `blue` 色相只差约 5°，远小于 15° 判据，
    #     同一张图上会被读成同一个颜色。
    #
    # 所以这两个**从 PAL 里删掉了**：留着等于暗示"可以用"。
    "orange": "#E69F00",
    "green": "#009E73",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "black": "#000000",
    # 语义别名 —— 脚本里用语义名，换配色时只改这里
    "primary": "#0072B2",     # 主序列（原 #2C7FB8）
    "highlight": "#D55E00",   # 阈值线 / 强调（原 #B2182B）
    "muted": "#999999",       # 次要参照（如随机基线）
    # 循环色扩到 10 色时补的（挑法见 PAL_CYCLE 的说明）
    "cyan": "#17BECF",
    "maroon": "#A50F15",
    "indigo": "#5B4FCF",
    "grey": "#A8A8A8",
}
# **分类循环色必须 >= 类别数，否则会撞色 —— 而撞色不报错。**
#
# 原来只有 5 色。前 5 个**保持原样**（blue / vermillion / green / purple /
# black）—— 这样 <=5 个类别的图外观不变，改动的影响面最小。
#
# 后 5 个是**按 `tools/check_palette.mjs` 的三条判据挑出来的，不是凭眼睛选的**：
# 两两色相 >=15°、三种色盲（protanopia / deuteranopia / tritanopia）下
# OKLab 距离 >=0.05、白底对比度 >=2.0。
#
# 挑的过程值得记下来，因为**它证明了色空间是饱和的**：试了 13 个候选色，
# 12 个都被拦下，而且每个只差一项 ——
#
# | 候选 | 被什么拦下 |
# |---|---|
# | `#3D3D00` 暗橄榄 | protanopia 下与 maroon 撞（Δ=0.018）|
# | `#7F7F7F` 中灰 | deuteranopia 下与 green 撞（Δ=0.033）|
# | `#1F5C3A` 暗绿 | 色相与 green 只差 9.9° |
# | `#4A2A00` 暗棕 | 色相与 orange 只差 10.1° |
# | `#008080` 青绿 | tritanopia 下与 blue 撞，且色相与 cyan 只差 11.7° |
# | `#9467BD` 紫罗兰 | protanopia 下与 blue 撞（Δ=0.019）|
# | `#8C564B` 棕 | 色相与 maroon 只差 5.2° |
# | `#D3D3D3` 浅灰 | 对比度 1.50:1 |
# | `#B4B4B4` / `#B0B0B0` / `#ADADAD` | protanopia 下与 cyan 撞 |
# | `#A0A0A0` | deuteranopia 下与 purple 撞 |
#
# **规律：红绿色盲把 20°–110° 的暖色区压成一条轴，只剩亮度能区分。**
# 暖色区已有三个亮度级（maroon 0.460 / vermillion 0.621 / orange 0.753），
# 第 4 个暖色无论放哪个亮度都会撞上其中之一。
# 所以最后补的是**无彩色**（grey）：它不受色相判据约束，
# 且亮度 0.72 与所有彩色都拉得开。**grey 排最后** ——
# 它最不显眼，让它承担第 10 个类别而不是第 8 个。
#
# **本仓库的域图不走这条循环色。** `03_spatial_domains.py` /
# `05_deconvolution.py` / `06_niche.py` 都用 `plt.get_cmap("tab20")`
# 显式取色（域有 13 个，远超 10）。所以扩这条循环色**不会改变域图的颜色**，
# 它管的是那些用默认循环的小类别图。
PAL_CYCLE = [PAL["blue"], PAL["vermillion"], PAL["green"], PAL["purple"],
             PAL["black"], PAL["orange"], PAL["cyan"], PAL["maroon"],
             PAL["indigo"], PAL["grey"]]

_STYLE_APPLIED = False


def apply_style(cfg: dict = None) -> None:
    """
    应用出版级样式。**幂等**，重复调用无副作用。

    样式文件在 assets/publication.mplstyle。找不到时退回手工设几个
    关键 rcParam 并警告 —— 静默用 matplotlib 默认样式会让图看起来
    "能出"但不符合任何投稿规范。
    """
    global _STYLE_APPLIED
    if _STYLE_APPLIED:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    style = Path(__file__).resolve().parent.parent.parent / "assets" / "publication.mplstyle"
    if style.exists():
        plt.style.use(str(style))
    else:
        log_warn(f"找不到样式文件 {style}，退回手工设置（图不符合投稿规范）")
        matplotlib.rcParams.update({
            "figure.constrained_layout.use": True,
            "savefig.bbox": "standard",
            "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
            "axes.spines.top": False, "axes.spines.right": False,
            "pdf.fonttype": 42, "ps.fonttype": 42,
        })
    matplotlib.rcParams["axes.prop_cycle"] = matplotlib.cycler(color=PAL_CYCLE)
    matplotlib.rcParams["image.cmap"] = "viridis"
    _STYLE_APPLIED = True


# ============================================================================
# 出图
# ============================================================================
def mm(*vals: float):
    """
    毫米 → 英寸。**投稿图的尺寸单位是毫米，不是英寸。**

    常用宽度（Nature 的规范，见参考 skill 的 journal_requirements.md）：
    89 mm 单栏，183 mm 双栏，120-136 mm 单栏半。
    直接用英寸写 figsize 会让"单栏图"到底多宽变成一个没人检查的猜测。
    """
    if len(vals) == 1:
        return vals[0] / 25.4
    return tuple(v / 25.4 for v in vals)


# 三种标准宽度，单位英寸。绘图脚本一律用它们，不写裸英寸。
W_SINGLE = mm(89)      # 单栏
W_ONE_HALF = mm(136)   # 单栏半（Nature 允许 120-136 mm）
W_DOUBLE = mm(183)     # 双栏（= 满版宽）


def _content_overflow(fig) -> dict:
    """
    检查内容有没有超出画布（= 被裁掉）。

    为什么需要这个检查：`savefig.bbox` 从 "tight" 改成 "standard" 之后，
    装不下的标签会被**直接裁掉**，而图文件照样生成、`check_figures.mjs`
    照样报"有墨迹" —— 只有打开图才看得出来。

    constrained layout 正常情况下会把内容塞进画布；这个检查兜住
    "某个图用了 add_axes / 手工 GridSpec，constrained layout 管不到"的情况。
    """
    try:
        fig.canvas.draw()
        tb = fig.get_tightbbox(fig.canvas.get_renderer())
        if tb is None:
            return {}
        w, h = fig.get_size_inches()
        # 留 2% 容差：constrained layout 会把 pad 也算进去，少量溢出是正常的
        ow = (tb.width - w) / w
        oh = (tb.height - h) / h
        bad = {}
        if ow > 0.02:
            bad["width_overflow_frac"] = round(float(ow), 4)
        if oh > 0.02:
            bad["height_overflow_frac"] = round(float(oh), 4)
        return bad
    except Exception:  # noqa: BLE001
        return {}


def plot_marker_dotplot(ax, frac, zmat, *, cmap="RdBu_r",
                        size_max=170, size_min=10, dot_edge=0.3):
    """
    手工画 marker dotplot：颜色 = 按基因做 z-score（跨簇可比偏离方向），
    点大小 = 表达该基因的细胞比例。

    **为什么不再用 `sc.pl.dotplot(..., standard_scale="var")`。**
    读 scanpy 源码（`_prepare_dot_data`）确认：`standard_scale="var"` 是
    **逐基因 min-max 归一化到 0–1**（减最小值、除最大值），**不是 z-score**。
    而副标题写"z-scored" —— 两者矛盾（用户 2026-09-24 指出的核心问题）。
    两个口径的科学含义不同：

    ================  ===============================================
    口径              含义
    ================  ===============================================
    min-max（旧）     0 = 该基因在所有簇里的最低表达，1 = 最高
                      —— 只有"相对排名"，没有偏离方向
    z-score（本版）   0 = 平均水平，正 = 高于均值，负 = 低于均值
                      —— 跨簇可比"哪个簇偏离更大"，且 0 有锚点含义
    ================  ===============================================

    自带的三个额外好处（都是用户逐条点名的）：
      * 标签不再被裁 —— 布局由本函数控制，不与 scanpy 的 grid 抢空间；
      * 色标对称 —— RdBu_r 以 0 为中点，正负偏离等权可见；
      * 图例间距 —— 点大小图例单独画、间距显式给定，不再粘连。

    :param ax: 目标 axes
    :param frac: DataFrame，index=分组（行），columns=基因（列），值 0..1
    :param zmat: DataFrame，同形状，值 = 按基因 z-score 后的值
    :param cmap: 色标（默认 RdBu_r：红=高表达，蓝=低表达，白=0）
    :param size_max: 100% 表达时点的面积（pt^2）
    :param size_min: 0% 表达时点的面积（仍画一个极小点，表示"测到但几乎不表达"）
    :param dot_edge: 点描边宽度
    :returns: (ScalarMappable, size_handles) —— 供调用方画共享色标与大小图例
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    genes = list(zmat.columns)
    groups = list(zmat.index)
    ng, ngrp = len(genes), len(groups)

    # **z 轴对称归一**：正负偏离等权，0 永远是白色 —— 这就是"色标含负值"
    # 且"白色 = 平均水平（不是 0 表达）"的来源。上限取全部 |z| 的 95 分位
    # 再放宽到 >=1.5，避免个别极端基因把色标撑得两极分化。
    vmax = max(1.5, float(np.nanpercentile(np.abs(zmat.to_numpy()), 95)))
    norm = Normalize(vmin=-vmax, vmax=vmax)
    cmap_obj = plt.get_cmap(cmap)

    xs = np.arange(ng)
    ys = np.arange(ngrp)
    for gi in range(ng):
        for ri in range(ngrp):
            f = float(frac.iloc[ri, gi])
            z = float(zmat.iloc[ri, gi])
            s = size_min + f * (size_max - size_min)
            ax.scatter(xs[gi], ys[ri], s=s,
                       color=cmap_obj(norm(z)),
                       edgecolor="black", linewidth=dot_edge, zorder=3)

    ax.set_xticks(xs)
    ax.set_xticklabels(genes, rotation=90, fontsize=7)
    ax.set_yticks(ys)
    ax.set_yticklabels(groups, fontsize=8)
    ax.set_xlim(-0.7, ng - 0.3)
    ax.set_ylim(ngrp - 0.5, -0.5)
    ax.tick_params(length=2, pad=2)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap_obj)
    # 大小图例句柄：25/50/75/100% 四档（用户第七轮：五档太挤）
    size_handles = [(f, size_min + f * (size_max - size_min)) for f in
                    (0.25, 0.5, 0.75, 1.0)]
    return sm, size_handles


def build_marker_dotplot_figure(frac_df, z_df, *, group_label, title,
                                subtitle, fig_width=None, fig_height_mm=126):
    """
    构建**整张** marker dotplot（主图 + 底部横带图例），返回 `(fig, size_handles)`。

    **为什么把它抽成函数（这是本图第八轮才做对的关键）。** 前七轮图例布局的
    代码是**内联在调用脚本里**的，而验证脚本又**照抄了一份** —— 同一段布局
    存在三份镜像（scrna 脚本 / spatial 脚本 / 本地验证脚本）。于是：

    * 改一处、另两处不同步 → 验证脚本测的不是生产代码，**"全 PASS"是假的**；
    * 每次微调都要改三处，而 matplotlib 的位置参数**看不到效果**，
      只能靠一轮轮烧 CI 去猜 —— 实测调了七轮、每轮都要你重新看图。

    抽成一个函数后，**脚本与验证脚本调用的是同一份代码**，改一次全体生效。

    **图例为什么放底部横带。** 用户提供的参考代码
    （`可参考代码/99.单细胞：自动注释`）用 `scCustomize::do_DotPlot(dot.scale=12)`
    —— Seurat 生态的标准范式：图例在主图**下方一条横带**，左半点大小、
    右半横向色标。前七轮把两块图例竖着塞进右侧 1/6 宽的窄列，两块图例加
    两个标题挤在 0.4 图高的竖条里，**没有不受挤的排法**（实测调四轮仍相撞）。
    改横带后主图横向延展到全宽，两个标题各在自己图例正上方，与点列物理分离。

    **必须 `set_layout_engine("none")`。** 本仓库全局开着
    `figure.constrained_layout.use: True`（AGENTS 规则 13），而它会**在每次
    draw 时重排手动设的坐标** —— 这正是前几轮"改了没效果、底部留白越调越多"
    的原因：`fig.add_axes([...])` 设的位置被 constrained layout 覆盖了。
    这张图的几何由本函数全权控制，所以显式关掉它。

    :param frac_df: DataFrame，index=分组，columns=基因，值 = 表达细胞比例
    :param z_df: DataFrame，同形状，值 = 按基因 z-score 后的表达
    :param group_label: y 轴标签（"Leiden cluster" / "Spatial domain"）
    :param title: suptitle（如 "Top markers per cluster"）
    :param subtitle: 主图标题第二行（轴语义说明）
    :param fig_width: 图宽（默认 W_DOUBLE）
    :param fig_height_mm: 图高（毫米）
    :returns: `(fig, size_handles)`
    """
    import matplotlib.pyplot as plt

    w = W_DOUBLE if fig_width is None else fig_width
    fig = plt.figure(figsize=(w, mm(fig_height_mm)))
    # **关掉 constrained layout** —— 见 docstring：它会覆盖下面所有手动坐标
    fig.set_layout_engine("none")

    # 主图占满上方（左右各留一点给 y 轴标签与右缘）
    ax = fig.add_axes([0.075, 0.315, 0.905, 0.535])
    sm, size_handles = plot_marker_dotplot(ax, frac_df, z_df)
    ax.set_xlabel("gene", fontsize=8)
    ax.set_ylabel(group_label, fontsize=8)
    fig.suptitle(title, fontsize=11, y=0.965)
    ax.set_title(subtitle, fontsize=8, pad=6, loc="left")

    # ---- 底部横带 · 左半：Percent Expressed (%)（四点横排）---------------
    lax = fig.add_axes([0.075, 0.045, 0.42, 0.13])
    lax.set_xlim(0, 1)
    lax.set_ylim(0, 1)
    lax.axis("off")
    lax.text(0.0, 0.88, "Percent Expressed (%)", ha="left", va="center",
             fontsize=7.5)
    for k, (f_, s_) in enumerate(size_handles):
        xx = 0.06 + k * 0.20
        lax.scatter([xx], [0.42], s=s_, color="gray",
                    edgecolor="black", linewidth=0.3)
        lax.text(xx, 0.02, f"{int(f_ * 100)}", ha="center", va="center",
                 fontsize=7.5)

    # ---- 底部横带 · 右半：Mean Expression（横向色标）--------------------
    cax = fig.add_axes([0.60, 0.085, 0.28, 0.035])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cax.text(0.0, 1.9, "Mean Expression", transform=cax.transAxes,
             ha="left", va="bottom", fontsize=7.5)
    cb.set_ticks([-1, 0, 1])
    cb.set_ticklabels(["Low", "Mid", "High"])
    cb.ax.tick_params(labelsize=7, top=False, bottom=True,
                      labeltop=False, labelbottom=True)
    return fig, size_handles



def fix_dotplot_legends(fig, size_title=None, cbar_title=None):
    """
    把 `sc.pl.dotplot` 的**整条图例列**整理成约定 v2 的样子：
    点大小图例竖排、色标竖排、两者上下排列互不重叠。

    **为什么必须后处理。** scanpy 的 `DotPlot` 没有暴露图例方向参数
    （`legend()` 只收 `width` / `show_size_legend` / `colorbar_title`）：

    * `_plot_size_legend()` 把示例点画在 **x 轴**上 —— 横排；
    * `_plot_colorbar()` 把色标**硬编码** `orientation="horizontal"` —— 横排。

    两个都要转竖排（用户约定 v2"纵向单列节约图幅"），而且**必须一起做**：
    只转一个，另一个还横着占着原来的宽度，两块会互相挤压或重叠
    （实测只转 size 图例时，colorbar 与它文字交叠 4 处）。

    **三处实测毛病，分别对应三个动作：**

    1. **size 图例的刻度标签被换行堆成两列**（`100/80/60` 挤在一起）——
       scanpy 给这块 axes 的宽度是按"横排一行点"算的，竖排后标签要单独占
       左侧一列，原宽度不够。→ 加宽 axes，并给刻度标签留出明确宽度。
    2. **colorbar 横向**。→ 找到 Colorbar 对象，竖向重建。
    3. **两块图例重叠**。→ 上下重新分区：size 在上、colorbar 在下，
       各自 `set_position` 不相交。

    :param fig: dotplot 所在的 figure
    :param size_title: 点大小图例标题（不传则保留原样）
    :param cbar_title: 色标标题（不传则保留原样）
    :returns: `dict(size=bool, colorbar=bool)` —— 各自是否成功转换
    """
    import numpy as np
    from matplotlib.axes import Axes
    from matplotlib.colorbar import Colorbar
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    done = {"size": False, "colorbar": False}

    # ---- 0. 定图例列的 x 位置 -------------------------------------------
    # **必须放在主图右侧、画布内**。主图（含 y 轴标签的那个）右缘一般在
    # x≈0.70（左边留给行标签）；图例列取主图右缘 + 一点间距。
    # 实测踩过：直接写 `1.0 - 宽度` 会把两块图例**推出画布右缘**——
    # size 刻度标签 6 个被裁、与色标文字重叠 5 处。
    main_ax = None
    for ax in fig.axes:
        if any(t.get_text() for t in ax.get_xticklabels()) \
                and ax.get_position().width > 0.4:
            main_ax = ax
            break
    leg_x = 0.80 if main_ax is None else min(0.94, main_ax.get_position().x1 + 0.055)

    # ---- 1. 点大小图例：横排 -> 纵排 ------------------------------------
    for ax in fig.axes:
        if done["size"]:
            break
        # 大小图例的判据：有 x 刻度标签、**没有** y 刻度标签、含散点
        if ax.get_yticklabels() and any(t.get_text() for t in ax.get_yticklabels()):
            continue
        colls = [c for c in ax.collections if hasattr(c, "get_sizes")]
        if not colls:
            continue
        labels = [t.get_text() for t in ax.get_xticklabels()]
        if not labels or not all(labels):
            continue
        # 只认"看起来像百分比数字"的刻度，避免误伤其它图
        try:
            [float(s) for s in labels]
        except ValueError:
            continue
        sizes = colls[0].get_sizes()
        if len(sizes) != len(labels):
            continue

        n = len(sizes)
        pos = ax.get_position()
        # **竖排后这块 axes 要"又高又窄"**：高度按点数给，宽度给刻度标签留
        # 足够列宽 —— 原宽度是按横排算的，标签会换行堆叠（实测 `100/80/60`
        # 挤成两列）。x 位置用上面算好的 `leg_x`（主图右侧、画布内）。
        need_h = min(max(0.030 * n + 0.06, 0.20), 0.50)
        need_w = 0.055
        ax.set_position([leg_x, pos.y1 - need_h, need_w, need_h])

        ax.clear()
        ys = np.arange(n)
        # 点贴近轴右缘、刻度标签在其左 —— 一行读作"标签 — 点"。
        # （标签不能放右侧：图例列在最右，再往右就出画布。）
        ax.scatter(np.full(n, 0.62), ys, s=sizes, color="gray",
                   edgecolor="black", linewidth=0.5, zorder=100)
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(-0.8, n - 0.2)
        ax.set_yticks(ys)
        ax.set_yticklabels(labels, fontsize="small")
        ax.set_xticks([])
        ax.tick_params(axis="x", bottom=False, labelbottom=False)
        ax.tick_params(axis="y", left=False, labelleft=True, pad=1)
        for sp in ax.spines.values():
            sp.set_visible(False)
        if size_title:
            ax.set_title(size_title, fontsize="small", pad=4, loc="left")
        done["size"] = True

    # ---- 2. 色标：横向 -> 纵向 -------------------------------------------
    # **先 draw 一次**（scanpy 是绘制时才把色标挂上去的，不 draw 找不到）。
    #
    # **色标轴的判据是"含 QuadMesh"，不是"含 Colorbar 对象"。** 实测踩过：
    # `Colorbar` 实例**不在** `ax.get_children()` 里（它挂在 figure 上），
    # 轴里能看到的只有色标本体的 `QuadMesh` 和边框 `_ColorbarSpine`。
    # 按类名找 Colorbar 永远返回 False —— helper 空转、图保持横向。
    fig.canvas.draw()
    src_cbar = None
    for ax in fig.axes:
        kinds = {type(c).__name__ for c in ax.get_children()}
        if "QuadMesh" in kinds and any(n.startswith("_ColorbarSpine") for n in kinds):
            src_cbar = ax
            break

    if src_cbar is not None:
        cax = src_cbar
        # 色标的 norm/cmap 要从**产生它的 mappable** 取。轴的 children 里只有
        # QuadMesh 本身，而 QuadMesh 带着创建时的 norm/cmap —— 从它读回。
        qm = next(c for c in cax.get_children() if type(c).__name__ == "QuadMesh")
        norm = getattr(qm, "norm", None) or Normalize()
        cmap = getattr(qm, "cmap", None) or plt.get_cmap("Reds")
        # 横向色标的刻度在 x 轴上；竖向要挂到 y 轴
        tick_vals = [t for t in cax.get_xticks()]
        old_title = cbar_title
        if old_title is None:
            old_title = cax.get_title()

        cax.clear()
        spos = cax.get_position()
        # **竖向色标要"窄而高"**，放在 size 图例正下方、同一列（`leg_x`）
        cax.set_position([leg_x + 0.012, spos.y0, 0.030, min(0.16, spos.y0)])
        sm = ScalarMappable(norm=norm, cmap=cmap)
        cb_new = Colorbar(cax, mappable=sm, orientation="vertical")
        # **刻度必须重新算，不能照抄横向时的 x 刻度。** 两个原因：
        # ① 横向时 x 轴的数值范围是 norm 的全程，竖向 y 轴也一样，但
        #    `get_xticks()` 返回的是"当时渲染出的位置"，直接 set_ticks 会
        #    和竖轴的实际范围错位；
        # ② 实测照抄会出现**镜像 + 叠字**（从上往下 1.0→0.0，且两位小数
        #    挤在一起）。正确做法：按竖轴范围取 3 个等距点，`set_yticks`。
        lo, hi = norm.vmin, norm.vmax
        new_ticks = list(np.linspace(lo, hi, 3)) if np.isfinite([lo, hi]).all() else []
        if new_ticks:
            cb_new.set_ticks(new_ticks)
            cb_new.ax.set_yticklabels([f"{v:.1f}" for v in new_ticks],
                                      fontsize="small")
        cax.tick_params(labelsize="small")
        if old_title:
            cax.set_title(old_title, fontsize="small", pad=4, loc="left")
        done["colorbar"] = True

    # ---- 3. 收口：两块上下排列，不相交 ------------------------------------
    # set_position 之后 constrained layout 会在下一帧重排，这里强制立即执行
    # 一次并做最终夹紧，保证两块不交叠（用户明确要求"变了后也不能重叠"）。
    # 判据同上：色标轴看 QuadMesh，size 图例轴看"y 刻度是数字"。
    fig.canvas.draw()
    size_ax = cbar_ax = None
    for ax in fig.axes:
        kinds = {type(c).__name__ for c in ax.get_children()}
        labs = [t for t in ax.get_yticklabels() if t.get_text()]
        if "QuadMesh" in kinds:
            cbar_ax = ax
        elif labs:
            try:
                [float(t.get_text()) for t in labs]
                size_ax = ax
            except ValueError:
                pass
    if size_ax is not None and cbar_ax is not None:
        sp, cp = size_ax.get_position(), cbar_ax.get_position()
        top = max(sp.y1, cp.y1)
        gap = 0.02
        h_size, h_cbar = sp.height, cp.height
        if h_size + h_cbar + gap > top:
            scale = (top - gap) / (h_size + h_cbar)
            h_size *= scale
            h_cbar *= scale
        size_ax.set_position([sp.x0, top - h_size, sp.width, h_size])
        cbar_ax.set_position([cp.x0, top - h_size - gap - h_cbar, cp.width, h_cbar])
        fig.canvas.draw()

    return done



def place_labels(ax, xs, ys, texts, fontsize=7, pad_px=2.0,
                 max_iters=400, seed=0):
    """
    在散点图上放**互不重叠**的文字标签（纯几何，不引第三方依赖）。

    **为什么需要它。** `adjustText` 不在本仓库依赖里（规则：不要临时引入），
    而"固定偏移 + 上下交替"这类纯参数法在**点挤成一条竖列**时必然失效 ——
    实测 `02-07-03-unit1-tf-specificity-scatter`：8 个 TF 的 x 都在 0–0.05，
    上下交替只把标签分成两层，同层内仍然互压（用户反馈"基因标签有重叠"）。

    做法：把每个标签候选位置按"离锚点由近及远"排序（右、左、上、下、
    四个对角共 8 个方向 × 多圈），**贪心**选第一个与已放标签不重叠的位置；
    都放不下就继续外扩。重叠判据用 matplotlib 的文本包围盒（渲染后实测，
    不是估算字宽）。

    **保证**：只要画布还有空间，返回的标签两两不重叠；实在放不下的会被
    推远（有引导线时仍可读）。返回实际使用的位置列表，便于日志记录。

    :param ax: 目标 axes（需已完成 scatter 且坐标范围已定）
    :param xs: 锚点 x（数据坐标）
    :param ys: 锚点 y（数据坐标）
    :param texts: 标签文字，与 xs/ys 等长
    :param fontsize: 标签字号（pt）
    :param pad_px: 标签之间要求的最小间隙（像素）
    :param max_iters: 每个标签最多尝试的候选位置数
    :param seed: 保留参数（本函数确定性，不用随机；留给调用方对齐口径）
    :returns: `[(x, y), ...]` 实际放置的**数据坐标**
    """
    import itertools

    fig = ax.figure
    # **必须先把 figure 的 dpi 对齐到"实际保存用的 dpi"再量包围盒。**
    #
    # 实测（2026-09-24）：本函数在创建 figure 时的默认 dpi（100）下量尺寸并
    # 排布，而 `save_fig()` 用 `dpi=300` 重新渲染 —— 文本的**字号是点、与 dpi
    # 无关，但包围盒是像素**，于是 300 dpi 下每个标签都放大 3 倍、而按点给的
    # 偏移没跟着放大，布局全散。实测 dpi=100 时 0 重叠、dpi=300 时
    # **28 对重叠 + 8 个越界** —— 正是 CI 图上"仍重叠"的原因。
    #
    # 所以：先设 dpi（与 `save_fig` 一致），再 draw、再量、再排。
    target_dpi = 300.0
    fig.set_dpi(target_dpi)
    fig.canvas.draw()  # 需要 renderer 才能量文本包围盒
    renderer = fig.canvas.get_renderer()

    # 候选方向：8 个方位，由近及远多圈外扩。
    # **优先向右/左上** —— 这类散点的点在左侧挤成竖列（x≈0），右侧是空的。
    dirs = [(1, 0), (1, 1), (1, -1), (0, 1), (0, -1),
            (-1, 1), (-1, -1), (-1, 0)]
    placed_boxes = []
    used = []
    # **标签必须留在坐标轴内。** 只判"标签之间不重叠"是不够的 ——
    # 实测（run 35974071689）：点在左上角挤成一团时，贪心把 8 个标签
    # 一路往右上推，**全部推出画布顶边**、还压住标题。
    # 所以加一条包含判据：标签包围盒必须完全落在 axes 内。
    ax_bb = ax.get_window_extent(renderer=renderer)

    def _mk(text, x, y, off_pt, dx, dy, arrow=False):
        """建一个标签。`off_pt` 是**以点为单位的偏移** —— 不是像素。"""
        kw = {}
        if arrow:
            kw["arrowprops"] = dict(arrowstyle="-", lw=0.4,
                                    color=PAL.get("muted", "#999999"))
        return ax.annotate(
            text, (x, y), xytext=off_pt, textcoords="offset points",
            fontsize=fontsize,
            ha="left" if dx > 0 else ("right" if dx < 0 else "center"),
            va="bottom" if dy > 0 else ("top" if dy < 0 else "center"),
            color=PAL.get("ink", "#1A1A1A"), **kw)

    for x, y, text in zip(xs, ys, texts):
        chosen = None
        best_outside = None   # 实在放不进时的次优（溢出最少的一个）
        for ring in itertools.count(1):
            if ring * len(dirs) > max_iters:
                break
            step_pt = 4.0 + (ring - 1) * 5.0   # 点（1/72 英寸），与 DPI 无关
            for dx, dy in dirs:
                # **`textcoords="offset points"` 要的是偏移量，不是绝对坐标。**
                # 实测踩过：这里原先传的是 `transData.transform()` 出来的
                # **像素绝对坐标**，于是标签被推到画布外几万像素处 ——
                # 包含判据全部拒绝，最后统统落进兜底分支挤在一起。
                off = (dx * step_pt, dy * step_pt)
                t = _mk(text, x, y, off, dx, dy)
                raw = t.get_window_extent(renderer=renderer)
                bb = raw.expanded(1 + pad_px / max(raw.width, 1),
                                  1 + pad_px / max(raw.height, 1))
                if any(bb.overlaps(o) for o in placed_boxes):
                    t.remove()
                    continue
                inside = (bb.x0 >= ax_bb.x0 and bb.x1 <= ax_bb.x1 and
                          bb.y0 >= ax_bb.y0 and bb.y1 <= ax_bb.y1)
                if not inside:
                    t.remove()
                    if best_outside is None:
                        over = (max(0, ax_bb.x0 - bb.x0) + max(0, bb.x1 - ax_bb.x1) +
                                max(0, ax_bb.y0 - bb.y0) + max(0, bb.y1 - ax_bb.y1))
                        best_outside = (over, off, dx, dy)
                    continue
                placed_boxes.append(bb)
                chosen = off
                used.append(off)
                break
            if chosen is not None:
                break
        if chosen is None:
            # **放不进就带引导线放最近处**，不静默丢弃、也不推到画布外。
            if best_outside is not None:
                _, off, dx, dy = best_outside
            else:
                off, dx, dy = (18.0, 0.0), 1, 0
            _mk(text, x, y, off, dx, dy, arrow=True)
            used.append(off)
    return used



def save_fig(cfg: dict, name: str, fig=None, tight: bool = False) -> list:
    """
    保存一张图为 PNG + PDF，并返回写出的路径。

    **PDF 与 PNG 都要出。** PDF 是矢量、可再编辑、字体以 Type 42 内嵌；
    PNG 用于快速查看与像素级非空白检查（check_figures 解 PNG）。

    **`bbox_inches` 默认不再是 "tight"。** 参考规范
    （scientific-visualization）明确写着 tight 会改变输出的物理尺寸 ——
    投稿要求"单栏 89 mm"时，tight 出来的就不是 89 mm。装不下由
    constrained layout 解决，另有 `_content_overflow()` 兜底并告警。
    """
    import matplotlib.pyplot as plt

    apply_style(cfg)
    figdir = Path(cfg["output"]["figures_dir"])
    figdir.mkdir(parents=True, exist_ok=True)
    f = fig if fig is not None else plt.gcf()

    bad = _content_overflow(f)
    if bad:
        log_warn(f"图 {name} 的内容超出画布（{bad}）—— 标签可能被裁掉。"
                 f"调大 figsize 或改用 layout='constrained'")

    written = []
    for ext in ("png", "pdf"):
        p = figdir / f"{name}.{ext}"
        f.savefig(p, dpi=cfg["analysis"]["figure_dpi"],
                  bbox_inches="tight" if tight else None)
        written.append(str(p))
    plt.close(f)
    return written


# ============================================================================
# 空间特有
# ============================================================================
def read_tissue_positions(path: Path):
    """
    读 Space Ranger 的 spot 位置文件，**兼容 v1 与 v2 两种格式**。

    v1: `tissue_positions_list.csv`，**无表头**，6 列：
        barcode, in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres
    v2: `tissue_positions.csv`，**有表头**，同名列

    踩过的坑：按无表头读 v2 文件，第一行数据会被当成表头吃掉 ——
    少一个 spot，而且**没有任何报错**。所以先探测第一行是不是表头。
    """
    import pandas as pd

    cols = ["barcode", "in_tissue", "array_row", "array_col",
            "pxl_row_in_fullres", "pxl_col_in_fullres"]
    if not path.exists():
        raise FileNotFoundError(f"找不到 spot 位置文件: {path}")

    with open(path, encoding="utf-8") as fh:
        first = fh.readline().strip()
    has_header = first.split(",")[0].strip().lower() in ("barcode", "barcodes")

    df = pd.read_csv(path, header=0 if has_header else None)
    if not has_header:
        if df.shape[1] != len(cols):
            raise ValueError(f"{path.name} 有 {df.shape[1]} 列，期望 {len(cols)} 列")
        df.columns = cols
    else:
        # v2 的列名可能有大小写差异
        df.columns = [str(c).strip().lower() for c in df.columns]
        missing = [c for c in cols if c not in df.columns]
        if missing:
            raise ValueError(f"{path.name} 缺少列: {missing}；实际列: {list(df.columns)}")
    return df, has_header


def load_scalefactors(path: Path) -> dict:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def spot_radius_plot_units(scalefactors: dict, which: str = "hires") -> float:
    """
    返回在**绘图坐标系**里 spot 的半径。

    这一步容易搞错：`spot_diameter_fullres` 是全分辨率图像上的直径，
    而绘图用的是 hires（或 lowres）图。必须乘对应的 scalef 才是绘图单位。
    不乘的话点的大小会离谱地大或小，而图看起来"只是有点怪"。
    """
    d = float(scalefactors["spot_diameter_fullres"])
    key = "tissue_hires_scalef" if which == "hires" else "tissue_lowres_scalef"
    return d * float(scalefactors[key]) / 2.0


def require_pkg(name: str, hint: str = "") -> None:
    try:
        __import__(name)
    except ImportError as e:
        extra = f"（{hint}）" if hint else ""
        raise RuntimeError(f"缺少依赖 {name}{extra}: {e}") from e


def df_to_records(df) -> list:
    import pandas as pd
    if df is None or len(df) == 0:
        return []
    return json.loads(df.where(pd.notna(df), None).to_json(orient="records"))


def load_registry(repo_root=None) -> dict:
    """
    读 assets/datasets.yml 的数据集注册表。

    **放在 common 里而不是 00_fetch 里。** 多数据集跑批时其他脚本也要
    能查注册表；定义在单个脚本内部，别人拿不到，只能复制一份。
    """
    import yaml

    root = (Path(repo_root) if repo_root
            else Path(__file__).resolve().parent.parent.parent)
    p = root / "assets" / "datasets.yml"
    if not p.exists():
        return {}
    with open(p, encoding="utf-8") as fh:
        return (yaml.safe_load(fh) or {}).get("datasets", {}) or {}


def spatial_xy(adata, scale: float = 1.0):
    """
    取 spot 的图像坐标，列序是 **(x, y) = (pxl_col_in_fullres, pxl_row_in_fullres)**。

    ---
    ## 为什么收敛成一个函数

    这个列序写反了**不会报错**：散点图仍然画出一个"看起来像组织形状"的
    点阵，只是转置了。只有把点叠到 H&E 图上、或者做定量的网格几何检查
    （见 `tools/verify_spatial_alignment.py`），才能发现不对。

    之前每个脚本各自写 `adata.obsm["spatial"]`，约定散落在 8 个文件里 ——
    任何一处写反都看不出来。收敛到这里之后：
      - 约定只在一个地方定义
      - `tools/check_py_syntax.mjs` 会拦住任何直接访问 `obsm["spatial"]` 的代码

    `scale` 用来把坐标从全分辨率像素换到某张图的尺度
    （hires 图要乘 `tissue_hires_scalef`）。

    **注意 matplotlib 的 imshow 用的是 (row, col) 顺序**，所以叠图时
    必须 `imshow(img)` 之后用 `scatter(x, y)`，而不是反过来。
    """
    import numpy as np

    xy = np.asarray(adata.obsm["spatial"])
    if xy.ndim != 2 or xy.shape[1] < 2:
        raise ValueError(f"obsm['spatial'] 形状异常: {xy.shape}")
    return xy[:, :2] * float(scale)


