"""
lib/common.py — 空间转录组流水线公共库

与 `scrna-pipeline-skill` / `geo-brca-microarray-skill` 同一套约定。
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
    ana.setdefault("figure_dpi", 150)
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
# §3.5 点名 CellChat，§3.6 点名 StPedf / SpaceFlow / ISORT / Stereopy-TGPI / stLearn。
# **本仓库一个都没用上。** 这一节把"为什么"逐条记下来，判据全部是实测的。
#
# ## 三种不同的"用不了"，不能混为一谈
#
#   r_package    R/Bioconductor 包，本仓库的 Python CI 里没有 rpy2 → 结构上装不了
#   not_on_pypi  真包不在 PyPI，只能从 GitHub 装（且通常还依赖 torch-geometric）
#   deps         PyPI 上有真包，但依赖链在当前 CI 上跑不动
#   name_taken   **PyPI 上那个名字是另一个不相干的包** ← 最危险的一类
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
    "BayesSpace": dict(
        kind="r_package", section="§3.2",
        reason="Bioconductor R 包，没有 Python 发行版；本仓库 CI 不装 R + rpy2",
    ),
    "STAGATE": dict(
        kind="not_on_pypi", section="§3.2",
        reason=("PyPI 上 STAGATE / STAGATE_pyG / stagate 三个名字**全部 404**；"
                "官方只发 GitHub，且依赖 torch + torch-geometric"),
    ),
    "SpaGCN": dict(
        kind="deps", section="§3.2",
        reason=("PyPI 有真包（1.2.7，作者 Jian Hu），但依赖 `louvain` —— "
                "**该包最新版 0.8.2 没有 py3.12 wheel，只有 sdist**，"
                "要从 2019 年的 C++/Cython 源码编译。CI 的 3.12 上未验证能编过，"
                "不敢押一轮 CI"),
    ),
    # ---- §3.3 解卷积 --------------------------------------------------------
    "RCTD": dict(
        kind="r_package", section="§3.3",
        reason="`spacexr` 是 R 包（Bioconductor/GitHub），PyPI 上无同名包",
    ),
    "cell2location": dict(
        kind="deps", section="§3.3",
        reason=("PyPI 有真包（0.1.5，BayraktarLab），但依赖 "
                "`scvi-tools>=1.3.0` + `torch>=1.9.0` + `pyro-ppl` + `opencv-python`；"
                "GPU 导向，CPU CI 上的时间与磁盘都吃不消"),
    ),
    # ---- §3.5 空间通讯 ------------------------------------------------------
    "CellChat": dict(
        kind="r_package", section="§3.5",
        reason="R 包（GitHub JinmiaoChenLab/CellChat），PyPI 上无同名包",
    ),
    # ---- §3.6 空间轨迹 ------------------------------------------------------
    "StPedf": dict(
        kind="not_on_pypi", section="§3.6",
        reason="PyPI 上无此包",
    ),
    "SpaceFlow": dict(
        kind="deps", section="§3.6",
        reason=("PyPI 有真包（1.0.4），但依赖 torch-geometric + torch-sparse + "
                "torch-scatter —— 后两者要按 torch 版本编译，是出了名的难装"),
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
        reason=("PyPI 有真包（1.4.1），但依赖 `numpy>=2.4.0` / `scipy>=1.17.0` / "
                "`scanpy>=1.12.0` / `zarr>=3.1` / spatialdata 全家桶 + torch + "
                "torchvision + geopandas + dask（20+ 个），与本仓库钉的 "
                "scanpy/numpy 区间冲突"),
    ),
    # ---- 细胞分割（§3.1 的平台分支）------------------------------------------
    "Bering": dict(
        kind="deps", section="§3.1",
        reason="PyPI 有真包（0.1.2，KANG-BIOINFO/Bering），但依赖 torch + torch-geometric",
    ),
    "BOMS": dict(
        kind="deps", section="§3.1",
        reason=("PyPI 有真包（1.1.0），依赖 `mkl` + `mkl-service` —— "
                "那是 conda 时代的 Intel MKL 绑定，pip 环境下不可靠"),
    ),
}


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
        "SpaGCN": "SpaGCN", "SpaceFlow": "spaceflow", "stLearn": "stlearn",
        "cell2location": "cell2location", "Bering": "Bering", "BOMS": "boms",
        "ISORT": "isort", "STAGATE": "STAGATE", "BayesSpace": None,
        "RCTD": None, "CellChat": None, "StPedf": None,
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
        if avail and log:
            # 登记说过不了、环境里却有 —— 登记过期了，必须报出来
            log_warn(f"工具 {tool} 登记为不可用，但环境里能 import —— 登记需要更新")
        out[tool] = {
            "available": bool(avail),
            "kind": meta["kind"],
            "section": meta["section"],
            "reason": meta["reason"],
        }
    return out


def named_tools_note() -> str:
    """一句话说明本仓库为什么一个 §3.2/§3.3 点名工具都没用上。"""
    return ("文档 §3.2 / §3.3 / §3.5 / §3.6 点名的工具本仓库**一个都没用上**："
            "要么是 R/Bioconductor 包（CI 无 rpy2），要么不在 PyPI，"
            "要么依赖链（torch / scvi-tools / torch-geometric / 旧版 mkl）"
            "在当前 CPU CI 上跑不动，要么 **PyPI 上那个名字是另一个不相干的包**。"
            "逐条理由见 domain_status.json 的 named_tools 字段。"
            "**这是缺口，不是「已覆盖」。**")



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
    # Okabe-Ito（Wong, Nature Methods 8:441, 2011）八个标准色
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "black": "#000000",
    # 语义别名 —— 脚本里用语义名，换配色时只改这里
    "primary": "#0072B2",     # 主序列（原 #2C7FB8）
    "highlight": "#D55E00",   # 阈值线 / 强调（原 #B2182B）
    "muted": "#999999",       # 次要参照（如随机基线）
}
# 对白底达到 3:1 对比度的五个，用作分类循环色
PAL_CYCLE = [PAL["blue"], PAL["vermillion"], PAL["green"], PAL["purple"],
             PAL["black"]]

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


