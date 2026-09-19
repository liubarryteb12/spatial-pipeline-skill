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
# 出图
# ============================================================================
def save_fig(cfg: dict, name: str, fig=None, tight: bool = True) -> list:
    import matplotlib.pyplot as plt

    figdir = Path(cfg["output"]["figures_dir"])
    figdir.mkdir(parents=True, exist_ok=True)
    f = fig if fig is not None else plt.gcf()
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


