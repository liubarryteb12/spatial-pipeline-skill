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


