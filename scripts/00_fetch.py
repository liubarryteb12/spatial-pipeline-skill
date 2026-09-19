#!/usr/bin/env python3
"""
00_fetch.py — 取 Visium 数据、组装 AnnData、硬校验

**Visium 的数据不是一个文件，是两半：**
  1. 表达矩阵（`*_filtered_feature_bc_matrix.h5`）
  2. 空间信息（`*_spatial.tar.gz`）—— 包含 spot 在组织切片上的像素坐标、
     缩放因子、以及 H&E 图像本身

**只有第一半的话，数据就退化成普通的单细胞数据了** —— 丢掉空间坐标
等于丢掉这个技术存在的理由。所以两半都必须拿到，缺任何一个都直接失败。

另外做硬校验：
  - 表达矩阵必须是原始整数计数（同 scRNA 流程的理由）
  - spot 数必须与位置文件里 in_tissue=1 的行数对得上
  - 空间坐标必须覆盖全部 spot（缺坐标的 spot 画不到图上）

输出：
  data/<id>/raw.h5ad          （含 obsm['spatial'] 与 uns['spatial']）
  data/<id>/dataset_info.json
  data/<id>/cache/            下载与解压缓存
"""

from __future__ import annotations

import json
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))

import numpy as np  # noqa: E402
from scipy import sparse  # noqa: E402

from alignment import write_alignment_check  # noqa: E402
from common import (ensure_dirs, load_config, load_registry,  # noqa: E402
                    load_scalefactors,
                    log_info, log_warn, parse_args, read_tissue_positions,
                    record_step, set_seed, write_json)

MIN_SPOTS = 100


def download(url: str, dest: Path, retries: int = 3) -> Path:
    """
    带重试的下载。**必须带 User-Agent。**

    10x 的 CDN（cf.10xgenomics.com）对 `Python-urllib/3.x` 直接返回
    403 Forbidden，而同样的 URL 用浏览器请求是 200。报错信息只有
    `HTTP Error 403: Forbidden`，看不出是 UA 的问题 —— 很容易误判成
    "链接失效了"然后去换数据集，方向完全错了。
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        log_info(f"已缓存，跳过下载: {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
        return dest
    last = None
    for attempt in range(1, retries + 1):
        try:
            log_info(f"下载 ({attempt}/{retries}): {url}")
            tmp = dest.with_suffix(dest.suffix + ".part")
            req = urllib.request.Request(url, headers={
                "User-Agent": ("Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
                               "(KHTML, like Gecko) Chrome/120.0 Safari/537.36"),
                "Accept": "*/*",
            })
            with urllib.request.urlopen(req, timeout=300) as r, open(tmp, "wb") as fh:
                total = int(r.headers.get("Content-Length") or 0)
                got = 0
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
                    got += len(chunk)
            if total and got != total:
                raise IOError(f"下载不完整: {got}/{total} 字节")
            if got == 0:
                raise IOError("下载到 0 字节")
            tmp.replace(dest)
            log_info(f"下载完成: {dest.name} ({got/1e6:.1f} MB)")
            return dest
        except Exception as e:  # noqa: BLE001
            last = e
            log_warn(f"下载失败 ({attempt}/{retries}): {e}")
    raise RuntimeError(f"下载 {url} 失败（重试 {retries} 次）: {last}")


def extract_spatial(tar_path: Path, cache_dir: Path) -> Path:
    """解压 spatial tarball，返回含 tissue_positions 的目录。"""
    out = cache_dir / "spatial"
    marker = out / ".done"
    if not marker.exists():
        out.mkdir(parents=True, exist_ok=True)
        log_info(f"解压 {tar_path.name}")
        with tarfile.open(tar_path, "r:gz") as tf:
            try:
                tf.extractall(out, filter="data")
            except TypeError:
                tf.extractall(out)
        marker.write_text("ok", encoding="utf-8")

    # tarball 里是 spatial/xxx，解压后变成 out/spatial/xxx
    cands = list(out.rglob("tissue_positions*.csv"))
    if not cands:
        raise RuntimeError(f"{tar_path.name} 里找不到 tissue_positions*.csv —— "
                           f"不是 Space Ranger 的 spatial 输出")
    return cands[0].parent


def validate_counts(adata) -> dict:
    """判断 X 是不是原始整数计数（同 scRNA 流程 —— 理由见那里的说明）。"""
    X = adata.X
    sample = X[:min(200, X.shape[0]), :]
    vals = sample.data if sparse.issparse(sample) else np.asarray(sample).ravel()
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {"is_counts": False, "reason": "矩阵里没有有限值"}
    is_int = bool(np.allclose(vals, np.round(vals), atol=1e-8))
    min_v = float(vals.min())
    return {
        "is_counts": bool(is_int and min_v >= 0),
        "min": min_v, "max": float(vals.max()),
        "frac_integer": round(float(np.mean(np.isclose(vals, np.round(vals),
                                                      atol=1e-8))), 6),
        "reason": ("非负整数 -> 判定为原始计数" if (is_int and min_v >= 0)
                   else f"**不是整数计数**（最小 {min_v:.4g}）—— 可能已 normalize/log 过"),
    }


def run_00_fetch(cfg: dict) -> dict:
    import anndata as ad
    import scanpy as sc

    ensure_dirs(cfg)
    set_seed(cfg)
    data_dir = Path(cfg["output"]["data_dir"])
    cache = data_dir / "cache"
    repo_root = Path(__file__).resolve().parent.parent

    src = dict(cfg.get("source") or {})
    if src.get("dataset"):
        reg = load_registry(repo_root)
        key = src["dataset"]
        if key not in reg:
            raise KeyError(f"source.dataset='{key}' 不在 assets/datasets.yml 里；"
                           f"可选: {', '.join(sorted(reg))}")
        entry = dict(reg[key])
        entry.update({k: v for k, v in src.items() if k != "dataset"})
        src = entry
        log_info(f"数据集 '{key}': {src.get('note', '')}")

    matrix_url = src.get("matrix_url")
    spatial_url = src.get("spatial_url")
    if not matrix_url or not spatial_url:
        raise ValueError(
            "配置的 source 段必须同时给 matrix_url 与 spatial_url"
            "（或给 source.dataset 指向 assets/datasets.yml 里的一项）。\n"
            "  **Visium 的数据是两半**：表达矩阵 + 空间信息（坐标/缩放/图像）。\n"
            "  只给矩阵的话数据会退化成普通单细胞数据 —— 丢掉空间坐标\n"
            "  等于丢掉这个技术存在的理由。")

    # ---- 1. 下载 ------------------------------------------------------------
    h5 = download(matrix_url, cache / Path(matrix_url).name)
    tar = download(spatial_url, cache / Path(spatial_url).name)

    # ---- 2. 读矩阵 ----------------------------------------------------------
    adata = sc.read_10x_h5(h5)
    adata.var_names_make_unique()
    log_info(f"表达矩阵: {adata.n_obs} spot x {adata.n_vars} 基因")

    # ---- 3. 读空间信息 ------------------------------------------------------
    spatial_dir = extract_spatial(tar, cache)
    pos_file = next(spatial_dir.glob("tissue_positions*.csv"))
    pos, had_header = read_tissue_positions(pos_file)
    log_info(f"位置文件: {pos_file.name}（{'有' if had_header else '无'}表头），"
             f"{len(pos)} 行")
    pos = pos.set_index("barcode")
    pos["in_tissue"] = pos["in_tissue"].astype(int)

    # ---- 4. 对齐 ------------------------------------------------------------
    common = adata.obs_names.intersection(pos.index)
    if len(common) == 0:
        raise RuntimeError(
            "表达矩阵的 barcode 与位置文件的 barcode 完全没有交集 —— "
            "这两个文件可能不是同一次 Space Ranger 运行的结果")
    n_missing_pos = int(adata.n_obs - len(common))
    if n_missing_pos:
        log_warn(f"{n_missing_pos} 个 spot 在位置文件里找不到坐标，将被剔除")

    adata = adata[common].copy()
    pos = pos.loc[common]
    adata.obs["in_tissue"] = pos["in_tissue"].values
    adata.obs["array_row"] = pos["array_row"].astype(int).values
    adata.obs["array_col"] = pos["array_col"].astype(int).values

    # obsm['spatial'] 的列序是 (x, y) = (pxl_col_in_fullres, pxl_row_in_fullres)
    # **顺序不能反。** 反了的话组织图会转置/镜像，而散点图看起来
    # "形状有点像" —— 只有叠到 H&E 图上才发现不对。
    adata.obsm["spatial"] = np.column_stack([
        pos["pxl_col_in_fullres"].astype(float).values,
        pos["pxl_row_in_fullres"].astype(float).values,
    ])

    n_in_tissue = int((adata.obs["in_tissue"] == 1).sum())
    log_info(f"in_tissue=1 的 spot: {n_in_tissue}/{adata.n_obs}")

    # ---- 5. 图像与缩放因子 --------------------------------------------------
    scale_file = next(spatial_dir.glob("scalefactors_json.json"))
    scalefactors = load_scalefactors(scale_file)
    log_info(f"缩放因子: hires={scalefactors['tissue_hires_scalef']:.6f}, "
             f"lowres={scalefactors['tissue_lowres_scalef']:.6f}, "
             f"spot_diameter_fullres={scalefactors['spot_diameter_fullres']:.2f}")

    images = {}
    for key, pat in (("hires", "tissue_hires_image.png"),
                     ("lowres", "tissue_lowres_image.png")):
        f = next(spatial_dir.glob(pat), None)
        if f is None:
            log_warn(f"缺少 {pat} —— 无法把 spot 叠到组织图上")
            continue
        from PIL import Image
        images[key] = np.asarray(Image.open(f).convert("RGB"))
        log_info(f"{key} 图像: {images[key].shape[1]}x{images[key].shape[0]}")

    if not images:
        raise RuntimeError("spatial tarball 里没有任何组织图像（hires/lowres）—— "
                           "没有图像就无法做空间可视化，Visium 的意义也没了")

    # scanpy/squidpy 期望的结构。
    # **注意读的是合并后的 `src`，不是 `cfg["source"]`。** library_id 通常
    # 只在 datasets.yml 里，config 里没有 —— 读 cfg 会静默拿到默认值
    # "library"，而 uns['spatial'] 的键跟着错，下游 squidpy 取图像时找不到。
    libname = str(src.get("library_id") or cfg["dataset_id"])
    adata.uns["spatial"] = {libname: {"images": images,
                                      "scalefactors": scalefactors}}
    adata.obs["library_id"] = libname

    # ---- 6. 硬校验 ----------------------------------------------------------
    counts = validate_counts(adata)
    log_info(f"计数性质: {counts['reason']}")
    if not counts["is_counts"]:
        raise RuntimeError(
            "输入矩阵不是原始整数计数，拒绝继续。\n"
            f"  依据: {counts['reason']}\n"
            "  为什么必须停: 空间域与 SVGs 依赖正确的计数分布；\n"
            "  解卷积的 NNLS 前提也会被破坏 —— 而图看起来完全正常。")

    if n_in_tissue < MIN_SPOTS:
        raise RuntimeError(f"in_tissue=1 的 spot 只有 {n_in_tissue} < {MIN_SPOTS}")
    if adata.n_vars < 50:
        raise RuntimeError(f"基因数 {adata.n_vars} < 50")

    # ---- 7. 落盘 ------------------------------------------------------------
    out = data_dir / "raw.h5ad"
    adata.write_h5ad(out)
    log_info(f"已写出 {out}（{adata.n_obs} spot x {adata.n_vars} 基因，"
             f"含 obsm['spatial'] 与 uns['spatial']）")

    # ---- 8. 坐标对齐验证 ----------------------------------------------------
    # **必须在这里做，不能只在 workflow 里单独跑。** 验收检查会读
    # `spatial_alignment_check.json`；如果这个文件在验收时还不存在，
    # 那一项会被**静默跳过** —— 实测 CI 里 25 项检查、本地 26 项，
    # 差的就是这一项。一个会因执行顺序而消失的检查等于没有。
    alignment = write_alignment_check(adata, data_dir, log_info, log_warn)
    if not alignment["current_is_best"]:
        raise RuntimeError(
            f"坐标对齐验证失败：当前坐标不是最佳假设，"
            f"最佳是「{alignment['best_hypothesis']}」。\n"
            f"  坐标可能转置或镜像了。散点图看起来仍然像组织形状，\n"
            f"  所以这个问题不检查就发现不了。\n"
            f"  详见 {data_dir / 'spatial_alignment_check.json'}")

    info = {
        "dataset_id": cfg["dataset_id"],
        "technology": "Visium (10x Genomics)",
        "matrix_url": matrix_url,
        "spatial_url": spatial_url,
        "library_id": libname,
        "source_note": src.get("note", ""),
        "organism": src.get("organism"),
        "tissue": src.get("tissue"),
        "condition": src.get("condition"),
        "n_spots_total": int(adata.n_obs),
        "n_spots_in_tissue": n_in_tissue,
        "n_genes": int(adata.n_vars),
        "positions_file": pos_file.name,
        "positions_had_header": bool(had_header),
        "scalefactors": scalefactors,
        "images": {k: [int(v.shape[1]), int(v.shape[0])] for k, v in images.items()},
        "counts_check": counts,
        "coordinate_note": ("obsm['spatial'] 列序为 (x, y) = "
                            "(pxl_col_in_fullres, pxl_row_in_fullres)；"
                            "顺序反了会让组织图转置/镜像"),
        "min_spots_gate": MIN_SPOTS,
    }
    write_json(data_dir / "dataset_info.json", info)
    return info


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        run_00_fetch(cfg)
        record_step(cfg, "fetch", "ok", time.time() - t0)
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "fetch", "failed", time.time() - t0, message=str(e))
        raise
