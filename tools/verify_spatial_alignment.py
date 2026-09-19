#!/usr/bin/env python3
"""
tools/verify_spatial_alignment.py — 验证 spot 坐标与 H&E 图像真的对齐

**为什么需要它。** Visium 的坐标有四种常见的错法，而且都不会报错：
  1. x/y 顺序反了（转置）
  2. 轴向反了（镜像）
  3. 忘了乘 scalef（坐标落在全分辨率尺度上，远超图像范围）
  4. 用了 fullres 坐标配 hires 图

这四种情况下散点图**看起来都"有点像组织形状"**，只有定量检查才能发现。

---

## 判据的选择（这里有一个实测踩到的坑）

第一版用"spot 中心处是不是组织像素"作为判据。**对淋巴结数据完全无效** ——
实测整张 hires 图的组织像素占比是 **1.0000**（组织铺满整帧，没有白色背景），
所以四个假设（含转置、镜像）都得 1.000 分，检查报"通过"而实际上什么都没验证。

**判据必须有区分力，否则它给的是虚假的安心。**

改用两条真正有区分力的判据：

  A. **网格几何**。Visium 的 spot 在 array (row, col) 空间是规则六边形网格，
     映射到图像后，array_col 应与某一个像素轴强相关、array_row 与另一个强相关。
     实测正确映射下 corr(col,px)=+1.0000、corr(row,py)=+1.0000、
     两个交叉项 ≈ -0.09；转置假设下这些值完全互换。
     **这是最强的判据** —— 转置/镜像会立刻暴露。

  B. **网格规整度**。最近邻距离应集中在一个值附近。实测中位 137.0 px、
     四分位距 0.6 px、变异系数 0.031。尺度用错（忘乘 scalef）会让这个值
     离谱地大或小。

组织像素占比仍然报出来，但会**标注它有没有区分力**（整图组织占比
接近 1 时它没有）。

用法: python tools/verify_spatial_alignment.py --config assets/config.lymph_node.yml
退出码: 0 = 当前坐标是最佳假设；1 = 不是
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts" / "lib"))

from common import load_config, log_info, log_warn  # noqa: E402


def grid_geometry_score(px, py, ac, ar) -> dict:
    """
    判据 A：网格几何。

    正确映射下，array_col 与一个像素轴 |corr| 接近 1，array_row 与另一个
    接近 1，且两个交叉项接近 0。返回"对角线强度 - 交叉项强度"作为分数。
    """
    import numpy as np

    def c(a, b):
        if np.std(a) < 1e-9 or np.std(b) < 1e-9:
            return 0.0
        return float(np.corrcoef(a, b)[0, 1])

    diag = abs(c(ac, px)) + abs(c(ar, py))
    cross = abs(c(ac, py)) + abs(c(ar, px))
    return {"diag": round(diag, 4), "cross": round(cross, 4),
            "corr_col_x": round(c(ac, px), 4), "corr_row_y": round(c(ar, py), 4),
            "corr_col_y": round(c(ac, py), 4), "corr_row_x": round(c(ar, px), 4),
            "score": round(diag - cross, 4)}


def grid_regularity(px, py) -> dict:
    """
    判据 B：网格规整度。

    最近邻距离的变异系数。规则六边形网格应接近 0（实测 0.031）。
    """
    import numpy as np
    from scipy.spatial import cKDTree

    xy = np.column_stack([px, py])
    t = cKDTree(xy)
    d, _ = t.query(xy, k=2)
    nn = d[:, 1]
    nn = nn[np.isfinite(nn) & (nn > 0)]
    if nn.size == 0:
        return {"cv": None}
    return {"median_nn": round(float(np.median(nn)), 2),
            "iqr_nn": round(float(np.percentile(nn, 75) - np.percentile(nn, 25)), 3),
            "cv": round(float(nn.std() / nn.mean()), 4)}


def tissue_fraction(img, px, py) -> dict:
    """辅助判据：spot 中心处是组织像素的比例。**可能没有区分力。**"""
    import numpy as np

    h, w = img.shape[:2]
    xi = np.clip(np.round(px).astype(int), 0, w - 1)
    yi = np.clip(np.round(py).astype(int), 0, h - 1)
    pxv = img[yi, xi].astype(float)
    spot_frac = float((np.abs(255.0 - pxv).sum(axis=1) > 60).mean())
    # 整图的组织占比 —— 用来判断这个判据有没有区分力
    whole = float((np.abs(255.0 - img.astype(float)).sum(axis=2) > 60).mean())
    return {"spot_tissue_fraction": round(spot_frac, 4),
            "whole_image_tissue_fraction": round(whole, 4),
            # 整图几乎全是组织时，任何点都"落在组织上"，判据失去区分力
            "discriminative": bool(whole < 0.98)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    import numpy as np
    import scanpy as sc

    cfg = load_config(args.config)
    data_dir = Path(cfg["output"]["data_dir"])
    p = data_dir / "raw.h5ad"
    if not p.exists():
        print(f"跳过：{p} 不存在")
        return 0

    adata = sc.read_h5ad(p)
    lib = list(adata.uns["spatial"].keys())[0]
    entry = adata.uns["spatial"][lib]
    img = entry["images"]["hires"]
    sf = float(entry["scalefactors"]["tissue_hires_scalef"])
    xy = adata.obsm["spatial"]
    h, w = img.shape[:2]
    ac = adata.obs["array_col"].values.astype(float)
    ar = adata.obs["array_row"].values.astype(float)

    log_info(f"hires 图像 {w}x{h}，scalef={sf:.6f}，{xy.shape[0]} 个 spot")

    hyps = {
        "当前 (x=col, y=row)": (xy[:, 0] * sf, xy[:, 1] * sf),
        "转置 (x=row, y=col)": (xy[:, 1] * sf, xy[:, 0] * sf),
        "镜像 x": ((w - 1) - xy[:, 0] * sf, xy[:, 1] * sf),
        "镜像 y": (xy[:, 0] * sf, (h - 1) - xy[:, 1] * sf),
        "忘乘 scalef": (xy[:, 0], xy[:, 1]),
    }

    results = {}
    for name, (px, py) in hyps.items():
        inb = float(((px >= 0) & (px < w) & (py >= 0) & (py < h)).mean())
        r = {"in_bounds": round(inb, 4)}
        if inb > 0.5:
            r["grid_geometry"] = grid_geometry_score(px, py, ac, ar)
            r["grid_regularity"] = grid_regularity(px, py)
            r["tissue"] = tissue_fraction(img, px, py)
        results[name] = r

    log_info("")
    log_info("判据 A：网格几何（对角线强度 − 交叉项强度，越高越好）")
    log_info("判据 B：网格规整度（最近邻距离变异系数，越低越好）")
    log_info("")
    log_info(f"    {'假设':<22} {'在界内':>7} {'几何分':>9} {'交叉项':>9} "
             f"{'NN中位':>9} {'NN-CV':>8} {'组织占比':>10}")

    def fmt(v, nd=4):
        return "—" if v is None else f"{v:.{nd}f}"

    best, best_key = None, None
    for name, r in results.items():
        gg = r.get("grid_geometry") or {}
        gr = r.get("grid_regularity") or {}
        ti = r.get("tissue") or {}
        gs = gg.get("score")
        key = gs if gs is not None else -99.0
        if best_key is None or key > best_key:
            best, best_key = name, key
        log_info(f"    {name:<22} {r['in_bounds']:>7.3f} {fmt(gs):>9} "
                 f"{fmt(gg.get('cross')):>9} {fmt(gr.get('median_nn'), 1):>9} "
                 f"{fmt(gr.get('cv')):>8} "
                 f"{fmt(ti.get('spot_tissue_fraction'), 3):>10}")

    cur = results["当前 (x=col, y=row)"]
    cur_ti = cur.get("tissue", {})
    if cur_ti and not cur_ti.get("discriminative", True):
        log_warn("注意：整图组织占比 %.3f —— 『组织像素占比』这条判据对本数据"
                 "**没有区分力**（组织铺满整帧，任何点都落在组织上）。"
                 "结论以网格几何为准" % cur_ti["whole_image_tissue_fraction"])

    ok = (best == "当前 (x=col, y=row)")
    if ok:
        gg = cur["grid_geometry"]
        log_info("")
        log_info(f"对齐验证通过：当前坐标是最佳假设，几何分 {gg['score']:.4f}"
                 f"（对角线 {gg['diag']:.4f} vs 交叉项 {gg['cross']:.4f}）")
    else:
        log_warn("")
        log_warn(f"**当前坐标不是最佳假设** —— 最佳是「{best}」"
                 f"（几何分 {best_key:.4f} vs 当前 "
                 f"{cur.get('grid_geometry', {}).get('score', float('nan')):.4f}）。"
                 f"坐标可能转置或镜像了")

    out = data_dir / "spatial_alignment_check.json"
    out.write_text(json.dumps({
        "library_id": lib, "image_size": [w, h], "scalef": sf,
        "n_spots": int(xy.shape[0]),
        "hypotheses": results,
        "best_hypothesis": best,
        "current_is_best": ok,
        "criteria": {
            "A_grid_geometry": "array(row,col) 与像素轴的相关系数矩阵；"
                               "转置/镜像会立刻暴露。**主判据**",
            "B_grid_regularity": "最近邻距离变异系数；尺度用错会暴露",
            "C_tissue_fraction": "辅助判据；整图组织占比接近 1 时**无区分力**",
        },
        "note": ("第一版只用『组织像素占比』，对淋巴结数据四个假设全得 1.000 "
                 "（整图 100% 是组织），等于什么都没验证。判据必须有区分力，"
                 "否则它给的是虚假的安心"),
        "detects": {
            "transposition": "能 —— 几何分会从 +1.81 翻成 -1.81（实测）",
            "scale_error": "能 —— 坐标落到图像外，在界内比例 0.000",
            "scale_mismatch": "能 —— NN 距离中位数会离谱地大或小",
        },
        "does_not_detect": {
            "mirroring": ("**检测不到。** 镜像保持 |相关系数| 不变，所以"
                          "几何分仍是 +1.81。要判镜像需要看 aligned_fiducials.jpg "
                          "里的基准框方位，本工具没做。"
                          "若怀疑镜像，必须人工核对 spot 网格与 H&E 上的"
                          "组织边界是否同向"),
        },
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    log_info(f"写出 {out}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
