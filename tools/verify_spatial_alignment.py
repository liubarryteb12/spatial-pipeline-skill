#!/usr/bin/env python3
"""
tools/verify_spatial_alignment.py — 验证 spot 坐标与 H&E 图像真的对齐

命令行封装。核心逻辑在 `scripts/lib/alignment.py`，因为
**流水线的 `00_fetch.py` 也要调用它** —— 对齐检查必须在
`main_analysis.py` 的验收检查之前产出 `spatial_alignment_check.json`，
否则那一项验收会因为文件不存在而**静默跳过**
（实测：CI 里 25 项检查、本地 26 项，差的就是这一项）。

判据的说明见 `scripts/lib/alignment.py` 的模块 docstring。

用法: python tools/verify_spatial_alignment.py --config assets/config.lymph_node.yml
退出码: 0 = 当前坐标是最佳假设；1 = 不是
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts" / "lib"))

from alignment import verify_alignment  # noqa: E402
from common import load_config, log_info, log_warn  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    import scanpy as sc

    cfg = load_config(args.config)
    p = Path(cfg["output"]["data_dir"]) / "raw.h5ad"
    if not p.exists():
        print(f"跳过：{p} 不存在")
        return 0

    res = verify_alignment(sc.read_h5ad(p), log_info=log_info, log_warn=log_warn)
    return 0 if res["current_is_best"] else 1


if __name__ == "__main__":
    sys.exit(main())
