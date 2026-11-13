#!/usr/bin/env python3
"""
gse192741_prepare.py — GSE192741 RAW tar（15 个 Visium library 的散装文件）
→ 按 library 组装成 Space Ranger 标准输出目录

GEO 补充包里的文件是平铺的 `GSM*_filtered_feature_bc_matrix_<LIB>.h5`、
`GSM*_tissue_positions_list_<LIB>.csv.gz` 等，而 spatial-pipeline-skill 的
00_fetch.py 吃的是「matrix_url（h5）+ spatial_url（tar.gz，含 tissue_positions/
scalefactors/hires/lowres 四件）」两半结构。本适配器把 tar 里的散文件按
library 重组：

  <outdir>/<LIB>/
    filtered_feature_bc_matrix.h5
    spatial/tissue_positions_list.csv
    spatial/scalefactors_json.json
    spatial/tissue_hires_image.png
    spatial/tissue_lowres_image.png

并落 <outdir>/libraries.json（library 清单 + GSM + organism/condition 元数据，
从随仓内嵌的 GSE192741 SOFT 映射来）。

用法（容器内）：
  python3 gse192741_prepare.py --tar <RAW.tar> --outdir <dir> [--library JBO001]
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import sys
import tarfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import log_info, log_warn  # noqa: E402

# GSM → (library, organism, condition)（GSE192741_family.soft 实测映射）
GSM_META = {
    'GSM5764414': ('JBO001', 'mouse', 'StSt_healthy'),
    'GSM5764415': ('JBO002', 'mouse', 'StSt_healthy'),
    'GSM5764416': ('JBO003', 'mouse', 'StSt_healthy'),
    'GSM5764417': ('JBO004', 'mouse', 'StSt_healthy'),
    'GSM5764418': ('CAP002', 'mouse', 'StSt_capsule'),
    'GSM5764419': ('JBO006', 'mouse', 'SD'),
    'GSM5764420': ('JBO009', 'mouse', 'WD'),
    'GSM5764421': ('JBO010', 'mouse', 'WD'),
    'GSM5764422': ('JBO012', 'mouse', 'WD'),
    'GSM5764423': ('BH001', 'mouse', 'StSt_healthy_withAbs'),
    'GSM5764424': ('JBO014', 'human', 'steatotic'),
    'GSM5764425': ('JBO015', 'human', 'steatotic'),
    'GSM5764426': ('JBO018', 'human', 'healthy'),
    'GSM5764427': ('JBO019', 'human', 'steatotic'),
    'GSM5764428': ('JBO022', 'human', 'healthy'),
}

# tar 内文件名模式 → 目标相对路径
PATTERNS = [
    ('_filtered_feature_bc_matrix_', '.h5', 'filtered_feature_bc_matrix.h5'),
    ('_tissue_positions_list_', '.csv.gz', 'spatial/tissue_positions_list.csv'),
    ('_scalefactors_json_', '.json.gz', 'spatial/scalefactors_json.json'),
    ('_tissue_hires_image_', '.png.gz', 'spatial/tissue_hires_image.png'),
    ('_tissue_lowres_image_', '.png.gz', 'spatial/tissue_lowres_image.png'),
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tar', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--library', default='', help='只处理指定 library（默认全部）')
    a = ap.parse_args()

    tar_path = Path(a.tar)
    out_root = Path(a.outdir)
    only = a.library or None

    made = {}
    with tarfile.open(tar_path, 'r') as tf:
        for m in tf.getmembers():
            if not m.isfile():
                continue
            name = Path(m.name).name
            gsm = name.split('_')[0]
            meta = GSM_META.get(gsm)
            if not meta:
                continue
            lib = meta[0]
            if only and lib != only:
                continue
            for key, suffix, dest in PATTERNS:
                if key in name and name.endswith(suffix):
                    lib_dir = out_root / lib
                    dest_path = lib_dir / dest
                    dest_path.parent.mkdir(parents=True, exist_ok=True)
                    fh = tf.extractfile(m)
                    data = fh.read()
                    if suffix.endswith('.gz'):
                        data = gzip.decompress(data)
                    dest_path.write_bytes(data)
                    made.setdefault(lib, []).append(dest)
                    break

    libs = {}
    for lib, files in sorted(made.items()):
        gsm = next(g for g, v in GSM_META.items() if v[0] == lib)
        need = {d for _, _, d in PATTERNS}
        missing = need - set(files)
        libs[lib] = {
            'gsm': gsm,
            'organism': GSM_META[gsm][1],
            'condition': GSM_META[gsm][2],
            'files': sorted(files),
            'complete': not missing,
            'missing': sorted(missing),
        }
        line = (f'{lib} (GSM{gsm[3:]}): {len(files)}/5 files'
                + (f' MISSING {missing}' if missing else ' OK'))
        (log_warn if missing else log_info)(line)

    (out_root / 'libraries.json').write_text(
        json.dumps(libs, ensure_ascii=False, indent=2), encoding='utf-8')
    log_info(f'共 {len(libs)} 个 library -> {out_root}')
    return 0 if all(v['complete'] for v in libs.values()) else 2


if __name__ == '__main__':
    sys.exit(main())
