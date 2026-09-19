#!/usr/bin/env python3
"""
main_analysis.py — 主流程：依次跑完所有步骤，产出验收报告

**每一步都是独立的可执行脚本**，这个主流程只是把它们串起来。
这样单步失败时可以单独重跑，而不用从头再来。

**验收检查（acceptance checks）在最后跑。** 它们检查的是"产物是否
真的存在且有内容"，不是"结果是否生物学正确" —— 后者需要人看。
一条检查失败会写进 `acceptance_report.json` 并以非零码退出，
但**前面的产物都已落盘**，便于诊断。
"""

from __future__ import annotations

import argparse
import importlib
import json
import sys
import time
import traceback
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts" / "lib"))
sys.path.insert(0, str(REPO / "scripts"))

from common import ensure_dirs, load_config, log_info, log_warn, write_json  # noqa: E402

STEPS = [
    ("fetch", "00_fetch", "run_00_fetch"),
    ("qc", "01_qc", "run_01_qc"),
    ("normalize", "02_normalize", "run_02_normalize"),
    ("spatial_domains", "03_spatial_domains", "run_03_spatial_domains"),
    ("svg", "04_svg", "run_04_svg"),
    ("deconvolution", "05_deconvolution", "run_05_deconvolution"),
    ("niche", "06_niche", "run_06_niche"),
    ("spatial_communication", "07_spatial_communication",
     "run_07_spatial_communication"),
]


def run_steps(cfg: dict, only: list = None) -> dict:
    ensure_dirs(cfg)
    results = {}
    for name, mod_name, fn_name in STEPS:
        if only and name not in only:
            continue
        log_info("")
        log_info("=" * 70)
        log_info(f"步骤 {name}  ({mod_name}.py)")
        log_info("=" * 70)
        t0 = time.time()
        try:
            mod = importlib.import_module(mod_name)
            fn = getattr(mod, fn_name)
            res = fn(cfg)
            results[name] = {"status": "ok", "seconds": round(time.time() - t0, 1),
                             "result_status": (res or {}).get("status", "ok")}
            log_info(f"步骤 {name} 完成（{time.time() - t0:.1f} 秒）")
        except Exception as e:  # noqa: BLE001
            results[name] = {"status": "failed", "seconds": round(time.time() - t0, 1),
                             "error": f"{type(e).__name__}: {e}",
                             "traceback": traceback.format_exc()[-3000:]}
            log_warn(f"步骤 {name} **失败**: {type(e).__name__}: {e}")
            log_warn("产物已落盘到失败前的位置，便于诊断")
            break
    return results


def run_acceptance(cfg: dict, step_results: dict) -> dict:
    """
    验收检查。

    分三类：
      - **required**：缺失就是流水线坏了（如没有 h5ad、没有图）
      - **content**：文件存在但可能是空的（如 CSV 只有表头）
      - **honesty**：**"如实记录"类检查** —— 检查每个步骤是否
        写明了它没做什么。这类检查是防止"静默跳过"的关键。

    `honesty` 类是本流水线的特色：如果某步没做（如没有参考就
    不解卷积），它必须写出 `status: not_done` 和原因，而不是
    悄悄不产出文件。
    """
    data_dir = Path(cfg["output"]["data_dir"])
    res_dir = Path(cfg["output"]["results_dir"])
    fig_dir = Path(cfg["output"]["figures_dir"])
    checks = []

    def chk(cid, kind, ok, detail, severity="required"):
        checks.append({"id": cid, "kind": kind, "passed": bool(ok),
                       "detail": detail, "severity": severity})

    # ---- required: 数据文件 ----
    for f in ("raw.h5ad", "qc_filtered.h5ad", "normalized.h5ad", "domains.h5ad"):
        p = data_dir / f
        chk(f"data:{f}", "required", p.exists() and p.stat().st_size > 1024,
            f"{p.name} {'存在' if p.exists() else '缺失'}"
            + (f" ({p.stat().st_size/1e6:.2f} MB)" if p.exists() else ""))

    # ---- required: 结果 JSON ----
    for f in ("qc_status.json", "normalize_status.json", "domain_status.json",
              "svg_status.json", "deconvolution_status.json", "niche_status.json",
              "communication_status.json"):
        p = res_dir / f
        chk(f"result:{f}", "required", p.exists(), f"{f} {'存在' if p.exists() else '缺失'}")

    # ---- required: 图 ----
    figs = sorted(fig_dir.glob("*.png")) if fig_dir.exists() else []
    chk("figures:count", "required", len(figs) >= 8,
        f"{len(figs)} 张图（要求 >=8）")

    # ---- content: 图不是空白 ----
    # 用像素标准差判断：全白/全黑的图标准差接近 0。
    # **这个检查是必要的**：matplotlib 在数据为空时会静默产出一张空白图，
    # 文件存在、大小正常，看不出问题。
    blank = []
    try:
        from PIL import Image
        import numpy as np
        for p in figs:
            try:
                a = np.asarray(Image.open(p).convert("L"), dtype=float)
                if a.std() < 3.0:
                    blank.append({"figure": p.name, "std": round(float(a.std()), 3)})
            except Exception as e:  # noqa: BLE001
                blank.append({"figure": p.name, "error": str(e)})
    except ImportError:
        blank = [{"note": "PIL 不可用，跳过空白图检查"}]
    chk("figures:non_blank", "content", not blank,
        "所有图都有内容" if not blank else f"疑似空白图: {blank}", severity="content")

    # ---- content: CSV 有数据行 ----
    for f in ("svg_results.csv", "domain_markers.csv",
              "communication_lr_scores.csv", "niche_enrichment_domains.csv"):
        p = res_dir / f
        if not p.exists():
            chk(f"content:{f}", "content", False, f"{f} 缺失", severity="content")
            continue
        n = sum(1 for _ in open(p, encoding="utf-8")) - 1
        chk(f"content:{f}", "content", n > 0, f"{f}: {n} 数据行", severity="content")

    # ---- honesty: 每步必须说明自己没做什么 ----
    # 这是防止"静默跳过"的检查。
    honest_map = {
        "qc_status.json": "threshold_note",
        "domain_status.json": "limitations",
        "svg_status.json": "limitations",
        "deconvolution_status.json": "limitations",
        "niche_status.json": "limitations",
        "communication_status.json": "limitations",
    }
    for f, key in honest_map.items():
        p = res_dir / f
        if not p.exists():
            continue
        try:
            d = json.loads(p.read_text(encoding="utf-8"))
        except Exception as e:  # noqa: BLE001
            chk(f"honesty:{f}", "honesty", False, f"无法解析: {e}", severity="honesty")
            continue
        if d.get("status") in ("disabled", "not_done", "not_configured",
                               "no_usable_pairs", "no_usable_signature",
                               "missing_domains"):
            # 没做的话，必须有 reason
            ok = bool(d.get("reason"))
            chk(f"honesty:{f}", "honesty", ok,
                f"状态 {d.get('status')}，原因: {str(d.get('reason'))[:120]}",
                severity="honesty")
        else:
            ok = key in d and d[key]
            chk(f"honesty:{f}", "honesty", ok,
                f"写明了 '{key}'" if ok else f"**缺少 '{key}'** —— 没说明局限",
                severity="honesty")

    # ---- honesty: 可用 LR 对数必须报 ----
    p = res_dir / "communication_status.json"
    if p.exists():
        d = json.loads(p.read_text(encoding="utf-8"))
        if d.get("status") == "ok":
            ok = "n_pairs_usable" in d and "n_pairs_total" in d
            frac = (d.get("usable_fraction") or 0)
            chk("honesty:lr_usability", "honesty", ok,
                f"可用 {d.get('n_pairs_usable')}/{d.get('n_pairs_total')} "
                f"({frac:.2f})" if ok else "**没有报可用 LR 对数**",
                severity="honesty")

    # ---- honesty: 解卷积必须标明是不是真解卷积 ----
    p = res_dir / "deconvolution_status.json"
    if p.exists():
        d = json.loads(p.read_text(encoding="utf-8"))
        if d.get("status") == "ok":
            ok = "is_deconvolution" in d and "method" in d
            chk("honesty:deconv_method", "honesty", ok,
                f"is_deconvolution={d.get('is_deconvolution')}；{str(d.get('method'))[:90]}"
                if ok else "**没有标明是真解卷积还是打分法**",
                severity="honesty")

    # ---- 空间对齐检查（如果有）----
    p = data_dir / "spatial_alignment_check.json"
    if p.exists():
        d = json.loads(p.read_text(encoding="utf-8"))
        chk("spatial:alignment", "content", d.get("current_is_best", False),
            f"最佳假设: {d.get('best_hypothesis')}；"
            f"当前是否最佳: {d.get('current_is_best')}", severity="content")

    n_req_fail = sum(1 for c in checks
                     if not c["passed"] and c["severity"] == "required")
    n_con_fail = sum(1 for c in checks
                     if not c["passed"] and c["severity"] == "content")
    n_hon_fail = sum(1 for c in checks
                     if not c["passed"] and c["severity"] == "honesty")

    report = {
        "dataset_id": cfg["dataset_id"],
        "steps": step_results,
        "n_checks": len(checks),
        "n_required_failed": n_req_fail,
        "n_content_failed": n_con_fail,
        "n_honesty_failed": n_hon_fail,
        "passed": n_req_fail == 0,
        "checks": checks,
    }
    write_json(res_dir / "acceptance_report.json", report)

    log_info("")
    log_info("=" * 70)
    log_info(f"验收: {len(checks)} 项检查，"
             f"required 失败 {n_req_fail}，content 失败 {n_con_fail}，"
             f"honesty 失败 {n_hon_fail}")
    log_info("=" * 70)
    for c in checks:
        if not c["passed"]:
            log_warn(f"  [未通过/{c['severity']}] {c['id']}: {c['detail']}")
    if n_req_fail == 0:
        log_info("所有 required 检查通过")
    return report


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--only", nargs="*", default=None,
                    help="只跑指定步骤（名字见 STEPS）")
    ap.add_argument("--skip-acceptance", action="store_true")
    args = ap.parse_args()

    cfg = load_config(args.config)
    t0 = time.time()
    log_info(f"数据集: {cfg['dataset_id']}")
    log_info(f"输出: {cfg['output']['data_dir']}")

    step_results = run_steps(cfg, args.only)
    report = None
    if not args.skip_acceptance:
        report = run_acceptance(cfg, step_results)

    log_info("")
    log_info(f"总耗时 {time.time() - t0:.1f} 秒")
    if any(v["status"] == "failed" for v in step_results.values()):
        return 1
    if report and report["n_required_failed"] > 0:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
