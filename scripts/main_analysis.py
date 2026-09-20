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

from common import (capture_versions, ensure_dirs, init_manifest,  # noqa: E402
                    load_config, log_info, log_warn, manifest_path,
                    manifest_summary, record_decision, record_human_review,
                    record_input, record_params, write_json)

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
    ("spatial_trajectory", "08_spatial_trajectory",
     "run_08_spatial_trajectory"),
]

# 每步会写的状态文件。**跑之前先删掉** —— 否则步骤崩溃时旧文件还在，
# 下游"产物存在"检查读的是**上一轮的**结果，给出虚假的通过。
# 姊妹项目 Part 2 实测踩过：一步因漏 import 崩了，而它的状态文件是上一轮的，
# 验收照样全绿。
STEP_STATUS_FILES = {
    "qc": "qc_status.json",
    "normalize": "normalize_status.json",
    "spatial_domains": "domain_status.json",
    "svg": "svg_status.json",
    "deconvolution": "deconvolution_status.json",
    "niche": "niche_status.json",
    "spatial_communication": "communication_status.json",
    "spatial_trajectory": "spatial_trajectory_status.json",
}

# 文档 §3「本部分人工复核节点」。**默认 pending，不是 confirmed** ——
# 自动化流水线不能替人签字，把未确认的节点记成已确认，等于把复核节点
# 变成摆设。验收里作为**可见但不阻断**的项列出。
HUMAN_REVIEW_NODES = [
    ("cell_segmentation",   "细胞分割参数调整",       False),
    ("domain_number",       "空间域数量确定",         True),
    ("deconv_reference",    "去卷积参考数据选择",     True),
    ("spatial_traj_direction", "空间拟时序轨迹方向验证", True),
]

# 需要登记哈希的输入（相对 data_dir）。(文件名, 中文说明, 是否必需)
INPUT_FILES = [
    ("raw.h5ad",          "原始计数矩阵", True),
    ("dataset_info.json", "数据集元信息",  True),
    ("qc_filtered.h5ad",  "QC 后矩阵",    True),
    ("normalized.h5ad",   "标准化后矩阵",  True),
    ("domains.h5ad",      "最终态矩阵",    True),
]


def run_steps(cfg: dict, only: list = None) -> dict:
    ensure_dirs(cfg)
    res_dir = Path(cfg["output"]["results_dir"])

    # ---- 模块零：建立本轮运行清单（§0.3 / §0.4）-----------------------------
    # **必须在任何步骤之前建，且先清掉上一轮** —— 清单描述的是本轮。
    # 清单和 state.json / acceptance_report.json 分开：那两个记"跑没跑成"，
    # 清单记"在什么条件下跑出来的"（证据，写入后不该再变）。
    init_manifest(cfg)
    capture_versions(cfg)
    record_params(cfg, {
        "seed": cfg.get("analysis", {}).get("seed"),
        "qc": cfg.get("qc", {}),
        "normalize": cfg.get("normalize", {}),
        "domains": cfg.get("domains", {}),
        "deconvolution": cfg.get("deconvolution", {}),
    })
    for node, label, req in HUMAN_REVIEW_NODES:
        record_human_review(cfg, node, required=req, status="pending",
                            note=f"{label} —— 需人工确认，本轮自动化未确认")
    log_info(f"运行清单：{manifest_path(cfg)}")

    results = {}
    for name, mod_name, fn_name in STEPS:
        if only and name not in only:
            continue
        log_info("")
        log_info("=" * 70)
        log_info(f"步骤 {name}  ({mod_name}.py)")
        log_info("=" * 70)
        # 先删本步的状态文件：崩溃时不留旧文件冒充本轮结果
        if name in STEP_STATUS_FILES:
            stale = res_dir / STEP_STATUS_FILES[name]
            if stale.exists():
                stale.unlink()
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

    # ---- 模块零：登记输入哈希（§0.4）----------------------------------------
    # 放在验收里 —— 所有步骤已经跑完，这轮到底产出了哪些文件此刻才确定。
    for fname, desc, req in INPUT_FILES:
        record_input(cfg, data_dir / fname, label=f"{desc} ({fname})",
                     required=req)
    msum = manifest_summary(cfg)

    checks = []

    def chk(cid, kind, ok, detail, severity="required"):
        checks.append({"id": cid, "kind": kind, "passed": bool(ok),
                       "detail": detail, "severity": severity})

    # ---- required: 每个步骤本身必须成功 ------------------------------------
    # 只检查"产物文件在不在"是不够的：步骤崩溃时旧文件还在，
    # 产物检查会通过而实际上本轮什么都没产出。
    for sid, res in step_results.items():
        ok = res.get("status") == "ok"
        chk(f"step:{sid}", "required", ok,
            "成功" if ok else f"**失败**: {str(res.get('error'))[:140]}")

    # ---- required: 数据文件 ----
    for f in ("raw.h5ad", "qc_filtered.h5ad", "normalized.h5ad", "domains.h5ad"):
        p = data_dir / f
        chk(f"data:{f}", "required", p.exists() and p.stat().st_size > 1024,
            f"{p.name} {'存在' if p.exists() else '缺失'}"
            + (f" ({p.stat().st_size/1e6:.2f} MB)" if p.exists() else ""))

    # ---- required: 结果 JSON ----
    for f in ("qc_status.json", "normalize_status.json", "domain_status.json",
              "svg_status.json", "deconvolution_status.json", "niche_status.json",
              "communication_status.json", "spatial_trajectory_status.json"):
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
              "communication_lr_scores.csv", "niche_enrichment_domains.csv",
              "spatial_pseudotime.csv"):
        p = res_dir / f
        if not p.exists():
            chk(f"content:{f}", "content", False, f"{f} 缺失", severity="content")
            continue
        n = sum(1 for _ in open(p, encoding="utf-8")) - 1
        chk(f"content:{f}", "content", n > 0, f"{f}: {n} 数据行", severity="content")

    # ---- honesty: §3.4 点名的 SpatialDE 是否真跑了 ----
    # **不判失败，但必须可见。** SpatialDE 在 scipy>=1.12 上需要垫片才能
    # 导入（见 04_svg.py），垫片失效时它会是 package_missing/failed ——
    # 那种情况必须红字显示，否则"Moran's I 跑了"会被当成"§3.4 做了"。
    sp = res_dir / "svg_status.json"
    if sp.exists():
        try:
            sd = json.loads(sp.read_text(encoding="utf-8")).get("spatialde") or {}
            st = sd.get("status")
            chk("svg:spatialde", "honesty", st == "ok",
                f"SpatialDE（§3.4）status={st}：{sd.get('reason') or ''}"
                + (f"，{sd.get('n_genes_run')} 个基因" if st == "ok" else ""))
            cmp_ = json.loads(sp.read_text(encoding="utf-8")) \
                .get("spatialde_vs_morans_i") or {}
            bg = cmp_.get("background_spearman_rho")
            chk("svg:spatialde_vs_morans", "honesty", bg is not None,
                (f"与 Moran's I 的一致性：背景组 rho={bg}"
                 f"（top 组 {cmp_.get('top_morans_spearman_rho')}）"
                 if bg is not None else
                 f"未量化：{cmp_.get('reason', '缺 spatialde_vs_morans_i')}"))
        except Exception as e:  # noqa: BLE001
            chk("svg:spatialde", "honesty", False, f"读取失败: {e}")

    # ---- honesty: §3.2 点名的 SpaGCN 是否真跑了、跑了是否可比 ----------------
    #
    # **和 SpatialDE 那条同构。** SpaGCN 这一轮才真正跑起来（`init="kmeans"`
    # 绕开了 louvain 的 py3.12 编译链），所以它的失败模式是新的：
    # torch 装不上、SpaGCN 装不上、kmeans 收敛不了 —— 任何一种都会让它退回
    # "内置实现 + 一句没做"。**那种情况必须红字显示**，否则一个漂亮的域划分
    # 会被当成"§3.2 的三个方法都跑了"。
    #
    # 跑了的情况：产物必须齐（CSV 有行 + 与主方法的一致性被量化），
    # 因为**孤立的一个划分说明不了任何事**。
    dp = res_dir / "domain_status.json"
    if dp.exists():
        try:
            dd = json.loads(dp.read_text(encoding="utf-8"))
            spg = (dd.get("domain_methods") or {}).get("SpaGCN") or {}
            st = spg.get("status")
            if st == "ok":
                csv_p = res_dir / "spagcn_domains.csv"
                n_row = (sum(1 for _ in open(csv_p, encoding="utf-8")) - 1
                         if csv_p.exists() else 0)
                ag = (dd.get("method_agreement") or {}).get("SpaGCN_vs_builtin") or {}
                chk("domain:spagcn", "honesty", n_row > 0 and bool(ag),
                    (f"SpaGCN（§3.2）status=ok：{spg.get('n_domains')} 域"
                     f"（kmeans init，n_clusters={spg.get('n_clusters')}，"
                     f"l={spg.get('length_scale_l')}），"
                     f"与主方法 ARI={ag.get('adjusted_rand_index')}、"
                     f"邻居同域率 {ag.get('neighbor_same_frac_a')} vs "
                     f"{ag.get('neighbor_same_frac_b')}，"
                     f"spagcn_domains.csv {n_row} 行"
                     if n_row > 0 and ag else
                     f"SpaGCN 报 status=ok 但产物不全："
                     f"spagcn_domains.csv {n_row} 行，method_agreement {bool(ag)}"))
            else:
                chk("domain:spagcn", "honesty", False,
                    f"SpaGCN（§3.2）**未跑成** status={st}：{spg.get('reason')} "
                    f"—— 上面的域划分是内置实现，不是 SpaGCN")
        except Exception as e:  # noqa: BLE001
            chk("domain:spagcn", "honesty", False, f"读取失败: {e}")

    # ---- honesty: §3.2/§3.3/§3.5/§3.6 点名的工具逐条登记 ----------------------
    #
    # **这是本仓库最容易被误读的地方。** 域划分、解卷积、通讯、轨迹四步都
    # 有产出、都有图，看起来"§3 做完了" —— 但文档点名的工具里，
    # **只有 SpaGCN 和 SpatialDE 真的跑了**，其余（BayesSpace / STAGATE /
    # RCTD / cell2location / CellChat / StPedf / SpaceFlow / ISORT /
    # Stereopy-TGPI / stLearn）都跑不了或没参考。
    #
    # 所以这几条检查的判据是"**每个点名工具都在登记里、且要么 `used: True`
    # 要么写了理由**"，不是"工具有没有跑"。工具将来能装了，这几条应该
    # **依然 PASS**（理由变成"已跑"），而不是因为"available=False"就变红 ——
    # 那会把"如实记录"惩罚成失败。
    #
    # **登记位置不唯一。** §3.2 的 SpaGCN 跑起来之后从 `named_tools` 挪到了
    # `domain_methods`（它不再是"用不了的工具"）。所以这里把两个表合并起来
    # 看 —— 只认 `named_tools` 会把"已经跑了的工具"报成"缺登记"，
    # 正好把好事判成坏事。
    tool_registry_files = {
        "domain_status.json": ("§3.2 空间域", ("BayesSpace", "STAGATE", "SpaGCN"),
                               (("domain_methods", None),)),
        "deconvolution_status.json": ("§3.3 解卷积", ("RCTD", "cell2location"),
                                      (("cell2location", "cell2location"),)),
        "svg_status.json": ("§3.4 SVG", ("SPARK-X", "SpatialDE2"), ()),
        "communication_status.json": ("§3.5 空间通讯", ("CellChat",), ()),
        "spatial_trajectory_status.json": ("§3.6 空间轨迹",
                                           ("StPedf", "SpaceFlow", "ISORT",
                                            "Stereopy-TGPI", "stLearn"), ()),
    }
    for f, (section, expect, extra_keys) in tool_registry_files.items():
        p = res_dir / f
        if not p.exists():
            continue
        try:
            d = json.loads(p.read_text(encoding="utf-8"))
            nt = d.get("named_tools")
        except Exception as e:  # noqa: BLE001
            chk(f"named_tools:{f}", "honesty", False, f"读取失败: {e}")
            continue
        if not isinstance(nt, dict) or not nt:
            chk(f"named_tools:{f}", "honesty", False,
                f"{section}：{f} 里没有 named_tools 登记 —— "
                "读者会以为文档点名的方法已经用上了")
            continue
        registry = dict(nt)
        for k, as_tool in extra_keys:
            v = d.get(k)
            if not isinstance(v, dict) or not v:
                continue
            if as_tool is None:
                # 本身就是 {工具名: 信息} 的表。**只接受值是 dict 的项** ——
                # 上游改形状时宁可少合并，也不要塞进字符串把验收炸掉
                registry.update({t: i for t, i in v.items()
                                 if isinstance(i, dict)})
            else:
                registry[as_tool] = v
        missing = [t for t in expect if t not in registry]
        # 未使用必须有 reason；使用了必须有 used=True —— 两者都没有 = 没说清
        unexplained = [t for t, i in registry.items()
                       if isinstance(i, dict)
                       and not i.get("reason") and not i.get("used")]
        n_used = sum(1 for i in registry.values()
                     if isinstance(i, dict) and i.get("used"))
        chk(f"named_tools:{f}", "honesty",
            not missing and not unexplained,
            (f"{section}：{len(registry)} 个登记项"
             f"（点名 + 各步自己的方法表），"
             f"其中 {n_used} 个实际产出了结果"
             if not missing and not unexplained else
             f"缺登记 {missing}；既没说 used 也没写理由 {unexplained}"))

    # ---- §0.4 Agent 决策链：点名工具到底落地了几个 --------------------------
    #
    # **决策链不是日志。** 日志会随 CI 滚动消失，而"这一轮文档点名的
    # §3.2–§3.6 工具，哪些真跑了、哪些是回退"是**结论的适用范围** ——
    # 半年后拿到 domain_status.json 的人必须能直接看到，
    # 而不是去翻几十个 JSON 自己数。
    #
    # 放在验收里而不是 run_steps 里：只有这里才把 named_tools 与各步的
    # domain_methods / cell2location / spatialde **合并**过，
    # 早写会漏掉"实际跑了"的工具。
    try:
        _tools, _used, _named = {}, [], set()
        for f, (_sec, _exp, extra_keys) in tool_registry_files.items():
            _p = res_dir / f
            if not _p.exists():
                continue
            _d = json.loads(_p.read_text(encoding="utf-8"))
            _reg = dict(_d.get("named_tools") or {})
            # **`named_tools` 里的才是"文档点名的工具"。** 各步自己的表
            # （`domain_methods` / `cell2location`）还包含**内置基线**
            # （`builtin_smooth_leiden`）—— 它当然产出了结果，但它不是
            # 文档点名的东西。混在一起数会把"点名 13 个"报成"14 个"。
            _named.update(_reg)
            for _k, _as_tool in extra_keys:
                _v = _d.get(_k)
                if not isinstance(_v, dict) or not _v:
                    continue
                # 和上面同一条规则：`None` 是"工具名 -> 信息"的表，
                # 否则是**一个**工具的信息（整个挂到它名下）。混着 update
                # 会把字段名当成工具名，值还是字符串。
                if _as_tool is None:
                    _reg.update({t: i for t, i in _v.items()
                                 if isinstance(i, dict)})
                else:
                    _reg[_as_tool] = _v
            for _t, _i in _reg.items():
                if not isinstance(_i, dict):
                    continue
                _tools[_t] = _i
                if _i.get("used"):
                    _used.append(f"{_t}（{_sec}）")
        _n = len(_tools)
        _n_named = len(_named)
        _n_extra = _n - _n_named
        record_decision(
            cfg, "named_tools",
            "文档 §3.2–§3.6 点名的空间方法，哪些真的产出了结果？",
            (f"点名工具登记 {_n_named} 个"
             + (f"（另有内置基线等 {_n_extra} 个）" if _n_extra else "")
             + f"；实际产出结果的是 {sorted(_used) or '无'}"
             if _used else
             f"点名工具登记 {_n_named} 个，"
             "**本轮没有一个点名工具产出结果** —— 全部回退到内置实现"),
            evidence=("逐条理由写在各步状态 JSON 的 named_tools / "
                      "domain_methods / cell2location 里；"
                      "SpaGCN 是 §3.2 唯一能真跑的点名工具（init=\"kmeans\"），"
                      "cell2location 属 needs_reference（缺参考，不是装不上）。"
                      "**`builtin_smooth_leiden` 不是点名工具** —— 它是本仓库的"
                      "基线，出现在产出结果名单里是因为它确实产出了结果，"
                      "不是因为文档点了它"))
        log_info(f"决策链已登记：点名工具 {_n_named} 个"
                 f"（另有内置基线等 {_n_extra} 个），"
                 f"实际产出结果的 {len(_used)} 个")
    except Exception as e:  # noqa: BLE001
        # 决策链写不进去不该让验收崩 —— 但必须可见
        log_warn(f"决策链登记失败（不影响其余验收）：{type(e).__name__}: {e}")

    # ---- honesty: §3.2 跑了点名工具就必须量化它与主方法的一致性 --------------
    #
    # **"跑通了"不是结论。** 一个点名工具跑出 13 个域，如果不说它和主方法
    # 的 ARI，读者只能看到一个孤立的划分 —— 而"两套划分是不是在说同一件事"
    # 才是这次交叉验证的全部意义。这条检查把"跑了但没比"变成红的。
    p = res_dir / "domain_status.json"
    if p.exists():
        try:
            d = json.loads(p.read_text(encoding="utf-8"))
            dm = d.get("domain_methods") or {}
            ag = d.get("method_agreement") or {}
            named_ran = [t for t, i in dm.items()
                         if (i or {}).get("used") and t != "builtin_smooth_leiden"]
            uncompared = [t for t in named_ran if f"{t}_vs_builtin" not in ag]
            chk("honesty:domain_method_agreement", "honesty",
                not uncompared,
                (f"§3.2 实际产出结果的点名方法 {named_ran}；"
                 f"与主方法的一致性对照 {sorted(ag)}"
                 if not uncompared else
                 f"{uncompared} 跑了但**没有量化与主方法的一致性** —— "
                 f"孤立的一个划分说明不了任何事"))
        except Exception as e:  # noqa: BLE001
            chk("honesty:domain_method_agreement", "honesty", False, f"读取失败: {e}")

    # ---- honesty: 每步必须说明自己没做什么 ----
    # 这是防止"静默跳过"的检查。
    honest_map = {
        "qc_status.json": "threshold_note",
        "domain_status.json": "limitations",
        "svg_status.json": "limitations",
        "deconvolution_status.json": "limitations",
        "niche_status.json": "limitations",
        "communication_status.json": "limitations",
        "spatial_trajectory_status.json": "limitations",
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

    # ---- honesty: §3.3 点名的 cell2location 落地状态必须自洽 ------------------
    #
    # **`needs_reference` 和"装不上"是两件事。** cell2location 的 PyPI 包
    # 装得上，缺的是带细胞类型标签的 scRNA 参考。把这两件事混成一句
    # "没做"，下一个人就会去折腾安装 —— 方向完全错。
    #
    # 判据：
    #   - status=ok        -> 必须有比例 CSV 且与 NNLS 的一致性被量化
    #   - needs_reference  -> 必须说清"缺的是什么参考"
    #   - 其它             -> 必须有 reason
    p = res_dir / "deconvolution_status.json"
    if p.exists():
        try:
            d = json.loads(p.read_text(encoding="utf-8"))
            c2l = d.get("cell2location") or {}
            st = c2l.get("status")
            if st == "ok":
                csv_p = res_dir / "deconvolution_proportions_cell2location.csv"
                n_row = (sum(1 for _ in open(csv_p, encoding="utf-8")) - 1
                         if csv_p.exists() else 0)
                cmp_ = c2l.get("vs_nnls") or {}
                chk("deconv:cell2location", "honesty",
                    n_row > 0 and (cmp_.get("mean_spearman") is not None
                                   or cmp_.get("compared") is False),
                    (f"cell2location（§3.3）status=ok：v{c2l.get('version')}，"
                     f"{c2l.get('n_types')} 种类型，{n_row} 行，"
                     f"与 NNLS 平均 Spearman {cmp_.get('mean_spearman')}"
                     if n_row > 0 else
                     f"cell2location 报 status=ok 但没有比例 CSV（{n_row} 行）"))
            elif st == "needs_reference":
                chk("deconv:cell2location", "honesty", bool(c2l.get("reason")),
                    f"cell2location（§3.3）status=needs_reference："
                    f"{str(c2l.get('reason'))[:150]}")
            else:
                chk("deconv:cell2location", "honesty", bool(c2l.get("reason")),
                    f"cell2location（§3.3）status={st}："
                    f"{str(c2l.get('reason'))[:150]}")
        except Exception as e:  # noqa: BLE001
            chk("deconv:cell2location", "honesty", False, f"读取失败: {e}")

    # ---- 空间对齐检查（如果有）----
    p = data_dir / "spatial_alignment_check.json"
    if p.exists():
        d = json.loads(p.read_text(encoding="utf-8"))
        chk("spatial:alignment", "content", d.get("current_is_best", False),
            f"最佳假设: {d.get('best_hypothesis')}；"
            f"当前是否最佳: {d.get('current_is_best')}", severity="content")

    # ---- 空间拟时序：读 status 的真实字段，不是只看文件在不在 ------------
    #
    # 姊妹项目 Part 1 的教训：TF 第一次跑时 tf_status.json 正常产出、
    # 验收全绿，而 figure_written 其实是 false（图一张没出）。
    p = res_dir / "spatial_trajectory_status.json"
    if p.exists():
        d = json.loads(p.read_text(encoding="utf-8"))
        if d.get("status") == "ok":
            has_both = (d.get("morans_I_expression_only") is not None
                        and d.get("morans_I_spatially_smoothed") is not None)
            chk("content:spatial_traj_comparison", "content", has_both,
                f"Moran's I 朴素 {d.get('morans_I_expression_only')} → "
                f"空间感知 {d.get('morans_I_spatially_smoothed')}"
                f"（增益 {d.get('morans_I_gain')}）" if has_both
                else "**没有记录朴素 vs 空间感知的对比数值** —— "
                     "那就只剩一个排序，说明不了空间感知有没有用",
                severity="content")
            ntu = d.get("named_tools_not_used") or {}
            chk("honesty:named_tools", "honesty", len(ntu) >= 3,
                f"写明了 {len(ntu)} 个具名工具未使用及原因：{','.join(list(ntu)[:5])}"
                if len(ntu) >= 3 else
                f"**只写了 {len(ntu)} 个** —— 文档点名的工具没用就要说明为什么",
                severity="honesty")
            chk("honesty:not_a_developmental_trajectory", "honesty",
                any("不是发育轨迹" in x or "空间排序" in x
                    for x in (d.get("limitations") or [])),
                "已写明这是空间排序而非发育轨迹"
                if any("不是发育轨迹" in x or "空间排序" in x
                       for x in (d.get("limitations") or []))
                else "**没有写明这不是发育轨迹** —— 一个空间排序会被读成分化过程",
                severity="honesty")
            chk("honesty:root_is_heuristic", "honesty",
                bool((d.get("root_selection") or {}).get("method")),
                f"选根方式: {(d.get('root_selection') or {}).get('method')}",
                severity="honesty")

    # ---- 模块零：运行清单（§0.3 / §0.4）-------------------------------------
    # 清单缺项不是"分析错了"，而是"这轮跑出来的东西没法追溯"。
    # **人工复核未确认不算失败** —— 默认就是 pending，那是设计如此；
    # 把它判成 FAIL 会让每个 job 都红，反而没人看。但必须可见。
    chk("manifest:present", "required", msum.get("present", False),
        (f"{msum.get('n_versions', 0)} 个包版本、{msum.get('n_inputs', 0)} 项输入、"
         f"{msum.get('n_decisions', 0)} 条决策"
         if msum.get("present") else "**缺 run_manifest.json**"))
    if msum.get("present"):
        chk("manifest:versions", "required", msum["n_versions"] >= 20,
            f"{msum['n_versions']} 个已安装包（pip freeze 全量）")
        chk("manifest:inputs", "required",
            # **只看 required 的缺失。** 可选项本来就可以不存在，
            # 算进来会让没有该文件的数据集全部误判失败。
            msum["n_inputs"] >= len(INPUT_FILES)
            and not msum["inputs_missing_required"],
            f"{msum['n_inputs']} 项"
            + (f"，必需缺失 {','.join(msum['inputs_missing_required'])}"
               if msum["inputs_missing_required"] else "，必需项齐全")
            + (f"（可选缺失 {','.join(msum['inputs_missing'])}）"
               if msum["inputs_missing"] else ""))
        pend = msum["human_review_pending"]
        chk("manifest:human_review", "honesty", True,
            f"{len(pend)} 个人工复核节点待确认（不阻断 job）："
            + (", ".join(pend) if pend else "全部已确认"))

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
        "manifest": msum,
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
