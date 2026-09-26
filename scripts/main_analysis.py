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
import re
import sys
import time
import traceback
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts" / "lib"))
sys.path.insert(0, str(REPO / "scripts"))

from common import (NAMED_TOOLS, capture_versions, ensure_dirs,  # noqa: E402
                    init_manifest, load_config, log_info, log_warn,
                    manifest_path, manifest_summary, probe_named_tools,
                    read_manifest, record_decision,
                    record_human_review, record_input, record_params,
                    write_json)

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

# 嵌套 status 里，哪些取值算"崩了"。其余（`needs_reference` /
# `package_missing` / `not_done` / `not_applicable` / `disabled` /
# `missing_domains` / `bad_root` …）都是**设计如此地没做**，只可见、不阻断。
NESTED_FAILED_VALUES = ("failed", "error", "fail")

# 本仓库的阶段号（geo=01 / scrna=02 / spatial=03）。
PART = "03"


def _strip_comments(src: str) -> str:
    """逐行剥注释，**引号内不剥** —— 与 `tools/check_fig_names.mjs` 的
    `stripComments` 同义。

    口径必须一致：门禁层用 JS 那份扫"声明了哪些图"，验收层用这份扫，
    两边算法不同就会对同一份源码给出不同的图名集合，而**没有任何东西
    能发现它们不一致**（门禁绿、验收也绿）。
    """
    out = []
    for line in src.split("\n"):
        q, cut = None, None
        for i, c in enumerate(line):
            if q:
                if c == q:
                    q = None
            elif c in "\"'":
                q = c
            elif c == "#":
                cut = i
                break
        out.append(line if cut is None else line[:cut])
    return "\n".join(out)


def _script_sources():
    for f in sorted((REPO / "scripts").glob("[0-9][0-9]_*.py")):
        # **必须排除 main_analysis.py**：本文件的验收层会引用图名，
        # 那是"检查对象"不是"出图声明"。不排除的话验收会要求自己
        # 引用过的每张图都存在，把口径搞反。
        if f.name == "main_analysis.py":
            continue
        yield f, _strip_comments(f.read_text(encoding="utf-8"))


_FIG_NAME_RE = re.compile(
    rf"^{PART}-\d{{2}}-\d{{2}}-unit\d+-[a-z0-9]+(?:-[a-z0-9]+)*$")


def declared_figures() -> list:
    """扫出**声明要出**的静态图名（字符串字面量）。

    Q-27 / E-48：原来图验收是 `len(figs) >= 8` —— 一个纯计数。**"该有的
    图没有"这一整类问题没有任何检查看得见**：姊妹仓库实测 5 张图因
    `figsize` 三元素元组从未产出过，而验收 70 项全绿。

    所以判据从"数量够不够"换成"**声明过的每张图在不在**"。声明从源码扫出来
    而不是手抄一张表 —— 手抄的表会漂移，而漂移的方向恰好是"新加的图
    不在表里"，也就是把 E-48 那个盲区原样再造一遍。

    含 `{ }` 的模板串跳过（运行时拼名，由 `DYNAMIC_FIG_BASES` 声明豁免）。
    """
    out = []
    pat = re.compile(rf'"{PART}-[^"]*"')
    for _f, src in _script_sources():
        for m in pat.finditer(src):
            nm = re.sub(r"\.(pdf|png)$", "", m.group(0)[1:-1])
            if "{" in nm or "}" in nm:
                continue
            if _FIG_NAME_RE.match(nm):
                out.append(nm)
    return sorted(set(out))


def dynamic_fig_bases() -> dict:
    """扫出 `DYNAMIC_FIG_BASES = {"<图号>": <张数>}` 声明。

    键补全成完整前缀 `03-<模块号>-<图号>`。这是**槽位上限**，不是精确值：
    实测 `03-03-04` 声明 3、实际只出 1（三对方法里只有一对算出了
    Jaccard 矩阵），少出是合法的，多出才说明有未声明的名字。
    """
    out = {}
    pat = re.compile(r"DYNAMIC_FIG_BASES\s*=\s*\{([^}]*)\}")
    for f, src in _script_sources():
        m = pat.search(src)
        if not m:
            continue
        for pair in m.group(1).split(","):
            if ":" not in pair:
                continue
            k, v = pair.split(":", 1)
            key = k.strip().strip("\"'")
            try:
                out[f"{PART}-{f.name[:2]}-{key}"] = int(v.strip())
            except ValueError:
                continue
    return out


# 条件产出的图：**图名 -> (状态文件, 判据路径, 能力就位时的取值, 说明)**。
#
# 为什么需要这张表：有些图**本来就可以合法地不存在**，因为它的代码路径
# 挂在一个能力探针后面。判据不能是"这张图必须在"，也不能是"把它从
# 声明里删掉"（删掉就等于永远不再检查它）。
#
# 语义：判据路径取到的值 == "能力就位"的值时，这张图**变成必需**；
# 否则允许缺失，但必须**可见**。
#
# **这不是"已知缺陷白名单"** —— 它写的是"为什么可以没有"，并且是
# **自愈**的：STAGATE 哪天装上了、`status` 变成 `ok`，这张图立刻
# 自动变成必需。白名单会把缺陷永久合法化，这张表不会。
CONDITIONAL_FIGURES = {
    "03-03-04-unit3-domains-stagate": (
        "domain_status.json", ("domain_methods", "STAGATE", "status"), "ok",
        "STAGATE_pyG 装不上时（`package_missing`）本就不该有这张图"),
    "03-05-03-unit1-deconvolution-error-map": (
        "deconvolution_status.json", ("is_deconvolution",), True,
        "没有带细胞类型标签的参考时走 marker 打分法，没有重建误差可画"),
}


def _dig(obj, path):
    """按路径取值；任一层缺失返回 `None`（**不抛异常**）。

    判据路径写错时必须表现为"探针读不到"而不是把验收整个炸掉 ——
    读不到时走的是"允许缺失、可见"分支，不会静默判成通过。
    """
    cur = obj
    for k in path:
        if isinstance(cur, dict):
            cur = cur.get(k)
        else:
            return None
        if cur is None:
            return None
    return cur


def _iter_nested_status(obj, path=()):
    """递归产出所有名为 `status` 的字段：`(路径, 值, 同级 reason)`。

    Q-27 / E-48：`status` 不只在顶层。`domain_status.json` 的
    `domain_methods.STAGATE.status`、`deconvolution_status.json` 的
    `cell2location.status`、`svg_status.json` 的 `spatialde.status` 都是
    **嵌套**的，而旧验收层只看顶层 —— 里面崩了、顶层还是 `ok`。

    结构上这和"`_content_overflow()` 打了 WARN 没人读"是同一个错误：
    `except` 把异常降级成了一个**没人读的字段**。这里把它读出来。
    """
    if isinstance(obj, dict):
        if "status" in obj:
            reason = obj.get("reason") or obj.get("message") or ""
            yield path, obj["status"], reason
        for k, v in obj.items():
            if isinstance(v, (dict, list)):
                yield from _iter_nested_status(v, path + (str(k),))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            if isinstance(v, (dict, list)):
                yield from _iter_nested_status(v, path + (str(i),))


def run_steps(cfg: dict, only: list = None) -> dict:
    ensure_dirs(cfg)
    res_dir = Path(cfg["output"]["results_dir"])

    # **溢出记录也要先清掉。** 和状态文件同理（规则 14）：图标题改短之后
    # 旧记录还在，验收会把**已经修好的**图判红 —— 那比没有检查更糟，
    # 因为它会训练人忽略红字。
    ovf_stale = res_dir / "figure_overflow.json"
    if ovf_stale.exists():
        ovf_stale.unlink()

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

    # ---- required: 步骤函数自己返回的 status 也要看（Q-27）------------------
    #
    # `run_steps` 把 `fn(cfg)` 的返回值记进了 `result_status`，但**全仓零
    # 消费者** —— 于是 `08_spatial_trajectory.py` 返回 `bad_root` /
    # `not_applicable` / `missing_pca` 时，验收照样记 `status="ok"`。
    #
    # **这不是"步骤崩了"，而是"步骤跑完了、但结果是『没做成』"** ——
    # 恰恰是最容易读成成功的一种。判据仍只把 `failed`/`error`/`fail` 判红：
    # `bad_root` 之类的取值是设计如此地没做成，只可见、不阻断。
    for sid, res in step_results.items():
        rs = res.get("result_status")
        if rs is None or rs == "ok":
            continue
        bad = str(rs).lower() in NESTED_FAILED_VALUES
        # 注意 `chk` 的第二个位置参数是 **kind**（归类），severity 是关键字参数 ——
        # 把 kind 当 severity 传会让"设计如此地没做成"也按 required 记账。
        chk(f"step_result:{sid}", "step_result", not bad,
            f"步骤 {sid} 返回 status={rs!r}"
            + ("—— **步骤跑完了但结果是失败**" if bad
               else "（设计如此地没做成，可见不阻断）"),
            severity="required" if bad else "honesty")

    # ---- required: 嵌套 status 的 failed 必须判红（Q-27 / E-48）-------------
    #
    # `status` 不只在顶层。`domain_status.json` 的
    # `domain_methods.STAGATE.status`、`deconvolution_status.json` 的
    # `cell2location.status`、`svg_status.json` 的 `spatialde.status` 都是
    # **嵌套**的，而旧验收层只看顶层 —— **里面崩了、顶层还是 `ok`**。
    #
    # 姊妹仓库实测（E-48）：顶层分布 `{ok: 7, not_configured: 1}`，
    # 嵌套分布里躺着唯一一条 `failed`（`figsize` 三元素元组让整段抛异常），
    # 而它所属文件的顶层 `status` 是 `ok`。**顶层全绿、里面已经崩了。**
    #
    # 判据必须**只**把 `failed`/`error`/`fail` 判红 —— 把
    # `needs_reference` / `package_missing` / `not_done` / `not_applicable`
    # / `disabled` / `missing_domains` / `bad_root` 判红会让每个 job 都红，
    # 反而没人看（规则 4 的同一条理由）。
    nested_bad, nested_seen = [], {}
    for sp in sorted(res_dir.glob("*status.json")):
        try:
            obj = json.loads(sp.read_text(encoding="utf-8"))
        except Exception as e:  # noqa: BLE001
            nested_bad.append(f"{sp.name}: 读取失败 {e}")
            continue
        for path, val, reason in _iter_nested_status(obj):
            key = ".".join(path) + ".status" if path else "status"
            nested_seen[str(val)] = nested_seen.get(str(val), 0) + 1
            if str(val).lower() in NESTED_FAILED_VALUES:
                where = f"{sp.name} → {key}"
                nested_bad.append(f"{where} = {val!r}"
                                  + (f"（{reason[:120]}）" if reason else ""))
    chk("status:nested", "required", not nested_bad,
        (f"全部 status 字段（含嵌套）无 failed，分布 {nested_seen}"
         if not nested_bad else
         f"**嵌套 status 里有失败**（顶层可能是 ok）：{nested_bad}；"
         f"全部取值分布 {nested_seen}"))

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

    # ---- required: 图 -------------------------------------------------------
    #
    # **从"数量够不够"换成"声明过的每张图在不在"。** 旧判据是
    # `len(figs) >= 8` —— 纯计数，于是"该有的图没有"这一整类问题没有任何
    # 检查看得见：姊妹仓库实测 5 张图因 `figsize` 三元素元组从未产出过，
    # 而验收 70 项全绿（E-48）。
    #
    # 声明从源码扫出来（`declared_figures()`），不手抄表 —— 手抄的表会漂移，
    # 漂移方向恰好是"新加的图不在表里"，等于把盲区原样再造一遍。
    figs = sorted(fig_dir.glob("*.png")) if fig_dir.exists() else []
    fig_names = {p.stem for p in figs}
    declared = declared_figures()
    dyn = dynamic_fig_bases()
    dyn_slots = sum(dyn.values())

    # 动态图名（`DYNAMIC_FIG_BASES` 声明的槽位）不逐名检查 —— 名字运行期才
    # 拼得出来。但它们**必须至少产出 1 张**：声明了槽位却一张都没有，
    # 说明那段循环整段没跑（E-48 的形态正是"整段没跑而验收不知道"）。
    #
    # **必须扣掉该前缀下的静态声明图。** `03-03-04` 这个前缀同时有静态图
    # （unit1/unit2/unit3）和动态图（unit8 的 Jaccard 矩阵）—— 不扣的话，
    # 动态循环整段没跑时计数仍是 2（全是静态的），判据**永远不会响**。
    # 这类"判据因为输入域重叠而永远为真"的错误在正向标定里看不出来，
    # 只有注入"整组动态图消失"才会暴露。
    declared_set = set(declared)
    dyn_missing = []
    for base, n_slots in dyn.items():
        got = [nm for nm in fig_names
               if nm.startswith(base + "-") and nm not in declared_set]
        if not got:
            dyn_missing.append(f"{base}（声明 {n_slots} 个槽位，实际 0 张）")

    # 条件产出的图（`CONDITIONAL_FIGURES`）：能力探针就位时变必需，
    # 否则允许缺失但**必须可见**。
    missing, waived = [], []
    for nm in declared:
        if nm in fig_names:
            continue
        cond = CONDITIONAL_FIGURES.get(nm)
        if cond:
            st_file, path, ready_val, why = cond
            p = res_dir / st_file
            got = None
            if p.exists():
                try:
                    got = _dig(json.loads(p.read_text(encoding="utf-8")), path)
                except Exception:  # noqa: BLE001
                    got = None
            if got != ready_val:
                waived.append(f"{nm}（{why}；{st_file} "
                              f"{'.'.join(str(x) for x in path)}={got!r}）")
                continue
        missing.append(nm)

    chk("figures:declared", "required", not missing,
        (f"声明 {len(declared)} 张静态图，全部产出"
         + (f"；另有 {len(waived)} 张条件图本轮不适用" if waived else "")
         if not missing else
         f"**声明了但没产出** {missing}"
         + (f"（条件图本轮不适用：{waived}）" if waived else "")
         + f" —— 实际产出 {len(figs)} 张"),)
    chk("figures:dynamic", "required", not dyn_missing,
        (f"{len(dyn)} 组动态图名共 {dyn_slots} 个槽位，各自至少产出 1 张"
         if not dyn_missing else
         f"**动态图名整组没产出** {dyn_missing} —— 槽位是上限不是精确值，"
         f"少出合法，但**一张都没有说明那段循环整段没跑**"),)
    # 保留计数作为**下限兜底**：声明扫描本身失效时（如源码结构大改导致
    # 一条字面量都扫不到）这条还能拦住"一张图都没有"。
    chk("figures:count", "required", len(figs) >= 8,
        f"{len(figs)} 张图（要求 >=8）")

    # ---- required: 溢出记录（E-49 的形态）------------------------------------
    #
    # `_content_overflow()` **正确检测到了**超宽、**正确打了 WARN**，
    # 然后**没有任何人读** —— 标题超宽 35%，在 `savefig.bbox: standard` 下
    # 被静默裁掉，图照样生成、门禁照样绿。
    #
    # **"检测到了"不等于"有人会知道"。** 所以溢出落盘成
    # `figure_overflow.json`，这里读它并判红。
    ovf_p = res_dir / "figure_overflow.json"
    ovf = {}
    if ovf_p.exists():
        try:
            ovf = json.loads(ovf_p.read_text(encoding="utf-8")).get("figures") or {}
        except Exception as e:  # noqa: BLE001
            ovf = {"<读取失败>": str(e)}
    chk("figures:overflow", "required", not ovf,
        "没有图内容超出画布" if not ovf else
        f"**{len(ovf)} 张图内容超出画布，会被静默裁掉**：{ovf}")

    # ---- content: 图不是空白 ----
    # 用像素标准差判断：全白/全黑的图标准差接近 0。
    # **这个检查是必要的**：matplotlib 在数据为空时会静默产出一张空白图，
    # 文件存在、大小正常，看不出问题。
    blank = []
    blank_skipped = None
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
    except ImportError as e:  # noqa: BLE001
        # **T4：「跳过」不等于「通过」。** 原先 `blank` 被塞进一条 note，
        # 于是 `not blank` 为 False → 记成失败，但 detail 写的是"跳过"，
        # 读起来像通过。现在拆成两个字段：真空白进 `blank` 判红，
        # PIL 缺失单独记 `blank_skipped`，**不伪装成通过、也不伪装成失败**。
        blank_skipped = f"PIL 不可用（{e}）—— 空白图检查没有执行"
    if blank_skipped:
        chk("figures:non_blank", "content", False,
            f"**空白图检查未执行**：{blank_skipped} —— 这不是通过",
            severity="content")
    else:
        chk("figures:non_blank", "content", not blank,
            "所有图都有内容" if not blank else f"疑似空白图: {blank}",
            severity="content")

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

    # ---- 全局：`NAMED_TOOLS` 里每一个都必须有交代 ---------------------------
    #
    # **上面那条是按步查的，查不出"整段没登记"。** 实测踩到：
    # `Bering` 与 `BOMS` 的 `section` 是 §3.1，而 §3.1 的两个步骤
    # （`01_qc` / `02_normalize`）**没有 `named_tools` 字段** ——
    # 于是这两个工具只存在于 `common.NAMED_TOOLS` 这个 Python 常量里，
    # **任何产物里都看不到它们**。按步检查全都通过，因为没有任何一步
    # "声明"过它们。
    #
    # 这正是"缺口没登记"最隐蔽的形态：不是登记错了，是**根本没登记**，
    # 而验收看起来是绿的。所以这里做一次全局清点。
    #
    # 判据：`NAMED_TOOLS` 的每个键，要么出现在某一步的 `named_tools` /
    # 方法表里（说明该步真的探测过），要么在下面的全局探测里有 `kind`
    # 和 `reason`（说明它是"文档点名但本环境用不了"，理由齐全）。
    try:
        _all_nt = dict(NAMED_TOOLS)
        _seen = set()
        for f, (_sec, _exp, extra_keys) in tool_registry_files.items():
            _p = res_dir / f
            if not _p.exists():
                continue
            _d = json.loads(_p.read_text(encoding="utf-8"))
            _seen.update((_d.get("named_tools") or {}).keys())
            for _k, _as_tool in extra_keys:
                _v = _d.get(_k)
                if not isinstance(_v, dict) or not _v:
                    continue
                if _as_tool is None:
                    _seen.update(t for t, i in _v.items()
                                 if isinstance(i, dict))
                else:
                    _seen.add(_as_tool)
        # `NAMED_TOOLS` 自己就是"文档点名 + 逐条理由"的登记表 ——
        # 所以缺的是"有没有被任何一步碰到过"这件事的可见性。
        _never_probed = sorted(t for t in _all_nt if t not in _seen)
        # 理由齐全性：`NAMED_TOOLS` 里每条都必须有 reason，且 kind 非空
        _no_reason = sorted(t for t, i in _all_nt.items()
                            if not (i or {}).get("reason"))
        _no_kind = sorted(t for t, i in _all_nt.items()
                          if not (i or {}).get("kind"))
        chk("honesty:named_tools_registry", "honesty",
            not _no_reason and not _no_kind,
            (f"`NAMED_TOOLS` 共 {len(_all_nt)} 条，全部有 kind 与 reason"
             if not _no_reason and not _no_kind else
             f"无理由 {_no_reason}；无 kind {_no_kind}"))
        # 没被任何一步探测过的：**可见但不判失败** ——
        # 它们是"文档点名、本环境用不了"，理由在 `NAMED_TOOLS` 里，
        # 只是没有哪一步负责把它们写进产物。判失败会让 job 红，
        # 但那不是分析错了。**必须可见**，否则下次又忘了。
        #
        # **T3：判据不能写死 True。** 原先第三个参数是字面量 `True`，
        # 于是这条检查在 `acceptance_report.json` 里**永远 `passed: true`**
        # —— `detail` 里明明列着"3 个工具没被探测过"，`passed` 却是真，
        # 任何按 `passed` 汇总的下游（含本文件的计数器）都会把它读成通过。
        # 现在判据是真实条件，但 severity 仍是 `honesty`（**不进 `passed`**），
        # 所以"不阻断 job"这个设计意图不变。
        chk("honesty:named_tools_never_probed", "honesty", not _never_probed,
            (f"全部 {len(_all_nt)} 个点名工具都至少被一步探测过"
             if not _never_probed else
             f"**{len(_never_probed)} 个点名工具没有任何一步探测过**"
             f"（只存在于 `common.NAMED_TOOLS` 常量里，产物中看不到）："
             f"{_never_probed} —— 逐条理由见清单 `named_tools` 决策"),
            severity="honesty")
        # **把全表探测结果写进清单。** 上面那条只报"哪些没被步骤碰到"，
        # 而这里把**全部 14 条连同 kind / section / reason** 落进
        # `run_manifest.json` —— 这样 §3.1 那种"没有任何步骤负责"的
        # 工具也有一条可查的理由，而不是只活在 Python 常量里。
        #
        # 不加 `only=` 过滤：这条就是全局清点，过滤掉谁都会漏。
        _probe_all = probe_named_tools(log=log_warn)
        record_params(cfg, {
            "named_tools_registry": {
                "n": len(_all_nt),
                "never_probed_by_any_step": _never_probed,
                "sections": sorted({(i or {}).get("section") or "?"
                                    for i in _all_nt.values()}),
            },
            "named_tools_probe": _probe_all,
        })
        log_info(f"点名工具全局清点：{len(_all_nt)} 条，"
                 f"其中 {len(_never_probed)} 条没有任何步骤探测过"
                 f"（理由已随清单落盘）")
    except Exception as e:  # noqa: BLE001
        log_warn(f"点名工具全局清点失败（不影响其余验收）："
                 f"{type(e).__name__}: {e}")

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

    # ---- §0.2 跨语言转换：本仓库没有，但要**说出来** ------------------------
    #
    # 姊妹项目 geo 是 R，它把 Part 1→Part 2 的 CSV 交接记进 `cross_language`
    # （`09_export_targets.R`）；scrna 也记了它读 Part 1 交接表那一段。
    #
    # **空间这边真的没有跨语言转换**：输入是 10x Cell Ranger 产出的
    # `filtered_feature_bc_matrix.h5` + `spatial.tar.gz`（Cell Ranger 是
    # 独立的 C++/Python 工具链，不是 R），Part 2 交接的参考也是 h5ad。
    #
    # 但"没有记录"和"忘了记"在清单里长得一样 —— 所以这里显式登记一条
    # 决策，把空数组解释掉。验收里 `manifest:cross_language` 会检查这一点：
    # **要么有记录，要么有解释。**
    record_decision(
        cfg, "cross_language",
        "§0.2 跨语言转换：本仓库有没有 R↔Python 的交接？",
        "**没有。** 全流程 Python，`cross_language` 为空是设计如此，不是遗漏",
        evidence=("输入是 10x Cell Ranger 的 `filtered_feature_bc_matrix.h5` "
                  "+ `spatial.tar.gz`（Cell Ranger 是独立工具链，不是 R）；"
                  "Part 2 交接的 scRNA 参考也是 `.h5ad`（同为 Python 生态）。"
                  "对比：geo 用 `record_cross_language()` 记 R→CSV→Python 的"
                  "靶基因表，scrna 记它读那张表 —— 两边都真有转换。"
                  "**所以这条不是「没做」，是「没有可做的」** —— "
                  "将来若加入 R 侧分析（如 BayesSpace），这条必须变成真记录"))

    # ---- content: 域标签打分有没有产出（S2）----------------------------------
    #
    # `03_spatial_domains.py` 的整段域标签打分包在裸 `try/except` 里，
    # 失败只 `log_warn`；而 `domain_annotation.json` 在本文件里
    # **原本没有任何消费者** —— 于是"域标签没算出来"和"算出来了"在
    # `acceptance_report.json` 里长得一模一样（都没有这条 id）。
    #
    # **"检测到了"不等于"有人会知道"。** 这条把"没产出"变成可见的红。
    p = res_dir / "domain_annotation.json"
    if not p.exists():
        chk("content:domain_annotation", "content", False,
            "**`domain_annotation.json` 不存在 —— 域标签打分整段没有执行**"
            "（该步的 try/except 会吞掉异常，所以缺文件就是唯一的信号）",
            severity="content")
    else:
        try:
            d = json.loads(p.read_text(encoding="utf-8"))
            st = d.get("status")
            n_dom = len(d.get("domains") or {})
            chk("content:domain_annotation", "content", st == "ok" and n_dom > 0,
                (f"域标签: {n_dom} 个域有标签（{d.get('signature_set')}）"
                 if st == "ok" and n_dom > 0 else
                 f"**域标签没有产出**：status={st}，{n_dom} 个域，"
                 f"原因 {str(d.get('reason'))[:140]}"),
                severity="content")
        except Exception as e:  # noqa: BLE001
            chk("content:domain_annotation", "content", False,
                f"读取失败: {e}", severity="content")

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
    #
    # **S4：没有 `else` 分支 = 文件缺失时这条检查根本不存在。**
    # 原先只在文件存在时 `chk` —— 于是"对齐检查从没跑过"和"对齐检查跑了
    # 且结论是最佳"在 `acceptance_report.json` 里看起来一样（都没有这条 id）。
    # `00_fetch.py` 的对齐检查会在 `current_is_best` 为假时 `raise`，
    # 所以正常情况下文件必然存在；**它不存在本身就是异常**。
    p = data_dir / "spatial_alignment_check.json"
    if p.exists():
        try:
            d = json.loads(p.read_text(encoding="utf-8"))
            chk("spatial:alignment", "content", d.get("current_is_best", False),
                f"最佳假设: {d.get('best_hypothesis')}；"
                f"当前是否最佳: {d.get('current_is_best')}", severity="content")
        except Exception as e:  # noqa: BLE001
            chk("spatial:alignment", "content", False,
                f"**对齐检查文件读不出来**: {e}", severity="content")
    else:
        chk("spatial:alignment", "content", False,
            "**`spatial_alignment_check.json` 不存在 —— 对齐检查没有执行过**。"
            "00_fetch 会在当前假设不是最佳时 raise，所以文件缺失说明这一步"
            "被跳过或产物被删了；坐标方向/镜像从未被验证过。",
            severity="content")

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
            # **S3：沿空间拟时序的基因分析也要有人看。**
            # 那段整块包在 `try/except` 里只 `log_warn`，而
            # `spatial_trajectory_genes.csv` 在本文件里没有任何消费者 ——
            # 于是"一个基因都没算出来"和"算出来 300 个"在验收层一样。
            gal = d.get("genes_along_pseudotime") or {}
            if gal:
                gal_ok = gal.get("status") == "ok" and (gal.get("n_genes_reported") or 0) > 0
                chk("content:spatial_traj_genes", "content", gal_ok,
                    (f"沿空间拟时序的基因 {gal.get('n_genes_reported')} 个"
                     if gal_ok else
                     f"**沿空间拟时序的基因分析失败**：status={gal.get('status')}，"
                     f"{gal.get('n_genes_reported')} 个基因，"
                     f"原因 {str(gal.get('reason'))[:140]}"),
                    severity="content")
            else:
                chk("content:spatial_traj_genes", "content", False,
                    "**status 里没有 `genes_along_pseudotime` 字段** —— "
                    "无法区分『这段没跑』和『跑了但没结果』",
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
        # ---- §0.2 跨语言转换：**空数组必须被解释** ------------------------
        #
        # 本仓库全程 Python（输入是 10x Cell Ranger 的 h5，Part 2 交接的
        # 参考也是 h5ad），**确实没有 R↔Python 转换**。但 `cross_language: []`
        # 和"这一步忘了做"长得一模一样 —— 这正是本仓库反复踩的坑
        # （规则 16 的旧状态文件、规则 20 的 `record_decision()` 没调用点）。
        #
        # 所以判据不是"必须有记录"，而是"**空的话必须有解释**"：
        # 要么 `cross_language` 非空，要么决策链里有一条说明为什么空不了。
        # 这样将来真的加了跨语言步骤而忘了记，这条会立刻变红。
        #
        # **但"有解释"这一半有个洞，这一轮补上了。** 原来的判据是
        # `_cl > 0 or "cross_language" in _dec_nodes` —— 只要决策链里
        # 有那条"本轮没有交接"的说明就通过。问题是：**配置成
        # `deconvolution.reference: h5ad` 时交接真的发生了**
        # （Part 2 的带标签 h5ad 是跨部分来的），而那条"没有交接"的
        # 说明如果还在，检查照样绿 —— **解释成了免检牌**。
        #
        # 所以先看配置说了什么：**配置要求了跨部分参考，就必须有真实记录**，
        # 不接受解释。只有没要求时才允许用解释说明"空是设计如此"。
        #
        # **读全量清单而不是 `msum`** —— `manifest_summary()` 只给
        # `n_decisions` 这个计数，看不到节点名（和 scrna 那边同一个写法）。
        #
        # 用 `read_manifest(cfg)` 而不是 `read_json(res_dir /
        # "run_manifest.json")`：两者等价，但**文件名只该有一处**。
        # `MANIFEST_NAME` 改了而这里写死字符串的话，这条检查会静默地
        # 读一个不存在的文件、`_dec_nodes` 变空 —— 然后报"没解释"，
        # 而真实原因是路径写错了。
        _cl = msum.get("n_cross_language", 0)
        _full = read_manifest(cfg)
        _dec_nodes = {d.get("node") for d in (_full.get("decisions") or [])}
        _dec_ref = str((cfg.get("deconvolution") or {}).get("reference", "builtin"))
        # 配置要求了 Part 2 的参考 → 必须有真实记录，解释不算
        _needs_real = _dec_ref == "h5ad"
        if _needs_real:
            _ok = _cl > 0
            _detail = (f"{_cl} 条交接记录（配置 reference=h5ad，**必须有真实记录**）"
                       if _ok else
                       "**配置了 `deconvolution.reference: h5ad` 却一条交接记录"
                       "都没有** —— Part 2 的参考确实跨了部分边界，"
                       "**不接受『本轮没有交接』的解释**")
        else:
            _ok = _cl > 0 or "cross_language" in _dec_nodes
            _detail = (f"{_cl} 条跨部分交接记录"
                       if _cl else
                       ("0 条，**但决策链里已说明本轮没有跨部分交接**"
                        f"（`reference: {_dec_ref}` 用的是本仓库的 marker 签名）"
                        if "cross_language" in _dec_nodes else
                        "**0 条且没有任何解释** —— 读者无法区分"
                        "『本来就没有』和『忘了记』"))
        chk("manifest:cross_language", "honesty", _ok, _detail)

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
