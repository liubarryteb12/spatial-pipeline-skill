"""静态验证：把 spatial workflow「汇总产物」里那段 python -c 抠出来，
喂一份合成 JSON 跑一遍，确认关键字段不再被截断。

**不碰流水线，只验证打印逻辑。**

**不能用 `import yaml`。** 这一步在 CI 里跑在 `pip install` **之前**
（静态检查就该在装依赖前），而 runner 的 Python 只有标准库 ——
第一版就是 `import yaml` 直接 ModuleNotFoundError 把 job 弄红了。
所以这里用纯文本扫描定位代码块，不引入任何依赖。
"""
import contextlib
import io
import json
import os
import pathlib
import re
import sys
import tempfile
import textwrap

WF = pathlib.Path(".github/workflows/spatial_analysis.yml")
if not WF.exists():
    print(f"找不到 {WF}", file=sys.stderr)
    raise SystemExit(1)

src = WF.read_text(encoding="utf-8")

# 纯文本扫描：找到**含 domain_methods 的**那段 `python -c "` ... `" 2>/dev/null`
code = None
for m in re.finditer(r'python -c "', src):
    start = m.end()
    end = src.find('" 2>/dev/null', start)
    if end < 0:
        continue
    body = src[start:end]
    if "domain_methods" in body:
        # **必须 dedent。** YAML 的块标量会把整段代码缩进 10 个空格；
        # 之前用 yaml.safe_load 时是它替我们剥掉的，改成纯文本扫描后
        # 就轮到我们自己剥 —— 不剥会 IndentationError。
        code = textwrap.dedent(body)
        break

if code is None:
    print("没找到打印 domain_methods 的那段 python -c", file=sys.stderr)
    raise SystemExit(1)
print("抠出的代码行数:", len(code.splitlines()))

fake = {
    "domain_methods": {
        "builtin_smooth_leiden": {
            "used": True, "role": "primary", "section": "§3.2",
            "status": "ok", "n_domains": 13, "note": "x" * 400,
        },
        "SpaGCN": {
            "attempted": True, "status": "ok", "reason": "", "section": "§3.2",
            "init": "kmeans", "why_kmeans": "y" * 400, "version": "1.2.7",
            "length_scale_l": 1.4142, "n_clusters": 13, "seed": 0,
            "torch_seeded": True,
            "seeded_before_train": ["random", "numpy", "torch"],
            "n_domains": 13, "mean_max_prob": 0.93,
        },
    },
    "method_agreement": {
        "SpaGCN_vs_builtin": {
            "adjusted_rand_index": 0.3708, "normalized_mutual_info": 0.5451,
            "neighbor_same_frac_a": 0.6682, "neighbor_same_frac_b": 0.6727,
        }
    },
    "reproducibility": {
        "stable": ["morans_I"], "unstable": ["ari"],
        "ari_observed_range": [0.3639, 0.4216],
    },
}

tmp = tempfile.mkdtemp(prefix="summary_check_")
res = pathlib.Path(tmp) / "results" / "lymph_node"
res.mkdir(parents=True)
(res / "domain_status.json").write_text(
    json.dumps(fake, ensure_ascii=False), encoding="utf-8")

buf = io.StringIO()
old = os.getcwd()
os.chdir(tmp)
try:
    with contextlib.redirect_stdout(buf):
        exec(compile(code, "<summary>", "exec"), {"__name__": "__main__"})
finally:
    os.chdir(old)
out = buf.getvalue()
print(out)

print("--- 断言 ---")
for key in ("seeded_before_train", "torch_seeded", "adjusted_rand_index",
            "ari_observed_range", "length_scale_l", "seed"):
    if key not in out:
        print(f"FAIL 关键字段 {key} 仍然没打出来", file=sys.stderr)
        raise SystemExit(1)
    print(f"OK  {key} 在日志里可见")
if "x" * 400 in out:
    print("FAIL 长散文没有截断，会淹没关键字段", file=sys.stderr)
    raise SystemExit(1)
print("OK  长散文（note/why_kmeans）被截断，没有淹没关键字段")

# ---- 反向检查：判据必须有区分力 -------------------------------------------
#
# AGENTS 规则 3：**判据必须有区分力，否则它给的是虚假的安心。**
# 上面那组断言只在"新逻辑正确"时通过，但它会不会在"退回旧逻辑"时
# 也照样通过？把旧的 `[:500]` 逻辑喂进去试一遍 —— 必须失败。
OLD = '''
import json, pathlib
for f, keys in (
        ('results/lymph_node/domain_status.json',
         ('domain_methods', 'method_agreement')),):
    p = pathlib.Path(f)
    if not p.exists():
        print(f, '（缺失）'); continue
    d = json.loads(p.read_text(encoding='utf-8'))
    for k in keys:
        v = d.get(k)
        print('  %s -> %s' % (k, json.dumps(v, ensure_ascii=False)[:500]))
'''
buf2 = io.StringIO()
os.chdir(tmp)
try:
    with contextlib.redirect_stdout(buf2):
        exec(compile(OLD, "<old>", "exec"), {"__name__": "__main__"})
finally:
    os.chdir(old)
old_out = buf2.getvalue()
# 注意**不要**把 adjusted_rand_index 放进来：旧的 [:500] 对
# `method_agreement` 这个短字典其实是够长的，那个字段当时**看得到**。
# 真正被藏掉的是 `domain_methods` 里的字段（它长）以及**根本没打印的**
# `reproducibility`。反向检查只针对确实丢失的那几个 ——
# 把没丢的也写进去会得到一个"永远失败"的假检查。
LOST_BEFORE = ("seeded_before_train", "torch_seeded", "ari_observed_range")
still_visible = [k for k in LOST_BEFORE if k in old_out]
if still_visible:
    print(f"FAIL 旧的 [:500] 逻辑下这些字段仍然可见：{still_visible} —— "
          f"断言没有区分力", file=sys.stderr)
    raise SystemExit(1)
print("OK  反向检查：退回旧的 [:500] 逻辑时，当时真正丢失的 3 个字段"
      "（seeded_before_train / torch_seeded / ari_observed_range）"
      "全部不可见 ——")
print("    所以上面的断言有区分力（回归就会红）")

print("\n全部通过")
