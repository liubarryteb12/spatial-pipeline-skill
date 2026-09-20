"""静态验证：把 spatial workflow「汇总产物」里那段 python -c 抠出来，
喂一份合成 JSON 跑一遍，确认关键字段不再被截断。

**不碰流水线，只验证打印逻辑。**
"""
import io
import json
import os
import pathlib
import tempfile

import yaml

WF = ".github/workflows/spatial_analysis.yml"
d = yaml.safe_load(io.open(WF, encoding="utf-8"))
step = [s for s in d["jobs"]["pipeline"]["steps"] if s.get("name") == "汇总产物"][0]
src = step["run"]

marker = 'python -c "'
code = None
pos = 0
while True:
    i = src.find(marker, pos)
    if i < 0:
        break
    i += len(marker)
    j = src.index('" 2>/dev/null', i)
    body = src[i:j]
    if "domain_methods" in body:
        code = body
        break
    pos = j
assert code is not None, "没找到打印 domain_methods 的那段"
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

old = os.getcwd()
os.chdir(tmp)
try:
    exec(compile(code, "<summary>", "exec"), {"__name__": "__main__"})
finally:
    os.chdir(old)

print("\n--- 断言 ---")
# 把输出重新捕获一次做断言
import contextlib
buf = io.StringIO()
os.chdir(tmp)
try:
    with contextlib.redirect_stdout(buf):
        exec(compile(code, "<summary>", "exec"), {"__name__": "__main__"})
finally:
    os.chdir(old)
out = buf.getvalue()
for key in ("seeded_before_train", "torch_seeded", "adjusted_rand_index",
            "ari_observed_range", "length_scale_l", "seed"):
    assert key in out, f"关键字段 {key} 仍然没打出来"
    print(f"OK  {key} 在日志里可见")
assert "x" * 400 not in out, "长散文没有截断"
print("OK  长散文（note/why_kmeans）被截断，没有淹没关键字段")
print("\n全部通过")

