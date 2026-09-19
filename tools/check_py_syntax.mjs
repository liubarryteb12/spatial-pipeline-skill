#!/usr/bin/env node
/**
 * tools/check_py_syntax.mjs — 所有 Python 脚本的语法检查
 *
 * **为什么需要它。** CI 里跑真实流程要装 scanpy 全家桶、下载数据、
 * 跑几分钟。语法错误在第 1 秒就能发现，却要等 5 分钟才暴露。
 * 这一步在依赖安装之前跑，把"手滑打错字"这类问题挡在最前面。
 *
 * 用 `python -m py_compile`（真正的编译，不是正则匹配）。它会写出
 * __pycache__，检查完清掉。
 *
 * 另外做两条本仓库特有的检查（见 AGENTS.md 规则）：
 *   1. 不允许在 scripts/ 里用 print() 直接输出日志 —— 必须走 common 的
 *      log_info/log_warn/log_error（否则 CI 日志没有时间戳和级别）
 *   2. 不允许硬编码 results/ 或 data/ 路径 —— 必须从 cfg 派生
 *
 * 用法: node tools/check_py_syntax.mjs
 */

import { execFileSync } from "node:child_process";
import { readdirSync, readFileSync, rmSync, existsSync, statSync } from "node:fs";
import { join, relative } from "node:path";

const REPO = new URL("..", import.meta.url).pathname.replace(/^\/([A-Za-z]:)/, "$1");
const SCAN_DIRS = ["scripts", "tools"];

function walk(dir, out = []) {
  if (!existsSync(dir)) return out;
  for (const e of readdirSync(dir, { withFileTypes: true })) {
    const p = join(dir, e.name);
    if (e.isDirectory()) {
      if (e.name === "__pycache__" || e.name.startsWith(".")) continue;
      walk(p, out);
    } else if (e.name.endsWith(".py")) {
      out.push(p);
    }
  }
  return out;
}

const files = SCAN_DIRS.flatMap((d) => walk(join(REPO, d)));
if (files.length === 0) {
  console.error("没有找到任何 .py 文件 —— 目录结构不对？");
  process.exit(1);
}

let failed = 0;

// ---- 1. 真正的语法编译 -----------------------------------------------------
console.log(`编译检查 ${files.length} 个 Python 文件`);
try {
  execFileSync("python", ["-m", "py_compile", ...files], { stdio: "pipe" });
} catch (e) {
  const out = `${e.stdout ?? ""}${e.stderr ?? ""}`;
  console.error("语法检查失败:\n" + out);
  failed++;
}

// 清掉 py_compile 产生的 __pycache__
for (const d of SCAN_DIRS) {
  const pc = join(REPO, d, "__pycache__");
  if (existsSync(pc)) rmSync(pc, { recursive: true, force: true });
}

// ---- 2. 本仓库特有的规则 ---------------------------------------------------
const RULES = [
  {
    name: "scripts/ 里不能用裸 print() 输出日志",
    // 允许 print 出现在 __main__ 的极少数地方？不允许 —— 一律走 common。
    test: (line) => /^\s*print\s*\(/.test(line),
    dirs: ["scripts"],
    hint: "改用 common 的 log_info / log_warn / log_error",
    // common.py 自己实现日志，内部用 print
    exempt: ["scripts/lib/common.py"],
  },
  {
    name: "不能硬编码 results/ 或 data/ 路径",
    test: (line) => /["'`](results|data)\//.test(line),
    dirs: ["scripts"],
    hint: "从 cfg['output']['results_dir'|'data_dir'|'figures_dir'] 派生",
    exempt: ["scripts/lib/common.py"],
  },
  {
    // **实测踩过。** 图的标题/轴标签里写中文时，matplotlib 用 DejaVu Sans，
    // 它没有 CJK 字形 —— 图上显示成一个个方框（豆腐块），
    // 而代码不报错、CI 是绿的、`check_figures.mjs` 也只看到"有墨迹"。
    // 只有打开图才发现标题不可读。
    //
    // 中文解释一律放注释和 JSON 产物里（那些地方中文完全没问题）；
    // 图上只用英文。
    name: "图标签里不能用中文（字体无 CJK 字形，会显示成豆腐块）",
    test: (line) =>
      /(set_title|suptitle|set_xlabel|set_ylabel|\.text|label\s*=|set_xticklabels|set_yticklabels)\s*\(/.test(line) &&
      /[\u4e00-\u9fff]/.test(line),
    dirs: ["scripts"],
    hint: "图上标签改用英文；中文说明放注释或 JSON 产物里",
    exempt: [],
  },
  {
    // `obsm['spatial']` 的列序是 (x, y) = (pxl_col_in_fullres, pxl_row_in_fullres)。
    // 写反了散点图看起来"也像组织形状"，只有叠到 H&E 上才发现是转置的。
    // 这条规则要求：凡是从 obsm 里取 spatial 的地方，附近必须出现 x/y 的
    // 明确注释 —— 无法静态判断对错，只能强制留痕。
    name: "取 obsm['spatial'] 必须注明列序 (x=col, y=row)",
    test: (line) => /obsm\s*\[\s*["']spatial["']\s*\]/.test(line),
    dirs: ["scripts"],
    hint: "同一行或注释里写明列序：(x, y) = (pxl_col_in_fullres, pxl_row_in_fullres)；" +
          "写反了只在叠 H&E 时才看得出来",
    exempt: [],
    // 只要文件里出现过说明就算通过（按文件级检查，不是行级）
    fileLevel: /pxl_col_in_fullres|array_col.*x|列序/,
  },
];

for (const f of files) {
  const rel = relative(REPO, f).replace(/\\/g, "/");
  if (!RULES.some((r) => r.dirs.some((d) => rel.startsWith(d + "/")))) continue;
  const src = readFileSync(f, "utf8");
  const lines = src.split(/\r?\n/);
  for (const rule of RULES) {
    if (!rule.dirs.some((d) => rel.startsWith(d + "/"))) continue;
    if (rule.exempt.includes(rel)) continue;
    // 文件级豁免：例如"必须注明列序"这类规则，只要文件里某处写明了就算通过
    if (rule.fileLevel && rule.fileLevel.test(src)) continue;
    lines.forEach((line, i) => {
      if (rule.test(line)) {
        console.error(`  ${rel}:${i + 1}  ${rule.name}`);
        console.error(`      ${line.trim()}`);
        console.error(`      -> ${rule.hint}`);
        failed++;
      }
    });
  }
}

if (failed > 0) {
  console.error(`\n${failed} 项检查失败`);
  process.exit(1);
}
console.log("语法与规则检查全部通过");
