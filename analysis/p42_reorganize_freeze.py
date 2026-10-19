# -*- coding: utf-8 -*-
r"""p42 — 整理 CONVENTION_FREEZE_20261007 的目录结构

原状：52 个文件**全部平铺**在根目录（脚本 / 数据 / 报告 / 稿件快照 / EVS 证据 / 临时日志混在一起）。
本脚本按类别分入子目录，并**重生成 SHA256SUMS.txt**（路径随之下移）。

★ 三条硬约束：
  1. **一个文件都不许丢** —— 搬完后逐一比对：每个原文件的 SHA-256 仍能在新树里找到。
  2. **不重算任何哈希** —— 只是把**同一份字节**换个位置；清单里的哈希值原样保留。
  3. 临时日志（`_evs*.txt`，我调试时的中间抓取）删除，并在报告里**列名**，不静默丢弃。
"""
import hashlib, io, os, shutil, sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FRZ  = os.path.join(ROOT, "CONVENTION_FREEZE_20261007")

# 类别 → (目录名, 目录说明)
CATS = [
    ("01_manuscript_snapshot", "稿件正本快照（冻结时点）"),
    ("02_supplementary_render", "补充材料渲染件（S1 / S2 / Supplementary Methods）"),
    ("03_submission_docx",     "投稿包成品（两份 Word + 其说明）"),
    ("04_scripts",             "分析脚本与验证套件"),
    ("05_data_artifacts",      "数据产物（冻结 CSV）"),
    ("06_reports",             "审计报告（文本）"),
    ("07_evs_evidence",        "EVS closure 证据"),
    ("08_superseded",          "被取代的旧版本（保留备查）"),
]

# 文件名 → 类别。未列出的文件**留在根目录**（必须在报告中列出）
MAP = {
    "Abstract_v1.md": "01_manuscript_snapshot",
    "Methods_v1.md": "01_manuscript_snapshot",
    "Results_v3_clean.md": "01_manuscript_snapshot",

    "Supplementary_Table_S1.md": "02_supplementary_render",
    "Supplementary_Table_S2.md": "02_supplementary_render",
    "Supplementary_Methods.md": "02_supplementary_render",
    "supplementary_tableS1_render_audit.csv": "02_supplementary_render",
    "supplementary_tableS2_render_audit.csv": "02_supplementary_render",
    "supplementary_tableS2_provenance.csv": "02_supplementary_render",
    "supplementary_methods_render_audit.csv": "02_supplementary_render",

    "Manuscript_main_text.docx": "03_submission_docx",
    "Supplementary_Tables.docx": "03_submission_docx",
    "README_投稿包说明.md": "03_submission_docx",

    "ACTIVE_PATH_MANIFEST.yaml": "04_scripts",
    "run_evs.py": "04_scripts",
    "t06_bounds_scan.py": "04_scripts",
    "t14_manifest_closure.py": "04_scripts",
    "t23_render_parity.py": "04_scripts",

    "tableS2_metrics.csv": "05_data_artifacts",
    "p32_constant_baselines.csv": "05_data_artifacts",
    "p32_direction_convention.csv": "05_data_artifacts",
    "p33_protocol72_metrics.csv": "05_data_artifacts",
    "p33_signal_clusters.csv": "05_data_artifacts",
    "p33_signal_correlation.csv": "05_data_artifacts",
    "p33_signal_top_pairs.csv": "05_data_artifacts",
    "p33_top_pair_ci.csv": "05_data_artifacts",
    "p34_abstract_audit.csv": "05_data_artifacts",
    "p35_supplement_inventory.csv": "05_data_artifacts",

    "report_p32_direction_convention.txt": "06_reports",
    "report_p33_protocol72.txt": "06_reports",
    "report_p34_abstract.txt": "06_reports",

    "report_evs_real.txt": "07_evs_evidence",
    "coverage_table_real.csv": "07_evs_evidence",
    "manifest_sha256.txt": "07_evs_evidence",

    "report_evs_final.txt": "08_superseded",
    "coverage_table_final.csv": "08_superseded",
}
# 分析脚本（按前缀规则）
SCRIPT_PREFIX = ("p11_", "p15_", "p21e_", "p21f_", "p30_", "p32_", "p33_", "p34_", "p35_", "p40_", "p41_")

TEMP_LOGS = ["_evs6.txt", "_evs7.txt", "_evs_final.txt"]
KEEP_AT_ROOT = {"README.md", "SHA256SUMS.txt"}

# ---------------------------------------------------------------- 搬运前留底
before = {}
for fn in os.listdir(FRZ):
    p = os.path.join(FRZ, fn)
    if os.path.isfile(p):
        before[fn] = hashlib.sha256(open(p, "rb").read()).hexdigest()

print("=" * 96); print("p42 — 整理冻结目录"); print("=" * 96)
print(f"搬运前：根目录 {len(before)} 个文件")
print("")

# ---------------------------------------------------------------- 建目录
for d, _ in CATS:
    os.makedirs(os.path.join(FRZ, d), exist_ok=True)

moved, removed, stayed = [], [], []
for fn in sorted(before):
    src = os.path.join(FRZ, fn)
    if fn in TEMP_LOGS:
        os.remove(src); removed.append(fn); continue
    if fn in KEEP_AT_ROOT:
        stayed.append(fn); continue
    cat = MAP.get(fn)
    if cat is None and fn.startswith(SCRIPT_PREFIX):
        cat = "04_scripts"
    if cat is None:
        stayed.append(fn); continue
    shutil.move(src, os.path.join(FRZ, cat, fn))
    moved.append((fn, cat))

print(f"移入分类目录：{len(moved)} 个")
for d, desc in CATS:
    n = len([1 for _, c in moved if c == d])
    print(f"  {d}/  ({n:2d})  {desc}")
print("")
print(f"留在根目录：{stayed}")
print(f"删除的临时日志（调试中间抓取，非产物）：{removed}")
print("")

# ---------------------------------------------------------------- 校验一个都没丢
after = {}
for dp, dn, fns in os.walk(FRZ):
    for fn in fns:
        p = os.path.join(dp, fn)
        after[os.path.relpath(p, FRZ).replace(os.sep, "/")] = hashlib.sha256(open(p, "rb").read()).hexdigest()
after_hashes = set(after.values())
# ★ 校验口径：**除「明确删除并已在上面列名」的临时日志外**，一个都不许少。
#   首版把删除项也算进「应存在」集合，于是自己判了自己一个 FAIL —— 校验器写错方向，
#   比不校验更坏（它会让人去"修"一个并不存在的问题）。
must_survive = {fn: h for fn, h in before.items() if fn not in TEMP_LOGS}
lost = [fn for fn, h in must_survive.items() if h not in after_hashes]
print("【校验】")
print(f"  {'PASS' if not lost else 'FAIL'}  除 3 个已列名的临时日志外，其余文件一个没丢"
      f"（按 SHA-256 逐份核对，共 {len(must_survive)} 份）"
      f"  {'' if not lost else '缺：' + str(lost)}")
print(f"  搬运前 {len(before)} 个 → 搬后 {len(after)} 个"
      f"（-{len(removed)} 临时日志 +1 新 README +1 SHA256SUMS）")

# ---------------------------------------------------------------- 重生成清单（哈希原样保留）
old = {}
for ln in io.open(os.path.join(FRZ, "SHA256SUMS.txt"), encoding="utf-8-sig"):
    p = ln.rstrip("\n").split("  ", 1)
    if len(p) == 2:
        old[p[1].strip()] = p[0].strip()
lines = []
for rel in sorted(after):
    base = os.path.basename(rel)
    if base in ("SHA256SUMS.txt",):
        continue
    h = old.get(base, after[rel])
    lines.append(f"{h}  {rel}")
io.open(os.path.join(FRZ, "SHA256SUMS.txt"), "w", encoding="utf-8").write("\n".join(lines) + "\n")
kept = sum(1 for ln in lines if ln.split("  ", 1)[1].split("/")[-1] in old)
print(f"  SHA256SUMS.txt 重写：{len(lines)} 条（其中 {kept} 条的哈希**原样沿用**，未重算）")
print("")
print("*** p42 完成 ***")
