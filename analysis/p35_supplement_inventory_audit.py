# -*- coding: utf-8 -*-
r"""p35 — 补充材料「产物 ↔ 引用」双向清点

================================ 为什么需要 ================================
Nature Communications 要求 Supplementary items **每一项都应在主文适当位置被引用**。
项目已有 `p18` 检查**补充图**的落点，但：

  · `p18` 的范围**只有图**（它读 `Supplementary_Citation_Map_v1.md`，该表只有 S1–S9）；
  · **补充表这一层既无引用地图，也无检查器**。

2026-10-07 实测发现的三个缺口（全部在本脚本建立前处于无人值守状态）：
  ① `Methods:139` 引用 "Supplementary Methods" —— **该产物不存在**（悬空引用）；
  ② `Supplementary_Tables_R1_R3.md` 已渲染 —— **正文从未引用**（孤儿表）；
  ③ `tableS1_tasks.csv` 只有数据件，无排版件，正文亦未引用。

本脚本做**双向**清点，两个方向都会漏，故都不能省：
  · **孤儿**：产物存在但正文从未引用 → 读者/审稿人会问「这张表为什么在这」；
  · **悬空**：正文引用了不存在的产物 → 引文完整性问题，比孤儿更硬。

边界：只读，不改任何文件，不重算任何统计量。
"""
import io, os, re, sys
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT  = os.path.join(ROOT, "分析结果")
SUPF = os.path.join(ROOT, "figures", "supplement")
REPF = os.path.join(ROOT, "figures", "data", "reporting_profiles")
REP  = open(os.path.join(OUT, "report_p35_supplement_inventory.txt"), "w", encoding="utf-8")
ROWS = []

def p(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); REP.write(s + "\n"); REP.flush()

# ---------------------------------------------------------------- 稿件正本
# MANUSCRIPT_MASTER_DECLARATION.md §2：正本分布在两个目录
MS_FILES = [
    os.path.join(ROOT, "Abstract_v1.md"),
    os.path.join(ROOT, "manuscript_renumbered", "Introduction_v1.md"),
    os.path.join(ROOT, "manuscript_renumbered", "Results_v3_clean.md"),
    os.path.join(ROOT, "manuscript_renumbered", "Discussion_v1.md"),
    os.path.join(ROOT, "manuscript_renumbered", "Methods_v1.md"),
    # 主图图注属于稿件本体，其中的 S 图引用算数
    os.path.join(ROOT, "figures", "final", "Figure_Captions.md"),
]
# ★★ 补充图例（Supplementary_Legends.md）**必须排除** ★★
#   它是**补充材料自身**，不是对补充材料的引用。首版把它算进来，于是
#   「某张 S 图从未被正文引用」这一整类缺陷**恒被漏掉** ——
#   因为图例里当然会出现 "Fig. S5"。
#   负向测试（删掉 Results 里唯一的 Fig. S5 引用）当场暴露了这一点：判据仍报 PASS。
MS_TXT = "\n".join(io.open(f, encoding="utf-8").read() for f in MS_FILES if os.path.exists(f))
FLAT = " ".join(MS_TXT.split())

# ★★ 记法展开（2026-10-07）★★
#   正文会用**简写**引用：`Fig. S1,S2`（逗号列举）、`Fig. S1–S3`（区间）。
#   只按单数形式匹配会把被列举/被区间覆盖的条目误判成孤儿。
#   首版即如此；本函数把简写**展开成单数形式**再参与匹配。
#   注：实测本稿当前**没有**用到这些简写（补充图/表全部是单数引用），
#   但判据不该依赖「作者恰好没这么写」。
def expand_citations(text):
    out = [text]
    # 区间：Fig. S1–S3 ／ Figure S1-S3 ／ Table S1–S2 ／ Tables R1–R3
    for m in re.finditer(r"(Fig\.|Figure|Table|Tables)\s*(S|R)(\d)\s*[–-]\s*(?:S|R)?(\d)", text):
        head, kind, a, b = m.group(1), m.group(2), int(m.group(3)), int(m.group(4))
        for i in range(a, b + 1):
            out.append(f"{head} {kind}{i}")
    # 逗号列举：Fig. S1,S2,S3 ／ Table S1,S2
    for m in re.finditer(r"(Fig\.|Figure|Table|Tables)\s*(S|R)(\d)((?:\s*,\s*(?:S|R)?\d)+)", text):
        head, kind, first, rest = m.group(1), m.group(2), int(m.group(3)), m.group(4)
        nums = [first] + [int(x) for x in re.findall(r"\d", rest)]
        for i in nums:
            out.append(f"{head} {kind}{i}")
    return " ".join(out)

FLAT_X = expand_citations(FLAT)

def cited(pat):
    """**孤儿方向**：该条目是否被正文以任何记法引用过 ⇒ 用**展开后**的文本。
       （`Fig. S1,S2` 要能证明 S2 被引用了。）"""
    return bool(re.search(pat, FLAT_X))

# ---------------------------------------------------------------- 产物清点
arts = []
# ★ 两类实体，引用要求不同：
#   · **编号条目**（Fig. Sx / Table Sx / Tables Rx）：NC 要求每一项都在主文被引用 ⇒ req_cite=True
#   · **结构性章节**（Supplementary Figure Legends、Supplementary Methods 之类）：
#     它们是补充材料的固定组成部分，**不要求**主文引用；但反过来，
#     若主文**指向**某个章节，该章节**必须存在**（悬空方向仍要查）。
for n in range(1, 10):
    f = os.path.join(SUPF, f"FigureS{n}.png")
    if os.path.exists(f):
        arts.append(dict(kind="figure", name=f"Fig. S{n}", file=f"figures/supplement/FigureS{n}.png",
                         cite_pat=rf"(?:Supplementary\s+)?(?:Fig\.|Figure)\s*S{n}\b", req_cite=True))
for sec_name, sec_pat, sec_file in [
    ("Supplementary Figure Legends", r"Supplementary Figure Legends",
     "figures/final/Supplementary_Legends.md"),
    ("Supplementary Methods", r"Supplementary Methods",
     "figures/data/reporting_profiles/Supplementary_Methods.md"),
]:
    if os.path.exists(os.path.join(ROOT, sec_file.replace("/", os.sep))):
        arts.append(dict(kind="section", name=sec_name, file=sec_file,
                         cite_pat=sec_pat, req_cite=False))
for fn in sorted(os.listdir(REPF)) if os.path.isdir(REPF) else []:
    if not fn.endswith(".md"):
        continue
    m = re.match(r"Supplementary_Table_S(\d+)\.md$", fn)
    if m:
        arts.append(dict(kind="table", name=f"Supplementary Table S{m.group(1)}",
                         file=f"figures/data/reporting_profiles/{fn}",
                         cite_pat=rf"Supplementary Table S{m.group(1)}\b", req_cite=True))
    m = re.match(r"Supplementary_Tables_R(\d+)_R(\d+)\.md$", fn)
    if m:
        lo, hi = int(m.group(1)), int(m.group(2))
        arts.append(dict(kind="table", name=f"Supplementary Tables R{lo}–R{hi}",
                         file=f"figures/data/reporting_profiles/{fn}",
                         cite_pat=rf"(?:Supplementary\s+)?Tables?\s+R{lo}\b", req_cite=True))

# ---------------------------------------------------------------- 方向一：孤儿
p("=" * 104); p("p35 — 补充材料「产物 ↔ 引用」双向清点"); p("=" * 104)
p(f"正本文件 {len(MS_FILES)} 份（含图注与补充图例）")
p("")
p("【方向一】孤儿：产物存在，正文从未引用")
bad_orphan = []
for a in arts:
    ok = cited(a["cite_pat"]) or not a["req_cite"]
    st = "PASS" if ok else "FAIL"
    note = a["file"] if a["req_cite"] else a["file"] + "（结构性章节，不要求主文引用）"
    ROWS.append(dict(direction="orphan", kind=a["kind"], item=a["name"], status=st, note=note))
    p("  [%-4s] %-34s %s" % (st, a["name"],
                            "" if ok else "**从未在正本出现**"))
    if not ok: bad_orphan.append(a["name"])
p("")

# ---------------------------------------------------------------- 方向二：悬空
p("【方向二】悬空：正文引用了，但产物不存在")
# 枚举正本里所有 "Supplementary X" / "Table Sx" / "Fig. Sx" 形式的引用 token
# ★ 交替**按最长优先**排列 —— Python 的 | 取第一个匹配成功的分支，
#   故带编号的具体形式必须排在裸前缀之前。首版把裸前缀排在前面，
#   于是 "Supplementary Fig. S8" 被截成 "Supplementary Fig."，产生 4 条假阳性。
CITE_RE = re.compile(
    r"Supplementary\s+Figure\s+S\d+"
    r"|Supplementary\s+Fig\.\s*S\d+"
    r"|Supplementary\s+Fig\s+S\d+"
    r"|Supplementary\s+Tables?\s+R\d+(?:\s*[–-]\s*R\d+)?"
    r"|Supplementary\s+Table\s+S\d+"
    r"|Tables?\s+R\d+"
    r"|Table\s+S\d+"
    r"|(?:Fig\.|Figure)\s*S\d+"
    r"|Supplementary\s+Figure\s+Legends"      # 必须排在裸前缀之前
    r"|Supplementary\s+[A-Z][A-Za-z]+"        # 章节级，如 Supplementary Methods
)
# ★ **悬空方向**用**原文**（FLAT），不用展开文本：
#   `Tables R1–R3` 展开成 R1/R2/R3 后会被误判成三个悬空条目。
#   两个方向要的东西不同 —— 孤儿问「有没有被引到」，悬空问「引用的是什么」。
found = {}
for m in CITE_RE.finditer(FLAT):
    t = " ".join(m.group(0).split())
    found[t] = found.get(t, 0) + 1
# 归一到「实体名」，便于与产物比对
def norm_token(t):
    """归一到产物名。「Supplementary Fig. S9」「Figure S9」「Fig. S9」必须归到同一个名字，
       否则同一个实体在引用侧与产物侧对不上，产生假阳性。"""
    t = " ".join(t.split())
    for pat, fmt in [
        (r"Supplementary\s+(?:Fig\.|Figure)\s*(S\d+)$", lambda m: f"Fig. {m.group(1)}"),
        (r"(?:Fig\.|Figure)\s*(S\d+)$",                   lambda m: f"Fig. {m.group(1)}"),
        (r"Supplementary\s+Tables?\s+(R\d+)\s*[–-]\s*(R\d+)$",
         lambda m: f"Supplementary Tables {m.group(1)}–{m.group(2)}"),
        (r"Supplementary\s+Tables?\s+(R\d+)$",
         lambda m: f"Supplementary Tables {m.group(1)}–R3"),
        (r"Tables?\s+(R\d+)$",                            lambda m: f"Supplementary Tables {m.group(1)}–R3"),
        (r"Supplementary\s+Table\s+(S\d+)$",              lambda m: f"Supplementary Table {m.group(1)}"),
        (r"Table\s+(S\d+)$",                              lambda m: f"Supplementary Table {m.group(1)}"),
    ]:
        m = re.match(pat, t)
        if m: return fmt(m)
    return t

known = {a["name"] for a in arts}
dangling = {}
for t, n in found.items():
    nt = norm_token(t)
    if nt in known:
        continue
    dangling.setdefault(nt, 0)
    dangling[nt] += n
# 允许的「章节级」引用需显式登记，不许默默放过
SECTION_LEVEL = {}      # 形如 {"Supplementary Methods": "…理由…"}；空 dict = 一律判悬空
for t, n in sorted(dangling.items(), key=lambda kv: -kv[1]):
    if t in SECTION_LEVEL:
        ROWS.append(dict(direction="dangling", kind="section", item=t, status="PASS",
                         note="已登记：" + SECTION_LEVEL[t]))
        p("  [PASS] %-34s 已登记（%d 次）" % (t, n))
    else:
        ROWS.append(dict(direction="dangling", kind="section", item=t, status="FAIL",
                         note=f"引用 {n} 次，无对应产物"))
        p("  [FAIL] %-34s **引用 %d 次，但无对应产物**" % (t, n))
p("")

# ---------------------------------------------------------------- 汇总
D = pd.DataFrame(ROWS); bad = D[D.status != "PASS"]
p("=" * 104)
p(f"p35 断言合计 {len(D)} 项：PASS {int((D.status=='PASS').sum())}｜**非 PASS {len(bad)}**")
if len(bad):
    p("")
    p(f"    孤儿 {sum(1 for _, r in bad.iterrows() if r['direction']=='orphan')} 项；"
      f"悬空 {sum(1 for _, r in bad.iterrows() if r['direction']=='dangling')} 项")
    for _, r in bad.iterrows():
        p(f"    [{r['direction']}] {r['item']}  {r['note']}")
    p(""); p("*** p35 未通过 ***")
else:
    p(""); p("*** p35 通过 —— 补充材料的产物与引用双向闭合 ***")
D.to_csv(os.path.join(OUT, "p35_supplement_inventory.csv"), index=False, encoding="utf-8-sig")
REP.close()
sys.exit(1 if len(bad) else 0)
