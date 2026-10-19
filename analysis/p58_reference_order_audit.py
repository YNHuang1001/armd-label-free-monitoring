# -*- coding: utf-8 -*-
r"""p58 — 参考文献**首现序**审计（References by order of first mention）

================================ 为什么需要 ================================
**INC-035（2026-10-12）**：`References_FINAL.md` 的表头写死了规则 ——

> 编号依据 Nature Communications 正文顺序 Introduction → Results → Discussion → Methods
> 的首次出现位置生成。

而实测**违反**该规则：`Data sources` 段被插到 Methods §1 开头后，
`[27][28][29][30]`（四个数据集/平台文献）出现在**第 13–16 行**，
而 `[21]` 在**第 45 行**、`[22]` 在 150 行、`[23]` 250 行、`[24][25]` 260 行、`[26]` 341 行。

⇒ 首现序实测为 `…20, 27, 28, 29, 30, 21, 22, …, 26`，**不是递增**。

**★ 为什么此前无人发现**：全库有 `p15/p16/p19/p22/p24/p35` 等一批引用相关审计，
它们查的是「**引用号是否存在、是否有对应条目、DOI 是否核实**」——
**没有一条查「出现的先后」**。而这条规则**明明白白写在被引文件的表头里**。
这就是 INC-021 族的老形状：**规则声明了，判据没落地。**

================================ 判据 ================================
  R-1  稿件正文（Intro→Results→Discussion→Methods）中每个 `[n]` 的**首现位置**递增
  R-2  `References_FINAL.md` 的条目数 == 正文引用的最大编号，且 1..N **无缺号**
  R-3  每个编号**至少被引用一次**（列了却不引 = 幽灵条目）
  R-4  `References_FINAL.md` 表头的规则句**仍在**（规则被删则本审计无意义，须判红）

**不查**（显式声明，避免"看着像全查过"）：DOI 有效性（`p22` 管）、
条目与引用号的一一对应（`p35` 管）、补充材料内的引用（不在 Nature 编号体系内）。

================================ 口径 ================================
像 `p19` 一样**先剥掉引用块**（`^> `）与参考区，否则文件头的修订记录会被当成正文
—— 这正是 **INC-034** 的教训。

`--selftest` 用三份手造样本验证判据本身：递增样本应绿、乱序样本应红、缺号样本应红。
"""
import io
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MR = os.path.join(ROOT, "manuscript_renumbered")
REFS = os.path.join(MR, "References_FINAL.md")

# 交付稿的**装订顺序**（p40 就是按这个顺序生成 docx 的）
ORDER = ["Introduction_v1.md", "Results_v3_clean.md", "Discussion_v1.md", "Methods_v1.md"]

REFS_MARKERS = ("## References cited", "## References used", "## References")
CITE_RE = re.compile(r"\[(\d+(?:\s*[,\u2013\u2014-]\s*\d+)*)\]")


def body_of(path):
    """正文 = 全文 － 文件头引用块（`^> `）－ 参考区。"""
    txt = io.open(path, encoding="utf-8").read()
    for mk in REFS_MARKERS:
        if mk in txt:
            txt = txt.split(mk)[0]
    return re.sub(r"^> .*$", "", txt, flags=re.M)


def expand(inner):
    """`21` / `4,5,12` / `17,18` / `27-29` / `27–29` -> [21] / [4,5,12] / ..."""
    out = []
    for part in re.split(r"[,\u2013\u2014]", inner):
        part = part.strip()
        if not part:
            continue
        m = re.fullmatch(r"(\d+)\s*-\s*(\d+)", part)
        if m:
            a, b = int(m.group(1)), int(m.group(2))
            out.extend(range(min(a, b), max(a, b) + 1))
        elif part.isdigit():
            out.append(int(part))
    return out


def first_appearance(files):
    """-> (seq, raw_counts)；seq = 按首现位置排序的编号列表。"""
    first, raw = {}, {}
    for f in files:
        p = os.path.join(MR, f)
        if not os.path.exists(p):
            raise IOError("缺文件：%s" % p)
        txt = body_of(p)
        for m in CITE_RE.finditer(txt):
            for n in expand(m.group(1)):
                raw[n] = raw.get(n, 0) + 1
                first.setdefault(n, (files.index(f), m.start()))
    seq = [n for n, _ in sorted(first.items(), key=lambda kv: kv[1])]
    return seq, raw


def ref_entry_numbers():
    """References_FINAL.md 的条目号（第一列）。"""
    txt = io.open(REFS, encoding="utf-8").read()
    return [int(m.group(1)) for m in re.finditer(r"^\|\s*(\d+)\s*\|", txt, flags=re.M)]


def audit(files, ref_nums, header_txt):
    """纯函数：给三份输入，回一组 (id, ok, expected, observed, note)。"""
    res = []
    try:
        seq, raw = first_appearance(files)
    except IOError as e:
        return [("R-1", False, "四份稿件文件齐备", str(e), "缺件即判红，不静默跳过")]

    # R-1 首现序递增
    bad = [(seq[i], seq[i + 1]) for i in range(len(seq) - 1) if seq[i + 1] < seq[i]]
    res.append(("R-1", not bad,
                "正文首现序严格递增（Intro → Results → Discussion → Methods）",
                "首现序 = %s%s" % (seq, "" if not bad else "；逆序对 = %s" % bad[:4]),
                "规则写在 References_FINAL.md 表头；编号必须按首次提及的先后"))

    # R-2 条目数 == 最大编号，且无缺号
    n_max = max(raw) if raw else 0
    missing = [n for n in range(1, n_max + 1) if n not in ref_nums]
    extra = [n for n in ref_nums if n > n_max]
    res.append(("R-2", (not missing) and (not extra) and (len(ref_nums) == n_max),
                "References 条目号 == 1..%d，无缺号、无超号" % n_max,
                "条目 %d 条；缺号 %s；超出正文最大引用号的条目 %s"
                % (len(ref_nums), missing or "无", extra or "无"),
                "条目集与引用集必须互相闭合"))

    # R-3 每条都被引用
    ghost = [n for n in ref_nums if n not in raw]
    res.append(("R-3", not ghost,
                "每个条目号在正文里至少被引用一次",
                "列了却未引 = %s" % (ghost or "无"),
                "幽灵条目：列了不引，读者无法定位其出现处"))

    # R-4 规则句仍在
    keep = "首次出现位置生成" in header_txt
    res.append(("R-4", keep,
                "References_FINAL.md 表头的编号规则句仍在",
                "规则句%s" % ("在场" if keep else "**已被删除**"),
                "规则被删则本审计失去判据来源（INC-021 族：规则声明与判据必须同处）"))
    return res


def selftest():
    """判据自测：三份手造样本 —— 递增应绿；乱序应红；缺号应红。"""
    global MR
    import tempfile
    ok = True
    cases = [
        ("递增样本", ["a.md"], [1, 2, 3],
         {"a.md": "x[1] y[2] z[3]\n"}, True),
        ("乱序样本", ["a.md"], [1, 2, 3],
         {"a.md": "x[3] y[1] z[2]\n"}, False),
        ("缺号样本", ["a.md"], [1, 2],
         {"a.md": "x[1] y[2] z[3]\n"}, False),
    ]
    d = tempfile.mkdtemp(prefix="p58st_")
    keep = MR
    try:
        for name, files, refnums, contents, want_ok in cases:
            MR = d
            for fn, c in contents.items():
                io.open(os.path.join(d, fn), "w", encoding="utf-8").write(c)
            r = audit(files, refnums, "首次出现位置生成")
            got = all(x[1] for x in r)
            flag = "OK " if got == want_ok else "自测不符"
            if got != want_ok:
                ok = False
            print("    [%s] %s：期望 %s，实测 %s（%s）"
                  % (flag, name, "全绿" if want_ok else "应判红",
                     "全绿" if got else "判红",
                     "; ".join("%s=%s" % (x[0], "PASS" if x[1] else "FAIL") for x in r)))
    finally:
        MR = keep
    return ok


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    import rephead

    REP = io.open(os.path.join(ROOT, "分析结果", "report_p58_reference_order.txt"),
                  "w", encoding="utf-8")

    def p(s=""):
        print(s)
        REP.write(s + "\n")

    files = [os.path.join(MR, f) for f in ORDER] + [REFS]
    for line in rephead.header(__file__, files, note="触发：INC-035 参考文献首现序违规"):
        p(line)
    p("=" * 100)
    p("p58 — 参考文献首现序审计")
    p("=" * 100)
    p("装订顺序：%s" % " → ".join(f.replace("_v1", "").replace("_v3_clean", "") for f in ORDER))
    p("")

    p("【自测】判据本身能否区分对错（先测判据，再测数据）")
    st = selftest()
    p("    自测结论：%s" % ("通过 —— 判据能区分递增/乱序/缺号" if st else "**不通过**"))
    p("")

    header_txt = io.open(REFS, encoding="utf-8").read().split("\n")[2] \
        if os.path.exists(REFS) else ""
    # 表头规则句在第二行（`共 30 条。编号依据 …`）
    hdr_all = io.open(REFS, encoding="utf-8").read()[:800] if os.path.exists(REFS) else ""

    ref_nums = ref_entry_numbers() if os.path.exists(REFS) else []
    res = audit(ORDER, ref_nums, hdr_all)

    seq, raw = None, None
    try:
        seq, raw = first_appearance(ORDER)
    except IOError:
        pass
    if seq:
        p("【实测首现序】")
        p("    %s" % seq)
        p("    引用出现次数：%s" % dict(sorted(raw.items())))
        p("")

    p("【判据】")
    nfail = 0
    for cid, ok, exp, obs, note in res:
        if not ok:
            nfail += 1
        p("    %-5s %-4s %s" % (cid, "PASS" if ok else "FAIL", exp))
        p("          实测：%s" % obs)
        p("          依据：%s" % note)
    p("")
    p("【汇总】FAIL %d ｜ PASS %d" % (nfail, len(res) - nfail))
    if nfail:
        p("")
        p("*** 首现序不符 —— 编号必须按首次提及的先后重排 ***")
        p("    重排方法：对**闭集**编号做两段式占位替换（先把 `[n]` 换成占位符，再换回新号），")
        p("    不要直接顺序替换 —— `[21]` 是 `[210]`/`[21,22]` 的子串。")
    else:
        p("")
        p("*** 首现序审计通过 ***")
    p("done -> %s" % os.path.join(ROOT, "分析结果", "report_p58_reference_order.txt"))
    REP.close()
    return 1 if nfail else 0


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(0 if selftest() else 1)
    sys.exit(main())
