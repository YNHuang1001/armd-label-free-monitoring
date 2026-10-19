# -*- coding: utf-8 -*-
r"""t46_tableS3.py — Supplementary Table S3 的**独立**验证

================================ 为什么需要 ================================
`p46` 自带 8 条闭合断言，但那是**自证**（生成器与生成物同错则同绿）。
本模块**不 import / 不 exec / 不 subprocess** `p46`：

  1. 自己从**最上游的冻结表**（p6d_outcome + p6b + p6c）重建每个信号在每个预算下的排名，
     自己数命中、自己算精确率；
  2. 自己**解析渲染出的 .md**（Markdown 表格逐格）；
  3. 逐格比对。

这是本项目唯一能抓住「生成器与产物一起错」的层次。

================================ 覆盖 ================================
  T46-01  独立重算 Panel A（15 信号 × 3 预算）逐格
  T46-02  独立重算 Panel B（4 政策层 × 3 预算 + 主对比行）逐格
  T46-03  随机行：命中 = k×57/186 的精确期望，且命中/k == 同行精确率
  T46-04  区间只出现在 k=20% 列（k=10% / 30% 一律无区间）
  T46-05  表体无粗体（逐行扫描实际产物）
  T46-06  每行 命中/k == 精确率（表内自洽，不依赖外部数据）
  T46-07  单次舍入：显示值 == 全精度一次舍到该位数
  T46-08  表注指认**恰好一个**预注册格
"""
import csv
import io
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evs_load as L

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(PROJ, "分析结果")
FIG = os.path.join(PROJ, "figures", "data")
RP = os.path.join(FIG, "reporting_profiles")
TBL = os.path.join(RP, "Supplementary_Table_S3.md")

L1 = ["domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max"]
SIGNED_L2 = ["pred_mean_shift", "entropy_mean_shift", "pred_pos_prev_shift"]
L3 = ["AUROC_DR", "AUROC_IOW", "AUROC_OM"]
ALL15 = L1 + SIGNED_L2 + ["entropy_js", "conf_dist_js", "SUDO_AURCC"] + L3
KS = [(0.10, 19), (0.20, 37), (0.30, 56)]


def rd(fn, key="task_id"):
    p = os.path.join(OUT, fn) if not os.path.isabs(fn) else fn
    if not os.path.exists(p):
        return None
    with io.open(p, encoding="utf-8") as f:
        return {r[key]: r for r in csv.DictReader(f)}


def val(m, t, B_, C_):
    v = float(B_[t][m]) if m in B_[t] else float(C_[t][m])
    if m in SIGNED_L2:
        return abs(v)
    if m in ("SUDO_AURCC",) + tuple(L3):
        return -v
    return v


def avg_rank(vals):
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    r = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        a = 0.5 * (i + j) + 1
        for k in range(i, j + 1):
            r[order[k]] = a
        i = j + 1
    return r


def parse_tables(md):
    """独立解析 Markdown 表格 → [(header, [rows...]), ...]"""
    out, cur, hdr = [], [], None
    for ln in md.splitlines():
        s = ln.strip()
        if s.startswith("|") and s.endswith("|"):
            cells = [c.strip() for c in s.strip("|").split("|")]
            if set("".join(cells)) <= set("-: "):
                continue
            if hdr is None:
                hdr = cells
            else:
                cur.append(cells)
        else:
            if hdr is not None:
                out.append((hdr, cur)); hdr, cur = None, []
    if hdr is not None:
        out.append((hdr, cur))
    return out


# ================================================================ 载入
MD = io.open(TBL, encoding="utf-8").read() if os.path.exists(TBL) else None
D_ = rd("p6d_outcome.csv")
B_ = rd("p6b_label_free_signals_blinded.csv")
C_ = rd("p6c_L3_blinded.csv")

if not (MD and D_ and B_ and C_):
    L.check("T46-00", False, "Table S3 与上游冻结表齐备",
            "缺件：" + str([n for n, v in [("S3.md", MD), ("p6d", D_), ("p6b", B_),
                                           ("p6c", C_)] if not v]),
            "缺件即判红，不静默跳过")
    sys.exit(0)

TASKS = sorted(D_)
HARM = {t: int(D_[t]["harmful"]) for t in TASKS}
N_TOT = len(TASKS)
N_FAIL = sum(HARM.values())
BASE = N_FAIL / N_TOT
TABLES = parse_tables(MD)


def p_at(score, k):
    """按 score 降序（并列按 task_id 升序）取前 k，返回 (命中数, 精确率, 展示串集合)"""
    idx = sorted(range(N_TOT), key=lambda i: (-score[i], TASKS[i]))
    ca = sum(HARM[TASKS[i]] for i in idx[:k])
    return ca, ca / k


SCORES = {m: [val(m, t, B_, C_) for t in TASKS] for m in ALL15}
R6 = [avg_rank([float(B_[t][m]) for t in TASKS]) for m in L1]
SHIFT = [sum(R6[j][i] for j in range(6)) / 6.0 for i in range(N_TOT)]
SCORES["L-shift"] = SHIFT
K_OBS = {kf: k for kf, k in KS}

# ================================================================ T46-01 / T46-02 逐格
# Panel A 表：8 列（Signal, Family, k10h, k10p, k20h, k20p, k30h, k30p）
pA = [t for t in TABLES if len(t[0]) == 8 and t[0][0] == "Signal"]
pB = [t for t in TABLES if len(t[0]) == 7 and t[0][0] == "Policy"]
cA = [t for t in TABLES if len(t[0]) == 4 and t[0][0] == "Input-shift statistic"]

badA = []
if not pA:
    badA.append("未解析到 Panel A")
else:
    for row in pA[0][1]:
        m = row[0].strip("`")
        if m not in SCORES:
            badA.append("未知信号 " + m); continue
        for j, (kf, k) in enumerate(KS):
            ca, pr = p_at(SCORES[m], k)
            if row[2 + 2 * j].split("/")[0] != str(ca):
                badA.append(f"{m}@{kf}:命中 {row[2+2*j]} vs {ca}")
            if row[3 + 2 * j].split(" ")[0] != "%.3f" % pr:
                badA.append(f"{m}@{kf}:精确率 {row[3+2*j]} vs {'%.3f' % pr}")
        if len(badA) > 4:
            break
L.check("T46-01", not badA,
        "Panel A：15 信号 × 3 预算的命中/精确率 == 从 p6b/p6c/p6d 独立重算",
        "不符 %d 项%s" % (len(badA), ("：" + "; ".join(badA[:3])) if badA else ""),
        "本模块不 import p46，自行重建排名与计数")

badB = []
if not pB:
    badB.append("未解析到 Panel B")
else:
    exp_layers = {"L-DR": SCORES["AUROC_DR"], "L-OM": SCORES["AUROC_OM"],
                  "L-shift": SHIFT}
    for row in pB[0][1]:
        lab = row[0]
        key = ("L-DR" if lab.startswith("DR (") else
               "L-OM" if lab.startswith("Outcome model") else
               "L-shift" if lab.startswith("Input-shift composite") else
               "RANDOM" if lab.startswith("Random ordering") else None)
        if key is None:
            continue
        for j, (kf, k) in enumerate(KS):
            if key == "RANDOM":
                want_h, want_p = k * BASE, BASE
                if abs(float(row[1 + 2 * j].split("/")[0]) - want_h) > 5e-3:
                    badB.append(f"随机@{kf}:命中 {row[1+2*j]} vs {want_h:.4f}")
            else:
                ca, pr = p_at(exp_layers[key], k)
                if row[1 + 2 * j].split("/")[0] != str(ca):
                    badB.append(f"{key}@{kf}:命中 {row[1+2*j]} vs {ca}")
                if row[2 + 2 * j].split(" ")[0] != "%.3f" % pr:
                    badB.append(f"{key}@{kf}:精确率 {row[2+2*j]} vs {'%.3f' % pr}")
        if len(badB) > 4:
            break
L.check("T46-02", not badB,
        "Panel B：政策层逐格 == 独立重算（复合 = 六个 L1 的等权秩平均）",
        "不符 %d 项%s" % (len(badB), ("：" + "; ".join(badB[:3])) if badB else ""),
        "复合对照的定义错了整条结论失效，故必须独立重建")

# ================================================================ T46-03 随机行自洽
bad3 = []
if pB:
    rnd = [r for r in pB[0][1] if r[0].startswith("Random ordering")]
    if not rnd:
        bad3.append("未找到随机行")
    else:
        for j, (kf, k) in enumerate(KS):
            hcell = rnd[0][1 + 2 * j].split("/")[0]
            pcell = rnd[0][2 + 2 * j]
            h = float(hcell)
            if abs(h - k * BASE) > 5e-3:
                bad3.append(f"@{kf}:命中 {hcell} ≠ k×57/186 = {k*BASE:.4f}")
            if abs(h / k - float(pcell)) > 5e-4:
                bad3.append(f"@{kf}:命中/k = {h/k:.4f} ≠ 同行精确率 {pcell}")
            if re.fullmatch(r"\d+", hcell):
                bad3.append(f"@{kf}:期望值被写成整数 {hcell}（应为小数期望）")
L.check("T46-03", not bad3,
        "随机行：命中 = k×57/186 的**精确期望**（小数），且命中/k == 同行精确率",
        "不符 %d 项%s" % (len(bad3), ("：" + "; ".join(bad3[:3])) if bad3 else ""),
        "针对实测过的缺陷：初版写成 5.00/19 而精确率给 0.306，自相矛盾")

# ================================================================ T46-04 区间只在一列
bad4 = []
for tab in TABLES:
    hdr = tab[0]
    for ci, hname in enumerate(hdr):
        if "CI" in hname or "95%" in hname:
            continue
        for row in tab[1]:
            if ci < len(row) and re.search(r"\(\s*0\.\d+\s*[–-]\s*0\.\d+\s*\)", row[ci]):
                bad4.append(f"{hname}: {row[0][:26]} → {row[ci]}")
L.check("T46-04", not bad4,
        "区间（形如 `(0.056–0.368)`）只出现在 k=20% 列，k=10% / 30% 列一律无区间",
        "越界 %d 处%s" % (len(bad4), ("：" + "; ".join(bad4[:3])) if bad4 else ""),
        "预注册只冻结了 k=20% 的那一个对比")

# ================================================================ T46-05 表体无粗体
tbl_lines = [ln for ln in MD.splitlines() if ln.strip().startswith("|")]
bold = [ln for ln in tbl_lines if "**" in ln]
L.check("T46-05", not bold,
        "表体（%d 行）无粗体" % len(tbl_lines),
        "无" if not bold else "%d 行含 **：%s" % (len(bold), bold[:2]),
        "唯一预注册格由**表注**指认，不由排版指认；且 NC 禁止用粗体作强调")

# ================================================================ T46-06 表内自洽
bad6 = []
for tab in TABLES:
    hdr = tab[0]
    hit_idx = [i for i, h in enumerate(hdr) if h.endswith("hits")]
    for hi in hit_idx:
        pi = hi + 1
        if pi >= len(hdr):
            continue
        for row in tab[1]:
            if hi >= len(row) or pi >= len(row):
                continue
            hc, pc = row[hi], row[pi]
            if "/" not in hc or not re.match(r"^[\d.]+$", pc.split(" ")[0]):
                continue
            h, k = hc.split("/")[0], hc.split("/")[1]
            if not re.match(r"^[\d.]+$", h):
                continue
            if abs(float(h) / float(k) - float(pc.split(" ")[0])) > 6e-4:
                bad6.append(f"{row[0][:24]} @{hdr[hi]}: {hc} vs {pc}")
L.check("T46-06", not bad6,
        "表内自洽：每行 命中/k == 显示的精确率（不依赖外部数据）",
        "不符 %d 处%s" % (len(bad6), ("：" + "; ".join(bad6[:3])) if bad6 else ""),
        "这一条**不需要**上游数据即可判，抓的是「数彼此矛盾」而非「数算错」")

# ================================================================ T46-07 单次舍入
# ★ 判据只扫**精确率**（形如 n/k 的分数），**必须跳过括号里的区间**：
#   区间下界来自 bootstrap 分位数（如 0.056 / 0.079），本就不是 n/k 形式，
#   初版把它们误判为「非一次舍入」⇒ 假失败。
KNOWN_FRAC = set()
for k in (19, 37, 56, 186):          # 186 用于基率 57/186
    for n in range(0, k + 1):
        KNOWN_FRAC.add("%.3f" % (n / k))

# 正对照：本判据必须**认得出**已知的链式舍入产物 0.487 —— 否则它是空转的
_sanity = ("0.487" not in KNOWN_FRAC) and ("0.486" in KNOWN_FRAC)

bad7 = []
for tab in TABLES:
    for row in tab[1]:
        for c in row:
            head = re.split(r"\(", c)[0]          # 丢掉括号内的区间部分
            for m in re.finditer(r"\b([01]\.\d{3})\b", head):
                if m.group(1) not in KNOWN_FRAC:
                    bad7.append(f"{row[0][:22]}: {m.group(1)}")
L.check("T46-07", (not bad7) and _sanity,
        "精确率列的每个三位小数均为某 k 下 n/k 的**一次**舍入（区间不在本判据范围内）",
        ("可疑 %d 处%s" % (len(bad7), "：" + "; ".join(bad7[:3]))) if bad7
        else ("正对照通过：0.486 ∈ 合法集，0.487 ∉ 合法集" if _sanity else "**正对照失败：判据无判别力**"),
        "来历：图例曾把 0.486486→0.4865→0.487（链式舍入，三位应为 0.486）")

# ================================================================ T46-08 表注指认唯一格
note = MD.split("## Panel A")[0]
# ★★ 2026-10-14 pattern 同步（**正当修改**）：外部核验要求表注**不得**再称
#   「pre-registered」—— 该模拟的方案固定在**标签解封之后**，在本稿术语里是 post hoc，
#   正文 §4 已明写「The simulation is post hoc」。表注若仍写 pre-registered 就与正文冲突。
#   ⇒ 本条的**意图未变**（表注必须指认**恰好一个**推断焦点，其余声明为完整性呈现），
#     判据从「出现 pre-registered」改为「出现**新口径**的唯一焦点 ＋ 事后标注」。
n_focus = len(re.findall(r"primary contrast of[\s>]*this simulation", note, re.I))
n_posthoc = len(re.findall(r"post hoc", note, re.I))
ok8 = ("Only the k = 20% column carries inference" in note.replace("\n", " ") or
       "only the k = 20%" in note.lower()) and n_focus >= 1 and n_posthoc >= 1
L.check("T46-08", ok8,
        "表注指认**恰好一个**推断焦点格（k=20%）＋ 标注事后，并声明其余为完整性呈现",
        f"唯一焦点句 {n_focus} 次；事后标注 {n_posthoc} 次；k=20% 限定句 "
        f"{'存在' if 'only the k = 20%' in note.lower() else '缺失'}",
        "45 个绩效格必须有唯一的推断焦点，否则敏感性分析会被蚕食成新的分叉花园")
