# -*- coding: utf-8 -*-
r"""t22_cs04_naming.py — CS-04（SUDO 保留 + 改名）的判据

依据
----
`FREEZE_pre_correction_20261005/CORRECTION_SPEC_v1.0.md` §1 CS-04
`FREEZE_pre_correction_20261005/CS04_SUDO_acceptance.md`（Phase A 规范）
`MANUSCRIPT_MASTER_DECLARATION.md`（逐文件正本认定）

CS-04 裁决：保留 prespecified 实现，正式改名 `SUDO-inspired single-branch AURCC`；
Methods 必须明确写出两处偏离（single branch / 自定义 aggregation）；删除"忠于原论文"。
验收：① 协议 §5.3.1 表头与 §5.3 引言的内部矛盾消除；② 不再出现 "followed the original
definition"；③ 全文命名统一。

★ 本模块**只读文本与冻结快照，不加载生产代码**（除 CS04-05 读两段块文本）。
   它读的是**已冻结的文档产物**，所以在补丁树上不会转绿 —— 它的生命周期是
   「文本修改前 RED → 文本修改后 GREEN」，与 t19（CS-03）不同，属于 artifact 级。

★ 基线（2026-10-06 实测，用于确认判据有靶子可打）
     protocol  `忠于原论文` = 1（:161）      `依原论文的` = 1（:172）
     manuscript `original definition` = 1（Methods:284）
     manuscript+图注 `SUDO-AURCC` = 8（Methods 2 / Results 2 / Discussion 1 / 图注 2 / 补注 1）
   若某条判据的基线本来就是 0，它就是**空转**的，本模块用 CS04-00 把这件事说出来。

设计上**刻意不做任何"按行首排除"**：实测协议里 `>` 开头的行不只在顶部修订区，
146/157/280/319/349 都在正文里，而 146 正是判据①指名的那一行；
稿件文件也每个都有 `>` 行（Supplementary_Legends 266 行里 91 行）。
按行首排除会在最关键的行上静默放行。因此本模块一律**全文扫描、零容忍**，
修订记录由作者在写入时避免复制被禁字符串（已核：`忠于原论文`/`依原论文的` 各只出现 1 次，
且 §14 的历史行 v1.2 里两者都没有）。
"""
import io
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import evs_load as L

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = L.PROJ
FREEZE = os.path.join(ROOT, "FREEZE_pre_correction_20261005")

PROTO = os.path.join(ROOT, "LabelFree_AMR_Deployment_Benchmark_Protocol_v1.0_FROZEN.md")

# 稿件正本 —— 逐文件，来自 MANUSCRIPT_MASTER_DECLARATION.md §2（不是"某个目录是正本"）
MASTERS = [
    "Title_v1.md",                                              # 根（renumbered 无对应）
    "Abstract_v1.md",                                           # 根（同）
    "manuscript_renumbered/Introduction_v1.md",
    "manuscript_renumbered/Methods_v1.md",
    "manuscript_renumbered/Results_v3_clean.md",
    "manuscript_renumbered/Discussion_v1.md",
]
CAPTIONS = [
    "figures/final/Figure_Captions.md",
    "figures/final/Supplementary_Legends.md",
]

# CS-04 验收①：协议**正文**里的两条假主张。刻意不用裸 `原论文` ——
# 实测 `原论文` 有 3 处，第 3 处在 §14 的 v1.2 历史行（:424），**必须保留**。
PROTO_BANNED = ["忠于原论文", "依原论文的"]

# CS-04 验收②：fidelity 主张。实测基线只有 `original definition` = 1，其余 = 0；
# 把 0 的那几条也放进来是为了防将来漂移，不是为了打靶。
FIDELITY_BANNED = [
    "original definition",
    "original SUDO",
    "the original method",
    "faithful to the original",
    "as published",
    "per the original",
]

# CS-04 验收③：旧名。首次定义之外不得出现。
OLD_NAME = "SUDO-AURCC"

# CS-04 验收④：协议 §5.3.1 必须写出**真实** aggregation 形式（不是"custom"一句话）。
# 依据 p6c:319-324 的实测代码，token 取代码里真正出现的算子。
AGG_TOKENS = ["argsort", "cumsum", "trapezoid"]

# CS-04-05：两处 SUDO 块不得因改名而变动。基线 = 冻结快照（真实文件，非硬编码哈希）。
# ★ 2026-10-06：这里原先硬写了两个块的行号（p6c 291-330 / p10b 116-148）。
#   实测 p10b 那一对**已陈旧** —— evs_load 因 p10a/p10b 行数漂移已重基线为 137-169。
#   同一个块在两个文件里各写一份行号 ⇒ 必然分叉（这正是 2026-10-06 那次事故的成因）。
#   现一律从 evs_load 取；那里有 BLOCK_FINGERPRINTS 做首行指纹校验。
BLOCK_P6C_SUDO = L.BLOCK_SUDO          # = (p6c_01_L3_blinded.py, 291, 330)


def _read(p):
    with io.open(p, encoding="utf-8") as f:
        return f.read()


def _count(txt, pat):
    return len(re.findall(re.escape(pat), txt))


def _targets():
    """返回 [(绝对路径, 标签)]；缺失的文件在 CS04-00 里报出来。"""
    out = [(PROTO, "protocol")]
    for m in MASTERS:
        out.append((os.path.join(ROOT, m.replace("/", os.sep)), "manuscript:" + m))
    for c in CAPTIONS:
        out.append((os.path.join(ROOT, c.replace("/", os.sep)), "caption:" + c))
    return out


# ------------------------------------------------------------------ CS04-00 非空转
def _c00():
    miss = [lab for p, lab in _targets() if not os.path.isfile(p)]
    if miss:
        return False, f"声明的文件缺失 {len(miss)} 个 :: {miss[:4]} —— 判据无靶子，直接 FAIL"
    tg = _targets()
    tot_lines = sum(len(_read(p).split("\n")) for p, _ in tg)
    tot_bytes = sum(os.path.getsize(p) for p, _ in tg)
    ok = len(tg) >= 8 and tot_lines >= 1500
    return ok, (f"扫到 {len(tg)} 个文件 / {tot_lines} 行 / {tot_bytes} 字节"
                f"；非空转判据：文件数>=8 且 行数>=1500 -> {ok}")


# ------------------------------------------------- CS04-01 协议正文不再主张"忠于原论文"
def _c01():
    txt = _read(PROTO)
    hits = []
    for pat in PROTO_BANNED:
        for m in re.finditer(re.escape(pat), txt):
            hits.append((pat, txt[:m.start()].count("\n") + 1))
    # 修订记录必须存在：删除主张而**不留记录**同样是缺陷
    has_v14 = bool(re.search(r"^>+\s*#*\s*v1\.4", txt, re.M))
    ok = (not hits) and has_v14
    return ok, (f"protocol 全文扫 {len(PROTO_BANNED)} 个禁用模式 -> 命中 {len(hits)}"
                + (f" :: {hits[:4]}" if hits else "")
                + f"；v1.4 修订记录存在 = {has_v14}"
                + f"（协议共 {len(txt.splitlines())} 行）")


# ------------------------------------------- CS04-02 稿件不再出现 fidelity 主张
def _c02():
    hits = []
    for p, lab in _targets():
        if lab.startswith("protocol"):
            continue
        txt = _read(p)
        for pat in FIDELITY_BANNED:
            for m in re.finditer(re.escape(pat), txt):
                hits.append(f"{lab}:{txt[:m.start()].count(chr(10)) + 1}:{pat}")
    ok = not hits
    return ok, (f"稿件正本+图注 全文扫 {len(FIDELITY_BANNED)} 个模式 -> 命中 {len(hits)}"
                + (f" :: {hits[:4]}" if hits else ""))


# ------------------------------------------------------ CS04-03 命名统一（旧名归零）
def _c03():
    hits = []
    for p, lab in _targets():
        if lab.startswith("protocol"):
            continue
        txt = _read(p)
        for m in re.finditer(re.escape(OLD_NAME), txt):
            hits.append(f"{lab}:{txt[:m.start()].count(chr(10)) + 1}")
    ok = not hits
    return ok, (f"`{OLD_NAME}` 在稿件正本+图注中命中 {len(hits)}"
                + (f" :: {hits[:6]}" if hits else ""))


# ------------------------------- CS04-04 协议 §5.3.1 写出真实 aggregation 形式
def _c04():
    lines = _read(PROTO).split("\n")
    try:
        i = next(k for k, s in enumerate(lines) if re.match(r"^#+\s*5\.3\.1", s))
    except StopIteration:
        return False, "找不到 §5.3.1 标题 —— 无法定位，判 FAIL"
    j = next((k for k in range(i + 1, len(lines)) if re.match(r"^#+\s*5\.3\.2", lines[k])), len(lines))
    sec = "\n".join(lines[i:j])
    found = [t for t in AGG_TOKENS if t in sec]
    ok = len(found) == len(AGG_TOKENS)
    return ok, (f"§5.3.1 区间 = 行 {i + 1}-{j}（{j - i} 行）；应出现的真实算子 "
                f"{AGG_TOKENS} -> 命中 {found}；缺 {[t for t in AGG_TOKENS if t not in found]}")


# ------------------------------ CS04-05 改名不得改动估计量（冻结快照对照）
def _c05():
    """★ 2026-10-06 语义收窄。原版有两半，(b) 会被合法改动永久误报：

      原 (b)「p10b 整文件 sha256 == 冻结记录」的前提**已被 CS-05/06/07/09 取代** ——
      那批补丁**合法地**改动了 p10b（12 处，重放验证全部命中，
      见 FREEZE.../APPLIED_CORRECTIONS_20261006.txt [5]）。
      整文件哈希钉此后只会一直红，而它想保护的东西（p10b 的 SUDO 与 p6c 是同一个估计量）
      已由**行为**判据覆盖：
        · t16_sudo_p10b.py     SUDO-P10B exact-identity（绿）
        · t07 FORK-05-WC-p10b  同一 fixture 上跨副本输出一致（绿）
      故 (b) **退役**，不是重新基线化 —— 重新基线化等于把「文件没被改过」当成不变量，
      而这个项目里「文件被改过」是常态；用行为判据代替才是对的。

      (a) 保留并改为**内容锚定**：按块首行在冻结备份中定位，不再依赖行号。
      本次 p6c 的行号恰好没变（291-330 两边一致），但依赖这一点本身就是隐患。"""
    probs = []

    # (a) p6c 的 SUDO 块 vs 冻结备份 —— 真实文件对照，非硬编码
    fn, a, b = BLOCK_P6C_SUDO
    base = os.path.join(FREEZE, "PRODUCTION_backup", fn)
    if not os.path.isfile(base):
        return False, f"冻结备份缺 {fn}"

    cur_lines = L.read_lines(os.path.join(L.OUT, fn))
    cur_seg = cur_lines[a - 1:b]
    fz_lines = L.read_lines(base)
    hits = [i for i, s in enumerate(fz_lines) if s == cur_seg[0]]
    if len(hits) != 1:
        return False, (f"{fn}:{a} 的首行在冻结备份中出现 {len(hits)} 次，无法按内容唯一锚定"
                       f"（{cur_seg[0][:60]!r}）—— 判据无靶子，判红")
    i0 = hits[0]
    fz_seg = fz_lines[i0:i0 + len(cur_seg)]
    if len(fz_seg) != len(cur_seg) or fz_seg != cur_seg:
        d = [i0 + k + 1 for k, (x, y) in enumerate(zip(fz_seg, cur_seg)) if x != y]
        return False, (f"{fn} SUDO 块与冻结备份不同（生产 {a}-{b} ↔ 冻结 {i0 + 1}-{i0 + len(cur_seg)}，"
                       f"差异行 {d[:5]}）")
    probs.append(f"{fn} SUDO 块 == 冻结备份**逐字节相同**"
                 f"（生产 {a}-{b} ↔ 冻结 {i0 + 1}-{i0 + len(cur_seg)}，{len(cur_seg)} 行，"
                 f"按块首行内容锚定；改名确为纯文本操作）")
    probs.append("(b) 整文件哈希钉**已退役**：其前提被 CS-05/06/07/09 的合法改动取代，"
                 "p10b SUDO 的同一性现由 t16_sudo_p10b（exact-identity）与 "
                 "t07 FORK-05-WC-p10b（行为对拍）承担，二者均绿")
    return True, "；".join(probs)


def _main():
    L.banner("CS-04 - SUDO 保留 + 改名（协议 / Methods / 全文命名）")
    L.guard("CS04-00", _c00, "CS-04 判据的非空转前提：声明的 8 个文本文件必须存在且非空")
    L.guard("CS04-01", _c01,
            "CS-04 验收①：协议正文不得再主张实现与原论文一致（`忠于原论文` / `依原论文的`），"
            "且必须存在 v1.4 修订记录（删主张而不留记录同样不合格）。"
            "★ 刻意不使用裸 `原论文`：它有 3 处，第 3 处在 §14 的 v1.2 历史行，必须保留")
    L.guard("CS04-02", _c02,
            "CS-04 验收②：稿件正本+图注不得再出现 any fidelity 主张（实测基线只有 "
            "`original definition` = 1 处，在 Methods:284；其余模式基线为 0，纳入以防漂移）")
    L.guard("CS04-03", _c03,
            f"CS-04 验收③：全文命名统一 —— 旧名 `{OLD_NAME}` 在稿件正本+图注中必须归零"
            f"（实测基线 8 处：Methods 2 / Results 2 / Discussion 1 / 图注 2 / 补注 1）")
    L.guard("CS04-04", _c04,
            "CS-04 验收④：协议 §5.3.1 必须写出**真实** aggregation 形式而不是一句 'custom'。"
            "token 取自 p6c:319-324 的实测代码（argsort / cumsum / trapezoid）")
    L.guard("CS04-05", _c05,
            "CS-04 不变式：改名是纯文本操作 —— p6c 的 SUDO 块（291-330）必须与冻结备份"
            "**逐字节相同**，按块首行内容锚定（不依赖行号）。"
            "★ p10b 整文件哈希钉已于 2026-10-06 退役（前提被 CS-05/06/07/09 合法改动取代；"
            "同一性改由 t16_sudo_p10b + t07 FORK-05-WC-p10b 行为覆盖），"
            "退役理由写在 observed 里，不是静默删除")
    return True


if __name__ == "__main__":
    try:
        _main()
    except Exception:
        import traceback
        L.emit("CS04-99", "FAIL", "模块不得抛出未捕获异常",
               traceback.format_exc().splitlines()[-1][:200], "")
