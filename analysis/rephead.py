# -*- coding: utf-8 -*-
r"""rephead — 审计报告的**自我声明头**

================================ 为什么存在 ================================
外部核验（GLM，2026-10-09）：**「报告自己就是工件」** —— 它被当作证据引用，
却**不声明自己的输入文件与生成时刻**。于是「这份报告是哪次跑的、读的是哪一版输入」
只能靠信里的叙述回答，而那正是 INC-021 族反复记录的形状：

    对象不声明自己  （INC-001 副本不声明代次 ／ INC-018 冻结目录不声明代次 ／
                     INC-021 判据不声明覆盖面 ／ **本条：报告不声明输入**）

药是同一种：**自声明字段** ＋ 普查。

================================ 用法 ================================
    import rephead
    for _l in rephead.header(__file__, [输入文件1, 输入文件2], note="触发：…"):
        p(_l)

输出块含：生成时刻(UTC) ｜ 脚本名+sha8 ｜ **上游宇宙指纹**（`genstamp.universe_id()`）
｜ 每个输入文件的 sha8 与**是否存在**（缺失即 `**MISSING**`，肉眼可见）。
"""
import hashlib
import os
import sys
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import genstamp  # noqa: E402


def _sha8(path):
    try:
        with open(path, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()[:8]
    except Exception:
        return "--------"


def _eligible():
    r"""`p6b_corrected_universe.csv` 是**候选格**（228 行）带一个 `eligible` 列（186 True）。
    只写「rows=228」会让人以为宇宙是 228 —— 故两个数都写，并把指纹的**覆盖面**说明白：
    sha 覆盖整个文件，故 eligible 集合一变，指纹必变。"""
    try:
        import pandas as pd
        d = pd.read_csv(genstamp.UNIVERSE)
        return int(d["eligible"].sum()), len(d)
    except Exception:
        return None, None


def header(script, inputs, note=""):
    u = genstamp.universe_id()
    nel, nrow = _eligible()
    L = ["=" * 100,
         "本报告**自我声明**（INC-021 族：对象必须声明自己）",
         "  生成时刻(UTC) : %s" % datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
         "  生成脚本      : %s  sha8=%s" % (os.path.basename(script), _sha8(script)),
         "  上游宇宙指纹  : %s（sha 覆盖**整份** p6b_corrected_universe.csv）" % u.get("universe_sha"),
         "                  候选 %s 行 ｜ **合格 %s** 行（eligible=True）"
         % (nrow, nel),
         "  输入（sha8 ｜ 存在）:"]
    for pth in inputs:
        L.append("    %-52s %s  %s" % (pth, _sha8(pth),
                                       "Y" if os.path.exists(pth) else "**MISSING**"))
    if note:
        L.append("  备注          : %s" % note)
    L.append("=" * 100)
    return L
