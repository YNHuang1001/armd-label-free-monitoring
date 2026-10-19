# -*- coding: utf-8 -*-
r"""冻结的 episode 定义 —— 单一实现来源（single source of truth）

规则（约束集 A-2，2026-09-29 冻结）：

    对每个 patient × organism × normalized specimen group，按日期升序排列；
    第一条记录建立 episode anchor；
    后续记录若 (date − current_anchor) < 14 天 → 并入当前 episode；
    否则建立新 episode，并把该记录日期作为新的 anchor。
    **anchor 在 episode 内不移动。**

明确后果：0 / 20 / 34 天 → **三个** episode（anchor 链）。
         这与"滚动 14 天窗"（rolling）不同：后者会把 0/5/12/19 并成一集。

实现说明（2026-09-29 审计更正）：
    `p2_10` 曾用**不动点迭代**近似该规则，在 0/20/34 这类序列上不收敛
    （在 [1,1,1] 与 [1,0,0] 之间振荡），4 轮上限后停在奇偶位置；
    MGB/EC 因此报 87,588 集而正确值为 86,344（1.8% 的组不同）。
    本模块是该规则的唯一权威实现。

运行自检：  python episode_definition.py
"""
import numpy as np
import pandas as pd

WINDOW_DAYS = 14
GROUP_KEYS = ("anon_id", "organism", "spec")
DATE_COL = "cdate"


def assign_episodes(df, group_keys=GROUP_KEYS, date_col=DATE_COL, window_days=WINDOW_DAYS):
    """返回与 **df 输入行序一致** 的 episode block id（全局递增整数，从 0 开始）。

    返回值对齐输入行序 —— 这是刻意的：早期实现返回的是内部排序后的顺序，
    导致调用方拿到的 block id 与自己的行序错位（静默错配）。
    """
    d = df[list(group_keys)].astype(str).agg("|".join, axis=1).values
    dt = pd.to_datetime(df[date_col], errors="coerce").values
    n = len(df)
    out = np.full(n, -1, dtype=int)
    order = np.lexsort((dt, d))          # 主键 = 分组键，次键 = 日期
    cur = None
    anchor = None
    b = -1
    for i in order:
        if d[i] != cur:                   # 新分组 -> 建立新 episode 与 anchor
            b += 1
            anchor = dt[i]
            cur = d[i]
        elif (dt[i] - anchor) / np.timedelta64(1, "D") >= window_days:
            b += 1                        # 超出窗口 -> 新 episode，anchor 前移
            anchor = dt[i]
        out[i] = b
    return out


def assign_epid(df, **kw):
    """在 assign_episodes 之上拼出字符串 epid = anon_id|spec|block。"""
    blk = assign_episodes(df, **kw)
    return (df["anon_id"].astype(str) + "|" + df["spec"].astype(str) + "|" + blk.astype(str)).values


# ============================== 永久 unit tests ==============================
def _mk(days, pid="P1", org="EC", spec="URINE"):
    return pd.DataFrame({
        "anon_id": pid, "organism": org, "spec": spec,
        "cdate": pd.to_datetime("2020-01-01") + pd.to_timedelta(days, unit="D"),
    })


def _n_ep(b):
    return int(b.max()) + 1


def run_tests():
    fails = []

    def chk(cond, name, got=None, want=None):
        if not cond:
            fails.append(f"{name}: got={got} want={want}")

    # --- 三条指定测试 ---
    b = assign_episodes(_mk([0, 13])); chk(_n_ep(b) == 1, "0,13 -> 1 episode", _n_ep(b), 1)
    b = assign_episodes(_mk([0, 14])); chk(_n_ep(b) == 2, "0,14 -> 2 episodes", _n_ep(b), 2)
    b = assign_episodes(_mk([0, 20, 34])); chk(_n_ep(b) == 3, "0,20,34 -> 3 episodes", _n_ep(b), 3)

    # --- anchor 不移动 vs 滚动窗（防止实现漂移的关键判据）---
    b = assign_episodes(_mk([0, 5, 12, 19]))
    chk(_n_ep(b) == 2, "0,5,12,19 -> 2 episodes (anchor 链，非滚动窗)", _n_ep(b), 2)
    chk(list(b) == [0, 0, 0, 1], "0,5,12,19 逐行 block", list(b), [0, 0, 0, 1])

    # --- 边界 ---
    b = assign_episodes(_mk([0])); chk(_n_ep(b) == 1, "单条 -> 1 episode", _n_ep(b), 1)
    b = assign_episodes(_mk([0, 14, 28, 42]))
    chk(_n_ep(b) == 4, "0,14,28,42 -> 4 episodes", _n_ep(b), 4)

    # --- 分组隔离：不同部位 / 不同菌 / 不同患者互不影响 ---
    df = pd.concat([_mk([0, 1], spec="URINE"), _mk([0, 1], spec="BLOOD")], ignore_index=True)
    b = assign_episodes(df)
    chk(len(set(b)), "跨部位分组数 = 2", len(set(b)), 2)
    df2 = pd.concat([_mk([0], pid="P1"), _mk([0], pid="P2")], ignore_index=True)
    chk(len(set(assign_episodes(df2))), "跨患者分组数 = 2", len(set(assign_episodes(df2))), 2)

    # --- 返回值必须对齐输入行序（乱序输入）---
    df3 = _mk([20, 0, 34])            # 故意乱序
    b3 = assign_episodes(df3)
    chk(list(b3) == [1, 0, 2], "乱序输入返回行序对齐", list(b3), [1, 0, 2])

    if fails:
        print("[FAIL] episode_definition unit tests")
        for f in fails:
            print("   -", f)
        return 1
    print("[OK] episode_definition unit tests 全部通过（10 项）")
    return 0


if __name__ == "__main__":
    raise SystemExit(run_tests())
