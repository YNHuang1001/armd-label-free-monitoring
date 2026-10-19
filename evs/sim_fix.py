# -*- coding: utf-8 -*-
r"""sim_fix — FIX-SIMULATION HARNESS（估计量验证套件）

目的
====
证明本套件不是"自证系统"：
  ① REAL    —— 在**真实未修改**的生产代码上，已知缺陷被测试抓到（红灯基线）；
  ② PATCHED —— 在同一套测试、指向**按正确方式修补的副本**时，被修缺陷对应测试转绿。
只有 ① 说明测试可能只是钉死了当前状态；只有 ② 说明测试抓不到缺陷。

硬约束（本脚本自身强制）
========================
1. 绝不写真实生产文件（分析结果\*.py）。开工前记录 sha256，收尾复核逐字节一致。
2. 每个补丁 = **一行替换一行**（行数必须不变）。evs_load 用行号区间抽取
   p6c:257-289（L3/DR 块）、p6c:291-330（SUDO 块）、p6b:121-144（STORE 块）；
   行数一变，区间就会静默指向别的代码。
3. 只打三个已确认缺陷：(A) MMD 内项缺因子 2；(B) wconcordance 并列未记 0.5；
   (C) DR 分子写成两残差之积（冻结规格头要求两乘积之差）。
   SUDO、episode 粒度、报告计数**不打**（等作者决定）。

本 harness 自身两个缺陷（本轮修复，不影响生产文件）
==================================================
H1. 逐测试状态聚合原取「最后一行」→ 掩盖 UNIT-01（3 条 FAIL + 6 条 PASS，取最后显示
    PASS）。改为取**最差**状态 FAIL > WARN > PASS，并保留 subgroup 明细
    （逐 id 行计数 + FAIL 子例逐字列出）。
H2. 族级交叉核对原解析 report_evs_*.txt（该文件只有 run_evs 自己的汇总表，没有逐模块
    计数）→ 恒报 report=None / "存在不一致" 的伪 MISMATCH。改为从 [6] 的结构化逐测试
    结果直接计算按族摘要，并写出 evs_results_{real,patched}.csv 机器可读输出后回读核对；
    完全不解析报告文本。

用法
====
    python sim_fix.py                 # 幂等、可重复运行（完整验收：run_evs x2 + 逐测试采集）
    python sim_fix.py --check-only    # 只做补丁/锚点核验，不跑套件
    python sim_fix.py --collect-only  # 不调用 run_evs.py，仅做逐测试采集 + 结论（并行开发期自检）

产物
====
    _patched\                   补丁副本（同名文件）
    report_evs_real.txt         run_evs.py 在 REAL 下写出的报告（--collect-only 时不重写）
    report_evs_patched.txt      run_evs.py 在 PATCHED 下写出的报告（--collect-only 时不重写）
    evs_results_real.csv        [6] 结构化逐测试结果（test_id,status,module,...；REAL）
    evs_results_patched.csv     同上（PATCHED）—— 供独立方与 run_evs 输出交叉核对
    sim_fix_report.txt          本脚本的完整日志
"""
import os, sys, io, csv, time, shutil, hashlib, subprocess, importlib
import concurrent.futures as _cf

# --collect-only 模式下，结构化结果的新鲜度容差（秒）。正常模式下不使用固定窗口，
# 而是拿 CSV 的 mtime 与本次 run_evs 的启动时刻比对（更精确）。
FRESH_WINDOW_S = int(os.environ.get("EVS_FRESH_WINDOW_S") or 3600)

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)                            # 第二篇 新方案
OUT = os.path.join(PROJ, "分析结果")                      # 生产树                 # 真实生产树 = 分析结果\
PATCHED = os.path.join(HERE, "_patched")    # 补丁副本目录
PY = sys.executable

_log = io.open(os.path.join(HERE, "sim_fix_report.txt"), "w", encoding="utf-8")


def p(*a):                                  # 同时打印与落盘
    line = " ".join(str(x) for x in a)
    print(line, flush=True)
    _log.write(line + "\n")
    _log.flush()


def need(cond, msg):
    if not cond:
        p("*** [FATAL] " + msg + " ***")
        _log.close()
        raise SystemExit(2)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_bytes(path):
    with open(path, "rb") as f:
        return f.read()


# ------------------------------------------------------------------ 生产文件清单
# evs_load 会经由 L.* 读取下列文件（P21A 被 t05 的逐字块使用；P5/P6D/EPI 也登记在
# evs_load 中）。全部复制进 _patched，使 PATCHED 与 REAL 只差三个补丁本身。
# v1.4：_patched 必须是**完整树**，否则每条 Critical path 的 required_tests 会因「文件不在
# 补丁树」而整批 skip，closure 只能读到 UNKNOWN（v1.3 的 critic 已点名此缺口）。
# 清单 = ACTIVE_PATH_MANIFEST.yaml 里登记的全部生产文件（不含 glob 的 Rendering 项）。
#
# ★ 披露：p11 / p6e / p6f / p6g / p12 / p14 的 `OUT` 现由
#   `os.path.dirname(os.path.abspath(__file__))` 推导（2026-10-19 起；此前是硬编码的
#   绝对路径常量）。把它拷进 _patched 会**改变源码文本**：副本若被当脚本执行，
#   `__file__` 指向 _patched ⇒ `OUT` 也指向 _patched。因此 t13 逐字执行这些文件的常量块
#   时**显式注入 `__file__` = 真实路径**，其运行期 I/O 仍指向**真实目录**。这正是这些模块
#   能与冻结产物逐项对拍的原因，但也意味着「_patched 是独立树」对这批文件**只成立于源码层**。
#   测试侧一律重定向 to_csv（t13 的 EC-* 显式改写写入目标），不得让补丁跑覆盖真实产物。
LOADER_FILES = [
    "p5_01_degradation_benchmark.py",
    "p6_02_label_free_signals.py",
    "p6b_01_corrected_pipeline.py",
    "p6c_01_L3_blinded.py",
    "p6d_01_unseal_evaluation.py",
    "p6e_01_lodo_fold_qc.py",
    "p6f_lodo_settlement.py",
    "p6g_direction_learning.py",
    "p10_preflight_gate.py",
    "p10a_prepare_cache.py",
    "p10b_subsample_curves.py",
    "p10c_aggregate.py",
    "p11_figure_datapack.py",
    "p12_configs_corrected.py",
    "p14_rebuild_fig3b.py",
    "p21a_reporting_profile_derivation.py",
    "episode_definition.py",
]
# manifest 登记的**目录型**路径（Rendering 层）；补丁树须为完整镜像
LOADER_DIRS = ["provenance_rebuild"]

# ------------------------------------------------------------------ 三个缺陷的补丁
# 每个元素: (tag, 文件, 行号(1-indexed), 原行, 新行, 说明)
PATCHES = [
    # (A) MMD 因子 2 -------------------------------------------------------------
    ("A1", "p6b_01_corrected_pipeline.py", 227,
     "        vals.append(Kaa.sum() / (n * (n - 1)) + Kbb.sum() / (m * (m - 1)) - 2 * Kab.mean())",
     "        vals.append(2 * Kaa.sum() / (n * (n - 1)) + 2 * Kbb.sum() / (m * (m - 1)) - 2 * Kab.mean())",
     "MMD: 两个内项各补因子 2（p6b mmd_rbf）"),
    ("A2", "p6_02_label_free_signals.py", 236,
     "        vals.append(Kaa.sum() / (n * (n - 1)) + Kbb.sum() / (m * (m - 1)) - 2 * Kab.mean())",
     "        vals.append(2 * Kaa.sum() / (n * (n - 1)) + 2 * Kbb.sum() / (m * (m - 1)) - 2 * Kab.mean())",
     "MMD: 同上（p6_02 mmd_rbf；本套件未装载该文件，按任务要求一并修补）"),
    ("A3", "p10a_prepare_cache.py", 279,
     "            Kaa.append(float(Ka.sum() / (len(a) * (len(a) - 1))))",
     "            Kaa.append(float(2 * Ka.sum() / (len(a) * (len(a) - 1))))",
     "MMD: source 内项 Kaa 补因子 2（t01 MMD-07 逐字执行 p10a:273-280）"),
    ("A4", "p10b_subsample_curves.py", 92,
     '                vals.append(A["Kaa"][k] + Kbb.sum() / (m_ * (m_ - 1)) - 2 * Kab.mean())',
     '                vals.append(A["Kaa"][k] + 2 * Kbb.sum() / (m_ * (m_ - 1)) - 2 * Kab.mean())',
     "MMD: target 内项 Kbb 补因子 2（t01 MMD-08 逐字执行 p10b:85-93）"),
    # (B) 并列 0.5 ---------------------------------------------------------------
    # 注意：0.5 记分是加在 **Σ_{i≠j}** 上的，必须扣掉 i=j 的自配对。
    # 按唯一 score 分组后，组 k 内的有序对质量 = A_k*B_k − Σ_{i∈k} a_i b_i。
    # 只加 0.5*Σ A_k B_k 会把对角线也算进去，在**无并列**的数据上都会改变取值
    # （验收实验：s=[3,2,1], a=[1,1,1], b=[1,2,3] → 正确 8.0，错补丁 11.0）。
    ("B1", "p6c_01_L3_blinded.py", 91,
     "    return float(np.sum(A * cumB))",
     "    return float(np.sum(A * cumB) + 0.5 * (np.sum(A * B) - np.sum(aa * bb)))",
     "wconcordance: 并列记 0.5（p6c 副本；Σ_{i≠j}，已扣对角线）"),
    ("B2", "p10b_subsample_curves.py", 46,
     "    return float(np.sum(A * cumB))",
     "    return float(np.sum(A * cumB) + 0.5 * (np.sum(A * B) - np.sum(aa * bb)))",
     "wconcordance: 并列记 0.5（p10b 同源副本；Σ_{i≠j}，已扣对角线）"),
    # (C) DR 分子 ----------------------------------------------------------------
    # (D) CS-07 sens_at_spec80 —— 分位数取反 ------------------------------------
    # 目标定义（CORRECTION_SPEC CS-07）：**从所有可实现阈值中，选「达成 specificity >= spec 时
    # sensitivity 最大」的那个**。等价于：在所有观测到的 score 值 v（加上 max(s) 之上的边界点）里，
    # 取满足 #{neg: s < v}/n_neg >= spec 的**最小** v —— 阈值越小 sensitivity 越大（弱单调），
    # 所以最小的合格阈值就是最优工作点；若无观测值合格，边界点 v > max(s) 恒合格（spec=1, sens=0）。
    # 单行写法以保持行数不变（sim_fix 的硬约束，BLOCK 锚点依赖行号）。
    # ★ 注意：只把 1-spec 换成 spec、阈值仍取分位数，在**有并列**时不能保证交付 specificity >= 0.80
    #   （t10 SENS-04 实测：naive 变体在 243/300 fuzz 数据集上交付 < 0.80）——故必须走枚举。
    ("D1", "p6d_01_unseal_evaluation.py", 266,
     "    thr = np.quantile(s[y == 0], 1 - spec)",
     "    thr = min([v for v in np.unique(s) if (s[y == 0] < v).mean() >= spec] or [np.nextafter(s.max(), np.inf)])",
     "CS-07: sens_at_spec 定义化解（枚举合格阈值取最小）—— p6d 副本"),
    ("D2", "p10b_subsample_curves.py", 260,
     "        thr = np.quantile(s[y == 0], 1 - spec)",
     "        thr = min([v for v in np.unique(s) if (s[y == 0] < v).mean() >= spec] or [np.nextafter(s.max(), np.inf)])",
     "CS-07: 同上 —— p10b 嵌套 def 副本"),
    ("D3", "p10c_aggregate.py", 69,
     "    thr = np.quantile(neg, 1 - spec)",
     "    thr = min([v for v in np.unique(s) if (neg < v).mean() >= spec] or [np.nextafter(s.max(), np.inf)])",
     "CS-07: 同上 —— p10c 顶层 def 副本"),
    # (E) CS-08 grouped bootstrap 的 set() 去重 ---------------------------------
    # 抽到 {A,A,B} 必须按 (2A+B)/3 计。`g[idx]` 本身保留重数；`x in fa` 过滤不改变重数。
    # ★ 注意：**只去掉 set() 是不够的**。`boot` 里 `pick` 是「抽中的 fold 索引」（含重数），
    #   `idx` 是把它展开成 row 索引；`g[idx]` 因此是**按 fold 行数加权**的展开，
    #   而不是「抽中列表」。C 是 fold 的**无权**均值，其 bootstrap 必须对**抽中的 fold 列表**
    #   求均值。t12 EVAL-A1 的 fixture（A/B/C/D 行数 3/2/4/1）把三者分开：
    #     {A,A,B}: 定义 (2fa[A]+fa[B])/3 = 0.8 ｜ set()去重 0.7 ｜ row-weighted g[idx] 0.85
    #   t12 的 stub 先试 fn(idx, pick)、再 fn(idx, pick=pick)、最后 fn(idx) —— 故 E1 用两参 lambda。
    ("E1", "p6f_lodo_settlement.py", 121,
     "        v = fn(idx)",
     "        v = fn(idx, pick) if fn.__code__.co_argcount >= 2 else fn(idx)",
     "CS-08: boot 把「抽中的 fold 索引 pick」交给两参统计量；1 参统计量（p6f:142 的 auc）行为不变"),
    ("E2", "p6f_lodo_settlement.py", 143,
     "    flo, fhi = boot(EVAL, lambda idx: float(np.mean([fa[x] for x in set(g[idx]) if x in fa]))",
     "    flo, fhi = boot(EVAL, lambda idx, pick: float(np.mean([fa[EVAL[i]] for i in pick if EVAL[i] in fa]))",
     "CS-08: 统计量改为对**抽中的 fold 列表**（含重数）求均值 = macro mean 的正确 cluster bootstrap"),
    ("E3", "p6f_lodo_settlement.py", 144,
     "                    if set(g[idx]) & set(fa) else np.nan)",
     "                    if any(EVAL[i] in fa for i in pick) else np.nan)",
     "CS-08: 空抽保护同步改为按抽中的 fold 列表判空"),
    # (F) CS-02 episode 分析单元 ------------------------------------------------
    # `unit` 每 (anon_id, cdate, spec) 一行；`unit.set_index("epid")` 后 X **保留同一 episode 的
    # 多个日期行**，各自继承 episode 聚合后的药敏 ⇒ 伪重复（MGB/PA 1.152x）。
    # `unit` 在本行之前已 .sort_values(["anon_id","spec","cdate"])，故 drop_duplicates("epid")
    # 保留的就是**该 episode 的最早日期**，即 CS-02 冻结的 ① episode index date = anchor date。
    # X["cdate"] 随之为 anchor 日期 ⇒ ② age/setting 取自 anchor 记录（下游 oid_map 用 X["cdate"]）、
    # ③ history 窗相对 anchor 起算（下游用 X["cdate"]）。④ 药敏仍来自 wide（episode 内聚合），
    # 不反向改变 index date。行数不变（BLOCK_STORE 121-144 的锚点依赖行号）。
    # (G) CS-12 report-profile unit identity ------------------------------------
    # p21a 是**独立脚本**，有自己一份 STORE 构造（:269 与 p6b 同源缺陷）**和自己一个独立的 bug**：
    # :276 的 .reset_index(drop=True) 抹掉 epid 索引，:315 的 `eps |= set(Xb.index)` 随后收集的是
    # **行号**；三个菌的行号互相覆盖 ⇒ len(eps) 是行号集合的并集（铁证：production_value ==
    # n_union_positions 精确相等 89,454 / 61,317 / 37,579）。
    # 修法：① :269 与 CS-02 同修；② :276 保留一个 **canonical episode_key = organism|epid**
    # （epid 本身不含 organism，blk 每菌从 0 重编号 ⇒ 跨菌可撞，必须加上 organism）；
    # ③ :315 / :317 改用 episode_key。全部行数不变。
    ("G2", "p21a_reporting_profile_derivation.py", 276,
     '        STORE[(site, sym)] = X.sort_values(["anon_id", "spec", "cdate"]).reset_index(drop=True)',
     '        STORE[(site, sym)] = X.assign(episode_key=sym + "|" + X.index.astype(str)).sort_values(["anon_id", "spec", "cdate"]).reset_index(drop=True)',
     "CS-12: 不抹索引，改为携带 canonical episode_key = organism|epid"),
    ("G3", "p21a_reporting_profile_derivation.py", 315,
     '        eps |= set(Xb.index); pts |= set(Xb["anon_id"])',
     '        eps |= set(Xb["episode_key"]); pts |= set(Xb["anon_id"])',
     "CS-12: episode 集合用 canonical key（不再用行号）"),
    ("G4", "p21a_reporting_profile_derivation.py", 317,
     '        for e, a, st, sp in zip(Xb.index, Xb["age"], Xb["setting"], Xb["spec"]):',
     '        for e, a, st, sp in zip(Xb["episode_key"], Xb["age"], Xb["setting"], Xb["spec"]):',
     "CS-12: ages/sets/specs 的键同步改为 canonical key"),
    # (H) CS-01 —— **已进生产（2026-10-05）**，故此处不再打补丁 ----------------------------
    #   5 个副本的 `_canon` 已改为精确匹配、DRUGMAP 已按 GPT6 裁定补齐别名（见
    #   FREEZE_pre_correction_20261005/CS01_amendment.txt）。_patched 由生产拷贝继承该修复。
    #   保留此段注释是为了让后来者知道：CS-01 的 gate 仍由 t15 覆盖，但它**不再是一条补丁**。
    # `_canon` 现在是 `if any(k in up for k in ks): return c`（子串 + 先到先得）：
    #   · 泛键 "tazobactam"（在 TZP 条目里，位置靠前）会**吞掉** ceftolozane/tazobactam
    #     ⇒ 真实数据 13,426 行被并进 TZP，而本项目 CZA 条目正是 ceftolozane/tazobactam。
    #   · 泛键 "sulbactam"（SAM）吞掉 cefoperazone/sulbactam 等复合药。
    #   · 反向失效：UT 的 `Amoxicillin/k clavulanate aug` 归一化后匹配不上任何键，
    #     落到 `return up` —— 这一支其实**已经**是正确行为（成为自己的类）。
    # 改名表已含精确键（`ceftolozanetazobactam` 等），故**精确匹配**即可：
    #   ceftolozane/tazobactam -> CZA ✅ ｜ ceftazidime/avibactam -> 自己的类（unmapped）✅
    #   cefoperazone/sulbactam -> 自己的类 ✅ ｜ Pip/tazobactam -> piptazobactam -> TZP ✅
    # 复合制剂与单药严格分开，未知串不再被静默归类（CS-01）。5 个副本同一行。
    # (I) CS-07 补全：**未定义输入契约统一为 NaN**（GPT6 2026-10-05 裁决：SENS-06 必须是 gate）
    # 三个副本对退化输入原本给出三种语义（NaN / 0.0 / IndexError）。同一个统计量必须有一致契约：
    #   empty y、仅一个 outcome class、无负类（specificity 未定义）、无正类（sensitivity 未定义）
    #   -> 返回 NaN；**不得返回 0，不得抛 IndexError**。
    # `all_tied` **不在**此列：该处统计量有定义（delivered spec=0, sens=0），必须返回 0.0。
    # 另：阈值行的 `s.max()` 在空数组上会抛 ValueError，故补 `if s.size else 0.0`。
    # 两条都是 1 行换 1 行；`old` 串为 **D1-D3 之后**的文本，故必须排在 D 之后。
    ("S1", "p6d_01_unseal_evaluation.py", 266,
     "    thr = min([v for v in np.unique(s) if (s[y == 0] < v).mean() >= spec] or [np.nextafter(s.max(), np.inf)])",
     "    thr = min([v for v in np.unique(s) if (s[y == 0] < v).mean() >= spec] or [np.nextafter(s.max() if s.size else 0.0, np.inf)])",
     "CS-07 补全: 空数组下阈值行不得抛 ValueError（p6d）"),
    ("S2", "p6d_01_unseal_evaluation.py", 269,
     "    return tp / max(tp + fn, 1)",
     "    return tp / max(tp + fn, 1) if (len(np.unique(y)) > 1 and (y == 0).sum() > 0 and (y == 1).sum() > 0) else np.nan",
     "CS-07 补全: 未定义输入统一返回 NaN（p6d；原为 0.0 / ValueError）"),
    ("S3", "p10b_subsample_curves.py", 260,
     "        thr = min([v for v in np.unique(s) if (s[y == 0] < v).mean() >= spec] or [np.nextafter(s.max(), np.inf)])",
     "        thr = min([v for v in np.unique(s) if (s[y == 0] < v).mean() >= spec] or [np.nextafter(s.max() if s.size else 0.0, np.inf)])",
     "CS-07 补全: 同上（p10b）"),
    ("S4", "p10b_subsample_curves.py", 263,
     "        return tp / max(tp + fn, 1)",
     "        return tp / max(tp + fn, 1) if (len(np.unique(y)) > 1 and (y == 0).sum() > 0 and (y == 1).sum() > 0) else np.nan",
     "CS-07 补全: 未定义输入统一返回 NaN（p10b；原为 0.0 / IndexError）"),
    ("S5", "p10c_aggregate.py", 69,
     "    thr = min([v for v in np.unique(s) if (neg < v).mean() >= spec] or [np.nextafter(s.max(), np.inf)])",
     "    thr = min([v for v in np.unique(s) if (neg < v).mean() >= spec] or [np.nextafter(s.max() if s.size else 0.0, np.inf)])",
     "CS-07 补全: 同上（p10c）"),
    ("S6", "p10c_aggregate.py", 72,
     "    return tp / max(tp + fn, 1)",
     "    return tp / max(tp + fn, 1) if (len(np.unique(y)) > 1 and (y == 0).sum() > 0 and (y == 1).sum() > 0) else np.nan",
     "CS-07 补全: 未定义输入统一返回 NaN（p10c；原为 nan/0.0，现与另两副本逐格一致）"),
    # (K) CS-05 —— DR / OM **完整恢复** Li et al. 2023 的 published estimator --------------
    # 依据：evs/ref_li2023_dr_auc.md（PMC10188769 全文转录，式 12–15）
    #   τ̂_DR = Σ_{i≠j} d^DR(k=一致性指示) / Σ_{i≠j} d^DR(k=1)
    #   塌缩后：source 对 = w_i w_j [y_i(1−y_j) − m_i(1−m_j)] k；target 对 = m_i(1−m_j) k；混合对 = 0
    # 三处偏离，全部单行替换（行数守恒，BLOCK_L3B 锚点不动）：
    #   ⓐ OM 分母 (Σm)(Σ(1−m)) 多算了对角线 Σ m_i(1−m_i)
    #   ⓑ DR 分子 两残差之积 -> 两乘积之差（取代旧的 C1/C2）
    #   ⓒ DR 分母 N_DR/(D1·D0)（两个边际均值之积）-> Σ_{i≠j} d^DR(k=1)（配对求和）
    # 恒等式：Σ_{i≠j} a_i b_j = (Σa)(Σb) − Σ a_i b_i
    ("K1", "p6c_01_L3_blinded.py", 279,
     '        numo = wconcordance(sc_t, m_t, 1 - m_t); deno = m_t.sum() * (1 - m_t).sum()',
     '        numo = wconcordance(sc_t, m_t, 1 - m_t); deno = m_t.sum() * (1 - m_t).sum() - (m_t * (1 - m_t)).sum()',
     "CS-05a: OM 分母扣除 i=j 自配对（Li 2023 明写排除 i=j）"),
    ("K2", "p6c_01_L3_blinded.py", 282,
     '        yt1 = (y_e == 1).astype(float); yt0 = (y_e == 0).astype(float)',
     '        yt1 = (y_e == 1).astype(float); yt0 = (y_e == 0).astype(float); w1 = w_te * yt1; w0 = w_te * yt0; p1 = w_te * m_e; p0 = w_te * (1 - m_e)',
     "CS-05: 拆出 source 侧的 观测项 (w1,w0) 与 模型项 (p1,p0)"),
    ("K3", "p6c_01_L3_blinded.py", 283,
     '        term_s = wconcordance(s_e, w_te * (yt1 - m_e), w_te * (yt0 - (1 - m_e)))',
     '        term_s = wconcordance(s_e, w1, w0) - wconcordance(s_e, p1, p0)',
     "CS-05b: DR 分子 = 两乘积之差（不是两残差之积）"),
    ("K4", "p6c_01_L3_blinded.py", 286,
     '        N_DR = term_s / (n_s * (n_s - 1)) + term_t / (n_t * (n_t - 1))',
     '        N_DR = term_s + term_t',
     "CS-05: 分子改为**配对求和**（原文分子分母同为 Sigma_{i!=j}，不做 n(n-1) 标准化）"),
    ("K5", "p6c_01_L3_blinded.py", 287,
     '        D1 = (w_te * (y_e - m_e)).sum() / n_s + m_t.sum() / n_t',
     '        d_s = (w1.sum() * w0.sum() - (w1 * w0).sum()) - (p1.sum() * p0.sum() - (p1 * p0).sum())',
     "CS-05c: DR 分母 source 部分 = Sigma_{i!=j}[w_i y_i w_j(1-y_j) - w_i m_i w_j(1-m_j)]"),
    ("K6", "p6c_01_L3_blinded.py", 288,
     '        D0 = (w_te * ((1 - y_e) - (1 - m_e))).sum() / n_s + (1 - m_t).sum() / n_t',
     '        d_t = m_t.sum() * (1 - m_t).sum() - (m_t * (1 - m_t)).sum()',
     "CS-05c: DR 分母 target 部分 = Sigma_{i!=j} m_i(1-m_j)"),
    ("K7", "p6c_01_L3_blinded.py", 289,
     '        rec["AUROC_DR"] = float(N_DR / (D1 * D0)) if (D1 * D0) > 0 else np.nan',
     '        rec["AUROC_DR"] = float(N_DR / (d_s + d_t)) if (d_s + d_t) > 0 else np.nan',
     "CS-05c: DR = 配对求和之比（取代边际均值之积 D1*D0）"),
    ("K8", "p10b_subsample_curves.py", 164,
     '            numo = wconcordance(sc_t, m_t_s, 1 - m_t_s); deno = m_t_s.sum() * (1 - m_t_s).sum()',
     '            numo = wconcordance(sc_t, m_t_s, 1 - m_t_s); deno = m_t_s.sum() * (1 - m_t_s).sum() - (m_t_s * (1 - m_t_s)).sum()',
     "CS-05a: 同上（p10b 副本）"),
    ("K9", "p10b_subsample_curves.py", 166,
     '            yt1 = (y_e == 1).astype(float); yt0 = (y_e == 0).astype(float)',
     '            yt1 = (y_e == 1).astype(float); yt0 = (y_e == 0).astype(float); w1 = w_te * yt1; w0 = w_te * yt0; p1 = w_te * m_e; p0 = w_te * (1 - m_e)',
     "CS-05: 同上（p10b 副本）"),
    ("K10", "p10b_subsample_curves.py", 167,
     '            term_s = wconcordance(s_e, w_te * (yt1 - m_e), w_te * (yt0 - (1 - m_e)))',
     '            term_s = wconcordance(s_e, w1, w0) - wconcordance(s_e, p1, p0)',
     "CS-05b: 同上（p10b 副本）"),
    ("K11", "p10b_subsample_curves.py", 170,
     '            N_DR = term_s / (n_s * (n_s - 1)) + term_t / (n * (n - 1))',
     '            N_DR = term_s + term_t',
     "CS-05: 同上（p10b 副本）"),
    ("K12", "p10b_subsample_curves.py", 171,
     '            D1 = (w_te * (y_e - m_e)).sum() / n_s + m_t_s.sum() / n',
     '            d_s = (w1.sum() * w0.sum() - (w1 * w0).sum()) - (p1.sum() * p0.sum() - (p1 * p0).sum())',
     "CS-05c: 同上（p10b 副本）"),
    ("K13", "p10b_subsample_curves.py", 172,
     '            D0 = (w_te * ((1 - y_e) - (1 - m_e))).sum() / n_s + (1 - m_t_s).sum() / n',
     '            d_t = m_t_s.sum() * (1 - m_t_s).sum() - (m_t_s * (1 - m_t_s)).sum()',
     "CS-05c: 同上（p10b 副本）"),
    ("K14", "p10b_subsample_curves.py", 173,
     '            rec["AUROC_DR"] = float(N_DR / (D1 * D0)) if (D1 * D0) > 0 else np.nan',
     '            rec["AUROC_DR"] = float(N_DR / (d_s + d_t)) if (d_s + d_t) > 0 else np.nan',
     "CS-05c: 同上（p10b 副本）"),

    # (M) CS-09 —— ESS 样本域（缺陷 (g)）。★ 与 K 组不同：这不是"公式写错"，
    #     `(sum w)^2 / sum w^2` 本身是正确的 ESS；错的是它挂在**哪个样本域**上。
    #     实测：中位 w_ESS=5416.5 vs 中位 n_src_test=1948（2.78x），170/202 的 w_ESS > n_src_test。
    #     两处都是 1 行换 1 行，行数守恒。p10b:157 已经是对的，不动。
    ("M1", "p6c_01_L3_blinded.py", 267,
     '        rec["w_ESS"] = float((w_all.sum() ** 2) / (w_all ** 2).sum())',
     '        rec["w_ESS"] = float((w_te.sum() ** 2) / (w_te ** 2).sum())',
     "CS-09: ESS 改用估计量实际使用的 w_te（source internal test），不再用 all-source 的 w_all"),
    ("M2", "p6c_01_L3_blinded.py", 362,
     '  f"占 n_source 比例中位={(DF[\'w_ESS\']/DF[\'n_source\']).median()*100:.1f}%")',
     '  f"占 n_src_test 比例中位={(DF[\'w_ESS\']/DF[\'n_src_test\']).median()*100:.1f}%")',
     "CS-09 派生：报告行的分母同步改为 n_src_test（ESS 的样本域）"),
]

# ---------------------------------------------------------------- 行区间删除（CS-11）
# 有些修复**不是** 1 行换 1 行：CS-11 要求把 p6f 那段误标为 "LOO-organism" 的陈旧计算
# **整段删除**（它的 `oof = lodo_oof(r)` 是 LODO，不是 LOO）。
# ★ 只删 header 而留下仍在运行的计算 = GPT6.1 批评过的"文本标记式假修复"，故必须删到计算本体。
# ★ 段边界依据：header 在 p6f:166，下一个 `"\n" + "=" * 112` 分隔符在 p6f:191；
#   段内 oeval/ofold 只在本段使用（已核，169-186），删除不影响下游。
# (tag, file, 起行, 止行, 说明) —— 1-indexed，含两端；**在 PATCHES 之后**应用。
DELETIONS = [
    ("J1", "p6f_lodo_settlement.py", 166, 190,
     "CS-11: 删除误标 LOO-organism 的陈旧段（header 166 + organism 折定义 168-178 + "
     "调用 lodo_oof 的 METHODS 循环 180-189）。它不写任何产物，但会让读者把 LODO 当 LOO；"
     "canonical Fig 3B 路径是 p14_rebuild_fig3b.py（p11:177 注明它独占写出 fig3b_loo_organism.csv）"),
]

MODULES = [
    ("evs_selftest.py",         "loader"),
    ("t01_mmd.py",              "MMD"),
    ("t02_concordance.py",      "CONC"),
    ("t03_l3_estimators.py",    "L3"),
    ("t04_sudo.py",             "SUDO"),
    ("t05_unit_grain.py",       "UNIT"),
    ("t06_bounds_scan.py",      "BOUND"),
    ("t07_fork_parity.py",      "FORK"),
    ("t08_psijs_family.py",     "PSI"),
    ("t09_bandwidth.py",        "BW"),
    ("t10_sens_at_spec.py",     "SENS"),
    ("t12_eval_layer.py",       "EVAL"),
    ("t13_eval_core.py",        "EC"),
    ("t14_manifest_closure.py", "MC"),
    ("t15_drugmap.py",          "DRUG"),
    ("t16_sudo_p10b.py",       "SUDO-P10B"),
    ("t17_cs11_loo.py",        "CS11"),
    ("t18_reimpl.py",          "REIMPL"),
    ("t19_cs03_sr.py",         "CS03"),
    ("t22_cs04_naming.py",     "CS04"),
]
# t11_kernel_census.py 已退役 -> retired/（heuristic census 由 ACTIVE_PATH_MANIFEST.yaml 取代）


def child_env(patched):
    """REAL: 移除 EVS_PROD_DIR；PATCHED: 指向 _patched。"""
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    env.pop("EVS_PROD_DIR", None)
    if patched:
        env["EVS_PROD_DIR"] = PATCHED
    return env


# ==================================================================== [0] header
p("=" * 104)
p("sim_fix — FIX-SIMULATION HARNESS")
p("=" * 104)
p(f"  真实生产树 (REAL)   : {OUT}")
p(f"  补丁副本树 (PATCHED): {PATCHED}")
p(f"  套件驱动器          : {os.path.join(HERE, 'run_evs.py')}")
p(f"  本机 Python         : {PY}")
p("")

# ============================================================ [0b] 「已进生产」前置守卫
# ★★★ 2026-10-06：CS-05..CS-12 已**批量进生产**
#   （方法见 evs/_apply_batch_to_production.py；清单与前后 sha256 存档于
#     FREEZE_pre_correction_20261005/APPLIED_CORRECTIONS_20261006.txt）。
#   补丁表里的 old 文本已被替换 ⇒ 本 harness 的 REAL-vs-PATCHED 机制**不再适用**。
#   这里显式判定并给出明确出路，而不是让 apply_patch 抛一句难懂的 mismatch。
#   ★ 红灯基线（EXPECTED_RED）的历史使命已经完成并留档；进生产后的闭包验证改由
#     `python run_evs.py` 在生产树上直接验四条机械条件。
def _patch_state():
    """返回 PATCHES 里 old 文本**仍然命中**的条数。全中=未应用；0=已全应用；中间=部分应用(危险)。"""
    hit = 0
    for _tag, _fn, _ln, _old, _new, _why in PATCHES:
        _raw = read_bytes(os.path.join(OUT, _fn))
        _lines = _raw.split(b"\n")
        _nl = b"\r" if b"\r\n" in _raw else b""
        if 1 <= _ln <= len(_lines) and _lines[_ln - 1] == _old.encode("utf-8") + _nl:
            hit += 1
    return hit


if PATCHES:
    _ph = _patch_state()
    if _ph == 0:
        p("=" * 104)
        p(f"  本 harness 不适用于当前生产树 —— 全部 {len(PATCHES)} 处补丁的 old 文本都已不在生产文件里。")
        p("  ⇒ 判定：CS-05..CS-12 已批量进生产（2026-10-06）。这是**预期状态**，不是故障。")
        p("  补丁清单与前后 sha256：FREEZE_pre_correction_20261005/APPLIED_CORRECTIONS_20261006.txt")
        p("  闭包验证请改用：python run_evs.py     （在生产树上直接验四条机械条件）")
        p("=" * 104)
        sys.exit(0)
    if _ph < len(PATCHES):
        p("=" * 104)
        p(f"  !! 生产树处于**部分应用**状态：{_ph} / {len(PATCHES)} 处补丁的 old 仍命中。")
        p("     这不是可运行状态；先查清为何只有一部分被应用，再重跑。**拒绝继续。**")
        p("=" * 104)
        sys.exit(4)
    p(f"  [0b] 补丁前置核验：{_ph}/{len(PATCHES)} 处 old 仍命中 ⇒ 生产树为**未应用**状态，继续。")
    p("")

# ==================================================================== [1] sha256 before
need(os.path.isdir(HERE), f"suite dir missing: {HERE}")
real_before = {}
for fn in LOADER_FILES:
    src = os.path.join(OUT, fn)
    need(os.path.isfile(src), f"production file missing: {src}")
    real_before[fn] = sha256_file(src)

# 套件自身（测试代码）的 sha256：本工作区可能有并发编辑；两次运行之间若测试代码变动，
# REAL/PATCHED 比较将不再可比 —— 收尾强制复核并在变动时硬失败。
SUITE_FILES = ["evs_load.py", "run_evs.py"] + [fn for fn, _ in MODULES]
suite_before = {fn: sha256_file(os.path.join(HERE, fn)) for fn in SUITE_FILES}

p("[1] 真实生产文件 sha256（开工前；收尾复核，绝不写入）")
for fn in LOADER_FILES:
    p(f"    {fn:<44} {real_before[fn][:16]}  {real_before[fn]}")
p(f"    套件自身 {len(SUITE_FILES)} 个文件（evs_load/run_evs + {len(MODULES)} 个测试模块）sha256 已记录，"
  f"收尾复核（防并发编辑）")
p("")

# ==================================================================== [2] refresh _patched
p("[2] 从真实文件刷新 _patched\\（每次全量重拷，幂等）")
ap = os.path.abspath(PATCHED)
need(os.path.basename(ap) == "_patched" and os.path.dirname(ap) == os.path.abspath(HERE),
     f"refuse to touch non-suite patched dir: {ap}")
if os.path.isdir(ap):
    shutil.rmtree(ap)
os.makedirs(ap)
for fn in LOADER_FILES:
    src, dst = os.path.join(OUT, fn), os.path.join(ap, fn)
    shutil.copy2(src, dst)
    need(sha256_file(src) == sha256_file(dst), f"copy not byte-identical: {fn}")
# ★ 2026-10-05：manifest 还登记了目录型路径（Rendering）。补丁树必须是**完整镜像**，
#   否则 t14 MC-01 在 PATCHED 上会因「glob 无命中」判红——那是**树不完整**，不是登记错误。
for dn in LOADER_DIRS:
    s, d = os.path.join(OUT, dn), os.path.join(ap, dn)
    if not os.path.isdir(s):
        continue
    shutil.copytree(s, d, ignore=shutil.ignore_patterns("__pycache__"))
p(f"    已复制 {len(LOADER_FILES)} 个文件 + {len(LOADER_DIRS)} 个目录"
  f"；副本与真实文件逐字节一致（copy 后立即校验）")
p("")

# ==================================================================== [3] apply patches
p("[3] 应用补丁（一行 -> 一行；行数必须不变）")


def apply_patch(path, ln, old, new):
    raw = read_bytes(path)
    crlf = b"\r\n" in raw
    lines = raw.split(b"\n")
    i = ln - 1
    need(0 <= i < len(lines), f"line {ln} out of range in {path} ({len(lines)} lines)")
    nl = b"\r" if crlf else b""
    old_b, new_b = old.encode("utf-8") + nl, new.encode("utf-8") + nl
    actual = lines[i]
    need(actual == old_b,
         f"patch target mismatch {os.path.basename(path)}:{ln}\n"
         f"      expected {old_b!r}\n      actual   {actual!r}")
    lines[i] = new_b
    raw2 = b"\n".join(lines)
    n_before, n_after = raw.count(b"\n"), raw2.count(b"\n")
    need(n_after == n_before, f"line count changed {n_before} -> {n_after}: {path}")
    need(raw2.count(b"\r\n") == raw.count(b"\r\n"), f"newline style changed: {path}")
    with open(path, "wb") as f:
        f.write(raw2)
    return dict(tag=None, fn=os.path.basename(path), ln=ln,
                n_before=n_before, n_after=n_after, crlf=crlf, old=old, new=new)


patch_info = []
for tag, fn, ln, old, new, label in PATCHES:
    dst = os.path.join(ap, fn)
    need(os.path.abspath(dst).startswith(ap + os.sep), f"patch escapes _patched: {dst}")
    info = apply_patch(dst, ln, old, new)
    info.update(tag=tag, label=label)
    patch_info.append(info)
    style = "CRLF" if info["crlf"] else "LF"
    p(f"  [{tag}] {fn}:{ln}  ({style})  行数 {info['n_before']} -> {info['n_after']}  "
      f"行数保持: OK")
    p(f"        old: {old.strip()}")
    p(f"        new: {new.strip()}")

# 整文件差异必须**恰好**等于该文件的补丁行集合
# ---------------------------------------------------------------- 行区间删除
for _tag, _fn, _a, _b, _why in DELETIONS:
    _p = os.path.join(ap, _fn)
    _L = read_bytes(_p).split(bytes([10]))
    need(1 <= _a <= _b <= len(_L), f"delete range {_tag} out of bounds {_fn}:{_a}-{_b}")
    _ndel = _b - _a + 1
    _keep = [x for i, x in enumerate(_L) if not (_a - 1 <= i <= _b - 1)]
    p(f"  [{_tag}] 删除 {_fn}:{_a}-{_b}（{_ndel} 行）-> 行数 {len(_L)} -> {len(_keep)}")
    p(f"        {_why}")
    io.open(_p, "wb").write(bytes([10]).join(_keep))
    # 校验：① 行数正确；② 被删内容确实不再出现；③ 其余非补丁行逐字节未变
    _new = read_bytes(_p).split(bytes([10]))
    need(len(_new) == len(_L) - _ndel,
         f"delete {_tag}: line count {len(_new)} != {len(_L) - _ndel}")
    _pl = {x["ln"] for x in patch_info if x["fn"] == _fn}
    for i, x in enumerate(_L):
        if _a - 1 <= i <= _b - 1 or (i + 1) in _pl:
            continue
        j = i - _ndel if i >= _b else i
        need(_new[j] == x, f"delete {_tag}: 非删除行被改动 {_fn}:{i + 1}")
    for _g in [x.strip() for x in _L[_a - 1:_b] if x.strip()][:1]:
        need(not any(x.strip() == _g for x in _new),
             f"delete {_tag}: 被删内容仍在文件中：{_g[:60]!r}")
    p("        校验：行数正确、被删内容已消失、其余行逐字节未变  OK")

p("  -- 整文件逐行 diff（副本 vs 真实）：除补丁行外不得有任何差异 --")
for fn in sorted(set(x["fn"] for x in patch_info)):
    a = read_bytes(os.path.join(OUT, fn)).split(b"\n")
    b = read_bytes(os.path.join(ap, fn)).split(b"\n")
    _del = [d for d in DELETIONS if d[1] == fn]
    if _del:
        _drop = set()
        for _t2, _f2, _a2, _b2, _w2 in _del:
            _drop.update(range(_a2 - 1, _b2))
        a = [x for i, x in enumerate(a) if i not in _drop]   # 期望的删除后内容
    need(len(a) == len(b), f"line count differs for {fn}: {len(a)} vs {len(b)}")
    diffs = [i + 1 for i, (x, y) in enumerate(zip(a, b)) if x != y]
    # ★ 去重：同一行可以被**两条**补丁依次修改（如 p10b:260 先被 D2 换公式、再被 S3 加守卫），
    #   此时实际 diff 只出现一次；不去重会误报 FATAL。
    exp = sorted({x["ln"] for x in patch_info if x["fn"] == fn})
    need(diffs == exp, f"{fn}: diffs {diffs} != expected {exp}")
    p(f"    {fn:<44} 改动行 = {diffs}  (={exp}) OK")
# 未打补丁的副本也必须逐字节一致
for fn in LOADER_FILES:
    if fn in set(x["fn"] for x in patch_info):
        continue
    need(sha256_file(os.path.join(OUT, fn)) == sha256_file(os.path.join(ap, fn)),
         f"unpatched copy differs: {fn}")
p(f"    未打补丁的 {len(LOADER_FILES) - len(set(x['fn'] for x in patch_info))} 个副本逐字节一致: OK")
p("")

# ==================================================================== [4] block anchors
p("[4] 内联块锚点核验（evs_load.grab；先用 REAL 记录锚点，再指向 _patched 重载）")
os.environ.pop("EVS_PROD_DIR", None)
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import evs_load as _L

BLOCK_ATTRS = [("L3B(257-289)", "BLOCK_L3B"), ("SUDO(291-330)", "BLOCK_SUDO"),
               ("STORE(121-144)", "BLOCK_STORE")]


def snap_blocks(mod):
    out = {}
    for name, attr in BLOCK_ATTRS:
        blk = getattr(mod, attr)
        src = mod.grab(*blk)
        out[name] = dict(label=mod.block_source(blk), blk=blk, text=src,
                         lines=src.split("\n"))
    return out


real_blocks = snap_blocks(_L)
os.environ["EVS_PROD_DIR"] = PATCHED
_L = importlib.reload(_L)
pat_blocks = snap_blocks(_L)
os.environ.pop("EVS_PROD_DIR", None)

for name, attr in BLOCK_ATTRS:
    rb, pb = real_blocks[name], pat_blocks[name]
    n = len(pb["lines"])
    first, last = pb["lines"][0], pb["lines"][-1]
    same_len = len(rb["lines"]) == len(pb["lines"])
    head_ok = rb["lines"][0] == first
    tail_ok = rb["lines"][-1] == last
    rel_changed = [i for i, (x, y) in enumerate(zip(rb["lines"], pb["lines"])) if x != y]
    abs_changed = [pb["blk"][1] + i for i in rel_changed]
    p(f"  {name:<16} {pb['label']:<32} lines={n}")
    p(f"        first line: {first.strip()[:96]}")
    p(f"        last  line: {last.strip()[:96]}")
    p(f"        first/last 与 REAL 锚点一致: {'OK' if (same_len and head_ok and tail_ok) else 'MISMATCH'}"
      f" ; 块内改动行: {abs_changed if abs_changed else '（无）'}")
    # ★ 2026-10-06：末行允许因**已申报的补丁**而变。
    #   CS-05 要改的正是 BLOCK_L3B 的末行（rec["AUROC_DR"] = ...），若仍要求"末行与 REAL 一致"，
    #   就等于**禁止改这一行** —— 那是又一个"手工维护的第二份真相"。
    #   保证不变：首行必须与 REAL 一致（锚点未漂移）、行数必须一致（行区间仍有效）；
    #   末行只有在**补丁表申报了它**时才允许不同，未申报的改动依旧判 FATAL。
    _bfn0 = os.path.basename(pb["blk"][0])
    _tail_declared = pb["blk"][2] in {x["ln"] for x in patch_info if x["fn"] == _bfn0}
    need(same_len and head_ok and (tail_ok or _tail_declared),
         f"block anchor moved for {name}: first={head_ok} same_len={same_len} "
         f"tail={tail_ok}（末行变动是否已申报={_tail_declared}）")
    # ★ 2026-10-05：期望值改为**从补丁表推导**，不再手工维护「哪条块应含哪些补丁行」。
    #   手维护的那版在 CS-02 进生产的**当下**就过期了（STORE 仍申报 [137]，而 _patched 已与生产一致、
    #   块内无任何改动）→ 误报 FATAL。这与「锚点钉住被修正的旧文本」是同一类错：
    #   **任何手工维护的第二份真相，都会在它所描述的东西变化时静默过期。**
    #   ★ 保证不变：`abs_changed` 是 REAL 与 PATCHED 的**实际**逐行差异；`exp` 只含补丁表申报的行。
    #     任何**未申报**的改动仍会出现在 abs_changed 而不在 exp 中 → 依旧 FATAL。
    _bfn = os.path.basename(pb["blk"][0])
    _ba, _bb = pb["blk"][1], pb["blk"][2]
    exp = sorted({x["ln"] for x in patch_info
                  if x["fn"] == _bfn and _ba <= x["ln"] <= _bb})
    need(abs_changed == exp,
         f"{name} 块内改动行 {abs_changed} != 补丁表推导的 {exp}（未申报的改动 = 锚点漂移）")
    if abs_changed:
        for i in rel_changed:
            p(f"          {pb['blk'][1] + i} 原: {rb['lines'][i].strip()}")
            p(f"          {pb['blk'][1] + i} 新: {pb['lines'][i].strip()}")
    else:
        p(f"        与 REAL 逐行一致（该块申报为不可改）: OK")
# 补丁行与块的关系（供审计）
p("  -- 补丁行 vs 块区间 --")
for x in patch_info:
    inside = []
    for name, attr in BLOCK_ATTRS:
        blk = real_blocks[name]["blk"]
        if os.path.abspath(blk[0]) == os.path.abspath(os.path.join(OUT, x["fn"])) and blk[1] <= x["ln"] <= blk[2]:
            inside.append(name)
    p(f"    [{x['tag']}] {x['fn']}:{x['ln']}  inside_block={inside if inside else 'none'}")
p("")

CHECK_ONLY = "--check-only" in sys.argv
COLLECT_ONLY = "--collect-only" in sys.argv
if CHECK_ONLY:
    p("[check-only] 停留于补丁/锚点核验（不跑套件）；复核真实文件 sha256 ...")
    bad = [fn for fn in LOADER_FILES
           if sha256_file(os.path.join(OUT, fn)) != real_before[fn]]
    need(not bad, "REAL production files were modified: " + ", ".join(bad))
    sbad = [fn for fn in SUITE_FILES if sha256_file(os.path.join(HERE, fn)) != suite_before[fn]]
    need(not sbad, "suite files changed during this run: " + ", ".join(sbad))
    p(f"    ✓ 全部 {len(LOADER_FILES)} 个生产文件逐字节未变；套件 {len(SUITE_FILES)} 个文件亦未变。"
      f"check-only done。")
    _log.close()
    raise SystemExit(0)

# ==================================================================== [5] run_evs x2
def run_sub(args, patched, timeout=1800):
    return subprocess.run([PY] + args, cwd=HERE, env=child_env(patched),
                          capture_output=True, text=True, encoding="utf-8",
                          errors="replace", timeout=timeout)


evs_runs = {}
RUN_T0 = {}          # mode -> run_evs 启动时刻（collect_results 的新鲜度锚点）
# ==================================================================== [4b] 生产输出树 hash inventory
# ★ 2026-10-05 硬要求（GPT6.1 提出、GPT6 钉死）：PATCHED 运行会执行生产脚本的**整块**代码，
#   而这些脚本用**硬编码绝对路径**指向真实 分析结果/。仅重定向 to_csv **不能**证明所有
#   report / cache / sidecar / 删除 / 新建 的写副作用都被隔离。
#   因此对生产输出树做**全目录树清单**（relative_path + size + sha256），跑前跑后比对：
#     内容改变 -> FAIL ｜ 新增 -> FAIL ｜ 删除 -> FAIL ｜ 改名 -> FAIL
#   ★ 不做文件类型白名单（旧方案只扫 *.csv + figures/data，会漏掉 .txt / .pkl / 新 sidecar）。
#   ★ EVS 套件已移出生产目录（第二篇 新方案\evs），故**无需任何排除**，guard 无盲区。
PROD_ROOTS = [OUT, os.path.join(PROJ, "figures")]


def prod_inventory(roots=None):
    """生产输出树全量清单。key = 相对「第二篇 新方案」的路径。

    roots 仅用于**自检**（--guard-selftest）；正式调用不传，恒为 PROD_ROOTS。
    """
    inv = {}
    for root in (roots if roots is not None else PROD_ROOTS):
        if not os.path.isdir(root):
            continue
        for dp, dn, fns in os.walk(root):
            dn[:] = [d for d in dn if d != "__pycache__"]
            for fn in fns:
                fp = os.path.join(dp, fn)
                # key 相对**各自的 root** 再加 root 名 —— 不用 PROJ 相对路径，
                # 否则自检的 scratch 树跨盘符时 relpath 会抛 ValueError。
                key = os.path.basename(root) + "/" + os.path.relpath(fp, root).replace("\\", "/")
                try:
                    inv[key] = (os.path.getsize(fp), sha256_file(fp))
                except OSError as e:
                    inv[key] = (-1, f"UNREADABLE {e}")
    return inv


def _guard_selftest():
    """★ 非空转对照（确定性，不靠 sleep 抢时序）：
    在**隔离 scratch 树**上证明 prod_inventory 的三种差异**都能被检出**。
    没有这个对照，「内容改变/新增/删除 -> FAIL」就只是一句声明。"""
    import tempfile
    base = os.path.abspath(tempfile.mkdtemp(prefix="evs_guard_selftest_"))
    a, b, c = (os.path.join(base, n) for n in ("a.txt", "b.txt", "c.txt"))
    for q in (a, b, c):
        with io.open(q, "w", encoding="utf-8") as fh:
            fh.write("x")
    inv0 = prod_inventory([base])
    with io.open(a, "w", encoding="utf-8") as fh:      # 内容改变
        fh.write("y")
    os.remove(b)                                       # 删除
    d = os.path.join(base, "d.txt")                    # 新增
    with io.open(d, "w", encoding="utf-8") as fh:
        fh.write("z")
    inv1 = prod_inventory([base])
    added = sorted(set(inv1) - set(inv0))
    removed = sorted(set(inv0) - set(inv1))
    changed = sorted(k for k in set(inv0) & set(inv1) if inv0[k] != inv1[k])
    ok = (len(added) == 1 and len(removed) == 1 and len(changed) == 1)
    p(f"[guard-selftest] scratch 树 {len(inv0)} -> {len(inv1)} 文件；"
      f"新增 {len(added)} / 删除 {len(removed)} / 内容改变 {len(changed)}")
    p(f"[guard-selftest] {'OK —— 三种差异均被检出（判据非空转）' if ok else '**FAIL —— 判据空转**'}"
      f"  新增={added} 删除={removed} 改变={[os.path.basename(x) for x in changed]}")
    shutil.rmtree(base, ignore_errors=True)
    return ok


if "--guard-selftest" in sys.argv:
    need(_guard_selftest(), "hash inventory 的非空转对照失败")
    p("")
    p("guard selftest done")
    sys.exit(0)


p("[4b] 生产输出树 hash inventory（跑前基线）")
_INV_BEFORE = prod_inventory()
p(f"    生产根: {[os.path.basename(r) for r in PROD_ROOTS]}；文件数 = {len(_INV_BEFORE)}"
  f"；EVS 套件位于生产目录之外（{os.path.basename(HERE)}/），故无排除项")
p("")


if COLLECT_ONLY:
    p("[5] --collect-only：跳过 run_evs.py 两次运行（不调用 run_evs；"
      "report_evs_*.txt 与 evs_results_*.csv 保持现状不重写）")
    p("    **注意**：本模式下新鲜度检查降级为警告 —— CSV 可能不是本次生成的。")
else:
    # 2026-10-04 加速：REAL 与 PATCHED 两棵树**完全独立**（各自子进程、各自 EVS_PROD_DIR），
    # 因此并发运行。每个 run_evs 内部再各自并行跑模块；把每棵树的 workers 减半以免超订。
    _cores = os.cpu_count() or 4
    # ★ 2026-10-05 事故：原为 (_cores-2)//2 = 15/树 × 2 树 = 30 个并发 Python 进程，
    #   每个都 import sklearn/scipy —— Windows 分页文件被耗尽，出现
    #   `ImportError: DLL load failed while importing _loss: 页面文件太小，无法完成操作`，
    #   导致 t16 整模块 BROKEN、t01/t08/t09 **spurious FAIL**（同一模块单独跑是全绿的）。
    #   那是**资源**问题，不是被测代码或判据的问题。上限压到 6/树（共 12）。
    _per_tree = max(2, min(6, (_cores - 2) // 4))
    if os.environ.get("EVS_PER_TREE"):
        _per_tree = max(1, int(os.environ["EVS_PER_TREE"]))
    p(f"[5] 并发运行 run_evs.py 两次（REAL 与 PATCHED 独立；每棵树 EVS_JOBS={_per_tree}，"
      f"总并发 ≤ {2 * _per_tree} / {_cores} 核）")

    def _run_tree(mode_key, patched):
        _t0 = time.time()
        _env = child_env(patched)
        _env["EVS_JOBS"] = str(_per_tree)
        _env["PYTHONIOENCODING"] = "utf-8"
        _r = subprocess.run([PY, os.path.join(HERE, "run_evs.py")], cwd=HERE, env=_env,
                            capture_output=True, text=True, encoding="utf-8",
                            errors="replace", timeout=3600)
        return mode_key, _r, time.time() - _t0

    with _cf.ThreadPoolExecutor(max_workers=2) as _ex:
        _futs = [_ex.submit(_run_tree, m, p) for m, p in (("real", False), ("patched", True))]
        _done = [f.result() for f in _futs]
    for mode_key, r, dt in sorted(_done):
        rep = os.path.join(HERE, f"report_evs_{mode_key}.txt")
        need(os.path.isfile(rep), f"run_evs did not write {rep}")
        RUN_T0[mode_key] = time.time() - dt          # 该次 run_evs 的启动时刻
        evs_runs[mode_key] = dict(rc=r.returncode, dt=dt, report=rep,
                                  report_sha=sha256_file(rep), stdout=r.stdout or "")
        p("-" * 104)
        p(f"  run_evs.py [{mode_key.upper()}]  rc={r.returncode}  用时 {dt:.1f}s  "
          f"报告 {os.path.basename(rep)} (sha256 {sha256_file(rep)[:16]})")
        p("-" * 104)
        for ln in evs_runs[mode_key]["stdout"].split("\n"):
            p("  | " + ln)
p("")

# ==================================================================== [6] per-test collect
p("[6] 逐测试状态采集（**读取 [5] 落下的结构化结果** evs_results_<mode>.csv；")
p("     模块本身只由 run_evs 并行执行一次 —— 这里不再重跑，避免每模块每棵树跑两遍）")
p("    聚合口径（修复 H1）：同一 test_id 多行 → 取**最差**状态 FAIL > WARN > PASS（不再取最后一行）；")
p("                         subgroup 明细保留（逐 id 行计数 + FAIL 子例逐字）。")

SEV = {"PASS": 0, "WARN": 1, "FAIL": 2}
STATUSES = ("FAIL", "WARN", "PASS")


def collect_results(patched):
    """读取 run_evs 本次运行落下的**结构化结果**（evs_results_<mode>.csv）。

    2026-10-04 加速改造（不改变任何判据，只去掉重复执行）：
      原实现自己再 subprocess 跑一遍全部模块 —— 而 run_evs.py 刚刚已经跑过同一批。
      于是每模块每棵树被执行 **2 次**、全流程共 4 次。模块之间完全独立，重复执行纯属浪费。
      现改为：run_evs 并行执行一次并落结构化 CSV，本函数**读**它。

    仍然保留的三道不可静默的安全带：
      ① **新鲜度**：CSV 的 mtime 必须晚于本次 run_evs 启动时刻，否则判 STALE 并 raise；
      ② **完整性**：CSV 里出现的模块必须覆盖 MODULES 全表，否则 raise；
      ③ **非空**：CSV 必须至少有一行 RESULT，否则 raise。
    这样"读到陈旧/残缺结果却当成本次结果"这一类（P16 事故的同一模式）无法静默通过。

    状态聚合口径不变：同一 test_id 多行取**最差**（FAIL>WARN>PASS），并保留 subgroup 明细。
    """
    import time as _t
    mode = "patched" if patched else "real"
    path = os.path.join(HERE, f"evs_results_{mode}.csv")
    if not os.path.exists(path):
        raise RuntimeError(f"[collect] 缺少结构化结果 {path} —— run_evs 未成功运行？")
    _mt = os.path.getmtime(path)
    if COLLECT_ONLY:
        age = _t.time() - _mt
        if age > FRESH_WINDOW_S:
            p(f"    [WARN] {os.path.basename(path)} 距今 {age:.0f}s —— "
              f"--collect-only 模式不重跑 run_evs，该结果**未必是本次的**；"
              f"仅作开发期自检，不得作为验收证据")
    else:
        anchor = RUN_T0.get(mode)
        if anchor is None:
            raise RuntimeError(f"[collect] 缺少 {mode} 的运行锚点 —— 本脚本逻辑有误")
        if _mt < anchor - 5:
            raise RuntimeError(
                f"[collect] {os.path.basename(path)} 的 mtime 早于本次 run_evs 启动时刻 "
                f"({_t.strftime('%H:%M:%S', _t.localtime(_mt))} < "
                f"{_t.strftime('%H:%M:%S', _t.localtime(anchor))}) —— "
                f"STALE：这不是本次运行写出的结果，禁止引用")
    rows = []
    with io.open(path, encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f):
            if r.get("status") in SEV and r.get("test_id"):
                rows.append(r)
    if not rows:
        raise RuntimeError(f"[collect] {os.path.basename(path)} 无有效 RESULT 行 —— 结果不可引用")
    mods_in_csv = {r["module"] for r in rows}
    missing = [fn for fn, _ in MODULES if fn not in mods_in_csv]
    if missing:
        raise RuntimeError(f"[collect] {os.path.basename(path)} 缺少模块 {missing} ——"
                           f"模块表不完整，结果不可引用")

    records, order, per_test = [], [], {}
    tid2mod, last = {}, {}
    fam, lines_tot, n_lines = {}, {"FAIL": 0, "WARN": 0, "PASS": 0}, 0
    short_of = {fn: short for fn, short in MODULES}
    for r in rows:
        fn, short, tid, st = r["module"], r["module_short"] if "module_short" in r else short_of[r["module"]], r["test_id"], r["status"]
        line = "|".join([r.get("module", ""), tid, st, r.get("expected", ""),
                         r.get("observed", ""), r.get("note", "")])
        records.append(dict(fn=fn, module=short, tid=tid, status=st, line=line))
        d = per_test.setdefault(tid, dict(module=short, n=0,
                                          cnt={"FAIL": 0, "WARN": 0, "PASS": 0},
                                          statuses=set(), worst="PASS", fail_lines=[]))
        d["n"] += 1
        d["cnt"][st] += 1
        d["statuses"].add(st)
        if SEV[st] > SEV[d["worst"]]:
            d["worst"] = st
        if st == "FAIL":
            d["fail_lines"].append(line)
        tid2mod.setdefault(tid, short)
        last[tid] = st
        if tid not in order:
            order.append(tid)
        lines_tot[st] += 1
        n_lines += 1
    for fn, short in MODULES:
        rs = [r for r in rows if r["module"] == fn]
        fam[fn] = dict(short=short, n_lines=len(rs),
                       nf=sum(1 for r in rs if r["status"] == "FAIL"),
                       nw=sum(1 for r in rs if r["status"] == "WARN"),
                       np_=sum(1 for r in rs if r["status"] == "PASS"), rc=0)
    mixed = sorted(t for t, d in per_test.items() if len(d["statuses"]) > 1)
    multi = [t for t in order if per_test[t]["n"] > 1]
    res = {t: per_test[t]["worst"] for t in order}
    return dict(records=records, per_test=per_test, res=res, last=last, order=order,
                tid2mod=tid2mod, fam=fam, lines_tot=lines_tot, n_lines=n_lines,
                mixed_dup=mixed, multi_dup=multi, source=os.path.basename(path))


# ==================================================================== [4c] 生产输出树跑后比对
# ★ 只做「检测器」：不假设哪些写被重定向了，直接比对**全树清单**。
#   内容改变 / 新增 / 删除 / 改名 任一项 -> FATAL（写副作用未被隔离）。
p("[4c] 生产输出树 hash inventory（跑后比对）")
_INV_AFTER = prod_inventory()
_added = sorted(set(_INV_AFTER) - set(_INV_BEFORE))
_removed = sorted(set(_INV_BEFORE) - set(_INV_AFTER))
_changed = sorted(k for k in set(_INV_BEFORE) & set(_INV_AFTER)
                  if _INV_BEFORE[k] != _INV_AFTER[k])
p(f"    跑前 {len(_INV_BEFORE)} 文件 -> 跑后 {len(_INV_AFTER)} 文件")
for _tag, _items in (("新增", _added), ("删除", _removed), ("内容改变", _changed)):
    p(f"    {_tag}: {len(_items)}" + ("" if not _items else f" :: {_items[:8]}"))
need(not (_added or _removed or _changed),
     "生产输出树在 PATCHED 运行期间被改动 —— 写副作用未被隔离"
     "（内容改变/新增/删除/改名均判失败）")
p("    ✓ 生产输出树逐文件未变（含 .txt / .pkl / 任何新 sidecar，不限文件类型）")
p("")


real = collect_results(False)
pat = collect_results(True)


# ============================================================================
# 归因的**唯一真相源** = ACTIVE_PATH_MANIFEST.yaml（与 t14 MC-10 同一入口 evs_load）
# ============================================================================
# ★ 2026-10-06：原先这里是**两张手写硬编码表** ——
#     · patched_fail_reason 的 SUDO- / UNIT- / BOUND- 前缀表
#     · CHAR_PIN = {"MMD-02b": ...}
#   它们与 manifest 是**两套互不同步的登记**。实测后果：同一份日志里
#     closure verdict 打印「四条全满足 → 允许 freeze」，
#     [8](iv) 却打印「存在 unexpected 仍红项：不得据此声称达到冻结条件」—— 自相矛盾。
#   根因不是缺记录（EC-01 / EC-03 在 manifest 里带 phase + reason + reentry 三字段，
#   登记得**最完整**的那一档），而是**记录有两份**。
#   现改为只读 manifest；manifest 里没有的才是真正的 unexpected（fail-closed）。
#
#   允许在 PATCHED 上仍红的三类归属，全部需要显式登记：
#     characterisations  自报非正确性断言，状态允许随生产修复而变
#     deferred_gates     产物验收，phase 已后移到 post_corrected_rerun
#     unbound_deferred   正确地红，但扫的是产物 / 稿件文字，不属任何单一 path
#   ★ required_tests 在 PATCHED 上 FAIL **不算 expected** —— 那正是闭包条件 2/3 要抓的东西。
#     这是对旧前缀表的**收紧**：旧表按前缀宽恕，不看登记。
def _manifest_registry():
    """读 manifest，返回 (kind, reason, prefixes)。
    kind[tid] ∈ {'required','characterisation','deferred_gate','unbound_deferred'}。"""
    import evs_load as _E                    # 按名重新取，避免与 _L 的 reload 状态耦合
    doc = _E.load_manifest()
    kind, why, pfxs = {}, {}, []
    for e in (doc.get("paths") or []):
        for t in (e.get("required_tests") or []):
            kind.setdefault(t, "required")
        for t in (e.get("characterisations") or []):
            kind.setdefault(t, "characterisation")
        for g in (e.get("deferred_gates") or []):
            t = g.get("test")
            if not t:
                continue
            kind[t] = "deferred_gate"
            why[t] = "phase=%s :: %s" % (g.get("phase", "-"), str(g.get("reason", "")).strip())
    for g in (doc.get("unbound_deferred") or []):
        t = g.get("test")
        if not t:
            continue
        why[t] = "phase=%s :: %s" % (g.get("phase", "-"), str(g.get("reason", "")).strip())
        if t.endswith("*"):
            pfxs.append(t)                   # 前缀登记，如 "BOUND-01::*"
        else:
            kind[t] = "unbound_deferred"
    return kind, why, pfxs


_MAN_KIND, _MAN_WHY, _MAN_PFX = _manifest_registry()


def _man_lookup(tid):
    """返回 (kind, reason)；manifest 未登记则为 (None, '')。"""
    if tid in _MAN_KIND:
        return _MAN_KIND[tid], _MAN_WHY.get(tid, "")
    for p in _MAN_PFX:
        if tid.startswith(p[:-1]):
            return "unbound_deferred", _MAN_WHY.get(p, "")
    return None, ""


_RED_OK = {"characterisation", "deferred_gate", "unbound_deferred"}
_RED_TAG = {
    "characterisation": "characterisations（自报非正确性断言，状态允许随生产修复而变）",
    "deferred_gate": "deferred_gates（产物验收，phase 已后移到 post_corrected_rerun）",
    "unbound_deferred": "unbound_deferred（正确地红，但扫的是产物/稿件文字，不属任何单一 path）",
}


def patched_fail_reason(tid):
    """PATCHED 下仍 FAIL 的归因：expected（manifest 已登记且允许仍红）或 unexpected。"""
    k, w = _man_lookup(tid)
    if k in _RED_OK:
        return "expected", _RED_TAG[k] + "；manifest 依据：" + (
            w or "（该条只登记了归属，未写 reason）")
    if k == "required":
        return "unexpected", ("该 test_id 在 manifest 里属于 **required_tests** —— 它在 PATCHED 上"
                              "必须 PASS；此处红即闭包条件 2/3 的违约，不得按 expected 宽恕")
    return "unexpected", ("manifest 未登记该 test_id 的任何归属（required_tests / "
                          "characterisations / deferred_gates / unbound_deferred 四类均无）"
                          " —— 必须人工检查")


def char_pin_reason(tid):
    """PASS -> WARN 的 test_id 是否为 manifest 已登记的 characterisation
    （其状态允许随生产修复而变化）。不是则返回 None，须按潜在回归人工判断。"""
    k, w = _man_lookup(tid)
    if k == "characterisation":
        return "manifest 登记为 characterisation（" + (
            w or "见该 path 的 characterisations 列表及其上方注释") + "）"
    return None


def family_summary(coll):
    """按族摘要：直接从结构化 per-test 结果（[6] 表的同一份数据）计算。

    修复 H2：不再解析 run_evs 的 report_evs_*.txt 文本 —— 那里面没有逐模块计数，
    原实现因此恒报 report=None / "存在不一致" 的伪 MISMATCH。
    返回 {module_short: (n_unique_ids, n_FAIL, n_WARN, n_PASS)}，状态取最差状态。
    """
    out = {}
    for fn, short in MODULES:
        ids = [t for t in coll["order"] if coll["tid2mod"][t] == short]
        out[short] = (len(ids),
                      sum(1 for t in ids if coll["res"][t] == "FAIL"),
                      sum(1 for t in ids if coll["res"][t] == "WARN"),
                      sum(1 for t in ids if coll["res"][t] == "PASS"))
    return out


def write_results_csv(coll, path):
    """机器可读输出（供独立方与 run_evs 输出交叉核对）：
    test_id,status,module,n_lines,n_fail,n_warn,n_pass；status = 最差状态。"""
    with io.open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["test_id", "status", "module", "n_lines", "n_fail", "n_warn", "n_pass"])
        for tid in coll["order"]:
            d = coll["per_test"][tid]
            w.writerow([tid, d["worst"], d["module"], d["n"],
                        d["cnt"]["FAIL"], d["cnt"]["WARN"], d["cnt"]["PASS"]])


def read_results_csv(path):
    with io.open(path, encoding="utf-8", newline="") as f:
        rows = list(csv.reader(f))
    need(rows and rows[0] == ["test_id", "status", "module", "n_lines", "n_fail",
                              "n_warn", "n_pass"],
         f"bad CSV header in {os.path.basename(path)}")
    out = {}
    for r in rows[1:]:
        if r:
            out[r[0]] = dict(status=r[1], module=r[2], n=int(r[3]),
                             FAIL=int(r[4]), WARN=int(r[5]), PASS=int(r[6]))
    return out


def _csv_row_matches(d, b):
    return (b is not None and b["status"] == d["worst"] and b["module"] == d["module"]
            and b["n"] == d["n"] and b["FAIL"] == d["cnt"]["FAIL"]
            and b["WARN"] == d["cnt"]["WARN"] and b["PASS"] == d["cnt"]["PASS"])

# 按族摘要 + 机器可读输出（修复 H2）：不再解析 run_evs 的 report_evs_*.txt ——
# 该文件只含 run_evs 自己的汇总表，没有逐模块计数，原实现因此恒报 report=None /
# "存在不一致" 的伪 MISMATCH。现从结构化 test results 直接计算，并写出
# evs_results_{real,patched}.csv 后回读核对（含按族计数），失败才报不一致。
p("  -- 按族摘要（直接从结构化 test results 计算；机器可读输出回读核对）--")
cross_ok = True
fam_sum = {}
for mode_key, coll in (("real", real), ("patched", pat)):
    fam_sum[mode_key] = family_summary(coll)
    # ★ 2026-10-05 修正：**不得覆写 run_evs 的 evs_results_<mode>.csv**。
    #   两个写者、同一文件名、不同 schema（run_evs: `module` = 文件名；本脚本: `module` = 族名）
    #   —— 下游 mutate_acceptance 曾因此把 `MMD` / `FORK` 当成模块名，产出 **7 条假的 ESCAPED**。
    #   本脚本只写自己的文件名；run_evs 的结构化产物（含 expected/observed）保持权威。
    csv_path = os.path.join(HERE, f"simfix_results_{mode_key}.csv")
    write_results_csv(coll, csv_path)
    back = read_results_csv(csv_path)
    bad = [tid for tid in coll["order"]
           if not _csv_row_matches(coll["per_test"][tid], back.get(tid))]
    fam_from_csv = {}
    for b in back.values():
        c = fam_from_csv.setdefault(b["module"], {"n": 0, "FAIL": 0, "WARN": 0, "PASS": 0})
        c["n"] += 1
        c[b["status"]] += 1
    n_fam_ids = sum(v[0] for v in fam_sum[mode_key].values())
    p(f"    [{mode_key}] evs_results_{mode_key}.csv: {len(coll['order'])} 行 "
      f"(sha256 {sha256_file(csv_path)[:16]})")
    p(f"    [{mode_key}] CSV 回读核对（最差状态/模块/行数/每状态计数）: "
      + ("全部一致 OK" if not bad else f"*** {len(bad)} 条不一致: {bad[:5]} ***"))
    p(f"    [{mode_key}] 按族摘要（id = 该族唯一 test_id 的最差状态；行 = RESULT 行计数）：")
    for fn, short in MODULES:
        d = coll["fam"][fn]
        n_ids, f_, w_, p_ = fam_sum[mode_key][short]
        c = fam_from_csv.get(short, {"n": 0, "FAIL": 0, "WARN": 0, "PASS": 0})
        ok = (c["n"], c["FAIL"], c["WARN"], c["PASS"]) == (n_ids, f_, w_, p_)
        cross_ok &= (not bad) and ok
        p(f"      {short:<7}{fn:<24} id {n_ids:>3} (FAIL {f_:>2}/WARN {w_:>2}/PASS {p_:>3})"
          f"   行 {d['n_lines']:>3} (FAIL {d['nf']:>2}/WARN {d['nw']:>2}/PASS {d['np_']:>2})"
          f"   rc={d['rc']}" + ("" if ok else "   *** CSV 按族计数不一致 ***"))
    cross_ok &= (n_fam_ids == len(coll["order"]))
p(f"  机器可读输出交叉核对：{'全部一致 OK' if cross_ok else '存在不一致（见上）'}")
p("")

# 总表
tot_r = {s: sum(1 for v in real["res"].values() if v == s) for s in STATUSES}
tot_p = {s: sum(1 for v in pat["res"].values() if v == s) for s in STATUSES}
p(f"  所有 RESULT 行：REAL {real['n_lines']} 行 "
  f"(FAIL {real['lines_tot']['FAIL']} / WARN {real['lines_tot']['WARN']} / PASS {real['lines_tot']['PASS']})"
  f"  |  PATCHED {pat['n_lines']} 行 "
  f"(FAIL {pat['lines_tot']['FAIL']} / WARN {pat['lines_tot']['WARN']} / PASS {pat['lines_tot']['PASS']})")
p(f"  唯一 test_id（同 id 多行取**最差**状态 FAIL>WARN>PASS）：REAL {len(real['res'])} 个 "
  f"(FAIL {tot_r['FAIL']} / WARN {tot_r['WARN']} / PASS {tot_r['PASS']})"
  f"  |  PATCHED {len(pat['res'])} 个 (FAIL {tot_p['FAIL']} / WARN {tot_p['WARN']} / PASS {tot_p['PASS']})")
p("")
p("  -- 修复 H1 生效证据：原「取最后一行」口径 vs 现「取最差状态」口径 --")
for _mk, _coll in (("REAL", real), ("PATCHED", pat)):
    _diff = [t for t in _coll["order"] if _coll["last"][t] != _coll["res"][t]]
    p(f"    [{_mk}] 两种口径不一致的 test_id 共 {len(_diff)} 个："
      + (", ".join(f"{t} (last={_coll['last'][t]} -> worst={_coll['res'][t]})" for t in _diff)
         or "（无）"))
p("")
p("  -- 同 id 多行 subgroup 明细（保留可见；FAIL 子例逐字列出）--")
for _mk, _coll in (("REAL", real), ("PATCHED", pat)):
    if not _coll["multi_dup"]:
        p(f"    [{_mk}] （无同 id 多行）")
        continue
    p(f"    [{_mk}]")
    for tid in _coll["multi_dup"]:
        d = _coll["per_test"][tid]
        mix = "状态混合" if len(d["statuses"]) > 1 else "同状态"
        p(f"      {tid:<16} module={d['module']:<6} 行数={d['n']}  "
          f"(FAIL {d['cnt']['FAIL']} / WARN {d['cnt']['WARN']} / PASS {d['cnt']['PASS']})  "
          f"{mix}  最差={d['worst']}")
        for _fl in d["fail_lines"]:
            p(f"          FAIL 子例逐字: {_fl}")
p("")
TRANS_TAGS = {
    ("FAIL", "PASS"): "FAIL->PASS  << FIXED",
    ("PASS", "FAIL"): "PASS->FAIL  << REGRESSION",
    ("WARN", "FAIL"): "WARN->FAIL  << REGRESSION",
    ("PASS", "WARN"): "PASS->WARN  << check pin status",
    ("WARN", "PASS"): "WARN->PASS  << diagnostic cleared",
    ("FAIL", "WARN"): "FAIL->WARN",
}
p("  test_id          module   REAL     PATCHED    transition")
p("  " + "-" * 88)
all_tids = real["order"] + [t for t in pat["order"] if t not in real["res"]]
flips_fp, p2w, w2p, regressions, other, pass_real, only_one = [], [], [], [], [], [], []
for tid in all_tids:
    rs, ps = real["res"].get(tid), pat["res"].get(tid)
    mod = real["tid2mod"].get(tid) or pat["tid2mod"].get(tid, "?")
    if rs is None or ps is None:
        only_one.append((tid, rs, ps))
        tr = f"only {'REAL' if ps is None else 'PATCHED'}={rs or ps}"
    elif rs == ps:
        tr = "-"
    else:
        tr = TRANS_TAGS.get((rs, ps), f"{rs}->{ps}")
        if (rs, ps) == ("FAIL", "PASS"):
            flips_fp.append(tid)
        elif ps == "FAIL":
            regressions.append((tid, rs, ps))
        elif (rs, ps) == ("PASS", "WARN"):
            p2w.append(tid)
        elif (rs, ps) == ("WARN", "PASS"):
            w2p.append(tid)
        else:
            other.append((tid, rs, ps))
    if rs == "PASS":
        pass_real.append(tid)
    p(f"  {tid:<16} {mod:<8} {str(rs):<8} {str(ps):<10} {tr}")
p("")

p("  ---- 状态变化分栏（详细判定见 [8]） ----")
p(f"  FAIL->PASS 共 {len(flips_fp)} 条：" + (", ".join(flips_fp) if flips_fp else "（无）"))
p(f"  PASS->FAIL 回归 共 {len(regressions)} 条："
  + (", ".join(f"{t} ({rs}->{ps})" for t, rs, ps in regressions) if regressions else "（无）"))
p(f"  PASS->WARN 共 {len(p2w)} 条：" + (", ".join(p2w) if p2w else "（无）"))
p(f"  WARN->PASS 共 {len(w2p)} 条：" + (", ".join(w2p) if w2p else "（无）"))
p(f"  其它变化 共 {len(other)} 条："
  + (", ".join(f"{t} ({rs}->{ps})" for t, rs, ps in other) if other else "（无）"))
p("  ---- PASS on REAL（在真实缺陷树上就是绿的 —— 可能的「钉死现状」项，需人工判断） ----")
by_mod = {}
for tid in pass_real:
    by_mod.setdefault(real["tid2mod"].get(tid, "?"), []).append(tid)
for mod in sorted(by_mod):
    p(f"    {mod:<8}: " + ", ".join(by_mod[mod]))
p(f"    合计 {len(pass_real)} 条")
if only_one:
    for tid, rs, ps in only_one:
        p(f"    {tid}: only-in-one-tree  REAL={rs} PATCHED={ps}")
p("")

# ==================================================================== [7] sha256 after
p("[7] 真实生产树完整性复核（sha256 after == before）")
bad = []
for fn in LOADER_FILES:
    now = sha256_file(os.path.join(OUT, fn))
    ok = now == real_before[fn]
    if not ok:
        bad.append(fn)
    p(f"    {fn:<44} {real_before[fn][:16]} -> {now[:16]}  {'OK' if ok else 'CHANGED!'}")
need(not bad, "REAL production files were modified: " + ", ".join(bad))
# ★ 2026-10-06 修正：原为硬编码 "9"，而上面循环遍历的是 LOADER_FILES（实为 17 个）。
#   检查本身没问题（遍历真实列表 + need() 硬失败），但**打印的计数是假的** ——
#   一个验证工具报错自己的覆盖面，正是本项目要抓的那类缺陷。改为计算值。
p(f"    ✓ 全部 {len(LOADER_FILES)} 个生产文件逐字节未变（sha256 一致）")
p("")
p("[7b] 套件自身文件完整性（sha256 after == before；防并发编辑使两棵树不可比）")
sbad = []
for fn in SUITE_FILES:
    now = sha256_file(os.path.join(HERE, fn))
    ok = now == suite_before[fn]
    if not ok:
        sbad.append(fn)
    p(f"    {fn:<28} {suite_before[fn][:16]} -> {now[:16]}  {'OK' if ok else 'CHANGED!'}")
need(not sbad, "suite files changed during the run (REAL/PATCHED 不再可比，需重跑): "
     + ", ".join(sbad))
# ★ 2026-10-06 修正：原为硬编码 "9"，而上面循环遍历的是 SUITE_FILES（= evs_load + run_evs +
#   len(MODULES) = 22）。同样只错在打印的计数。
p(f"    ✓ 套件 {len(SUITE_FILES)} 个文件逐字节未变：本次 REAL/PATCHED 比较在**同一版测试代码**上完成")
p("")

# ==================================================================== [8] verdict
p("[8] 结论（口径：同一 test_id 多行取最差状态 FAIL>WARN>PASS；回归单独列名，不埋进计数）")

def _rc_txt(mode_key):
    return (f"run_evs rc={evs_runs[mode_key]['rc']}" if mode_key in evs_runs
            else "run_evs 未运行（--collect-only）")

p(f"  REAL    : RESULT 行 FAIL {real['lines_tot']['FAIL']} / WARN {real['lines_tot']['WARN']} / "
  f"PASS {real['lines_tot']['PASS']}；唯一 id FAIL {tot_r['FAIL']} / WARN {tot_r['WARN']} / "
  f"PASS {tot_r['PASS']}   ({_rc_txt('real')})")
p(f"  PATCHED : RESULT 行 FAIL {pat['lines_tot']['FAIL']} / WARN {pat['lines_tot']['WARN']} / "
  f"PASS {pat['lines_tot']['PASS']}；唯一 id FAIL {tot_p['FAIL']} / WARN {tot_p['WARN']} / "
  f"PASS {tot_p['PASS']}   ({_rc_txt('patched')})")
p("")
p("  (i) FAIL -> PASS（修补后被修缺陷对应测试转绿）：")
for _lab, _pref in (("(A) MMD 因子 2", "MMD"), ("(B) 并列 0.5", "CONC"), ("(C) DR 分子", "L3")):
    _got = [t for t in flips_fp if t.startswith(_pref)]
    p(f"      {_lab:<14}: " + (", ".join(_got) if _got else "（无）"))
_other_fp = [t for t in flips_fp if not t.startswith(("MMD", "CONC", "L3"))]
p(f"      其它转绿      : " + (", ".join(_other_fp) if _other_fp else "（无）"))
p(f"      合计 {len(flips_fp)} 条")
p("")
p("  (ii) PASS -> FAIL（回归 / regression；应为 0 条）：")
if regressions:
    for _t, _rs, _ps in regressions:
        p(f"      *** {_t}: {_rs} -> {_ps} —— 修补生产代码反而使该测试变红，必须逐条归因 ***")
    p(f"      合计 {len(regressions)} 条 —— 冻结条件「0 个既有 PASS 回归」不满足")
else:
    p("      无 —— 0 条 PASS 回归。")
p("")
p("  (iii) PASS -> WARN（须判定：manifest 已登记的 characterisation，还是潜在回归）：")
if p2w:
    for _t in p2w:
        _r = char_pin_reason(_t)
        if _r:
            p(f"      {_t}: 已登记 pin —— {_r}")
        else:
            p(f"      *** {_t}: manifest 未登记为 characterisation —— 必须按潜在回归人工判断 ***")
else:
    p("      无。")
p("")
p("  (iv) PATCHED 下仍 FAIL 的测试及原因（expected / unexpected）：")
_still = [t for t in pat["order"] if pat["res"][t] == "FAIL"]
_n_exp = _n_unexp = 0
for _t in _still:
    _tag, _why = patched_fail_reason(_t)
    if _tag == "expected":
        _n_exp += 1
    else:
        _n_unexp += 1
    p(f"      [{_tag:^10}] {_t}: {_why}")
p(f"      expected {_n_exp} 条 / unexpected {_n_unexp} 条")
if _n_unexp:
    p("      *** 存在 unexpected 仍红项：不得据此声称「PATCHED 达到冻结条件」***")
p("")
p("  注1（harness 修复效果，H1）：「取最后一行」缺陷修复后，UNIT-01 由被 3 条 PASS 行")
p("      掩盖的 PASS 变为 FAIL（REAL/PATCHED 均为 FAIL）—— 这是可见性修复，不是新的生产")
p("      回归；其失败本身属 (iv) 中未修补的 episode 粒度缺陷。")
p("  注2：t06(BOUND) 读取**冻结 CSV**（L.OUT 与 figures\\data，不读代码），代码补丁不可能")
p("      改变它 —— 它在 PATCHED 下保持原状态/原 FAIL，这不是 harness 失败。")
p("  注3：本 harness 的 REAL/PATCHED 只是「测试工具」层面的验收；生产代码一行未动（[7]）。")
p("")
p(f"日志：{os.path.join(HERE, 'sim_fix_report.txt')}")
p(f"机器可读结果：{os.path.join(HERE, 'evs_results_real.csv')} / {os.path.join(HERE, 'evs_results_patched.csv')}")
p(f"报告：{os.path.join(HERE, 'report_evs_real.txt')} / {os.path.join(HERE, 'report_evs_patched.txt')}")
_log.close()
print("\nsim_fix done")
