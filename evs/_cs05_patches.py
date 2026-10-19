# -*- coding: utf-8 -*-
r"""CS-05 补丁集：DR / OM **完整恢复** Li et al. 2023 的 published estimator。
依据 evs/ref_li2023_dr_auc.md（式 12–15）。
移除被取代的 C1/C2（按**行区间**删除），加入 K1–K14（全部单行替换、行数守恒）。"""
import io, os, re, sys, ast

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
F = "sim_fix.py"
L = io.open(F, encoding="utf-8").read().split("\n")


def find(sub, start=0):
    for i in range(start, len(L)):
        if sub in L[i]:
            return i
    raise LookupError(sub)


# ---- 1) 按行区间移除 C1..C2 ----
i0 = find('("C1", "p6c_01_L3_blinded.py", 283,')
i1 = find('"DR 分子: 两乘积之差（p10b:167 同源副本）"),')
assert i0 < i1, (i0, i1)
print("移除 C1..C2：行", i0 + 1, "-", i1 + 1)
del L[i0:i1 + 1]

NEW = '''    # (K) CS-05 —— DR / OM **完整恢复** Li et al. 2023 的 published estimator --------------
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
]'''.split("\n")
# ★ NEW 的字面量末尾自带一个 `]`，而原有列表的 `]` 会被保留 —— 必须去掉这一个，否则重复。
if NEW[-1].strip() == "]":
    NEW.pop()

# 找到 PATCHES 列表的收尾 ']'：**从最后一条补丁之后**找第一个纯 ']' 行。
# （不能"从后往前找紧邻 MODULES 的 ]"——DELETIONS 列表也在 MODULES 前面，会插错位置。）
k0 = find('("S6", "p10c_aggregate.py", 72,')
j = None
for i in range(k0, len(L)):
    if L[i].strip() == "]":
        j = i
        break
assert j is not None, "找不到 PATCHES 的收尾 ]"
L[j:j] = NEW
out = "\n".join(L)
ast.parse(out)
io.open(F, "w", encoding="utf-8", newline="\n").write(out)

ns = {}
m = re.search(r"^PATCHES = \[(.*?)^\]", out, re.S | re.M)
exec("PATCHES = [" + m.group(1) + "]", ns)
print("补丁表现在", len(ns["PATCHES"]), "条：", [x[0] for x in ns["PATCHES"]])

# ---- dry-check：按序模拟应用 ----
OUT = os.path.join(os.path.dirname(HERE), "分析结果")
cur, bad = {}, 0
for tag, fn, ln, old, new, why in ns["PATCHES"]:
    if fn not in cur:
        cur[fn] = io.open(os.path.join(OUT, fn), encoding="utf-8").read().split("\n")
    A = cur[fn]
    if not (1 <= ln <= len(A)) or A[ln - 1].rstrip("\r") != old:
        print(f"  !! {tag} {fn}:{ln} 不命中")
        print(f"     old={old[:95]!r}")
        print(f"     act={(A[ln-1].rstrip(chr(13))[:95] if 1<=ln<=len(A) else 'OOB')!r}")
        bad += 1
        continue
    A[ln - 1] = new
print("按序模拟应用：不一致 =", bad)
