# -*- coding: utf-8 -*-
"""金丝雀：证明 L3-02 的非空转守卫**真的会因退化 fixture 而判红**，不是装饰。

不改任何文件；只 import t03 的 helper 与 fixture。四段：
  A 真 fixture（连续 m）—— 守卫成立，判据的鉴别力来自 9.76e-03 的间距
  B 把 m_t 二值化 —— 证明两种分母形式在二值 m 上给出**同一个数**（间距恰为 0）
  C 二值 fixture 直接跑 _test_02 —— 守卫报 nonvacuous=False
  D 隔离实验 —— 人为令 code == ref（absdiff=0），判红的**唯一**来源就是守卫
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import t03_l3_estimators as T

fx = T.build_fixture()

print("=== A. true fixture (continuous m) ===")
ok, obs = T._test_02(fx)
print("  ok =", ok)
print(" ", obs)

print()
print("=== B. binarise m_t -- the degenerate fixture family ===")
m = np.asarray(fx["m_t"], float)
m_bin = (m >= 0.5).astype(float)
print("  m_bin value set =", sorted(set(m_bin.tolist())))
den_pairs = T._bf_om_den_pairs(m_bin)
den_prod = T._loop_sum(m_bin) * T._loop_sum(1.0 - m_bin)
gap = abs(den_prod - den_pairs)
print("  den_pairs = %.10f" % den_pairs)
print("  den_prod  = %.10f" % den_prod)
print("  gap       = %.3e   <-- identically zero for binary m" % gap)
print("  guard  gap > 1e-9  ->", gap > 1e-9)

print()
print("=== C. run _test_02 on the binarised fixture ===")
ok_c, obs_c = T._test_02(dict(fx, m_t=m_bin))
print("  ok =", ok_c)
print(" ", obs_c)

print()
print("=== D. isolation: force code == degenerate reference (absdiff = 0) ===")
fx_d = dict(fx, m_t=m_bin)
ref_d = T._bf_pair(fx_d["sc_t"], m_bin, 1.0 - m_bin, tie=False) / T._bf_om_den_pairs(m_bin)
rec_d = dict(fx["rec"])
rec_d["AUROC_OM"] = ref_d
fx_d["rec"] = rec_d
ok_d, obs_d = T._test_02(fx_d)
print("  forced code = ref = %.12f  (so absdiff is 0)" % ref_d)
print("  ok =", ok_d)
print(" ", obs_d)
print()
print("  >>> absdiff==0 yet ok=%s -- the ONLY source of the red is the non-vacuity guard."
      % ok_d)
