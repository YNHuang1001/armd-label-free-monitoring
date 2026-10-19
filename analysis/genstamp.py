# -*- coding: utf-8 -*-
r"""genstamp — 内部产物的**代次头**（A3.3）

================================ 为什么存在（INC-001 的病历）================================
`task_profile.csv` 陈旧了三代（202 宇宙时代），**没有任何东西发现它** ——
因为**文件不声明自己的一代**，陈旧只能靠内容法医识别。
202 能活三轮，就是这个缺口。

================================ 做什么 ================================
每个内部产物旁写一个 sidecar `<artifact>.gen.json`，记录：
  · `producer`       —— 产出它的脚本
  · `universe_sha`   —— **上游宇宙（`p6b_corrected_universe.csv`）在产出时点的 sha256**
  · `universe_rows`  —— 该宇宙的行数（人可读的快速判据）
  · `written_at`     —— UTC 时间戳
  · `artifact_sha`   —— 产物自身的 sha256

⇒ 宇宙一变，`universe_sha` 立刻对不上，**陈旧从「法医问题」变成「比对一个字符串」**。

**不写进 CSV 本身**：加列会改 schema，波及 p21c / p21f / BOUND 规则；
sidecar 零侵入。

用法（在生产者脚本里）：
    import genstamp
    genstamp.stamp(path, "p21c")
"""
import hashlib
import io
import json
import os
import shutil
from datetime import datetime, timezone

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UNIVERSE = os.path.join(ROOT, "分析结果", "p6b_corrected_universe.csv")


def _sha(path):
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def universe_id():
    """上游宇宙的指纹（sha256 前 16 位 ＋ 行数）。"""
    h = _sha(UNIVERSE)
    n = None
    if h:
        with io.open(UNIVERSE, encoding="utf-8") as f:
            n = sum(1 for _ in f) - 1
    return {"universe_sha": (h[:16] if h else None), "universe_rows": n}


def stamp(path, producer, extra=None):
    """给产物 `path` 写代次头 sidecar。返回 sidecar 路径。"""
    d = {"producer": producer,
         "artifact": os.path.basename(path),
         "artifact_sha": (_sha(path) or "")[:16],
         "written_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
    d.update(universe_id())
    if extra:
        d.update(extra)
    out = path + ".gen.json"
    with io.open(out, "w", encoding="utf-8") as f:
        json.dump(d, f, ensure_ascii=False, indent=2, sort_keys=True)
    return out


def read(path):
    p = path + ".gen.json"
    if not os.path.exists(p):
        return None
    with io.open(p, encoding="utf-8") as f:
        return json.load(f)


def is_current(path):
    """产物是否与**当前**宇宙同代。无 sidecar ⇒ 判为「未声明代次」（既非当前也非陈旧）。"""
    g = read(path)
    if g is None:
        return None
    cur = universe_id()
    return (g.get("universe_sha") == cur["universe_sha"]
            and g.get("universe_rows") == cur["universe_rows"])


# ================================================================ 覆盖前自动归档（A3.6）
# GLM 原话：「归档该是**跑批器的属性** —— 任何要覆盖注册产物的脚本，覆盖前自动归档前一版。
#            与『防陈旧靠失效传播不靠警觉』是**同一个动作的两次出现**。」
#
# 病历：INC-002 —— 重跑 `p21f` 前未归档其输出，两个文件被覆盖，只能靠存档输入**按构造重建**。
#       当时的补救之所以成立，靠的是「旧输入恰好被存过」这个运气，不是机制。
#
# 设计：**写侧助手**（本项目无统一跑批器，链是手工按序跑的）。
#   `backup(path, producer)` 在**覆盖之前**把现存版本复制到 `_artifact_archive/`，
#   文件名带 UTC 时间戳与内容 sha 前 8 位 ⇒ 每一版都可指认、可比对。
ARCHIVE = os.path.join(ROOT, "分析结果", "_artifact_archive")


def backup(path, producer):
    """覆盖前归档前一版。返回归档路径；文件不存在则返回 None（首次生产，无需归档）。"""
    if not os.path.exists(path):
        return None
    h = (_sha(path) or "")[:8]
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    os.makedirs(ARCHIVE, exist_ok=True)
    dst = os.path.join(ARCHIVE, "%s.%s.%s" % (os.path.basename(path), ts, h))
    shutil.copy2(path, dst)
    return dst
