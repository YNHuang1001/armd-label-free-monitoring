# -*- coding: utf-8 -*-
r"""t04_sudo.py — family code SUDO. Estimator validation suite, module 4.

PRODUCTION UNDER TEST (verbatim, read-only):
    L.BLOCK_SUDO = p6c_01_L3_blinded.py lines 291-330  ("L3-A: SUDO-AURCC")
    executed byte-for-byte via L.exec_block(L.BLOCK_SUDO, env, names=["rec"]).
Nothing outside estimator_validation/ is written.

DEFINITION THE BRUTE-FORCE REFERENCE IS WRITTEN FROM
----------------------------------------------------
Kiyasseh et al., Nat Commun 15, 2888 (2024), doi:10.1038/s41467-024-46000-9, the
SUDO / AURCC construction, as stated in the task definition handed to this module:

  * For EACH predicted-probability interval (bin) of the in-the-wild (target) points,
    the SAME sampled in-the-wild points are pseudo-labelled class 0 AND separately
    class 1 — i.e. TWO auxiliary branches per bin, drawn with the same random seed.
  * One auxiliary classifier is trained per branch, against ground-truth points of
    the OPPOSITE class (branch pseudo-label c  <->  ground-truth class 1-c).
  * SUDO is the DISCREPANCY BETWEEN THE TWO auxiliary classifiers' performances
    (a signed difference), NOT a single classifier's 2*(AUC-0.5).
  * AURCC sweeps symmetric low/high probability-threshold PAIRS (A,B); its
    reliability is the mean |SUDO| over thresholds, completeness the fraction of
    in-the-wild points covered.

INDEPENDENCE OF THE REFERENCE — stated honestly (an earlier version claimed the two
implementations "share no lines of code"; that was false): `binning()` below is a
verbatim copy of the production's frozen bin rule (p6c:292-293, only nbin/NBIN
renamed); `independent_valid_bins()` re-implements the production gates (p6c:298
len(idx)>=10, p6c:301 len(pool)>=10) and shares that copied bin rule, so the expected
denominator is derived from the same bin/gate logic as the code under test and
bin-rule or gate-threshold defects are OUTSIDE this module's scope; the reference's
per-round scaffolding (nb formula, vstack, label vector, LogisticRegression,
evaluation labels) is a renamed copy of p6c:304-315. What is genuinely independent is
the part under test: the TWO-branch loop and the cross-branch DISCREPANCY statistic,
which the production block does not have. The one interface convention deliberately
shared with production (disclosed, not hidden) is the evaluation rule "the auxiliary
classifier is scored on source-validation points with positive class = ground-truth
label == the branch pseudo-label" (quoted from the p6c frozen spec header); it is
applied identically to both branches of the reference.

PRODUCTION DEFECT EXPOSED (already diagnosed): the production block trains ONE
pseudo-label branch per bin (yb = 1 if bin mean pt >= 0.5 else 0), picks the negative
pool as the OPPOSITE class, and reports 2*(AUC-0.5) as the "signed SUDO". The
original method requires TWO branches and a DISCREPANCY between them.

EXPECTED RED-LIGHT PATTERN on the current (known-buggy) code:
  SUDO-01  FAIL  count 50 != 100 AND every training round uses a single pool class
                 (no complementary second branch): structural, not merely a count
  SUDO-01F PASS  structural count control on the two-branch reference (its rounds carry
                 complementary pools). It does NOT prove that this module's instrumentation
                 accepts a conformant production block — that is what the PATCHED run is for.
  SUDO-01P PASS  control: the proxies do not alter the production block's own numbers
  SUDO-02  PASS  characterisation assertion (warn_only=True) — the code DOES use one label
                 per bin, which is exactly the defect; green BY CONSTRUCTION, not a
                 correctness check, and NOT an acceptance criterion after a correct fix
                 (it degrades to WARN there by design)
  SUDO-03  FAIL  the per-bin signed statistic does not negate under the branch-role swap:
                 the single-branch block has no branch datasets to exchange (0 replayed pool
                 draws, 50 fresh), so the swap degenerates into a fresh-draw re-run and
                 |S + S_inverted| reaches ~1.96; guarded by a magnitude floor max|S| >= 0.25
                 so a degenerate zero-information statistic cannot pass vacuously
  SUDO-03C PASS  structural antisymmetry control (reference standalone: invert=True is a pure
                 relabelling of the reference's own two branches, exact by construction)
  SUDO-04  WARN  structural divergence of the AURCC aggregation (diagnostic only)
  SUDO-04X PASS  consistency guard: AURCC recomputed from per-bin attribution == block rec
  SUDO-05  FAIL  conformity verdict = NOT CONFORMANT (same structural criterion as SUDO-01)

ACCEPTANCE NOTE: SUDO-01's expected count 2*SUDO_REP*n_valid assumes no round is skipped;
if a conformant implementation skips rounds, the per-bin round counts reported in the
observed field show where the accounting differs.

DATA: synthetic only (no 202-task cache is touched).
OUTPUT: one RESULT|... line per assertion, via L.emit/L.check/L.guard.
The module always exits 0, even when assertions fail.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import evs_load as L
from sklearn.linear_model import LogisticRegression as RealLR
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score as real_auc

L.banner("t04_sudo — SUDO family (SUDO-AURCC, p6c lines 291-330)")

# ------------------------------------------------------------------ constants
SEED = 20261005
NBIN, SUDO_REP, SUDO_CAP, SUDO_FRAC = 10, 5, 200, 0.5
MIN_N = 10                       # production per-bin / per-pool floor (p6c:298,301)
TOL_FLIP = 1e-6                  # exact-antisymmetry tolerance for the sign test
MAG_FLOOR = 0.25                 # SUDO-03 power floor: max_b |S_b| must clear it, else the
                                 # antisymmetry test is vacuous (a zero statistic would pass)
BLOCK_REF = "p6c_01_L3_blinded.py:291-330"

# ------------------------------------------------------- synthetic data (only)
def make_data():
    g = np.random.default_rng(SEED)
    NF = 6
    Dtr = g.normal(size=(600, NF))
    Dva = g.normal(size=(300, NF))
    Dt = g.normal(size=(600, NF))
    w = g.normal(size=NF)
    # balanced by construction (median split) so both pseudo-label pools have equal size:
    # the inverted run's pool draws then have the same length as the main run's, which the
    # SUDO-03 branch-exchange replay requires (size-matched replay).
    ztr = Dtr @ w - 0.3
    ytr = (ztr > np.median(ztr)).astype(int)
    zva = Dva @ w
    yva = (zva > np.median(zva)).astype(int)
    sc = StandardScaler().fit(Dtr)
    mdl = RealLR(max_iter=2000, random_state=1).fit(sc.transform(Dtr), ytr)
    Ds_std = sc.transform(Dtr)
    Dt_std = sc.transform(Dt)
    pt = mdl.predict_proba(Dt_std)[:, 1]
    M = dict(sc=sc, ytr=ytr, yva=yva, Dva=Dva)
    return dict(NF=NF, Dtr=Dtr, Dva=Dva, Dt=Dt, ytr=ytr, yva=yva, sc=sc, mdl=mdl,
                Ds_std=Ds_std, Dt_std=Dt_std, pt=pt, M=M)


D = make_data()
PT, DT_STD, DS_STD = D["pt"], D["Dt_std"], D["Ds_std"]
YTR, YVA, SC, DVA = D["ytr"], D["yva"], D["sc"], D["Dva"]


def binning(pt, nbin=NBIN):
    """VERBATIM copy of the production's frozen bin rule (p6c:292-293; only nbin/NBIN
    renamed). Used ONLY to attribute observed training calls to bins; every assertion on
    the pseudo-label rule is made on instrumented observations, not on this helper. A
    defect in the bin rule itself would be invisible here (outside this module's scope)."""
    qs = np.quantile(pt, np.linspace(0, 1, nbin + 1))
    qs[0], qs[-1] = -np.inf, np.inf
    return np.digitize(pt, qs[1:-1])


BID = binning(PT)


def independent_valid_bins(pt, ytr, nbin=NBIN, min_n=MIN_N):
    """Bins that reach the training stage under the production's own gates
    (p6c:298 len(idx)>=10 and p6c:301 len(pool)>=10) with the 0.5 bin-mean rule.
    Denominator of the expected branch count. NOTE: this re-implements p6c:296-301 and
    shares the copied bin rule with production, so it is not independent of the code
    under test for gate/threshold defects (declared scope limitation)."""
    bid = binning(pt, nbin)
    out = []
    for b in range(nbin):
        idx = np.where(bid == b)[0]
        if len(idx) < min_n:
            continue
        yb = 1 if pt[idx].mean() >= 0.5 else 0
        pool = np.where(ytr == 1 - yb)[0]
        if len(pool) < min_n:
            continue
        out.append(dict(bin=b, yb=yb, n_idx=len(idx), n_pool=len(pool),
                        mean_pt=float(pt[idx].mean())))
    return out


# ------------------------------------------------------------- instrumentation
# Current training round = the most recent TARGET draw seen by RecordingRng. Every fit
# and every AUC call is stamped with it, so a round may contain k pool draws, k fits and
# k AUC calls (k = 2 for a conformant two-branch block; k = 1 for the buggy block).
ROUND = {"i": -1, "bin": None}


class CountingLR(object):
    """Faithful proxy for LogisticRegression: delegates every call to the real sklearn
    estimator and records (X, y, round) for every .fit(...).  Counts auxiliary classifiers
    actually trained and ties each one to the target-draw round it belongs to."""
    FIT = []

    def __init__(self, **kw):
        self._in = RealLR(**kw)
        self._kw = dict(kw)

    def fit(self, X, y, **kw):
        CountingLR.FIT.append((np.asarray(X), np.asarray(y), ROUND["i"]))
        self._in.fit(X, y, **kw)
        return self

    def decision_function(self, X):
        return self._in.decision_function(X)

    def predict_proba(self, X):
        return self._in.predict_proba(X)

    def predict(self, X):
        return self._in.predict(X)

    def get_params(self, deep=True):
        return self._in.get_params(deep)

    def __getattr__(self, name):
        return getattr(object.__getattribute__(self, "_in"), name)


class RecordingRng(object):
    """Proxy for np.random.Generator recording every .choice(...) call as a dict
    (target?, round index, bin, size, drawn values, replayed?).  Everything else is
    delegated.

    A draw is a TARGET draw when its population lies inside ONE score-quantile bin (the
    block samples a bin's target points from such a population); every other draw is a
    pool draw.  Each target draw opens a new round, which is how fits and AUC calls are
    attributed to the bin even when a round issues several pool draws / fits / AUC calls.

    REPLAY (used only for the SUDO-03 role-swap run): if `replay` maps round -> {ground-
    truth class: [drawn arrays in call order]} — recorded from the main run — then a pool
    draw whose population is a single ground-truth class is answered with the main run's
    recorded draw OF THAT CLASS (same round, same size).  The underlying generator is
    still advanced so the call schedule is unchanged.  This ENFORCES the role swap for
    any two-branch block (each branch receives the dataset the main run's opposite branch
    trained on) instead of hoping an independent re-draw reproduces it; a single-branch
    block requests only the complementary class, finds no recorded draw of it, and falls
    back to a genuine fresh draw (counted in n_fresh)."""
    def __init__(self, gen, replay=None):
        self._g = gen
        self.calls = []
        self.replay = replay or {}
        self.n_replay = 0
        self.n_fresh = 0

    def choice(self, a, size=None, replace=True, p=None):
        out = self._g.choice(a, size=size, replace=replace, p=p)
        pop = np.asarray(a).ravel()
        try:
            bins = np.unique(BID[pop]) if pop.size else np.array([], dtype=int)
            cls = np.unique(YTR[pop]) if pop.size else np.array([], dtype=int)
        except IndexError:
            bins, cls = np.array([], dtype=int), np.array([], dtype=int)
        is_target = bool(pop.size and bins.size == 1)
        if is_target:
            ROUND["i"] += 1
            ROUND["bin"] = int(bins[0])
        drawn = np.asarray(out).ravel().copy()
        replayed = False
        if (not is_target) and cls.size == 1 and size is not None:
            slot = self.replay.get(ROUND["i"], {}).get(int(cls[0]))
            if slot:
                cand = slot.pop(0)
                if cand.size == drawn.size:
                    drawn, replayed = cand.copy(), True
        if not is_target:
            if replayed:
                self.n_replay += 1
            else:
                self.n_fresh += 1
        self.calls.append(dict(target=is_target, round=ROUND["i"],
                               bin=(int(bins[0]) if is_target else None),
                               size=(None if size is None else int(size)),
                               drawn=drawn, replayed=replayed))
        return drawn if replayed else out

    def __getattr__(self, name):
        return getattr(object.__getattribute__(self, "_g"), name)


def build_replay(calls):
    """round -> {ground-truth class of the drawn pool: [drawn arrays, call order]} from
    the main run's recorded calls (used to enforce the branch exchange in SUDO-03)."""
    out = {}
    for c in calls:
        if c["target"] or c["drawn"].size == 0:
            continue
        try:
            cls = np.unique(YTR[c["drawn"]])
        except IndexError:
            continue
        if cls.size == 1:
            out.setdefault(c["round"], {}).setdefault(int(cls[0]), []).append(c["drawn"])
    return out


AUC_LOG = []          # (round index, y_true used for the evaluation, auc value)


def recording_auc(y_true, y_score, **kw):
    v = real_auc(y_true, y_score, **kw)
    AUC_LOG.append((ROUND["i"], np.asarray(y_true).copy(), float(v)))
    return v


def branch_of(y_true, yva):
    """Branch pseudo-label of an evaluation call: the class c whose ground-truth
    validation labels are exactly what the block scored (the disclosed shared interface
    convention).  None if neither class matches."""
    yt = np.asarray(y_true).astype(int)
    for c in (1, 0):
        if np.array_equal(yt, (np.asarray(yva) == c).astype(int)):
            return c
    return None


def match_pool(X, pools):
    """Index into `pools` of the pool draw whose rows form one half of the fit matrix X
    (exact content match), else None.  Used to attribute each fit to the pool it actually
    trained against; the caller applies an order fallback when content matching fails."""
    Xa = np.asarray(X)
    for j, (drawn, _cls) in enumerate(pools):
        rows = DS_STD[np.asarray(drawn)]
        n = rows.shape[0]
        if n == 0 or Xa.shape[0] < n:
            continue
        if np.array_equal(Xa[:n], rows) or np.array_equal(Xa[Xa.shape[0] - n:], rows):
            return j
    return None


def rows_idx(rows):
    """Indices of the DS_STD rows that are exact copies of `rows` (same order); None if
    any row is not found there."""
    rows = np.asarray(rows)
    if rows.ndim != 2 or rows.shape[0] == 0:
        return None
    eq = np.all(DS_STD[None, :, :] == rows[:, None, :], axis=2)      # (n_rows, n_src)
    if not eq.any(axis=1).all():
        return None
    return np.array([np.where(r)[0][0] for r in eq])


def rows_class(rows):
    """Ground-truth source-train class of a block of DS_STD rows under the GLOBAL labels
    (None if not found or spanning both classes).  Used by the SUDO-01F control."""
    idx = rows_idx(rows)
    if idx is None:
        return None
    cls = np.unique(YTR[idx])
    return int(cls[0]) if cls.size == 1 else None


def fit_pool_class(X, ytr):
    """Class of the source rows inside a fit matrix, by exact content match against DS_STD
    (tries both halves).  Fallback attribution for blocks that draw their pools without
    the injected rng (so no pool draw was recorded by RecordingRng)."""
    Xa = np.asarray(X)
    nb = Xa.shape[0] // 2
    for half in (Xa[nb:], Xa[:nb]):
        idx = rows_idx(half)
        if idx is not None:
            cls = np.unique(np.asarray(ytr)[idx])
            if cls.size == 1:
                return int(cls[0])
    return None


def run_production(M, LRclass, rng, aucfn):
    """Exec the verbatim production block once. Returns (rec, fit_log, auc_log, rng)."""
    CountingLR.FIT = []
    del AUC_LOG[:]
    ROUND["i"], ROUND["bin"] = -1, None      # round indices are run-local (replay matching)
    env = dict(np=np, LogisticRegression=LRclass, roc_auc_score=aucfn, M=M,
               Ds_tr_std=DS_STD, Dt_std=DT_STD, pt=PT, NBIN=NBIN, SUDO_REP=SUDO_REP,
               SUDO_CAP=SUDO_CAP, SUDO_FRAC=SUDO_FRAC, SEED=SEED, rng=rng, rec={})
    out, _src = L.exec_block(L.BLOCK_SUDO, env, names=["rec"])
    return out["rec"], list(CountingLR.FIT), list(AUC_LOG), rng


def attribute_per_bin(fit_log, auc_log, rng, ytr=None, yva=None):
    """Attribute every fit and every AUC call to the most recent TARGET draw (a round =
    one target draw plus the pool draws / fits / AUC calls that follow it; rounds and
    bins come from RecordingRng, stamps from CountingLR.fit / recording_auc).  This is
    exact for k pool draws / k fits / k AUC calls per round, so a conformant two-branch
    block (1 target draw + 2 pool draws + 2 fits + 2 AUC calls) is handled; the previous
    positional (target, pool) pairing is only used as an order fallback for the fit->pool
    match, never for the fit/AUC->round match.

    For each fit the pool class of the source points it actually trained against is
    recorded (content match against the round's pool draws, order fallback), so "the two
    fits of a round use complementary pool classes {0,1}" is a per-fit statement.

    `ytr` / `yva` are the run's own labels (pass 1-YTR / 1-YVA for the inverted run).
    Returns per_bin: bin -> dict(n_rounds, n_fits, round_classes, aucs [(branch, value)],
    round_auc, pool_classes(set[int]), mixed_pool(bool), labels_ok(bool = label vector is
    [1]*nb+[0]*nb), n_target, n_pool, two_branch(bool), signed(float)).
    """
    ytr = YTR if ytr is None else np.asarray(ytr)
    yva = YVA if yva is None else np.asarray(yva)
    rounds = {}
    for c in rng.calls:
        r = rounds.setdefault(c["round"], dict(bin=None, pools=[], n_target=0, n_pool=0))
        if c["target"]:
            r["bin"] = c["bin"]
            r["n_target"] = c["size"]
        else:
            cls = np.unique(ytr[c["drawn"]]) if c["drawn"].size else np.array([])
            r["pools"].append((c["drawn"], int(cls[0]) if cls.size == 1 else None))
            r["n_pool"] = c["size"]
    fits_by_round = {}
    for X, y, r in fit_log:
        fits_by_round.setdefault(r, []).append((np.asarray(X), np.asarray(y)))
    if not fits_by_round:
        raise RuntimeError("instrumentation not applicable: no auxiliary fit was recorded")
    if -1 in fits_by_round:
        raise RuntimeError("instrumentation not applicable: %d fit(s) recorded before any "
                           "target draw" % len(fits_by_round[-1]))
    for r, rinfo in rounds.items():
        if rinfo["bin"] is not None and r not in fits_by_round:
            raise RuntimeError("instrumentation not applicable: target draw in round %d is "
                               "followed by zero fits" % r)

    def _entry(b):
        return per_bin.setdefault(b, dict(n_rounds=0, n_fits=0, round_classes=[], aucs=[],
                                          round_auc={}, pool_classes=set(), mixed_pool=False,
                                          labels_ok=True, n_target=0, n_pool=0,
                                          two_branch=False, signed=float("nan")))

    per_bin = {}
    for r in sorted(fits_by_round):
        rinfo = rounds.get(r)
        if rinfo is None or rinfo["bin"] is None:
            raise RuntimeError("instrumentation not applicable: %d fit(s) in round %d cannot "
                               "be attributed to a target draw"
                               % (len(fits_by_round[r]), r))
        rr = _entry(rinfo["bin"])
        rr["n_rounds"] += 1
        rr["n_fits"] += len(fits_by_round[r])
        rr["n_target"] = rinfo["n_target"]
        rr["n_pool"] = rinfo["n_pool"]
        classes = []
        for k, (X, y) in enumerate(fits_by_round[r]):
            j = match_pool(X, rinfo["pools"])
            if j is None and k < len(rinfo["pools"]):
                j = k                     # order fallback: k-th fit <-> k-th pool draw
            cls = rinfo["pools"][j][1] if j is not None else None
            if cls is None:
                cls = fit_pool_class(X, ytr)   # content fallback (pools drawn without rng)
            if cls is None:
                rr["mixed_pool"] = True
            else:
                classes.append(cls)
                rr["pool_classes"].add(cls)
            nb = X.shape[0] // 2
            rr["labels_ok"] &= bool(np.array_equal(y, np.r_[np.ones(nb), np.zeros(nb)]))
        rr["round_classes"].append(frozenset(classes))
    for r, y_true, v in auc_log:
        rinfo = rounds.get(r)
        if rinfo is None or rinfo["bin"] is None:
            raise RuntimeError("instrumentation not applicable: AUC call in round %d cannot "
                               "be attributed to a target draw" % r)
        rr = _entry(rinfo["bin"])
        br = branch_of(y_true, yva)
        rr["aucs"].append((br, float(v)))
        rr["round_auc"].setdefault(r, {})[br] = float(v)
    for rr in per_bin.values():
        rr["two_branch"] = bool(rr["round_classes"]) and all(
            s == frozenset({0, 1}) for s in rr["round_classes"])
        if rr["two_branch"]:
            bad = [r for r, d in rr["round_auc"].items() if set(d) != {0, 1}]
            if bad:
                raise RuntimeError("instrumentation not applicable: round(s) %s have no AUC "
                                   "for both branch pseudo-labels (evaluation y_true does "
                                   "not match the frozen convention)" % sorted(bad)[:3])
            rr["signed"] = float(np.mean([d[1] - d[0] for d in rr["round_auc"].values()]))
        else:
            vals = [v for d in rr["round_auc"].values() for v in d.values()]
            rr["signed"] = (2.0 * (float(np.mean(vals)) - 0.5)) if vals else float("nan")
    return per_bin


# --------------------------------------------------------- reference (original)
def two_branch_sudo(LRclass, aucfn, invert=False):
    """Two-branch SUDO per the stated original definition (see module docstring).

    For every bin and repeat: ONE target sample ib is drawn and used by BOTH branches
    (the same sampled in-the-wild points, c=0 and c=1); each branch trains its own
    auxiliary classifier against source-train points of the OPPOSITE class (1-c) and is
    scored with positive class = ground-truth validation points of class c.  The per-bin
    statistic is the DISCREPANCY mean_r(AUC_{c=1}) - mean_r(AUC_{c=0}).

    invert=True applies the full role swap: pseudo-label rule inverted (c -> 1-c),
    opposite-class pool swapped, evaluation labels complemented.  Because both branches
    live in one run, this must be the exact negation of the non-inverted statistic.
    Returns (per_bin_signed: dict, fit_log: list).
    """
    src = {0: np.where(YTR == 0)[0], 1: np.where(YTR == 1)[0]}
    sv = SC.transform(DVA)
    CountingLR.FIT = []
    per_bin = {}
    for b in range(NBIN):
        idx = np.where(BID == b)[0]
        if len(idx) < MIN_N:
            continue
        nb = min(int(SUDO_FRAC * len(idx)), SUDO_CAP, len(idx), len(src[0]), len(src[1]))
        if nb < MIN_N:
            continue
        dd = []
        for r in range(SUDO_REP):
            # one target draw per (bin, repeat), shared by both branches
            ib = np.random.default_rng([SEED, b, r]).choice(idx, nb, replace=False)
            aucs = {}
            for c in (0, 1):
                cc = (1 - c) if invert else c
                # pool draw keyed by the PSEUDO-LABEL cc (not by the branch index): the
                # branch pseudo-labelled cc always draws the same source points, so
                # invert=True is a pure relabeling of the two branches.
                sb = np.random.default_rng([SEED, 991, b, r, cc]).choice(
                    src[1 - cc], nb, replace=False)
                Xa = np.vstack([DT_STD[ib], DS_STD[sb]])
                ya = np.r_[np.ones(nb), np.zeros(nb)]
                clf = LRclass(max_iter=2000, random_state=SEED + r).fit(Xa, ya)
                yv = (YVA == cc).astype(int)
                aucs[c] = aucfn(yv, clf.decision_function(sv))
            dd.append(aucs[1] - aucs[0])
        per_bin[b] = float(np.mean(dd))
    return per_bin, list(CountingLR.FIT)


# =============================================================== shared runs ==
RUNS = {}
RUN_ERR = {}


def _do_run(tag, fn):
    try:
        RUNS[tag] = fn()
    except Exception as e:
        RUN_ERR[tag] = "%s: %s" % (type(e).__name__, e)


_do_run("main", lambda: run_production(dict(D["M"]), CountingLR,
                                       RecordingRng(np.random.default_rng(SEED + 11)),
                                       recording_auc))
_do_run("plain", lambda: run_production(dict(D["M"]), RealLR,
                                        np.random.default_rng(SEED + 11), real_auc))


def _replay_map():
    """Main run's pool draws by round and ground-truth class, for SUDO-03's enforced swap."""
    return build_replay(RUNS["main"][3].calls) if "main" in RUNS else {}


_do_run("inverted", lambda: run_production(dict(sc=SC, ytr=1 - YTR, yva=1 - YVA, Dva=DVA),
                                           CountingLR,
                                           RecordingRng(np.random.default_rng(SEED + 11),
                                                        _replay_map()),
                                           recording_auc))
_do_run("ref", lambda: two_branch_sudo(CountingLR, recording_auc, invert=False))
_do_run("ref_inv", lambda: two_branch_sudo(CountingLR, recording_auc, invert=True))


# =================================================================== SUDO-01 ==
def sudo_conformity():
    """Shared SUDO-01 / SUDO-05 criterion, re-used verbatim by SUDO-05.

    Conformant iff (a) the number of auxiliary classifiers trained equals
    2 * SUDO_REP * n_valid_bins, (b) EVERY training round (one target draw) contains fits
    whose source pools carry the complementary ground-truth classes {0,1} — one per
    pseudo-label branch — and (c) every definition-valid bin receives SUDO_REP rounds.
    (a) alone is insufficient: a single-branch block that fits each classifier twice also
    reaches 2*SUDO_REP*n_valid_bins, which is why (b) is the decisive structural clause.
    Returns (ok, observed)."""
    if "main" in RUN_ERR:
        return False, "production run failed: " + RUN_ERR["main"]
    rec, fit_log, auc_log, rng = RUNS["main"]
    per_bin = attribute_per_bin(fit_log, auc_log, rng)
    valid = independent_valid_bins(PT, YTR)
    n_valid = len(valid)
    observed = len(fit_log)
    expected = 2 * SUDO_REP * n_valid
    counts = ",".join("b%d:%d" % (b, per_bin[b]["n_rounds"]) for b in sorted(per_bin))
    pats = {}
    for b in sorted(per_bin):
        for s in per_bin[b]["round_classes"]:
            pats[tuple(sorted(s))] = pats.get(tuple(sorted(s)), 0) + 1
    n_rounds = sum(per_bin[b]["n_rounds"] for b in per_bin)
    n_two = sum(1 for b in per_bin for s in per_bin[b]["round_classes"]
                if s == frozenset({0, 1}))
    valid_bins = sorted(v["bin"] for v in valid)
    rounds_per_bin_ok = (sorted(per_bin) == valid_bins and
                         all(per_bin[b]["n_rounds"] == SUDO_REP for b in per_bin))
    count_ok = (observed == expected)
    struct_ok = (n_rounds > 0 and n_two == n_rounds)
    ok = count_ok and struct_ok and rounds_per_bin_ok
    obs = ("aux_classifiers_fit=%d; expected=2*%d_reps*%d_valid_bins=%d (count_ok=%s); "
           "per-bin rounds=[%s]; bins_trained=%d; every definition-valid bin has %d rounds=%s; "
           "rounds with complementary pool classes {0,1}=%d/%d (fit-pool-class patterns: %s); "
           "code rec SUDO_nbin=%s" %
           (observed, SUDO_REP, n_valid, expected, count_ok, counts, len(per_bin), SUDO_REP,
            rounds_per_bin_ok, n_two, n_rounds,
            ";".join("%s:%d" % (list(p), n) for p, n in
                     sorted(pats.items(), key=lambda kv: -kv[1])[:4]),
            rec.get("SUDO_nbin")))
    return ok, obs


L.guard("SUDO-01", sudo_conformity,
        "two-branch SUDO (Kiyasseh 2024 doi:10.1038/s41467-024-46000-9): EVERY training round "
        "(one target draw, per repeat) trains TWO auxiliary classifiers whose source pools "
        "have the complementary ground-truth classes {0,1} — one per pseudo-label branch — so "
        "the number of classifiers is 2 * SUDO_REP * n_valid_bins; n_valid_bins from the "
        "definition (bin n>=10 and opposite-class pool n>=10)",
        "per-fit instrumentation of p6c:303-310 with a CountingLR/RecordingRng proxy around "
        "the verbatim block; the complementary-pool clause is what a double-fit single-branch "
        "block cannot satisfy")

# SUDO-01F: the same structural criterion must be met by a genuine two-branch implementation
# — evidence that the criterion is satisfiable (it is a structural count control).
def t01f():
    if "ref" in RUN_ERR:
        return False, "reference run failed: " + RUN_ERR["ref"]
    ref, fit_log = RUNS["ref"]
    expected = 2 * SUDO_REP * len(ref)
    bad, pats = 0, {}
    for k in range(0, len(fit_log) - len(fit_log) % 2, 2):
        cs = []
        for X, y, _r in fit_log[k:k + 2]:
            nb = len(y) // 2
            c = rows_class(np.asarray(X)[nb:])          # reference stacks pool LAST
            if c is None:
                c = rows_class(np.asarray(X)[:nb])
            cs.append(c)
        key = tuple(sorted(-1 if c is None else c for c in cs))
        pats[key] = pats.get(key, 0) + 1
        if set(cs) != {0, 1}:
            bad += 1
    return ((len(fit_log) == expected) and bad == 0), (
        "structural count control on the two-branch reference: aux_classifiers_fit=%d; "
        "expected=2*%d_reps*%d_bins=%d; rounds (consecutive fit pairs)=%d; rounds whose two "
        "fits trained against complementary source classes=%d; patterns=%s"
        % (len(fit_log), SUDO_REP, len(ref), expected, len(fit_log) // 2,
           len(fit_log) // 2 - bad, pats))


L.guard("SUDO-01F", t01f,
        "structural count control: the two-branch reference raises the same counter to "
        "2*SUDO_REP*n_bins AND every one of its rounds carries complementary pool classes, so "
        "SUDO-01's structural criterion is satisfiable. This does NOT prove that the module's "
        "instrumentation accepts a conformant production block — the PATCHED (EVS_PROD_DIR) "
        "run is the proof of that",
        "independent two-branch reference written from the Nat Commun 2024 definition, run "
        "through the SAME CountingLR proxy; pool classes read by exact content match of each "
        "fit's negative half against the source-train matrix")

# SUDO-01P: proxy faithfulness — the instrumented run must reproduce the uninstrumented
# sklearn run exactly.
def t01p():
    if "main" in RUN_ERR or "plain" in RUN_ERR:
        return False, "run failed: %s %s" % (RUN_ERR.get("main"), RUN_ERR.get("plain"))
    r1 = RUNS["main"][0]
    r2 = RUNS["plain"][0]
    keys = ["SUDO_AURCC", "SUDO_nbin", "SUDO_signed_min", "SUDO_signed_max", "SUDO_absmean"]
    same, diffs = True, []
    for k in keys:
        if k not in r1 or k not in r2:
            same = False
            diffs.append("%s MISSING (proxy=%s plain=%s)" % (k, k in r1, k in r2))
            continue
        a, b = r1.get(k), r2.get(k)
        if a is None or b is None:
            same = False
            diffs.append("%s None (proxy=%r plain=%r)" % (k, a, b))
            continue
        if isinstance(a, float) and isinstance(b, float):
            if not (np.isfinite(a) and np.isfinite(b) and abs(a - b) <= 1e-12):
                same = False
                diffs.append("%s %.12g!=%.12g" % (k, a, b))
        elif a != b:
            same = False
            diffs.append("%s %r!=%r" % (k, a, b))
    obs = ("proxy rec == real-sklearn rec (same seed/data): %s; AURCC proxy=%.10f plain=%.10f"
           % (same, r1.get("SUDO_AURCC", np.nan), r2.get("SUDO_AURCC", np.nan)))
    if diffs:
        obs += "; DIFFS " + ";".join(diffs)
    return same, obs


L.guard("SUDO-01P", t01p,
        "instrumentation faithfulness: identical rec values with and without the "
        "CountingLR/RecordingRng proxies (same seed, same data) — the proxy must not alter "
        "what the production block computes",
        "main run (proxies) vs plain run (real LogisticRegression, real Generator)")


# =================================================================== SUDO-02 ==
def t02():
    if "main" in RUN_ERR:
        return False, "production run failed: " + RUN_ERR["main"]
    rec, fit_log, auc_log, rng = RUNS["main"]
    per_bin = attribute_per_bin(fit_log, auc_log, rng)
    valid = {v["bin"]: v for v in independent_valid_bins(PT, YTR)}
    rule_ok, rule_txt = True, []
    for b in sorted(per_bin):
        yb_rule = valid[b]["yb"] if b in valid else None
        pool_obs = sorted(per_bin[b]["pool_classes"])
        pool_expected = [1 - yb_rule] if yb_rule is not None else []
        ok_b = (pool_obs == pool_expected)
        rule_ok &= ok_b
        rule_txt.append("b%d(mean_pt=%.3f->yb=%s,obs_pool_class=%s)%s" %
                        (b, valid[b]["mean_pt"] if b in valid else float("nan"),
                         yb_rule, pool_obs, "" if ok_b else "!MISMATCH"))
    one_label = all((len(per_bin[b]["pool_classes"]) == 1) and (not per_bin[b]["mixed_pool"])
                    for b in per_bin)
    target_always_pos = all(per_bin[b]["labels_ok"] for b in per_bin)
    ok = rule_ok and one_label and target_always_pos
    obs = ("threshold rule yb=1[bin_mean_pt>=0.5] with pool=(1-yb): matches for %d/%d bins "
           "[%s]; exactly one pool class per bin (no complementary second branch)=%s; every "
           "fit label vector == [1]*nb+[0]*nb (sampled target points always positive, never "
           "relabelled 0)=%s" %
           (sum(1 for t in rule_txt if "!MISMATCH" not in t), len(per_bin),
            ";".join(rule_txt), one_label, target_always_pos))
    return ok, obs


# CHARACTERISATION ONLY (warn_only=True): on the buggy code this line is GREEN by
# construction (it asserts the observed fact — one label per bin — which IS the defect);
# on a correctly fixed block it degrades to WARN instead of emitting a false FAIL. It is
# deliberately not an acceptance criterion and must never be read as evidence of
# conformance either way.
try:
    _ok02, _obs02 = t02()
except Exception as _e:                      # noqa: BLE001 - mirror L.guard behaviour
    _ok02, _obs02 = False, "EXCEPTION %s: %s" % (type(_e).__name__, _e)

L.check("SUDO-02", _ok02,
        "characterisation (current-state observation; NOT a correctness bound): p6c frozen "
        "spec header line 31 and code p6c:299-300 state ONE label per bin, yb = 1 if "
        "pt[idx].mean() >= 0.5 else 0, pool = src_by_class[1-yb], no complementary branch. "
        "The two-branch original (Kiyasseh 2024) is the contrast that makes this a defect. "
        "GREEN on the current buggy code by design; WARN on a conformant block",
        _obs02,
        "pseudo-label rule recovered from instrumented rng.choice populations + fit label "
        "vectors; warn_only=True (characterisation line, excluded from acceptance)",
        warn_only=True)


# =================================================================== SUDO-03 ==
def t03():
    if "main" in RUN_ERR or "inverted" in RUN_ERR:
        return False, "run failed: %s %s" % (RUN_ERR.get("main"), RUN_ERR.get("inverted"))
    pb1 = attribute_per_bin(RUNS["main"][1], RUNS["main"][2], RUNS["main"][3], YTR, YVA)
    pb2 = attribute_per_bin(RUNS["inverted"][1], RUNS["inverted"][2], RUNS["inverted"][3],
                            1 - YTR, 1 - YVA)
    rg = RUNS["inverted"][3]
    n_poolcalls = sum(1 for c in rg.calls if not c["target"])
    if n_poolcalls == 0:
        swapped = ("branch datasets NOT exchangeable: no pool draw was issued through the "
                   "injected rng, so the role swap could not be enforced and the two runs "
                   "sampled independently (exact negation is then not expected)")
    else:
        swapped = ("branch datasets exchanged: pool draws replayed from the main run of the "
                   "requested ground-truth class=%d, fresh draws (no recorded draw of the "
                   "requested class in that round - single-branch signature)=%d"
                   % (getattr(rg, "n_replay", 0), getattr(rg, "n_fresh", 0)))
    shared = sorted(set(pb1) & set(pb2))
    if not shared:
        return False, "no shared bins between the two runs"
    devs, same_sign, txt = [], 0, []
    for b in shared:
        s1, s2 = pb1[b]["signed"], pb2[b]["signed"]
        devs.append(abs(s1 + s2))
        if np.sign(s1) == np.sign(s2):
            same_sign += 1
        txt.append("b%d:%+.3f/%+.3f" % (b, s1, s2))
    maxdev = float(np.max(devs))
    maxmag = float(np.max([abs(pb1[b]["signed"]) for b in shared]))
    if not np.isfinite(maxmag) or maxmag < MAG_FLOOR:
        return False, ("degenerate statistic - antisymmetry test has no power: max|S| = %s < "
                       "%g floor (a zero-information statistic makes S' = -S vacuously true); "
                       "bins=%d; per-bin S=[%s]"
                       % ("%.6f" % maxmag if np.isfinite(maxmag) else str(maxmag),
                          MAG_FLOOR, len(shared), ";".join(txt)))
    ok = maxdev <= TOL_FLIP
    return ok, ("max|S + S_inverted| = %.6f over %d bins (tol %g); max|S| = %.3f (power floor "
                "%g); bins whose SIGN did NOT flip = %d/%d; %s; per-bin S(as-is)/S(inverted)"
                "=[%s]" %
                (maxdev, len(shared), TOL_FLIP, maxmag, MAG_FLOOR, same_sign, len(shared),
                 swapped, ";".join(txt)))


L.guard("SUDO-03", t03,
        "role-swap antisymmetry of a two-branch DISCREPANCY statistic: swapping the two "
        "branch roles (each branch receives the OTHER branch's sampled auxiliary dataset - "
        "same target draw, complementary-class pool - and the evaluation labels complement) "
        "must negate the per-bin signed statistic exactly (S' = -S); guarded by the power "
        "floor max|S| >= %g so a degenerate zero-information statistic cannot pass vacuously. "
        "A single-branch block has no second branch dataset to exchange: every pool draw "
        "falls back to a fresh draw and its statistic does not negate" % MAG_FLOOR,
        "production run vs production run on M with ytr->1-ytr and yva->1-yva, fresh rng with "
        "the identical seed (evaluation labels complemented AND the source pools actually "
        "swapped); the swap is enforced by replaying the main run's recorded pool draws of "
        "the requested ground-truth class, so the identity holds for any two-branch block "
        "regardless of its random-draw schedule; the per-bin statistic is the branch-aware "
        "discrepancy when the round carries complementary branches and the block's own "
        "2*(AUC-0.5) when it does not")


# SUDO-03C: control — the same swap on the genuine two-branch reference is exact.
def t03c():
    if "ref" in RUN_ERR or "ref_inv" in RUN_ERR:
        return False, "reference run failed: %s %s" % (RUN_ERR.get("ref"), RUN_ERR.get("ref_inv"))
    d1, d2 = RUNS["ref"][0], RUNS["ref_inv"][0]
    shared = sorted(set(d1) & set(d2))
    if not shared:
        return False, "no shared bins in the reference runs"
    maxdev = float(np.max([abs(d1[b] + d2[b]) for b in shared]))
    return (maxdev <= TOL_FLIP), (
        "two-branch reference max|D + D_inverted| = %.3e over %d bins; sample D values: %s" %
        (maxdev, len(shared), ";".join("b%d:%+.3f" % (b, d1[b]) for b in shared[:5])))


L.guard("SUDO-03C", t03c,
        "structural antisymmetry control (reference standalone): on the genuine two-branch "
        "reference the role swap negates its per-bin statistic, so SUDO-03's criterion is "
        "satisfiable; exactness here is by construction (invert=True is a pure relabelling of "
        "the reference's own two branches, which key the pool draw on the pseudo-label), and "
        "this control does NOT exercise the module's instrumentation path",
        "reference run invert=False vs invert=True (identical sampled points, same seeds)")


# =================================================================== SUDO-04 ==
def t04():
    rec, fit_log, auc_log, rng = RUNS["main"]
    per_bin = attribute_per_bin(fit_log, auc_log, rng)
    valid = {v["bin"]: v for v in independent_valid_bins(PT, YTR)}
    # production construction: sort bins by |bin mean - 0.5| descending and accumulate
    # (the block additionally drops bins with <3 valid repeats, p6c:316 — mirrored here)
    rows = [(b, valid[b]["n_idx"] / len(PT), abs(valid[b]["mean_pt"] - 0.5),
             abs(per_bin[b]["signed"])) for b in sorted(per_bin)
            if b in valid and len(per_bin[b]["aucs"]) >= 3]
    rows = sorted(rows, key=lambda r: -r[2])
    fr = np.array([r[1] for r in rows])
    rb = np.array([r[3] for r in rows])
    comp = np.cumsum(fr) / np.sum(fr)
    rel = np.cumsum(rb * fr) / np.cumsum(fr)
    aurcc_prod = float(np.trapezoid(rel, comp))
    # original construction: sweep symmetric threshold pairs (A,B) = (0.5-d, 0.5+d)
    ds = np.linspace(0.0, 0.5, 26)
    cgrid, rgrid = [], []
    for d in ds:
        A, B = 0.5 - d, 0.5 + d
        covered = np.where((PT >= A) & (PT <= B))[0]
        cgrid.append(len(covered) / len(PT))
        vals = [abs(per_bin[b]["signed"]) for b in per_bin
                if b in valid and A <= valid[b]["mean_pt"] <= B]
        rgrid.append(float(np.mean(vals)) if vals else 0.0)
    cgrid, rgrid = np.array(cgrid), np.array(rgrid)
    order = np.argsort(cgrid, kind="mergesort")
    aurcc_sweep = float(np.trapezoid(rgrid[order], cgrid[order]))
    rec_aurcc = rec.get("SUDO_AURCC")
    got = float(rec_aurcc) if rec_aurcc is not None else float("nan")
    L.check("SUDO-04X", (rec_aurcc is not None) and abs(aurcc_prod - got) <= 1e-12,
            "instrumentation consistency (exact identity): the AURCC recomputed from this "
            "module's own per-bin attribution of the verbatim block under the block's frozen "
            "aggregation (bins ordered by |bin mean-0.5| DESC; completeness=cumsum(frac)/"
            "sum(frac); reliability=cumsum(|S|*frac)/cumsum(frac); trapezoid) must equal the "
            "block's own rec['SUDO_AURCC']; expectation = exact float identity, no threshold",
            "recomputed=%.17g; block rec=%.17g; |diff|=%.3g; bins_attributed=%d" %
            (aurcc_prod, got, abs(aurcc_prod - got), len(rows)),
            "free regression guard on the AURCC aggregation: verified bit-identical on the "
            "current code (0.89029877954144621); the per-bin statistics are measured by the "
            "same instrumentation SUDO-01/SUDO-03 use")
    obs = ("production: bins sorted by |binmean-0.5| DESC, completeness=cumsum(frac)/sum(frac), "
           "reliability=cumsum(|S|*frac)/cumsum(frac), AURCC_prod=%.4f (block rec SUDO_AURCC=%.4f, "
           "nbin=%s); original: sweep of %d symmetric threshold pairs (0.5-d,0.5+d), "
           "completeness=frac(pt in [A,B]), mean|SUDO| over covered bins, AURCC_sweep=%.4f" %
           (aurcc_prod, rec_aurcc if rec_aurcc is not None else float("nan"),
            rec.get("SUDO_nbin"), len(ds), aurcc_sweep))
    L.emit("SUDO-04", "WARN",
           "original AURCC: symmetric probability-threshold-pair sweep (A,B); completeness = "
           "fraction of in-the-wild points covered; reliability = mean |SUDO| over thresholds "
           "(Kiyasseh 2024). Production: one curated ordering of bins by |bin mean-0.5|",
           obs,
           "STRUCTURAL DIVERGENCE (diagnostic, not a numeric bound): production ranks bins by "
           "distance of the bin mean to 0.5; the original sweeps symmetric threshold pairs. "
           "The two curves are different objects; sweep value is this module's reading of the "
           "stated definition, so this line is WARN-only")


try:
    t04()
except Exception as e:
    L.check("SUDO-04X", False,
            "instrumentation consistency (exact identity): recomputed AURCC must equal the "
            "block's own rec['SUDO_AURCC']; expectation = exact float identity",
            "EXCEPTION %s: %s" % (type(e).__name__, e),
            "consistency guard could not be computed")
    L.emit("SUDO-04", "WARN",
           "AURCC aggregation construction (structural comparison, diagnostic)",
           "EXCEPTION %s: %s" % (type(e).__name__, e),
           "comparison could not be computed; recorded as WARN (no numeric bound is asserted)")


# =================================================================== SUDO-05 ==
def t05():
    ok, det = sudo_conformity()          # same structural criterion as SUDO-01, not a re-count
    obs = ("CONFORMANT" if ok else "NOT CONFORMANT") + (" — " + det) + (
        "; clause (a) alone (aux classifier count) is satisfiable by a single-branch block "
        "that fits each classifier twice — the decisive clause is that every round's fits "
        "span both complementary pool classes {0,1}, which the current block cannot satisfy")
    return ok, obs


L.guard("SUDO-05", t05,
        "original two-branch SUDO definition: two pseudo-label branches per interval, SUDO = "
        "discrepancy between the two auxiliary classifiers' performances (Kiyasseh et al., "
        "Nat Commun 2024, doi:10.1038/s41467-024-46000-9); conformant = count 2*SUDO_REP*"
        "n_valid_bins AND every round's fits carry complementary pool classes {0,1} AND every "
        "definition-valid bin receives SUDO_REP rounds",
        "conformity verdict on the verbatim block " + BLOCK_REF + " using the SAME structural "
        "criterion as SUDO-01 (count + complementary pool classes + rounds per valid bin); "
        "expected NOT CONFORMANT on the current single-branch code")

print("\n[t04_sudo done]")
