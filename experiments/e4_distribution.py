"""E4 — distributional incidence, by income quintile and by wealth decile.

The paper's stated contribution over Bi–Foerster–Traum is heterogeneous
households, and until now nothing in the model exercised that margin: the only
distributional statistics were Gini coefficients, which docs/STATE.md (DIST-1)
flags as specifically the wrong object for the Greek crisis — measured inequality
barely moved, because Greece was already highly unequal pre-crisis and the
worst-affected households dropped out of the surveys. No Gini is computed here.

TWO CUTS, AND ONLY ONE OF THEM IS BEHAVIOURAL.

*Income quintiles* bin on the exogenous productivity state. The marginal
distribution over that state is the stationary distribution of the exogenous
Markov chain, so the mass in each bin is INVARIANT to the shock: households move
between income states, but the share in each state does not move. The per-capita
consumption response of a quintile is therefore a clean behavioural object. **This
is the cut to report.**

*Wealth deciles* bin on the steady-state deposit distribution with the boundaries
then held fixed. Membership churns hard as the whole deposit distribution shifts
across those fixed thresholds, and the net per-capita number turns out to be a
small residue of two large, nearly-cancelling terms — measured for the bottom
decile: −41.6 (consumption) against −44.4 (mass), netting +2.8. The arithmetic is
exact and the object is well defined, but it is overwhelmingly composition, not
behaviour, and must not be described as "how poor households responded". Retained
for completeness and reported with that caveat attached.

Both are defined on the SOLVED steady state, so this is a two-pass procedure:
solve the SS, cut the distribution, bake the masks into the hetoutputs, rebuild
the model with the augmented household block, then solve the Jacobian.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (os.path.join(ROOT, "code"), os.path.join(ROOT, "diagnostics", "regimes")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

CACHE = os.path.join(HERE, "cache_e4_deciles.npz")
N_DEC = 10


def decile_masks(D_ss, dep_pol):
    """Indicator arrays over the (z, dep) grid for each SS wealth decile.

    D_ss is the SS distribution, dep_pol the SS deposit policy — both shaped
    (nZ, nDep). Households are ranked by deposits and cut at population deciles,
    so each bin holds ~10% of households by construction (exactly 10% is not
    attainable on a discrete grid; the realised masses are returned and used as
    the denominator rather than assumed).
    """
    flat_dep = dep_pol.ravel()
    flat_D = D_ss.ravel()
    order = np.argsort(flat_dep, kind="stable")
    cum = np.cumsum(flat_D[order]) / flat_D.sum()

    masks, edges = [], []
    for k in range(N_DEC):
        lo, hi = k / N_DEC, (k + 1) / N_DEC
        sel = (cum > lo - 1e-12) & (cum <= hi + 1e-12) if k else (cum <= hi + 1e-12)
        m = np.zeros_like(flat_dep)
        m[order[sel]] = 1.0
        masks.append(m.reshape(dep_pol.shape))
        # Upper deposit edge of this bin, for the figure's axis annotation.
        edges.append(float(flat_dep[order[sel]].max()) if sel.any() else np.nan)
    return masks, np.array(edges)


# Masks are injected here before the household block is built. A module global is
# used because SSJ parses a hetoutput's SOURCE with inspect.getsource to discover
# its output names — an exec-generated or functools-wrapped closure has no source
# and fails with "could not get source code". The function below must therefore be
# written out literally, with its return statement naming every output.
_MASKS_D = None


def decile_out_D(c_D):
    """Mass-weighted consumption and mass, per wealth decile AND income quintile.

    SSJ aggregates a hetoutput as sum(D * output), so `c_D * mask_k` aggregates to
    the total consumption OF bin k and `mask_k` alone to its mass. Emitting both
    rather than a pre-divided per-capita number keeps the composition and
    behavioural parts separable downstream.
    """
    m = _MASKS_D
    cdec01_D = c_D * m[0]
    cdec02_D = c_D * m[1]
    cdec03_D = c_D * m[2]
    cdec04_D = c_D * m[3]
    cdec05_D = c_D * m[4]
    cdec06_D = c_D * m[5]
    cdec07_D = c_D * m[6]
    cdec08_D = c_D * m[7]
    cdec09_D = c_D * m[8]
    cdec10_D = c_D * m[9]
    o = np.ones_like(c_D)
    mdec01_D = o * m[0]
    mdec02_D = o * m[1]
    mdec03_D = o * m[2]
    mdec04_D = o * m[3]
    mdec05_D = o * m[4]
    mdec06_D = o * m[5]
    mdec07_D = o * m[6]
    mdec08_D = o * m[7]
    mdec09_D = o * m[8]
    mdec10_D = o * m[9]
    # ONE line, NO wrapping parentheses, and NO trailing comment. SSJ discovers
    # output names by splitting this statement's raw source text on commas, so a
    # parenthesised tuple yields "(cdec01_D" and a trailing "# noqa" corrupts the
    # LAST name (both were hit while building this). Keep the line bare.
    q = _QMASKS_D
    cqnt1_D = c_D * q[0]
    cqnt2_D = c_D * q[1]
    cqnt3_D = c_D * q[2]
    cqnt4_D = c_D * q[3]
    cqnt5_D = c_D * q[4]
    mqnt1_D = o * q[0]
    mqnt2_D = o * q[1]
    mqnt3_D = o * q[2]
    mqnt4_D = o * q[3]
    mqnt5_D = o * q[4]
    # ONE line, NO wrapping parentheses, NO trailing comment. SSJ discovers output
    # names by splitting this statement's raw source on commas, so a parenthesised
    # tuple yields "(cdec01_D" and a trailing "# noqa" corrupts the LAST name. Both
    # cuts live in ONE function because chaining two hetoutputs drops the first
    # one's outputs from the steady state (KeyError 'mdec01_D' in _steady_state).
    return cdec01_D, cdec02_D, cdec03_D, cdec04_D, cdec05_D, cdec06_D, cdec07_D, cdec08_D, cdec09_D, cdec10_D, mdec01_D, mdec02_D, mdec03_D, mdec04_D, mdec05_D, mdec06_D, mdec07_D, mdec08_D, mdec09_D, mdec10_D, cqnt1_D, cqnt2_D, cqnt3_D, cqnt4_D, cqnt5_D, mqnt1_D, mqnt2_D, mqnt3_D, mqnt4_D, mqnt5_D


_QMASKS_D = None
N_QNT = 5


def income_masks(D_ss, nZ):
    """Exact-quintile weights over the exogenous income states.

    FRACTIONAL, not 0/1. There are only 15 Rouwenhorst states and their stationary
    distribution is bell-shaped, so whole states cannot be partitioned into
    equal-mass quintiles — an all-or-nothing assignment gives masses like
    [0.090, 0.306, 0.000, 0.393, 0.212], with one bin empty. A state straddling a
    quintile boundary is therefore split, with the overlapping mass fraction going
    to each side.

    This is exact rather than an approximation *for an income cut*: every household
    in a given state has the same income, so any fraction of that state is a
    representative sample of it. The weight is applied uniformly across the deposit
    dimension, which is what makes that true. The same trick would NOT be
    legitimate for a wealth cut, where households within a bin differ in the very
    variable being binned on.
    """
    pz = D_ss.sum(axis=1)
    pz = pz / pz.sum()
    cum = np.concatenate([[0.0], np.cumsum(pz)])
    masks = []
    for k in range(N_QNT):
        lo, hi = k / N_QNT, (k + 1) / N_QNT
        w = np.zeros(len(pz))
        for j in range(len(pz)):
            overlap = max(0.0, min(hi, cum[j + 1]) - max(lo, cum[j]))
            w[j] = overlap / pz[j] if pz[j] > 1e-15 else 0.0
        masks.append(np.repeat(w[:, None], D_ss.shape[1], axis=1))
    return masks


INCOME_NAMES_D = ([f"cqnt{k + 1}_D" for k in range(N_QNT)]
                  + [f"mqnt{k + 1}_D" for k in range(N_QNT)])
INCOME_AGG_D = [n.upper() for n in INCOME_NAMES_D]
CQNT_AGG = INCOME_AGG_D[:N_QNT]
MQNT_AGG = INCOME_AGG_D[N_QNT:]


# Per-gridpoint hetoutput names (what decile_out_D returns; these live in the het
# block's internals).
DECILE_NAMES_D = ([f"cdec{k + 1:02d}_D" for k in range(N_DEC)]
                  + [f"mdec{k + 1:02d}_D" for k in range(N_DEC)])
# SSJ UPPERCASES a het block's output names to form the aggregate it exposes to the
# rest of the model (c_D -> C_D, dep_D -> DEP_D). These are the names that appear in
# the steady state's toplevel and in G.outputs, and therefore in the cache.
DECILE_AGG_D = [n.upper() for n in DECILE_NAMES_D]
CDEC_AGG = DECILE_AGG_D[:N_DEC]
MDEC_AGG = DECILE_AGG_D[N_DEC:]


def build():
    """Solve the SS, cut deciles, re-solve the Jacobian with decile outputs."""
    import sequence_jacobian as sj
    import tpi
    from calibration import get_calibration
    from depreciation_calibration import calibrate_depreciation
    from full_model import build_and_solve, solve_jacobian_padded
    from ic_delta_calibration import calibrate_ic_delta
    from regime_model import _ss_tpi, build_tpi_model_main
    from steady_state import solve_steady_state

    cal = get_calibration()
    ssr = calibrate_depreciation(calibrate_ic_delta(solve_steady_state(cal)))
    res = build_and_solve(ssr)
    ss = res["ss_final"]

    inner = ss.internals["hh_D"]
    D_ss, dep_ss, c_ss = inner["D"], inner["dep_D"], inner["c_D"]
    masks, dep_edges = decile_masks(D_ss, dep_ss)

    realised = np.array([float((D_ss * m).sum()) for m in masks])
    c_bin = np.array([float((D_ss * m * c_ss).sum()) for m in masks])
    print("SS decile masses:", np.round(realised, 4))
    assert abs(realised.sum() - 1.0) < 1e-8, f"decile masses sum to {realised.sum()}"
    assert realised.min() > 0.01, f"an SS decile bin is nearly empty: {realised}"

    global _MASKS_D, _QMASKS_D
    _MASKS_D = masks
    _QMASKS_D = income_masks(D_ss, D_ss.shape[0])
    q_mass = np.array([float((D_ss * m).sum()) for m in _QMASKS_D])
    q_c = np.array([float((D_ss * m * c_ss).sum()) for m in _QMASKS_D])
    print("SS income-quintile masses:", np.round(q_mass, 4))
    assert np.allclose(q_mass, 1.0 / N_QNT, atol=1e-9), (
        f"income quintiles are not equal-mass: {q_mass}. The fractional-weight "
        f"construction in income_masks should make each exactly 0.2.")
    hh_D_aug = tpi.hh_extended_D.add_hetoutputs([decile_out_D])

    # The model's steady state was solved with the PLAIN household block, so it has
    # no decile arrays — and a het block's Jacobian looks up ss[k] for every one of
    # its outputs (KeyError 'cdec01_D' otherwise). Re-evaluate the augmented block
    # at the solved aggregates and merge its internals back in. The masks are pure
    # functions of the SS policy and distribution, so this cannot move the SS; it
    # only materialises the extra outputs. Asserted below.
    hh_ss = hh_D_aug.steady_state(ss.toplevel)
    for k, v in hh_ss.internals["hh_D"].items():
        ss.internals["hh_D"][k] = v
    for nm in DECILE_AGG_D + INCOME_AGG_D:
        ss.toplevel[nm] = float(hh_ss[nm])
    # The decile bins partition the state space, so their consumption must sum to
    # aggregate consumption. This is the check that the masks are a partition and
    # that the hetoutput did not perturb the steady state.
    c_agg = float((D_ss * c_ss).sum())
    agg = sum(float(ss.toplevel[n]) for n in CDEC_AGG)
    assert abs(agg - c_agg) < 1e-9, \
        f"decile consumption sums to {agg}, not aggregate {c_agg} — masks are not a partition"
    m_agg = sum(float(ss.toplevel[n]) for n in MDEC_AGG)
    assert abs(m_agg - 1.0) < 1e-9, f"decile masses sum to {m_agg}, not 1"
    q_agg = sum(float(ss.toplevel[n]) for n in CQNT_AGG)
    assert abs(q_agg - c_agg) < 1e-9, \
        f"income-quintile consumption sums to {q_agg}, not aggregate {c_agg}"

    model = build_tpi_model_main(tpi, res["financial_solved_D"],
                                 res["financial_solved_F"], hh_D=hh_D_aug)
    ss_tpi = _ss_tpi(ss, float(cal["kappa_cb_F"]))
    T = res["T"]
    print(f"Solving Jacobian with {len(DECILE_AGG_D)} extra decile outputs (T={T}) ...",
          flush=True)
    G = solve_jacobian_padded(model, ss_tpi, res["unknowns_tp"],
                              res["targets_tp"],
                              ["Z_D", "shock_def_D", "Z_F", "shock_def_F",
                               "cb_buy_D"], T)

    out = {"T": np.array(T), "dShock_def_D": np.asarray(res["dShock_def_D"]),
           "decile_mass_ss": realised, "decile_c_ss": c_bin, "dep_edges": dep_edges,
           "qnt_mass_ss": q_mass, "qnt_c_ss": q_c}
    for o in ["spread_rb", "Y_D", "C_D", "n_inter_D"] + DECILE_AGG_D + INCOME_AGG_D:
        for i in ("shock_def_D", "cb_buy_D"):
            try:
                out[f"{o}__{i}"] = np.asarray(G[o][i], dtype=float)
            except KeyError:
                out[f"{o}__{i}"] = np.zeros((T, T))
    for key in ("Y_D", "C_D", "n_inter_D", "beta_D"):
        out[f"{key}_ss"] = np.array(float(ss[key]))
    np.savez_compressed(CACHE, **out)
    print(f"wrote {CACHE}")
    return CACHE


if __name__ == "__main__":
    build()
