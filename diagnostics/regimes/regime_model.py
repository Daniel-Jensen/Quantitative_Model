"""
regime_model.py — solve MAIN's pipeline once, compute G_tpi at psi_lambda_B in
{1.1793 (main's calibrated value), 0}, and cache per-output response matrices
M[o][i] (i in {shock_def_D, cb_buy_D}) plus SS meta.

REBUILT FOR MAIN (2026-07-23). Differs from the ms-regime version:
  * Main has no build_tpi_model helper — the TPI model is assembled inline here
    from main's block list (referenced via the `tpi` module namespace, which
    imports every block), matching code/tpi.py's run_tpi exactly. Correctness is
    guarded by main's own G_tpi[cb=0] vs baseline-Jacobian sanity check (<1e-8).
  * Main's ECB capital-key conduit needs cb_flow_D and kappa_cb_F in the SS.
  * recovery_rate=0.30 on main, so EL_price is NOT 0.102491 — it is logged, not
    asserted against the ms-regime anchor.
  * psi_lambda_B is read live from calibration.py (the 150bp-target value; 8.5
    under BANK_SCOPE="broad" since 2026-07-31). The breakdown ceiling is
    calibration-dependent and re-derived per calibration — see PSILAM_BREAKDOWN.
  * On main the CB backstop COMPRESSES the spread (A_cb<0): the ms-regime SA-1
    pathology is absent (capital-key conduit socialises the funding to F).
"""
import os, sys, copy, datetime, hashlib
import numpy as np
import sequence_jacobian as sj

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), "code"))
LOG = os.path.join(HERE, "regimes_log.md")

def _live_psilam():
    """psi_lambda_B from the live calibration — this is BOTH the cache key and the
    provenance anchor, so it must track calibration.py rather than be hardcoded.
    Was hardcoded 1.1793 (the EBA-era 150bp-target value); after the 2026-07-30
    pre-EBA revert the live value is 3.0 and a hardcoded anchor silently pointed at
    a stale cache built under a different model. (code/ is already on sys.path from
    the module-level insert above.)"""
    from calibration import get_calibration  # noqa: PLC0415 - dynamic sys.path
    return float(get_calibration()["psi_lambda_B_D"])


PSILAM_MAIN = _live_psilam()

# Linear-approximation-breakdown threshold. CALIBRATION-DEPENDENT — it scales with
# bank net worth, so it must be re-derived whenever n_inter moves:
#   n_inter_D = 0.408 (EBA 2011, CT1 scope) -> breakdown from ~1.5-2.0
#   n_inter_D = 3.0   (pre-EBA placeholder) -> breakdown from ~4-5
#   n_inter_D = 2.138 (EBA 2011, BANK_SCOPE="broad", LIVE) -> see below
#
# RE-DERIVED 2026-07-31 for the broad scope by diagnostics/psilam_breakdown_sweep.py
# (16-point grid, steady state solved once because psi_spread is exactly linear in
# psi_lambda_B, Jacobian re-solved per point with BOTH dials moved together):
#
#   psi_lambda_B   8.5    14     15     18     20     25      26      27      28
#   peak spread   150.3  223.8  225.8  234.6  273.6  625.3  1034.5  8903.8   41.0(flip)
#   n_inter_D[0]  -3.38  -4.82  -4.63  -4.33  -5.29 -15.47  -27.72 -264.92  +29.87
#
# There is a POLE between 27 and 28 — the amplification denominator crosses zero, the
# response runs up superlinearly and comes back sign-flipped. The A7 >1000bp flag
# first fires at 26. But the FIRST pathology is milder and earlier: over [14,18]
# n_inter_D[0] *shrinks* while the peak spread keeps rising, which is already an
# internally inconsistent doom loop.
#
# The guard is set from that first pathology, not from the pole, so it has real
# margin (15.0 is 55% of the pole and 1.76x the live psi_lambda_B=8.5) rather than
# sitting on the edge of a singularity. run_regimes.py's A7 flag remains the
# empirical backstop — it catches breakdown by measurement, not by threshold.
PSILAM_BREAKDOWN = 15.0

# Bump whenever REQUIRED / OPTIONAL / SS_META change. It goes in the cache
# FILENAME: the calibration fingerprint alone cannot detect a schema change, so
# without this an old cache would reload under the same name missing the new
# keys — silently, because irf_all discovers outputs by scanning cache keys.
CACHE_SCHEMA = 3

REQUIRED = ["spread_rb", "rb_D", "rb_F", "q_b_D", "q_b_F", "Y_D", "C_D", "I_D",
            "NX_D", "K_D", "n_inter_D", "b_D_D", "b_D_F", "b_gov_D", "U_D", "U_F",
            "TAX_D", "P_CES_D",
            # Added for the experiments package (schema 2):
            #   Phi_D, def_rate_D — Phi_D closes the market_clearing_D identity
            #   for E2; def_rate_D is the off-path expected-loss leg for E1's
            #   A5-1 reporting (cb_pnl reads it).
            "Phi_D", "def_rate_D"]
# OPTIONAL: logged loudly if missing, never silently dropped.
# (cb_flow_D excluded — it's the CB inter-block conduit flow, unused downstream, and
#  SSJ returns its cb_buy_D Jacobian as a non-array object; keep the cache clean.)
#
# T_D is OPTIONAL, not REQUIRED, on purpose: T0=T1=0 so the macroprudential bond
# tax is identically zero and SSJ may omit it from G.outputs entirely. Zero-filling
# is the CORRECT value here rather than a silent hole — and E2's closure assertion
# catches it either way if that ever stops being true.
OPTIONAL = ["G_D", "ra_D", "lambda_gk_D", "theta_D", "GINI_WEALTH", "GINI_C",
            "div_fund_D", "T_D"]
SS_META  = ["q_b_D_ss:q_b_D", "b_D_D_ss:b_D_D", "b_gov_D_ss:b_gov_D", "Y_D_ss:Y_D",
            "C_D_ss:C_D", "I_D_ss:I_D", "NX_D_ss:NX_D", "n_inter_D_ss:n_inter_D",
            "K_D_ss:K_D", "TAX_D_ss:TAX_D", "P_CES_D_ss:P_CES_D",
            "beta_D:beta_D", "beta_F:beta_F", "EL_price_D:EL_price_D",
            # schema 2: needed by E1's cb_pnl port and E2's identity
            # schema 3: delta_b_F is NOT delta_b_D (0.0568 vs 0.0777 — the two
            # countries' bank books have different measured durations). E1's
            # cb_pnl computes the SS yield on each leg as delta_b*(1/q_b_ss - 1),
            # so using D's duration on the F leg puts the SS spread at -9.2e-04
            # instead of its true -5.0e-08 and silently contaminates carry_ss_pv.
            "delta_b_D_ss:delta_b_D", "delta_b_F_ss:delta_b_F", "q_b_F_ss:q_b_F",
            "Phi_D_ss:Phi_D"]


def log(m=""):
    print(m, flush=True)
    with open(LOG, "a") as f:
        f.write(m + "\n")


def _calibration_fingerprint():
    """Short hash of the WHOLE live calibration, used in the cache filename.

    Without this the filename keys only on psi_lambda_B, so the psi_lambda_B=0
    cache built under one calibration is silently reused under another — exactly
    the failure mode that bit PSILAM_MAIN when it was hardcoded, and that the
    2026-07-31 EBA rebuild would have hit (its psilam=0 baseline is a different
    model from the pre-EBA one bearing the same filename). Any calibration change
    now yields a new filename, so a stale cache can never be picked up.
    """
    from calibration import get_calibration  # noqa: PLC0415 - dynamic sys.path
    cal = get_calibration()
    payload = ";".join(f"{k}={float(v):.12g}" for k, v in sorted(cal.items())
                       if isinstance(v, (int, float)) and not isinstance(v, bool))
    return hashlib.sha256(payload.encode()).hexdigest()[:8]


# Kept as a provenance snapshot of the calibration at import. NOT used for cache
# filenames — see cache_path, which must read the live calibration so a
# calibration_override (experiments/common.py) mints its own filename instead of
# clobbering the baseline cache.
CAL_FINGERPRINT = _calibration_fingerprint()


def cache_path(psilam, fingerprint=None):
    """Cache filename for a given psi_lambda_B at the LIVE calibration.

    The fingerprint is computed at CALL time, not import time. An override applied
    after this module was imported must produce a different filename.
    """
    fp = fingerprint or _calibration_fingerprint()
    tag = f"{psilam:.2f}".replace(".", "p")
    return os.path.join(HERE, f"cache_G_main_v{CACHE_SCHEMA}_psilam{tag}_cal{fp}.npz")


def build_tpi_model_main(tpi, financial_solved_D, financial_solved_F,
                         hh_D=None, hh_F=None):
    """Assemble main's TPI-extended model — identical block list to code/tpi.py's
    run_tpi (blocks referenced via the tpi module, which imports them all). The
    two financial_solved blocks are runtime-constructed, passed in.

    hh_D / hh_F optionally REPLACE the household blocks with versions carrying
    extra hetoutputs (experiments/e4_distribution.py adds per-decile consumption).
    Substituting here rather than assembling a second block list keeps this the
    single place the model is defined — a second copy is how the retired
    audit_artifacts/ harness drifted into testing a different model."""
    t = tpi
    hh_D = t.hh_extended_D if hh_D is None else hh_D
    hh_F = t.hh_extended_F if hh_F is None else hh_F
    return sj.create_model([
        t.deposit_return_D, t.tax_rule_D, hh_D, t.ghh_composite_D,
        t.sdf_D, t.sdf_banker_D, t.government_default_D, financial_solved_D,
        t.bond_return_D, t.bank_return_D, t.capital_fund_D, t.cap_adj_cost_inter_D, t.macro_pru_tax_D,
        t.intermediation_P2_D, t.intermediation_P3_D, t.k_balance_sheet_D,
        t.capital_adj_D, t.capital_producer_profit_D, t.budget_residual_D_tpi,
        t.labor_D, t.labor_market_D, t.labor_demand_D, t.banker_div_res_D,
        t.market_clearing_D, t.welfare_agg_D,
        t.deposit_return_F, t.tax_rule_F, hh_F, t.ghh_composite_F,
        t.sdf_F, t.sdf_banker_F, t.government_default_F, financial_solved_F,
        t.bond_return_F, t.bank_return_F, t.capital_fund_F, t.cap_adj_cost_inter_F, t.macro_pru_tax_F,
        t.intermediation_P2_F, t.intermediation_P3_F, t.k_balance_sheet_F,
        t.capital_adj_F, t.capital_producer_profit_F, t.budget_residual_F_tpi,
        t.labor_F, t.labor_market_F, t.labor_demand_F, t.banker_div_res_F,
        t.market_clearing_F, t.welfare_agg_F,
        t.ces_price_D, t.import_demand_D, t.ces_price_F, t.import_demand_F,
        t.trade_balance, t.external_account_D_tpi, t.domestic_bond_clearing_tpi,
        t.bond_yield, t.portfolio_level_anchors, t.divert_portfolio_adj,
        t.divert_bond_foc_D, t.divert_bond_foc_F, t.global_goods_mkt,
    ], name="Full 2-Country MU HANK — TPI Extension (regimes cache, main)")


def _ss_tpi(ss_final, kappa_cb_F):
    ss = copy.deepcopy(ss_final)
    ss.toplevel["cb_buy_D"] = 0.0
    ss.toplevel["cb_flow_D"] = 0.0
    ss.toplevel["kappa_cb_F"] = kappa_cb_F
    return ss


def _solve_G(model, ss_tpi, unk, tgt, T, label):
    log(f"- {datetime.datetime.now():%Y-%m-%d %H:%M:%S} solving G_tpi at psi_lambda_B={label} ...")
    return model.solve_jacobian(ss_tpi, unknowns=unk, targets=tgt,
                                inputs=["Z_D", "shock_def_D", "Z_F", "shock_def_F", "cb_buy_D"], T=T)


def _col(G, o, i, T):
    """Response matrix of output o to input i. SSJ omits zero-response columns from
    an output's Jacobian dict, so an absent (o, i) means o does not respond to i at
    this calibration — a genuine zero T x T matrix, filled and logged (not silent)."""
    try:
        return np.asarray(G[o][i], dtype=float)   # coerce: reject non-float (object) Jacobians
    except KeyError:
        log(f"  note: output `{o}` has no Jacobian column for `{i}` at this calibration "
            f"-> zero response (filled 0 T x T); economically = o does not respond to i.")
        return np.zeros((T, T))


def _extract(G, ss, T, dshock, psilam):
    out = {"T": np.array(T), "dShock_def_D": np.asarray(dshock), "psi_lambda_B": np.array(psilam),
           "cal_fingerprint": np.array(_calibration_fingerprint())}
    missing_req = [o for o in REQUIRED if o not in G.outputs]
    if missing_req:
        raise RuntimeError(f"REQUIRED outputs missing from main's G_tpi: {missing_req}. "
                           f"Available: {sorted(G.outputs)}")
    for o in REQUIRED:
        out[f"{o}__shock_def_D"] = _col(G, o, "shock_def_D", T)
        out[f"{o}__cb_buy_D"]    = _col(G, o, "cb_buy_D", T)
    for o in OPTIONAL:
        if o in G.outputs:
            out[f"{o}__shock_def_D"] = _col(G, o, "shock_def_D", T)
            out[f"{o}__cb_buy_D"]    = _col(G, o, "cb_buy_D", T)
        else:
            log(f"  **MISSING OPTIONAL OUTPUT `{o}`** — not in main's G_tpi.outputs; "
                f"panel zero-filled/omitted WITH a caption note, never silently.")
    for spec in SS_META:
        name, key = spec.split(":")
        out[name] = np.array(float(ss[key]))
    return out


def build_caches(force=False):
    # Read live, not from the import-time PSILAM_MAIN: under a calibration_override
    # the two can differ, and the override must win.
    psilam_live = _live_psilam()
    paths = {psilam_live: cache_path(psilam_live), 0.0: cache_path(0.0)}
    if not force and all(os.path.exists(p) for p in paths.values()):
        return paths
    from calibration import get_calibration
    from steady_state import solve_steady_state
    from ic_delta_calibration import calibrate_ic_delta
    from depreciation_calibration import calibrate_depreciation
    from full_model import build_and_solve
    import tpi

    cal = get_calibration()
    assert cal["psi_lambda_B_D"] < PSILAM_BREAKDOWN, (
        f"HARD GUARD: psi_lambda_B >= {PSILAM_BREAKDOWN} — linear-approximation-"
        "breakdown region for this net-worth calibration (see PSILAM_BREAKDOWN note)")
    # Consistency, not provenance: psilam_live is read from the same calibration, so
    # this only fires if calibration.py changed under a live interpreter session.
    assert abs(cal["psi_lambda_B_D"] - psilam_live) < 1e-9, (
        f"calibration drifted mid-run: live={cal['psi_lambda_B_D']} vs cache key {psilam_live}")
    kappa_cb_F = float(cal["kappa_cb_F"])

    ssr = calibrate_depreciation(calibrate_ic_delta(solve_steady_state(cal)))
    res = build_and_solve(ssr)
    ss, T, dshock = res["ss_final"], res["T"], res["dShock_def_D"]
    unk, tgt = res["unknowns_tp"], res["targets_tp"]
    irfs_def_D = res["irfs_def_D"]
    log(f"\n## Cache build (main model) — {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
    log(f"- calibration: psi_lambda_B={cal['psi_lambda_B_D']}, mv_rule={cal['mv_rule_D']}, "
        f"recovery_rate={cal['recovery_rate_D']}, kappa_cb_F={kappa_cb_F}")
    log(f"- EL_price_D = {float(ss['EL_price_D']):.6f} (main recovery=0.30; NOT the ms-regime 0.102491 anchor)")

    model = build_tpi_model_main(tpi, res["financial_solved_D"], res["financial_solved_F"])

    ss28 = _ss_tpi(ss, kappa_cb_F)
    G28 = _solve_G(model, ss28, unk, tgt, T, f"{psilam_live}")
    # main's own correctness guard: G_tpi[cb=0] spread must match the baseline Jacobian
    _chk = G28 @ {"Z_D": np.zeros(T), "Z_F": np.zeros(T), "shock_def_D": dshock,
                  "shock_def_F": np.zeros(T), "cb_buy_D": np.zeros(T)}
    _err = float(np.max(np.abs(_chk["spread_rb"][:50] - irfs_def_D["spread_rb"][:50])))
    log(f"- model-build sanity: G_tpi[cb=0] vs baseline spread_rb max|err| = {_err:.2e} (expect <1e-8)")
    assert _err < 1e-8, f"model build mismatch (err={_err:.2e}) — block list wrong"
    _acb0 = float(np.array(G28["spread_rb"]["cb_buy_D"])[0, 0])
    log(f"- cross-check vs SA-1 probe: d(spread_rb)/d(cb_buy)[0,0] = {_acb0:+.5e} "
        f"(probe found -1.9455e-2 → expect match; A_cb<0 = backstop COMPRESSES on main)")
    np.savez_compressed(paths[psilam_live], **_extract(G28, ss28, T, dshock, psilam_live))

    ss0 = _ss_tpi(ss, kappa_cb_F)
    ss0.toplevel["psi_lambda_B_D"] = 0.0; ss0.toplevel["psi_lambda_B_F"] = 0.0
    ss0.toplevel["psi_spread_D"]   = 0.0; ss0.toplevel["psi_spread_F"]   = 0.0
    G0 = _solve_G(model, ss0, unk, tgt, T, "0.0")
    np.savez_compressed(paths[0.0], **_extract(G0, ss0, T, dshock, 0.0))

    log(f"- caches written: {[os.path.basename(p) for p in paths.values()]}")
    return paths


def load_cache(psilam):
    """Load a cache, asserting it was built under the live calibration.

    cache_path already keys on the live fingerprint, so a mismatch normally shows
    up as a missing file. This second check catches the case where a file was
    hand-copied or renamed — it fails loudly instead of returning another model.
    """
    live = _calibration_fingerprint()
    path = cache_path(psilam, fingerprint=live)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"No cache at {os.path.basename(path)}. The live calibration has no cache "
            f"built for it — run:  /opt/anaconda3/envs/ssj/bin/python "
            f"diagnostics/regimes/regime_model.py")
    # allow_pickle for backward-compat with caches that stored a stray object entry
    # (cb_flow_D, now excluded); the matrices this module reads are all plain float.
    with np.load(path, allow_pickle=True) as d:
        cache = {k: d[k] for k in d.files if not d[k].dtype == object}
    stored = str(cache["cal_fingerprint"])
    assert stored == live, (
        f"cache fingerprint {stored} != live calibration {live} — stale or "
        f"hand-renamed cache; rebuild with regime_model.py --force")
    return cache


if __name__ == "__main__":
    build_caches(force="--force" in sys.argv)
