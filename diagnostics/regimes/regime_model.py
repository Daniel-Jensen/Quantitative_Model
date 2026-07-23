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
  * psi_lambda_B=1.1793 (main's 150bp-target value). Main warns psi_lambda_B must
    not exceed ~1.5-2.0 (linear-approximation-breakdown region) — guarded here.
  * On main the CB backstop COMPRESSES the spread (A_cb<0): the ms-regime SA-1
    pathology is absent (capital-key conduit socialises the funding to F).
"""
import os, sys, copy, datetime
import numpy as np
import sequence_jacobian as sj

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), "code"))
LOG = os.path.join(HERE, "regimes_log.md")

PSILAM_MAIN = 1.1793

# Spec §8.1 output set, resolved against main's model. REQUIRED: hard error if absent.
REQUIRED = ["spread_rb", "rb_D", "rb_F", "q_b_D", "q_b_F", "Y_D", "C_D", "I_D",
            "NX_D", "K_D", "n_inter_D", "b_D_D", "b_D_F", "b_gov_D", "U_D", "U_F",
            "TAX_D", "P_CES_D"]
# OPTIONAL: logged loudly if missing, never silently dropped.
# (cb_flow_D excluded — it's the CB inter-block conduit flow, unused downstream, and
#  SSJ returns its cb_buy_D Jacobian as a non-array object; keep the cache clean.)
OPTIONAL = ["G_D", "ra_D", "lambda_gk_D", "theta_D", "GINI_WEALTH", "GINI_C", "div_fund_D"]
SS_META  = ["q_b_D_ss:q_b_D", "b_D_D_ss:b_D_D", "b_gov_D_ss:b_gov_D", "Y_D_ss:Y_D",
            "C_D_ss:C_D", "I_D_ss:I_D", "NX_D_ss:NX_D", "n_inter_D_ss:n_inter_D",
            "K_D_ss:K_D", "TAX_D_ss:TAX_D", "P_CES_D_ss:P_CES_D",
            "beta_D:beta_D", "beta_F:beta_F", "EL_price_D:EL_price_D"]


def log(m=""):
    print(m, flush=True)
    with open(LOG, "a") as f:
        f.write(m + "\n")


def cache_path(psilam):
    return os.path.join(HERE, f"cache_G_main_psilam{f'{psilam:.2f}'.replace('.', 'p')}.npz")


def build_tpi_model_main(tpi, financial_solved_D, financial_solved_F):
    """Assemble main's TPI-extended model — identical block list to code/tpi.py's
    run_tpi (blocks referenced via the tpi module, which imports them all). The
    two financial_solved blocks are runtime-constructed, passed in."""
    t = tpi
    return sj.create_model([
        t.deposit_return_D, t.tax_rule_D, t.hh_extended_D, t.ghh_composite_D,
        t.sdf_D, t.sdf_banker_D, t.government_default_D, financial_solved_D,
        t.bond_return_D, t.bank_return_D, t.capital_fund_D, t.cap_adj_cost_inter_D, t.macro_pru_tax_D,
        t.intermediation_P2_D, t.intermediation_P3_D, t.k_balance_sheet_D,
        t.capital_adj_D, t.capital_producer_profit_D, t.budget_residual_D_tpi,
        t.labor_D, t.labor_market_D, t.labor_demand_D, t.banker_div_res_D,
        t.market_clearing_D, t.welfare_agg_D,
        t.deposit_return_F, t.tax_rule_F, t.hh_extended_F, t.ghh_composite_F,
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
    out = {"T": np.array(T), "dShock_def_D": np.asarray(dshock), "psi_lambda_B": np.array(psilam)}
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
    paths = {PSILAM_MAIN: cache_path(PSILAM_MAIN), 0.0: cache_path(0.0)}
    if not force and all(os.path.exists(p) for p in paths.values()):
        return paths
    from calibration import get_calibration
    from steady_state import solve_steady_state
    from ic_delta_calibration import calibrate_ic_delta
    from depreciation_calibration import calibrate_depreciation
    from full_model import build_and_solve
    import tpi

    cal = get_calibration()
    assert cal["psi_lambda_B_D"] < 1.5, ("HARD GUARD: psi_lambda_B >= 1.5 — main's "
                                         "linear-approximation-breakdown region (STATE.md)")
    assert abs(cal["psi_lambda_B_D"] - PSILAM_MAIN) < 1e-9, f"expected main default {PSILAM_MAIN}"
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
    G28 = _solve_G(model, ss28, unk, tgt, T, f"{PSILAM_MAIN}")
    # main's own correctness guard: G_tpi[cb=0] spread must match the baseline Jacobian
    _chk = G28 @ {"Z_D": np.zeros(T), "Z_F": np.zeros(T), "shock_def_D": dshock,
                  "shock_def_F": np.zeros(T), "cb_buy_D": np.zeros(T)}
    _err = float(np.max(np.abs(_chk["spread_rb"][:50] - irfs_def_D["spread_rb"][:50])))
    log(f"- model-build sanity: G_tpi[cb=0] vs baseline spread_rb max|err| = {_err:.2e} (expect <1e-8)")
    assert _err < 1e-8, f"model build mismatch (err={_err:.2e}) — block list wrong"
    _acb0 = float(np.array(G28["spread_rb"]["cb_buy_D"])[0, 0])
    log(f"- cross-check vs SA-1 probe: d(spread_rb)/d(cb_buy)[0,0] = {_acb0:+.5e} "
        f"(probe found -1.9455e-2 → expect match; A_cb<0 = backstop COMPRESSES on main)")
    np.savez_compressed(paths[PSILAM_MAIN], **_extract(G28, ss28, T, dshock, PSILAM_MAIN))

    ss0 = _ss_tpi(ss, kappa_cb_F)
    ss0.toplevel["psi_lambda_B_D"] = 0.0; ss0.toplevel["psi_lambda_B_F"] = 0.0
    ss0.toplevel["psi_spread_D"]   = 0.0; ss0.toplevel["psi_spread_F"]   = 0.0
    G0 = _solve_G(model, ss0, unk, tgt, T, "0.0")
    np.savez_compressed(paths[0.0], **_extract(G0, ss0, T, dshock, 0.0))

    log(f"- caches written: {[os.path.basename(p) for p in paths.values()]}")
    return paths


def load_cache(psilam):
    # allow_pickle for backward-compat with caches that stored a stray object entry
    # (cb_flow_D, now excluded); the matrices this module reads are all plain float.
    with np.load(cache_path(psilam), allow_pickle=True) as d:
        return {k: d[k] for k in d.files if not d[k].dtype == object}


if __name__ == "__main__":
    build_caches(force="--force" in sys.argv)
