"""CB-block audit — live pipeline probe.

Runs main.py's stages 1-5 plus run_tpi (NO figure generation, no source edits),
then dumps everything the audit's Steps 0/4/5 need to
diagnostics/cb_audit/probe_pipeline.npz + .json.

Audit-only. Reads the production modules; writes nothing outside cb_audit/.
"""
import os, sys, json, datetime
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "code"))
sys.path.insert(0, HERE)


def ts():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(m):
    print(f"[{ts()}] {m}", flush=True)


def main():
    from calibration import get_calibration, EBA_CALIBRATION, BANK_SCOPE
    from steady_state import solve_steady_state, gk_feasibility_margin
    from ic_delta_calibration import calibrate_ic_delta
    from depreciation_calibration import calibrate_depreciation
    from full_model import build_and_solve
    import tpi as tpi_mod

    out = {"timestamp": ts(), "git_head": os.popen("git -C %s rev-parse HEAD" % ROOT).read().strip()}
    arrays = {}

    # ── Step 0: which base? ───────────────────────────────────────────────────
    log("STEP 0: calibration flags")
    cal = get_calibration()
    out["EBA_CALIBRATION"] = bool(EBA_CALIBRATION)
    out["BANK_SCOPE"] = BANK_SCOPE
    for k in ("psi_lambda_B_D", "psi_lambda_B_F", "Delta_bD_D", "Delta_bF_D",
              "Delta_bD_F", "Delta_bF_F", "zeta_writeoff_D", "zeta_writeoff_F",
              "writeoff_enabled_D", "writeoff_enabled_F", "recovery_rate_D",
              "kappa_cb_F", "size_F", "phi_bD_D_ss", "phi_bF_D_ss",
              "phi_bD_F_ss", "phi_bF_F_ss", "f_D", "f_F", "mv_rule_D",
              "phi_lamb_D", "psi_bF_D", "psi_bD_F"):
        if k in cal:
            out[f"cal.{k}"] = float(cal[k])
    log(f"  EBA_CALIBRATION={EBA_CALIBRATION}  BANK_SCOPE={BANK_SCOPE}")

    log("solving steady state (stages 2-4)...")
    ssr = solve_steady_state(cal)
    ssr = calibrate_ic_delta(ssr)
    ssr = calibrate_depreciation(ssr)
    ss = ssr["ss_final"] if "ss_final" in ssr else ssr
    log("steady state solved.")

    def g(k, default=None):
        try:
            return float(ss[k])
        except Exception:
            try:
                return float(ss.toplevel[k])
            except Exception:
                return default

    for k in ("lambda_gk_D", "lambda_gk_F", "Omega_D", "Omega_F",
              "nu_K_D", "nu_K_F", "nu_bD_D", "nu_bF_D", "nu_bD_F", "nu_bF_F",
              "eta_D", "eta_F", "theta_D", "theta_F",
              "n_inter_D", "n_inter_F", "q_b_D", "q_b_F", "b_D_D", "b_D_F",
              "b_F_D", "b_F_F", "b_gov_D", "b_gov_F", "K_D", "K_F", "Q_D", "Q_F",
              "Delta_bD_eff_D", "Delta_bF_eff_D", "Delta_bD_eff_F", "Delta_bF_eff_F",
              "def_rate_D", "def_rate_F", "EL_load_D", "EL_load_F",
              "rk_D", "rk_F", "rdep_D", "rdep_F", "rb_exp_D", "rb_exp_F",
              "rb_actual_D", "rb_actual_F", "delta_b_D", "delta_b_F",
              "Y_D", "Y_F", "beta_D", "beta_F", "beta_inter_D",
              "gk_wedge_F_D_ss", "gk_wedge_D_F_ss", "kappa_cb_F", "cb_buy_D",
              "cb_flow_D", "rem_cb_D", "rem_cb_F",
              "goods_mkt_D", "goods_mkt_F", "ca_res_D", "ic_res_D", "ic_res_F"):
        out[f"ss.{k}"] = g(k)

    # concentration ratios actually in force
    for c, own, cross in (("D", "b_D_D", "b_F_D"), ("F", "b_F_F", "b_D_F")):
        qo = "q_b_D" if c == "D" else "q_b_F"
        qc = "q_b_F" if c == "D" else "q_b_D"
        n = g(f"n_inter_{c}")
        out[f"ss.phi_own_{c}"] = g(qo) * g(own) / n
        out[f"ss.phi_cross_{c}"] = g(qc) * g(cross) / n
    out["ss.gk_margin_D"] = gk_feasibility_margin(
        g("theta_D"), out["cal.f_D"], out["ss.phi_own_D"], out["ss.phi_cross_D"],
        g("Delta_bD_eff_D"), g("Delta_bF_eff_D"))
    out["ss.gk_margin_F"] = gk_feasibility_margin(
        g("theta_F"), out["cal.f_F"], out["ss.phi_own_F"], out["ss.phi_cross_F"],
        g("Delta_bF_eff_F"), g("Delta_bD_eff_F"))

    log("STEP 0 key numbers:")
    for k in ("ss.lambda_gk_D", "ss.lambda_gk_F", "ss.Omega_D", "ss.Omega_F",
              "ss.nu_K_D", "ss.nu_K_F", "ss.phi_own_D", "ss.phi_own_F",
              "ss.gk_margin_D", "ss.gk_margin_F"):
        log(f"  {k:22s} = {out[k]}")

    # ── Step 5 (SS part): CB block SS-neutrality, evaluated directly ─────────
    log("STEP 5: CB block SS neutrality")
    import copy
    ss_tpi = copy.deepcopy(ss)
    ss_tpi.toplevel["cb_buy_D"] = 0.0
    ss_tpi.toplevel["cb_flow_D"] = 0.0
    ss_tpi.toplevel["kappa_cb_F"] = float(cal["kappa_cb_F"])

    def _blk_ss(blk, ssobj):
        """Evaluate a @simple block at a steady state, return its outputs."""
        args = {k: ssobj[k] for k in blk.inputs}
        return blk.steady_state(args)

    step5 = {}
    for name, blk in (("budget_residual_D_tpi", tpi_mod.budget_residual_D_tpi),
                      ("budget_residual_F_tpi", tpi_mod.budget_residual_F_tpi),
                      ("external_account_D_tpi", tpi_mod.external_account_D_tpi),
                      ("domestic_bond_clearing_tpi", tpi_mod.domestic_bond_clearing_tpi)):
        try:
            r = _blk_ss(blk, ss_tpi)
            step5[name] = {k: float(v) for k, v in r.items()}
        except Exception as e:                       # noqa: BLE001 - audit probe
            step5[name] = {"ERROR": repr(e)}
        log(f"  {name}: {step5[name]}")
    # and the NON-TPI counterparts, for the difference
    from equations_D import budget_residual_D
    from equations_global import external_account_D, domestic_bond_clearing
    for name, blk in (("budget_residual_D", budget_residual_D),
                      ("external_account_D", external_account_D),
                      ("domestic_bond_clearing", domestic_bond_clearing)):
        try:
            r = _blk_ss(blk, ss_tpi)
            step5[name] = {k: float(v) for k, v in r.items()}
        except Exception as e:                       # noqa: BLE001
            step5[name] = {"ERROR": repr(e)}
        log(f"  {name}: {step5[name]}")
    out["step5_ss_blocks"] = step5

    # ── Steps 2/4: full TPI solve ────────────────────────────────────────────
    log("building dynamic model + baseline Jacobian (this is the ~6 min leg)...")
    mr = build_and_solve(ssr)
    log("baseline done. running TPI...")
    tr = tpi_mod.run_tpi(mr)
    log("TPI done.")

    G = tr["G_tpi"]
    T = tr["T"]
    A_cb = np.array(G["spread_rb"]["cb_buy_D"])
    A_def = np.array(G["spread_rb"]["shock_def_D"])
    arrays["A_cb"] = A_cb
    arrays["A_def"] = A_def
    for o in ("q_b_D", "n_inter_D", "b_D_D", "K_D", "Y_D", "theta_D",
              "b_gov_D", "rem_cb_D", "rem_cb_F", "cb_flow_D", "TAX_D", "TAX_F"):
        if o in G.outputs:
            for inp in ("cb_buy_D", "shock_def_D"):
                # A pure-CB object (rem_cb_D, cb_flow_D) has no shock_def_D column:
                # it responds to nothing but the CB's own purchases. Recorded, not
                # assumed away.
                try:
                    arrays[f"G_{o}__{inp}"] = np.array(G[o][inp])
                except KeyError:
                    out.setdefault("G_missing_columns", []).append(f"{o}__{inp}")
    out["G_outputs_has_rem_cb_D"] = "rem_cb_D" in G.outputs
    out["G_outputs_has_rem_cb_F"] = "rem_cb_F" in G.outputs
    out["G_outputs_has_cb_flow_D"] = "cb_flow_D" in G.outputs

    out["A_cb_00"] = float(A_cb[0, 0])
    out["A_cb_diag_first10"] = [float(A_cb[i, i]) for i in range(10)]
    out["A_cb_col0_first10"] = [float(A_cb[i, 0]) for i in range(10)]
    ev = np.linalg.eigvals(A_cb)
    out["A_cb_max_real_eig"] = float(np.max(ev.real))
    out["A_cb_spectral_radius"] = float(np.max(np.abs(ev)))
    # closed-loop pole
    from lottery_math_shim import pole_scan
    out["closed_loop_pole"] = pole_scan(A_cb, T)

    for g_ in tr["gamma_values"]:
        irf = tr["irfs_tpi"][g_]
        out[f"gamma{g_}.peak_spread_bp"] = float(irf["spread_rb"][:100].max() * 4e4)
        out[f"gamma{g_}.b_gov_D_499"] = float(irf["b_gov_D"][min(499, T - 1)])
        out[f"gamma{g_}.n_inter_D_0"] = float(irf["n_inter_D"][0])
        out[f"gamma{g_}.Y_D_0"] = float(irf["Y_D"][0])
        out[f"gamma{g_}.max_ca_res_D"] = float(np.max(np.abs(irf["ca_res_D"])))
        out[f"gamma{g_}.max_goods_mkt_F"] = float(np.max(np.abs(irf["goods_mkt_F"])))
        out[f"gamma{g_}.max_goods_mkt_D"] = float(np.max(np.abs(irf["goods_mkt_D"])))
        for o in ("spread_rb", "n_inter_D", "q_b_D", "Y_D", "b_gov_D", "b_D_D",
                  "K_D", "cb_buy_D", "TAX_D", "TAX_F", "C_D", "C_F"):
            if o in irf:
                arrays[f"irf_g{g_}_{o}"] = np.array(irf[o])
        for o in ("rem_cb_D", "rem_cb_F", "cb_flow_D"):
            if o in irf:
                arrays[f"irf_g{g_}_{o}"] = np.array(irf[o])
    out["pnl_by_gamma"] = {str(k): v for k, v in tr["pnl_by_gamma"].items()}
    out["kappa_cb_F_used"] = float(tr["kappa_cb_F"])
    arrays["gammas_fine"] = tr["gammas_fine"]
    arrays["peak_arr"] = tr["peak_arr"]
    arrays["loading_arr"] = tr["loading_arr"]

    np.savez_compressed(os.path.join(HERE, "probe_pipeline.npz"), **arrays)
    with open(os.path.join(HERE, "probe_pipeline.json"), "w") as f:
        json.dump(out, f, indent=2, default=str)
    log("WROTE probe_pipeline.npz / .json")


if __name__ == "__main__":
    main()
