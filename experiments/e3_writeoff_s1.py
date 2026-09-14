"""E3 — S-1: does the sovereign default produce realised bank losses?

Two nested variants, because the two switches do different things:

                     writeoff_enabled  zeta_writeoff   steady state
  baseline                  0               0.0        --
  E3a  coupon-only          1               0.0        STRICTLY INVARIANT
  E3b  full writeoff        1               1.0        MOVES, via EL_price

In bond_return_D / government_ss_D / bond_price_ss_D / budget_residual_D:

    current_payoff = delta_b * (1 - def_rate*haircut*writeoff_enabled)
    continuation   = (1-delta_b)*q_b * (1 - zeta*def_rate*haircut*writeoff_enabled)

Both legs carry def_rate, which is 0 at SS, so writeoff_enabled is SS-neutral. But
zeta_writeoff ALSO appears in the EL_price anchor (code/steady_state.py:107-112),
and there it is NOT gated by writeoff_enabled:

    EL_price = (1-recovery) * [delta_b + zeta*(1-delta_b)*q_b] / q_b

At the live calibration (recovery=0.30, delta_b_D=0.0777006, q_b_D=0.968941) that
takes EL_price_D from 0.056134 to ~0.7017 — about 12.5x. EL_price is the loading's
DENOMINATOR, so this lands directly on SPEC Live Claim 1. Reported, never re-tuned
away.

Recovery stays at 0.30 (EL-1's resolved Greek-PSI NPV value) rather than
docs/STATE.md's older recovery=0.40 suggestion, which predates EL-1 and would move
two dials at once.
"""
import numpy as np

from common import calibration_override, load_cache, provenance, write_results

VARIANTS = {
    "e3a_coupon_only": {"writeoff_enabled_D": 1.0, "writeoff_enabled_F": 1.0,
                        "zeta_writeoff_D": 0.0, "zeta_writeoff_F": 0.0},
    "e3b_full": {"writeoff_enabled_D": 1.0, "writeoff_enabled_F": 1.0,
                 "zeta_writeoff_D": 1.0, "zeta_writeoff_F": 1.0},
}

SS_INVARIANT_KEYS = ["q_b_D_ss", "Y_D_ss", "C_D_ss", "I_D_ss", "NX_D_ss",
                     "n_inter_D_ss", "K_D_ss", "TAX_D_ss", "P_CES_D_ss",
                     "b_gov_D_ss", "b_D_D_ss"]


def expected_EL_price(cal, q_b_D):
    """The closed form from code/steady_state.py:107-109.

    Note this must be evaluated at the VARIANT's own solved q_b_D, not the
    baseline's: EL_price feeds the bond FOC, so under E3b the price moves and the
    closed form has to be checked against the price the model actually settled on.
    """
    return ((1.0 - cal["recovery_rate_D"])
            * (cal["delta_b_D"] + cal["zeta_writeoff_D"] * (1.0 - cal["delta_b_D"]) * q_b_D)
            / q_b_D)


def build_variant(overrides):
    """Full SS + Jacobian re-solve under an overridden calibration.

    A full re-solve is REQUIRED, not merely safer. The cheap route — patch
    ss.toplevel and re-solve only the Jacobian, as regime_model.build_caches does
    for its psi_lambda_B=0 cache — PRESUMES the SS invariance that E3a exists to
    test, which would make the check circular. E3b genuinely moves the SS, so it
    has no cheap route either.

    load_cache must also run inside the override: cache_path keys on the live
    calibration fingerprint, so outside the context it would look for (and fail to
    find, or worse, find) the baseline file.
    """
    from regime_model import build_caches
    with calibration_override(**overrides):
        build_caches()
        return load_cache()


def summarise(cache, gammas):
    """The E1 regime table under one variant, at a GIVEN set of gammas.

    gammas is passed in rather than re-solved per variant, and that is deliberate.
    Re-solving the compression targets under each variant would change the policy
    and the model at the same time, so a difference in the table could not be
    attributed to the writeoff switch. Holding gamma fixed compares the SAME
    policy across settings.

    It is also a necessity: under E3b the peak spread is NOT monotone in gamma, so
    `gamma_for_compression`'s bisection is invalid there and the named regimes are
    undefined. See `compression_feasible`.
    """
    import e1_backstop_schedule as e1
    from common import BP_ANN, irf_from_cache
    from lottery_math import closed_loop

    A_def, A_cb = cache["spread_rb__shock_def_D"], cache["spread_rb__cb_buy_D"]
    eps = np.asarray(cache["dShock_def_D"])
    Y_ss, n_ss = float(cache["Y_D_ss"]), float(cache["n_inter_D_ss"])

    out = {"EL_price_D": float(cache["EL_price_D"]),
           "A_cb_impact": float(A_cb[0, 0]),
           "regimes": {}}
    for name, gamma in gammas.items():
        _spread, cb = closed_loop(A_def, A_cb, eps, float(gamma))
        irf = irf_from_cache(cache, cb, eps)
        pnl = e1.cb_pnl(irf, cache)
        out["regimes"][name] = {
            "gamma": float(gamma),
            "peak_spread_bp_ann": float(np.max(np.asarray(irf["spread_rb"])[:100]) * BP_ANN),
            "Y_D_impact_pct_ss": float(np.asarray(irf["Y_D"])[0] * 100.0 / Y_ss),
            "n_inter_D_impact_pct_ss": float(np.asarray(irf["n_inter_D"])[0] * 100.0 / n_ss),
            "expected_loss_pv_pct_Y": 100.0 * pnl["el_pv"] / Y_ss,
            "premium_pv_pct_Y": 100.0 * pnl["prem_pv"] / Y_ss,
            "loading": (pnl["prem_pv"] / pnl["el_pv"]) if pnl["el_pv"] > 1e-16 else None,
        }
    return out


def compression_feasible(cache, gamma_max=15.0, n=40):
    """Is peak spread monotone decreasing in gamma, so compression targeting works?

    `gamma_for_compression` bisects on peak spread and raises if it is not
    monotone. Under E3b it is not, which means the named regimes (defined as
    25%/50% peak-spread compression) are UNDEFINED under that setting — there is
    no unique gamma delivering a given compression. That is a result about the
    policy experiment's own parameterisation, not a numerical nuisance, so it is
    measured and reported rather than worked around.
    """
    from common import BP_ANN
    from lottery_math import closed_loop

    A_def, A_cb = cache["spread_rb__shock_def_D"], cache["spread_rb__cb_buy_D"]
    eps = np.asarray(cache["dShock_def_D"])
    gammas = np.linspace(0.0, gamma_max, n)
    peaks = np.array([float(np.max(closed_loop(A_def, A_cb, eps, float(g))[0][:100])) * BP_ANN
                      for g in gammas])
    d = np.diff(peaks)
    monotone = bool(np.all(d <= 1e-9))
    viol = gammas[1:][d > 1e-9]
    return {
        "monotone_decreasing": monotone,
        "gamma_grid": gammas.tolist(),
        "peak_spread_bp_ann": peaks.tolist(),
        "first_violation_gamma": float(viol[0]) if viol.size else None,
        "n_violations": int(viol.size),
        "peak_at_gamma0_bp": float(peaks[0]),
        "min_peak_bp": float(peaks.min()),
        "gamma_at_min_peak": float(gammas[int(peaks.argmin())]),
    }


def run():
    # Import the MODULE, not the function. `from calibration import
    # get_calibration` here would bind the original function object before the
    # override context is ever entered, so `expected_EL_price` would be handed the
    # baseline zeta and the closed-form check would compare E3b's solved EL_price
    # against the zeta=0 value. This is the exact footgun documented in
    # common.calibration_override's docstring; it was written anyway and the
    # assertion caught it.
    import calibration
    from common import named_regime_gammas

    baseline_cache = load_cache()
    baseline_ss = {k: float(baseline_cache[k]) for k in SS_INVARIANT_KEYS}
    baseline_EL = float(baseline_cache["EL_price_D"])

    # Solved ONCE, on the baseline, then held fixed across every variant so the
    # comparison changes the model without also changing the policy.
    gammas = named_regime_gammas(baseline_cache)

    payload = {"provenance": provenance(),
               "gammas": gammas,
               "gamma_note": "Solved on the BASELINE (0/25/50% peak-spread "
                             "compression) and held fixed across variants, so a "
                             "difference in the table is attributable to the writeoff "
                             "switch alone. Under e3b_full the peak spread is not "
                             "monotone in gamma, so these targets are not even "
                             "well-defined there — see compression.",
               "baseline": summarise(baseline_cache, gammas),
               "compression": {"baseline": compression_feasible(baseline_cache)},
               "variants": {}, "checks": {}}

    for name, overrides in VARIANTS.items():
        cache = build_variant(overrides)
        payload["variants"][name] = summarise(cache, gammas)
        payload["compression"][name] = compression_feasible(cache)

        with calibration_override(**overrides):
            cal = calibration.get_calibration()   # resolved at USE time, inside the override
        q_b_D = float(cache["q_b_D_ss"])
        el_expected = expected_EL_price(cal, q_b_D)
        el_actual = float(cache["EL_price_D"])
        assert abs(el_actual - el_expected) < 1e-12, (
            f"{name}: EL_price_D={el_actual:.9f} != closed form {el_expected:.9f}. "
            f"code/steady_state.py:107-109 no longer matches this experiment's model "
            f"of it — reconcile before reporting.")

        drift = {k: float(cache[k]) - baseline_ss[k] for k in SS_INVARIANT_KEYS}
        max_drift = max(abs(v) for v in drift.values())
        payload["checks"][name] = {"EL_price_expected": el_expected,
                                   "EL_price_actual": el_actual,
                                   "EL_price_vs_baseline_ratio": el_actual / baseline_EL,
                                   "max_ss_drift": max_drift,
                                   "ss_drift": drift}

        if name == "e3a_coupon_only":
            # writeoff_enabled multiplies terms that already carry def_rate_ss = 0,
            # and zeta is unchanged, so the SS must be bit-identical. Drift is a bug.
            assert max_drift < 1e-10, (
                f"E3a moved the steady state (max drift {max_drift:.3e}). "
                f"writeoff_enabled is supposed to be SS-neutral — every writeoff term "
                f"is multiplied by def_rate_ss=0. Investigate before reporting.")
        else:
            # E3b DOES move the SS, through EL_price. Asserting invariance here would
            # be wrong; the closed-form check above is the check that applies. What we
            # do assert is that the override reached the solve at all.
            assert abs(el_actual - baseline_EL) > 1e-6, (
                f"E3b's EL_price ({el_actual:.9f}) is indistinguishable from baseline "
                f"({baseline_EL:.9f}) — the zeta override did not reach the SS solve.")

    write_results("e3_writeoff_s1", payload)
    return payload


if __name__ == "__main__":
    res = run()
    print(f"\n{'setting':>18} {'EL_price':>10} {'peak bp (passive)':>18} "
          f"{'loading (medium)':>17} {'loading (aggr.)':>16}")
    print("-" * 84)
    rows = [("baseline", res["baseline"])] + list(res["variants"].items())
    for name, r in rows:
        def ld(reg):
            v = r["regimes"][reg]["loading"]
            return "n/a" if v is None else f"{v:.2f}"
        print(f"{name:>18} {r['EL_price_D']:>10.4f} "
              f"{r['regimes']['passive']['peak_spread_bp_ann']:>18.1f} "
              f"{ld('medium'):>17} {ld('aggressive'):>16}")
    print("-" * 84)
    for name, c in res["checks"].items():
        print(f"{name}: EL_price {c['EL_price_actual']:.6f} (closed form "
              f"{c['EL_price_expected']:.6f}, {c['EL_price_vs_baseline_ratio']:.2f}x "
              f"baseline), max SS drift {c['max_ss_drift']:.3e}")

    print(f"\n{'setting':>18} {'compression targeting':>22} {'peak@g=0':>10} "
          f"{'min peak':>10} {'argmin g':>9}")
    print("-" * 74)
    for name, c in res["compression"].items():
        verdict = "OK (monotone)" if c["monotone_decreasing"] else \
                  f"INFEASIBLE (g>={c['first_violation_gamma']:.2f})"
        print(f"{name:>18} {verdict:>22} {c['peak_at_gamma0_bp']:>10.1f} "
              f"{c['min_peak_bp']:>10.1f} {c['gamma_at_min_peak']:>9.2f}")
    print("-" * 74)
    print("Where compression targeting is INFEASIBLE the named regimes are undefined "
          "(no unique gamma delivers a given compression), so all rows above are "
          "evaluated at the BASELINE's gammas held fixed.")
    print("\npsi_lambda_B was tuned to 150bp with realised losses OFF. Any overshoot "
          "here is a REPORTABLE FACT about whether the target survives S-1, not a "
          "number to re-tune away.")
