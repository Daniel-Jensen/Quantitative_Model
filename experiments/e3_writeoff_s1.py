"""E3 — S-1 and the payoff specification.

REBASED 2026-08-18 on the structural GK refactor. The baseline now PRICES the full
default loss (`zeta_writeoff = 1`), so the old E3 question — "what if the principal
were written down too?" — is the baseline, not a variant. Two questions remain, and
the two switches answer one each.

                              writeoff_enabled  zeta_writeoff   steady state
  baseline                           0                1.0       --
  e3a_realised_writeoff              1                1.0       STRICTLY INVARIANT
  e3b_coupon_only_pricing            0                0.0       STRICTLY INVARIANT

`writeoff_enabled` selects which BRANCH the impulse response traces; `zeta_writeoff`
governs what is PRICED. In bond_return_D / government_ss_D / budget_residual_D:

    realised coupon = delta_b * (1 - def_rate*haircut*writeoff_enabled)
    realised cont.  = (1-delta_b)*q_b * (1 - zeta*def_rate*haircut*writeoff_enabled)

Both legs carry def_rate, which is 0 at the steady state, so **e3a is exactly
SS-neutral** — it is the clean S-1 test of whether realised (as opposed to priced)
losses change the transmission.

**e3b is the §12 Arm-3 diagnostic**: the pre-refactor coupon-only payoff. Measured
2026-08-18, it is ALSO exactly SS-neutral, and the reason is worth stating because it
is not obvious. zeta multiplies the continuation-value haircut `zeta*def_rate*h`, and
`def_rate_ss = 0` kills that term inside `rb_exp` just as it kills the realised one —
so the priced and realised payoffs coincide at the steady state whatever zeta is, and
`q_b_D = 0.974906` in every arm. zeta is allocation-neutral while still changing the
LINEARISED bond pricing equation, and hence every dynamic result. Expected loss per
unit of default probability falls from

    (1-rec) * [delta_b + (1-delta_b) q_b] / q_b  ~ 0.7014     (zeta = 1)
to  (1-rec) *  delta_b                   / q_b  ~ 0.0558     (zeta = 0),

a factor of ~12.6. The gap is the principal/continuation loss on a 12.9-quarter
claim, and it is what the deleted free parameter `psi_spread_D` used to stand in for.
e3b is reported to size that contribution, NOT as an economic specification.

Recovery stays at 0.30 (EL-1's resolved Greek-PSI NPV value).
"""
import numpy as np

from common import calibration_override, load_cache, provenance, write_results

VARIANTS = {
    "e3a_realised_writeoff": {"writeoff_enabled_D": 1.0, "writeoff_enabled_F": 1.0,
                              "zeta_writeoff_D": 1.0, "zeta_writeoff_F": 1.0},
    "e3b_coupon_only_pricing": {"writeoff_enabled_D": 0.0, "writeoff_enabled_F": 0.0,
                                "zeta_writeoff_D": 0.0, "zeta_writeoff_F": 0.0},
}
SS_INVARIANT_VARIANT = "e3a_realised_writeoff"

SS_INVARIANT_KEYS = ["q_b_D_ss", "Y_D_ss", "C_D_ss", "I_D_ss", "NX_D_ss",
                     "n_inter_D_ss", "K_D_ss", "TAX_D_ss", "P_CES_D_ss",
                     "b_gov_D_ss", "b_D_D_ss"]


def expected_EL_load(cal, q_b_D):
    """The closed form behind bond_return_D's EL_load_D output.

    EL_load is ENDOGENOUS now — an output of the payoff block, not the deleted
    EL_price_D anchor — so this is a cross-check of the block against its algebra.
    Evaluate at the VARIANT's own solved q_b_D, not the baseline's: the priced loss
    feeds the GK portfolio FOC, so a variant that changes zeta moves the price too.
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

    out = {"EL_load_D": float(cache["EL_load_D"]),
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
    baseline_EL = float(baseline_cache["EL_load_D"])

    # Solved ONCE, on the baseline, then held fixed across every variant so the
    # comparison changes the model without also changing the policy.
    gammas = named_regime_gammas(baseline_cache)

    payload = {"provenance": provenance(),
               "gammas": gammas,
               "gamma_note": "Solved on the BASELINE and held fixed across variants, "
                             "so a difference in the table is attributable to the "
                             "switch alone. medium = 25% peak-spread compression; "
                             "aggressive is NOT 50% -- that target sits beyond a "
                             "closed-loop pole at gamma ~ 27.3, so it falls back to the "
                             "maximum feasible intervention (~46.6%). See "
                             "common.named_regime_gammas.",
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
        el_expected = expected_EL_load(cal, q_b_D)
        el_actual = float(cache["EL_load_D"])
        assert abs(el_actual - el_expected) < 1e-12, (
            f"{name}: EL_load_D={el_actual:.9f} != closed form {el_expected:.9f}. "
            f"code/equations_D.py bond_return_D no longer matches this experiment's "
            f"model of it — reconcile before reporting.")

        drift = {k: float(cache[k]) - baseline_ss[k] for k in SS_INVARIANT_KEYS}
        max_drift = max(abs(v) for v in drift.values())
        payload["checks"][name] = {"EL_load_expected": el_expected,
                                   "EL_load_actual": el_actual,
                                   "EL_load_vs_baseline_ratio": el_actual / baseline_EL,
                                   "max_ss_drift": max_drift,
                                   "ss_drift": drift}

        if name == SS_INVARIANT_VARIANT:
            # writeoff_enabled multiplies terms that already carry def_rate_ss = 0,
            # and zeta is unchanged from baseline, so the SS must be bit-identical.
            assert max_drift < 1e-10, (
                f"{name} moved the steady state (max drift {max_drift:.3e}). "
                f"writeoff_enabled is supposed to be SS-neutral — every writeoff term "
                f"is multiplied by def_rate_ss=0. Investigate before reporting.")
        else:
            # e3b is SS-neutral too (zeta multiplies def_rate_ss = 0 inside rb_exp), so
            # assert BOTH: the allocation does not move, and the priced loading does.
            # The second half is what proves the override actually reached the solve.
            assert max_drift < 1e-10, (
                f"{name} moved the steady state (max drift {max_drift:.3e}). zeta is "
                f"allocation-neutral: it multiplies def_rate_ss = 0 in both rb_exp and "
                f"rb_actual. Investigate before reporting.")
            assert abs(el_actual - baseline_EL) > 1e-6, (
                f"{name}'s EL_load ({el_actual:.9f}) is indistinguishable from baseline "
                f"({baseline_EL:.9f}) — the zeta override did not reach the SS solve.")

    write_results("e3_writeoff_s1", payload)
    return payload


if __name__ == "__main__":
    res = run()
    print(f"\n{'setting':>24} {'EL_load':>10} {'peak bp (passive)':>18} "
          f"{'loading (medium)':>17} {'loading (aggr.)':>16}")
    print("-" * 84)
    rows = [("baseline", res["baseline"])] + list(res["variants"].items())
    for name, r in rows:
        def ld(reg):
            v = r["regimes"][reg]["loading"]
            return "n/a" if v is None else f"{v:.2f}"
        print(f"{name:>24} {r['EL_load_D']:>10.4f} "
              f"{r['regimes']['passive']['peak_spread_bp_ann']:>18.1f} "
              f"{ld('medium'):>17} {ld('aggressive'):>16}")
    print("-" * 84)
    for name, c in res["checks"].items():
        print(f"{name}: EL_load {c['EL_load_actual']:.6f} (closed form "
              f"{c['EL_load_expected']:.6f}, {c['EL_load_vs_baseline_ratio']:.2f}x "
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
    print("\nThe baseline no longer tunes any parameter to a 150bp moment — psi_lambda_B "
          "is 0 and psi_spread is deleted. Whatever peak spread these variants produce "
          "is a REPORTABLE FACT about the payoff specification, not a number to re-tune "
          "away.")
