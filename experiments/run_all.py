"""Run every experiment and render docs/experiments_results.md.

E1 and E2 are seconds (post-Jacobian numpy on one cache). E3 is ~30-40 min because
it re-solves the model twice under overridden calibrations. Pass --skip-e3 to
regenerate the document from existing results without paying for that.

The rendered document is GENERATED. Do not hand-edit it — edit the experiment and
re-run, or the next run silently discards the edit.
"""
import json
import os
import sys

from common import RESULTS_DIR, ROOT, provenance

DOC = os.path.join(ROOT, "docs", "experiments_results.md")


def _load(name):
    path = os.path.join(RESULTS_DIR, f"{name}.json")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def _stamp(p):
    """Provenance line. Everything read live from the calibration, never hardcoded —
    run_regimes.py once shipped a hardcoded 'market-value rule' caption while
    actually running the par rule."""
    dirty = p.get("git_dirty")
    dirty_note = "" if dirty is False else (
        " · **working tree DIRTY**" if dirty else " · tree-state unknown")
    return (f"*Generated {p['generated']} from `{p['git_sha']}`{dirty_note} · "
            f"calibration `{p['cal_fingerprint']}` · scope `{p['BANK_SCOPE']}` · "
            f"`psi_lambda_B={p['psi_lambda_B_D']}` · `mv_rule={p['mv_rule_D']:g}` · "
            f"`recovery_rate={p['recovery_rate_D']}` · "
            f"`writeoff_enabled={p['writeoff_enabled_D']:g}` · "
            f"`zeta_writeoff={p['zeta_writeoff_D']:g}`*")


def _fmt_loading(v):
    return "n/a" if v is None else f"{v:.2f}"


def _render_e1(e1, L):
    L += ["## E1 — Backstop schedule", "", _stamp(e1["provenance"]), "",
          "γ selection: γ is **solved** for peak-spread compression, not chosen. "
          "`medium` = 25% (spec section 7). **`aggressive` is NOT 50%**: since the "
          "2026-08-18 GK structural refactor that target lies beyond a closed-loop pole "
          "at γ ≈ 27.3, so it falls back to the strongest intervention the model can "
          "represent — γ just below the pole, achieving **≈46.6%**. See "
          "`common.named_regime_gammas` / `lottery_math.closed_loop_pole`.", "",
          f"<sub>Rule recorded in the results file at run time: {e1['gamma_selection_rule']}</sub>", "",
          "| regime | γ | peak spread (bp ann) | Y_D[0] (% SS) | C_D[0] (% SS) "
          "| I_D[0] (% SS) | n_inter_D[0] (% SS) | loading |",
          "|---|---|---|---|---|---|---|---|"]
    for name, r in e1["regimes"].items():
        i = r["impact"]
        L.append(f"| {name} | {r['gamma']:.4f} | {r['peak_spread_bp_ann']:.1f} | "
                 f"{i['Y_D_pct_ss']:+.4f} | {i['C_D_pct_ss']:+.4f} | "
                 f"{i['I_D_pct_ss']:+.4f} | {i['n_inter_D_pct_ss']:+.3f} | "
                 f"{_fmt_loading(r['loading'])} |")

    L += ["", "### A5-1 — three separate objects", "",
          "**Do not sum these into an 'implicit transfer'.** Expected loss is computed "
          "off-path; reading it off the realised path shows the CB mechanically "
          "profiting, because the excess-return flow is a first-order deviation times a "
          "steady-state level and so does not vanish along the no-default branch.", "",
          "| regime | exposure PV (% Y) | expected loss PV (% Y) | `pd_D` differential PV |",
          "|---|---|---|---|"]
    for name, r in e1["regimes"].items():
        L.append(f"| {name} | {r['a5_1_exposure_pv_pct_Y']:.4f} | "
                 f"{r['a5_1_expected_loss_pv_pct_Y']:.5f} | "
                 f"{r['a5_1_greek_fiscal_saving_pv']:+.6f} |")
    L += ["", "> **The third column is misnamed in the code and its sign needs an author "
          "decision.** It reports `Σ β^t (pd_passive − pd_intervention)`, which is "
          "negative because the backstop lets Greece run a *larger* primary deficit — it "
          "relaxes required austerity. So negative means Greece is better off, the "
          "opposite of what \"fiscal saving\" implies. Flip the sign or rename it "
          "(\"austerity relief, PV\") before quoting it. Magnitudes are unaffected.", ""]

    sched = e1.get("loading_schedule")
    if sched:
        pts = [(g, l) for g, l in zip(sched["gamma"], sched["loading"]) if l is not None]
        mono = all(pts[i][1] >= pts[i + 1][1] for i in range(len(pts) - 1))
        L += [f"### Loading schedule (Live Claim 5)", "",
              f"{len(pts)} finite grid points over γ ∈ [{pts[0][0]:.2f}, {pts[-1][0]:.2f}]: "
              f"loading falls from **{pts[0][1]:.2f}** to **{pts[-1][1]:.2f}**. "
              f"Monotone decreasing: **{'YES' if mono else 'NO'}**.", "",
              "The self-extinguishing premium is the *decline*, so the schedule — not any "
              "single point — is the object.", "",
              "![loading schedule](../experiments/figures/fig_e1_loading_schedule.png)", ""]

    L += ["### Welfare (secondary)", "", f"{e1['welfare_caveat']}", "",
          "| regime | W_D | W_F |", "|---|---|---|"]
    for name, r in e1["regimes"].items():
        L.append(f"| {name} | {r['welfare_W_D_secondary']:+.4f} | "
                 f"{r['welfare_W_F_secondary']:+.4f} |")
    L.append("")
    return L


def _render_e2(e2, L):
    L += ["## E2 — ΔY decomposition", "", _stamp(e2["provenance"]), "",
          "Identity (`market_clearing_D`): "
          "`dY = P_ss·dC + C_ss·dP_CES + dI + dG + dΦ + dT + dNX`. `goods_mkt_D` is a "
          "targeted residual (≤1e−14), so this closes to solver tolerance — the "
          "decomposition is **self-verifying**, and the runner asserts closure at "
          f"{e2['closure_tol']:.0e} rather than warning.", "",
          "| regime | Y_D[0] (% SS) | Y_D trough (% SS) | dI PV | dNX PV | dC(qty) PV "
          "| dC(price) PV | max\\|residual\\| |",
          "|---|---|---|---|---|---|---|---|"]
    for name, r in e2["regimes"].items():
        c = r["components_pv"]
        L.append(f"| {name} | {r['dY_impact_pct_ss']:+.4f} | {r['dY_trough_pct_ss']:+.4f} | "
                 f"{c['investment']:+.3e} | {c['net_exports']:+.3e} | "
                 f"{c['consumption_quantity']:+.3e} | {c['consumption_price']:+.3e} | "
                 f"{r['max_abs_residual']:.2e} |")

    L += ["", "### Impact (t=0) decomposition, level deviations", "",
          "| component | " + " | ".join(e2["regimes"]) + " |",
          "|---" * (len(e2["regimes"]) + 1) + "|"]
    comps = ["consumption_quantity", "consumption_price", "investment",
             "government", "portfolio_cost", "macropru_tax", "net_exports"]
    for c in comps:
        row = " | ".join(f"{e2['regimes'][r]['components_impact'][c]:+.3e}"
                         for r in e2["regimes"])
        L.append(f"| {c} | {row} |")
    L.append("| **dY[0] total** | " + " | ".join(
        f"**{e2['regimes'][r]['dY_path'][0]:+.3e}**" for r in e2["regimes"]) + " |")

    # Computed live, never hardcoded. This caption asserted "+4.9e-04 / +2.2e-03 /
    # -1.9e-03, each roughly 4x the headline" until 2026-08-06 — flex-price numbers that
    # outlived the sticky-price re-tune. Under sticky prices the ordering REVERSES (the
    # headline becomes the larger object), so a hardcoded caption here does not merely go
    # stale, it states the opposite of what the table above it shows.
    _names = list(e2["regimes"])
    _lo, _hi = e2["regimes"][_names[0]], e2["regimes"][_names[-1]]
    _dY = _hi["dY_path"][0] - _lo["dY_path"][0]
    _dI = _hi["components_impact"]["investment"] - _lo["components_impact"]["investment"]
    _dNX = (_hi["components_impact"]["net_exports"]
            - _lo["components_impact"]["net_exports"])
    _big = max(abs(_dI), abs(_dNX))
    _residue = abs(_dY) < _big
    _lead = ("**The headline output number is the residue of larger offsetting "
             "channels.**" if _residue else
             "**The headline output number is no longer a residue of larger offsetting "
             "channels — it now exceeds each of them.**")
    _rel = (f"each roughly {_big / abs(_dY):.1f}x the headline and opposite in sign"
            if _residue else
            f"the largest single channel is {_big / abs(_dY):.2f}x the headline")
    L += ["", f"> {_lead} {_names[0]} → {_names[-1]}, `dY[0]` moves by "
          f"{_dY:+.2e} while investment moves {_dI:+.2e} and net exports {_dNX:+.2e} — "
          f"{_rel}. **Report the decomposition, not the headline ΔY** — the channels "
          "still offset, and `docs/SPEC.md`'s standing caution is about their "
          "cancellation, not about which term happens to be largest.", "",
          "> `government`, `portfolio_cost` and `macropru_tax` are **verified** zero, not "
          "merely uncached: `G_D` is constant and absent from the Jacobian, `Phi_D` has no "
          "Jacobian column (the portfolio adjustment cost is quadratic about its anchor, "
          "so its level deviation is second-order), and `T_D` is identically zero at "
          "`T0=T1=0`. The identity closes *because* all three are genuinely zero.", ""]
    return L


def _render_e3(e3, L):
    L += ["## E3 — S-1 writeoff", "", _stamp(e3["provenance"]), "",
          "The two switches answer different questions. `writeoff_enabled` selects which "
          "BRANCH the impulse response traces and is steady-state-neutral: every realised "
          "writeoff term is multiplied by `def_rate_ss = 0`. `zeta_writeoff` governs what "
          "is PRICED — whether a default writes down the perpetuity's continuation value "
          "alongside its coupon. Both are SS-neutral, for the same reason: every writeoff "
          "term is multiplied by `def_rate_ss = 0`, inside `rb_exp` as well as inside "
          "`rb_actual`. So `zeta_writeoff` is allocation-neutral while still changing the "
          "linearised pricing equation, and hence every dynamic result. Since the "
          "2026-08-18 refactor the "
          "baseline is `zeta_writeoff = 1`; `e3b_coupon_only_pricing` is the §12 Arm-3 "
          "diagnostic showing what the pre-refactor coupon-only payoff was worth.", "",
          "| setting | `writeoff_enabled` | `zeta_writeoff` | EL_load_D | peak spread, "
          "passive (bp ann) | loading (medium) | loading (aggressive) |",
          "|---|---|---|---|---|---|---|"]
    flags = {"baseline": ("0", "1.0"), "e3a_realised_writeoff": ("1", "1.0"),
             "e3b_coupon_only_pricing": ("0", "0.0")}
    for name, r in [("baseline", e3["baseline"])] + list(e3["variants"].items()):
        we, ze = flags.get(name, ("?", "?"))
        L.append(f"| {name} | {we} | {ze} | {r['EL_load_D']:.6f} | "
                 f"{r['regimes']['passive']['peak_spread_bp_ann']:.1f} | "
                 f"{_fmt_loading(r['regimes']['medium']['loading'])} | "
                 f"{_fmt_loading(r['regimes']['aggressive']['loading'])} |")

    L += ["", "γ note: solved on the BASELINE and held fixed across variants, so a "
          "difference in the table is attributable to the switch alone. `medium` = 25% "
          "peak-spread compression; **`aggressive` is ≈46.6%, not 50%** — see the E1 note "
          "above.", "",
          f"<sub>Note recorded in the results file at run time: {e3['gamma_note']}</sub>", "",
          "### What the payoff specification is worth", "",
          "`EL_load_D` is the expected loss per unit of default probability implied by "
          "the bond contract. Coupon-only pricing (`zeta_writeoff = 0`) puts it at "
          "`(1-rec)·delta_b/q_b`; full pricing puts it at "
          "`(1-rec)·[delta_b + (1-delta_b)q_b]/q_b`, larger by roughly "
          "`[delta_b + (1-delta_b)q_b]/delta_b ≈ 12.6` on a 12.9-quarter claim. Read the "
          "loading column with that denominator in mind: it is premium income per unit of "
          "expected loss ABSORBED, so a bigger, better-specified loss lowers it "
          "mechanically without the central bank earning any less.", "",
          "| variant | regime | EL PV (% Y) | premium PV (% Y) | loading |",
          "|---|---|---|---|---|"]
    for name, r in [("baseline", e3["baseline"])] + list(e3["variants"].items()):
        for reg, v in r["regimes"].items():
            if v["loading"] is None:
                continue
            L.append(f"| {name} | {reg} | {v['expected_loss_pv_pct_Y']:.5f} | "
                     f"{v['premium_pv_pct_Y']:.5f} | {_fmt_loading(v['loading'])} |")

    L += ["", "### Verification", "",
          "| variant | EL_load (closed form) | EL_load (solved) | ×baseline | max SS drift |",
          "|---|---|---|---|---|"]
    for name, c in e3["checks"].items():
        L.append(f"| {name} | {c['EL_load_expected']:.6f} | {c['EL_load_actual']:.6f} | "
                 f"{c['EL_load_vs_baseline_ratio']:.2f}× | {c['max_ss_drift']:.3e} |")
    L += ["", "Both variants require a full SS + Jacobian re-solve: patching the solved "
          "SS and re-solving only the Jacobian would presume the very invariance these "
          "variants exist to test.", "",
          "> **Both variants must show zero SS drift.** Measured 2026-08-18: `q_b_D = "
          "0.974906` in all three arms. Every writeoff term — priced or realised — is "
          "multiplied by `def_rate_ss = 0`, so neither switch moves an allocation. What "
          "`zeta_writeoff` does move is `EL_load_D` (0.0558 -> 0.7014, 12.6x) and hence "
          "the linearised pricing equation: peak spread on a 1pp shock goes 12.3bp -> "
          "205.9bp. Allocation-neutral, dynamically decisive.", ""]

    comp = e3.get("compression")
    if comp:
        L += ["### Is compression targeting even well-defined?", "",
              "The named regimes are *defined* as 25%/50% peak-spread compression, found "
              "by bisection — which requires peak spread to be monotone in γ.", "",
              "| setting | monotone in γ? | peak @ γ=0 (bp) | min peak (bp) | γ at min | violations |",
              "|---|---|---|---|---|---|"]
        for name, c in comp.items():
            verdict = "yes" if c["monotone_decreasing"] else \
                      f"**NO** (from γ≈{c['first_violation_gamma']:.2f})"
            L.append(f"| {name} | {verdict} | {c['peak_at_gamma0_bp']:.1f} | "
                     f"{c['min_peak_bp']:.1f} | {c['gamma_at_min_peak']:.2f} | "
                     f"{c['n_violations']} |")
        L += ["", "> **Under full writeoff the named-regime construction itself breaks.** "
              "Peak spread stops being monotone in γ, so \"25% compression\" no longer "
              "identifies a unique γ. The violations are two of 39 grid steps: a trivial "
              "one at γ≈0.39 and a large spike at γ≈3.46 (144.4 → 166.6 bp), after which "
              "the curve resumes falling. That isolated spike sits where `I − γ·A_cb` is "
              "plausibly near-singular, so read it as a linear-algebra pathology rather "
              "than economics until confirmed — but compression-targeted regimes cannot "
              "be defined under this setting, which is why every row above is evaluated "
              "at the baseline's γ held fixed.", ""]

    # Read the amplification live: this caption hardcoded "8.5" until 2026-08-06 and
    # silently outlived the sticky-price re-tune to 7.85, asserting a stale number in
    # a GENERATED document. Every parameter quoted in prose must come from provenance.
    L += [f"> `psi_lambda_B = {e3['provenance']['psi_lambda_B_D']:g}` was tuned to 150 bp "
          "with realised losses **off**. The "
          "overshoot above is a reportable fact about whether that target survives S-1 — "
          "**not** a number to re-tune away. Whether to re-tune is a separate author "
          "decision this result informs.", ""]
    return L


def render():
    e1, e2, e3 = (_load("e1_backstop_schedule"), _load("e2_dy_decomposition"),
                  _load("e3_writeoff_s1"))
    L = ["# Policy experiments — standard results", "",
         _stamp(provenance()), "",
         "Generated by `experiments/run_all.py`. **Do not edit by hand** — edit the "
         "experiment and re-run, or the next run discards the edit. Design spec: "
         "`docs/superpowers/specs/2026-08-01-policy-experiments-design.md`.", ""]

    missing = [n for n, v in (("E1", e1), ("E2", e2), ("E3", e3)) if v is None]
    if missing:
        L += [f"> **Incomplete:** no results on disk for {', '.join(missing)}. "
              f"Run `experiments/run_all.py` to populate.", ""]

    if e1:
        L = _render_e1(e1, L)
    if e2:
        L = _render_e2(e2, L)
    if e3:
        L = _render_e3(e3, L)

    with open(DOC, "w") as fh:
        fh.write("\n".join(L) + "\n")
    return DOC


def main():
    import e1_backstop_schedule
    import e2_dy_decomposition

    if "--render-only" not in sys.argv:
        # E2 first: it is self-verifying, so it validates the cache before anything
        # else reports numbers off it.
        e2_dy_decomposition.run()
        e1_backstop_schedule.run()
        if "--skip-e3" not in sys.argv:
            import e3_writeoff_s1
            e3_writeoff_s1.run()
    print(f"Wrote {render()}")


if __name__ == "__main__":
    main()
