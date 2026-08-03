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
          f"γ selection: {e1['gamma_selection_rule']}.", "",
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

    L += ["", "> **The headline output number is the residue of two much larger "
          "offsetting channels.** Passive → aggressive, `dY[0]` moves by ~+4.9e−04 while "
          "investment moves ~+2.2e−03 and net exports ~−1.9e−03 — each roughly 4× the "
          "headline and opposite in sign. This confirms `docs/SPEC.md`'s standing caution "
          "as a measured property of this calibration. **Report the decomposition, not "
          "the headline ΔY.**", "",
          "> `government`, `portfolio_cost` and `macropru_tax` are **verified** zero, not "
          "merely uncached: `G_D` is constant and absent from the Jacobian, `Phi_D` has no "
          "Jacobian column (the portfolio adjustment cost is quadratic about its anchor, "
          "so its level deviation is second-order), and `T_D` is identically zero at "
          "`T0=T1=0`. The identity closes *because* all three are genuinely zero.", ""]
    return L


def _render_e3(e3, L):
    L += ["## E3 — S-1 writeoff", "", _stamp(e3["provenance"]), "",
          "`writeoff_enabled` is steady-state-neutral: every writeoff term is multiplied "
          "by `def_rate_ss = 0`. `zeta_writeoff` is **not** — it enters the `EL_price` "
          "anchor at `code/steady_state.py:107-112` *ungated by* `writeoff_enabled`, and "
          "`EL_price` is the loading's denominator. S-1 is therefore two nested variants, "
          "not one switch.", "",
          "| setting | `writeoff_enabled` | `zeta_writeoff` | EL_price_D | peak spread, "
          "passive (bp ann) | loading (medium) | loading (aggressive) |",
          "|---|---|---|---|---|---|---|"]
    flags = {"baseline": ("0", "0.0"), "e3a_coupon_only": ("1", "0.0"),
             "e3b_full": ("1", "1.0")}
    for name, r in [("baseline", e3["baseline"])] + list(e3["variants"].items()):
        we, ze = flags.get(name, ("?", "?"))
        L.append(f"| {name} | {we} | {ze} | {r['EL_price_D']:.6f} | "
                 f"{r['regimes']['passive']['peak_spread_bp_ann']:.1f} | "
                 f"{_fmt_loading(r['regimes']['medium']['loading'])} | "
                 f"{_fmt_loading(r['regimes']['aggressive']['loading'])} |")

    L += ["", "### Verification", "",
          "| variant | EL_price (closed form) | EL_price (solved) | ×baseline | max SS drift |",
          "|---|---|---|---|---|"]
    for name, c in e3["checks"].items():
        L.append(f"| {name} | {c['EL_price_expected']:.6f} | {c['EL_price_actual']:.6f} | "
                 f"{c['EL_price_vs_baseline_ratio']:.2f}× | {c['max_ss_drift']:.3e} |")
    L += ["", "E3a asserts the steady state is **strictly invariant** (drift < 1e−10). "
          "E3b legitimately moves it, so invariance is *not* asserted there — the "
          "closed-form `EL_price` check applies instead. Both variants require a full "
          "SS + Jacobian re-solve: patching the solved SS and re-solving only the "
          "Jacobian would presume the very invariance E3a exists to test.", "",
          "> `psi_lambda_B = 8.5` was tuned to 150 bp with realised losses **off**. Any "
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
