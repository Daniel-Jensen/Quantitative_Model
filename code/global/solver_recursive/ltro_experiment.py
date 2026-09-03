# THE LTRO BACKSTOP: A STOCHASTIC FACILITY THAT STABILISES BY BEING BELIEVED.
# With per-period probability phi (cal["phi_ltro"]) the central bank offers
# collateralised credit of size cal["ltro_D"]. It is Bocola's own instrument
# (residual_model_ltro_firstperiod.m): CB funding both LEAVES the divertable base and
# COUNTS as equity in the incentive constraint,
#     mu_ratio = N'/(lambda*A')  ->  (N' + m)/(lambda*(A' - m)),
# which to first order is (1 + leverage) = 6x the constraint relief of a bond purchase of
# the same size. Lent at the deposit rate it moves NO budget identity -- it changes the
# COMPOSITION of the bank's funding, not its size or its cost -- so the whole effect is
# in mu (test_recursive_nesting N4 asserts exactly that).
#
# THE HEADLINE READ IS THE NEVER-FIRED PATH. Agents price the facility whether or not it
# is drawn, so the regime-(d=0, m=0) realisation -- announced, never used -- is where the
# OMT fact lives: no euros spent, spreads compressed anyway. run() reports that path
# against phi = 0, and alongside it the two channels that decide the sign:
#   RISK PREMIUM (stabilising, indirect): the facility lowers Om' most in the states
#     where the constraint is tightest, which are the states where payD is lowest, so
#     cov(Om, payD) shrinks and today's price rises EVERYWHERE, CB present or not.
#   FRANCHISE VALUE (destabilising, indirect): a looser future lowers alpha', hence
#     E[Om], which RAISES today's mu = max(1 - E[Om]R n/lev, 0). The bank's charter value
#     IS its collateral. This channel is not second-order -- mu is a small difference of
#     numbers near one -- and whether stabilisation beats it is the experiment.
# See docs/ltro_backstop_plan.md. Serial pointwise solves; __main__ guard kept regardless.
import numpy as np

from config.calibration import get_calibration
from config.steady_state import solve_steady_state
from solver_recursive.state_grid import s_process_params
from solver_recursive.recursive_main import calibrate_household_anchors
from solver_recursive.recursive_experiment import (solve_recursive, dynamic_irf,
                                                   stochastic_rest_point, read_at,
                                                   s_from_pd, liquidity_ceiling_report)
from solver_recursive.output_decomposition import (simulate, s_decay_path,
                                                   decompose_bond_price, BOND_CHANNELS)
from solver_recursive.state_grid import IS
from solver_recursive.accuracy import accuracy_report
from reporting.plots import ACTIVATION_RAMP
from reporting.prints import bp_ann, print_sovereign_spread
from solver_recursive.output_decomposition import sovereign_spread_legs, SOVEREIGN_LEGS

DEFAULT_ACTIVATIONS = (0.0, 0.5, 1.0)      # per-period activation probabilities phi
# backstop strength is an ORDERED variable, so it takes the sequential ramp plots.py
# defines for it (one hue light->dark), NOT a categorical palette
_LTRO_COLORS = ACTIVATION_RAMP
# The shock is stated as a TARGET one-quarter-ahead default probability, the same units
# main.py's RISK_SHOCK_PD uses, so this overlay and the headline IRF are the SAME shock.
PD_SHOCK, T_IRF, DECOMP_T = 0.0198, 21, 25


def _specs(activations):
    # MAP ACTIVATION PROBABILITIES TO (value, label, color) TRIPLES FOR SOLVE + PLOT.
    # Backstop strength is an ORDERED variable, so the colour must be ordered too. The
    # three-stop ramp is INTERPOLATED rather than cycled: at eight activations cycling
    # would give three activations the same colour and the overlay would read as three
    # groups instead of one gradient.
    import matplotlib.colors as mcolors
    cmap = mcolors.LinearSegmentedColormap.from_list("ltro", _LTRO_COLORS)
    n = max(len(activations) - 1, 1)
    out = []
    for i, a in enumerate(activations):
        lab = "0% — no backstop" if a == 0.0 else f"{round(a * 100):.0f}% chance"
        out.append((a, lab, cmap(i / n)))
    return out


def irf_series(rules, cal, ss, sproc, pd_shock=PD_SHOCK, T=T_IRF):
    # THE NEVER-FIRED IRF, on the same footing as the headline one.
    # A thin adapter over dynamic_irf rather than its own loop: dynamic_irf starts at the
    # model's own stochastic rest point, advances an unshocked reference path with the
    # shared law of motion, differences quarter by quarter and reports box escapes --
    # Bocola's generate_irf.m, both halves. It reads regime 0, which under the four-regime
    # table is (d=0, m=0): no default and NO FACILITY DRAWN. That is exactly the path this
    # experiment is about, and it is why no new simulation machinery is needed.
    P = dynamic_irf(rules, cal, ss, sproc, pd_shock=pd_shock, T=T, rest_verbose=False)
    return dict(pd=P["pd"], Y_D=P["Y"], C_D=P["C"], I_D=P["I"], K_D=P["K"], n_D=P["n"],
                Q_bD=P["dQ_bD"], spread=P["d_spread"], sov_bp=P["sov_bp"],
                m_ltro=P["m_ltro"], mu=P["mu"], E_Om=P["E_Om"])


def rest_point_diagnostics(rules, cal, ss, sproc):
    # E3: THE FRANCHISE-VALUE COUNTER-TEST, read at the model's own rest point.
    # If the facility stabilises, mu falls with phi. If the charter-value channel wins,
    # mu RISES with phi even though no facility is ever drawn on this path -- and that is
    # a result about standing liquidity backstops, not a failed run. Reporting E[Om] and
    # alpha beside mu is what separates the two: both channels move mu, but only the
    # franchise one moves it by moving E[Om].
    S = stochastic_rest_point(rules, cal, ss, sproc, verbose=False)
    o = read_at(rules, cal, ss, sproc, S.copy())
    Sk = S.copy(); Sk[IS] = s_from_pd(PD_SHOCK)
    ok_ = read_at(rules, cal, ss, sproc, Sk)
    # lev_IC is lambda*(assets - m), so this ratio is DIVERTABLE assets per unit of net
    # worth -- the quantity the constraint is actually about, not the accounting leverage
    yD = cal["delta_b_D"] * (1.0 - o["Q_bD"]) / o["Q_bD"]
    yF = cal["delta_b_F"] * (1.0 - o["Q_bF"]) / o["Q_bF"]
    return dict(mu=o["mu_D"], E_Om=o["E_Om_D"], alpha=o["alpha_D"],
                slack=o["slack_D"], Q_bD=o["Q_bD"], n=o["n_D"],
                lev_div=o["lev_IC_D"] / (cal["lambda_K_D"] * o["n_D"]),
                spread_bp=bp_ann(cal["lambda_K_D"] * o["mu_D"] / o["E_Om_D"]),
                r_wc_bp=bp_ann(o["r_wc_D"]), sov_bp=bp_ann(yD - yF),
                Y=o["Y_D"], C=o["C_D"], I=o["I_D"],
                legs_rest=sovereign_spread_legs(o, cal),
                legs_shock=sovereign_spread_legs(ok_, cal))


def run(cal, ss, sproc, mu=1, activations=None, mu_vec=None, pd_shock=PD_SHOCK,
        accuracy=True, accuracy_T=600, s_refine=None, decompose=True):
    # SOLVE THE ACTIVATION SCENARIOS AND WRITE THE OVERLAY FIGURE (shared SS).
    from reporting.plots import (plot_activation_irf, plot_certainty_curve,
                                 OUTDIR)
    import os, time
    os.makedirs(OUTDIR, exist_ok=True)
    specs = _specs(DEFAULT_ACTIVATIONS if activations is None else activations)
    print("=== LTRO backstop — the same shock under several activation probabilities ===",
          flush=True)
    print(f"    facility {100 * cal['ltro_D'] / ss['ss_firm_D']['Y_ss']:.1f}% of quarterly"
          f" GDP to D, {100 * cal['ltro_F'] / ss['ss_firm_D']['Y_ss']:.1f}% to F;"
          f" read along the NEVER-FIRED path", flush=True)
    phi0 = cal.get("phi_ltro", 0.0)
    t0 = time.perf_counter()
    scenarios, rest, decs = [], [], []
    # THE BASELINE IS SOLVED ONCE AND REUSED. d0, the haircut homotopy and the
    # no-facility joint solve are phi-INDEPENDENT; re-solving them per activation costs
    # ~7 min a point on the coarse grid and buys nothing. Every activation, including
    # phi = 0, then runs on the SAME grid and the SAME regime count, so a kink at the
    # left-hand end of the curve is economics and not a change of configuration.
    base_out, solved_ok = [], []
    for a, label, color in specs:
        cal["phi_ltro"] = a
        print(f"  --- solving [{label}] ---", flush=True)
        kw = {} if s_refine is None else dict(s_refine=s_refine)
        rules = solve_recursive(cal, ss, sproc, mu=mu, verbose=True, mu_vec=mu_vec,
                                with_cb=True,
                                base=(base_out[0] if base_out else None),
                                base_out=(None if base_out else base_out), **kw)
        solved_ok.append(bool(getattr(rules, "solve_ok", True)))
        P = irf_series(rules, cal, ss, sproc, pd_shock=pd_shock)
        rest.append(rest_point_diagnostics(rules, cal, ss, sproc))
        if decompose:
            S0 = stochastic_rest_point(rules, cal, ss, sproc, verbose=False)
            s_path = s_decay_path(sproc, s_from_pd(pd_shock), DECOMP_T)
            sim = simulate(rules, cal, ss, sproc, s_path, S_init=S0)
            ref = simulate(rules, cal, ss, sproc,
                           np.full(DECOMP_T, sproc["s_star"]), S_init=S0)
            decs.append(decompose_bond_price(sim, ref, cal))
        if accuracy:
            accuracy_report(rules, cal, ss, sproc,
                            stochastic_rest_point(rules, cal, ss, sproc, verbose=False),
                            T=accuracy_T, label=f"LTRO {label}")
        scenarios.append((label, P, color))
        print(f"  [{label}] solved ({time.perf_counter() - t0:.0f}s)  "
              f"impact Y_D={P['Y_D'][0]:+.4f}%  Q_bD={P['Q_bD'][0]:+.3f}%  "
              f"sov spread={P['sov_bp'][0]:+.0f}bp  drawn={P['m_ltro'][0]:.2f}% of GDP",
              flush=True)

    plot_activation_irf(scenarios)

    # E1 -- the headline table
    print("\n  E1  IMPACT (t=0) ON THE NEVER-FIRED PATH, by activation probability",
          flush=True)
    print("      (the facility is announced and NOT drawn: 'drawn' must be 0.00)",
          flush=True)
    print("   scenario                  Y_D%     C_D%    I_D%   Q_bD%   sov_bp  drawn%",
          flush=True)
    for label, P, _ in scenarios:
        print(f"   {label:22s} {P['Y_D'][0]:+8.4f} {P['C_D'][0]:+8.4f} "
              f"{P['I_D'][0]:+7.4f} {P['Q_bD'][0]:+7.3f} {P['sov_bp'][0]:+8.0f} "
              f"{P['m_ltro'][0]:7.2f}", flush=True)

    # E3 -- the two channels, at the rest point
    print("\n  E3  REST POINT vs phi: does the constraint LOOSEN or TIGHTEN when nothing"
          " is drawn?", flush=True)
    print("   scenario                    mu_D      E[Om_D]    alpha_D    lev_div"
          "  spread_bp", flush=True)
    for (label, _, _), r, okf in zip(scenarios, rest, solved_ok):
        print(f"   {label:22s} {r['mu']:10.6f} {r['E_Om']:11.6f} {r['alpha']:10.5f} "
              f"{r['lev_div']:9.3f} {r['spread_bp']:10.1f}"
              f"{'' if okf else '   <- DID NOT REACH THE ACCEPTANCE FLOOR'}", flush=True)
    if len(rest) > 1:
        d_mu = rest[-1]["mu"] - rest[0]["mu"]
        verdict = ("the CHARTER-VALUE channel dominates: a looser future TIGHTENS the "
                   "constraint today" if d_mu > 0 else
                   "the RELIEF channels dominate: the constraint is looser even with "
                   "nothing drawn")
        print(f"   -> d(mu)/d(phi) = {d_mu:+.6f} over the range; {verdict}.", flush=True)
        # THE LEVEL SHIFT, WHICH THE IRF CANNOT SHOW. dynamic_irf differences each phi's
        # shocked path against ITS OWN unshocked path, so it reports the response
        # CONDITIONAL on the policy regime and silently removes the shift in the ergodic
        # point -- which is most of what "the economy is stabilised" means here. Reading
        # every rest point against the phi = 0 one puts that shift back.
        b = rest[0]
        print("\n  E3b WHERE THE ECONOMY RESTS, against the no-backstop rest point"
              " (nothing ever drawn)", flush=True)
        print("   scenario                   Y_D%     C_D%     I_D%    Q_bD%"
              "   sov_bp   credit_bp", flush=True)
        for (label, _, _), r in zip(scenarios, rest):
            print(f"   {label:22s} {100*(r['Y']/b['Y']-1):+8.4f} "
                  f"{100*(r['C']/b['C']-1):+8.4f} {100*(r['I']/b['I']-1):+8.4f} "
                  f"{100*(r['Q_bD']/b['Q_bD']-1):+8.3f} "
                  f"{r['sov_bp']-b['sov_bp']:+8.1f} {r['spread_bp']-b['spread_bp']:+11.1f}",
                  flush=True)

    # THE SPREAD, LEG BY LEG. The policy targets the sovereign spread, so the table
    # that matters is which part of that spread each instrument can reach.
    for (label, _, _), r in zip(scenarios, rest):
        print_sovereign_spread(r["legs_shock"], f"{label}, at the shock")
    print("\n  WHAT THE FACILITY MOVES, leg by leg (annualised bp, vs the first scenario)")
    print(f"   {'scenario':<22s}" + "".join(f"{nm[:13]:>15s}" for _, nm in SOVEREIGN_LEGS)
          + f"{'TOTAL':>10s}")
    b = rest[0]["legs_shock"]
    for (label, _, _), r in zip(scenarios, rest):
        L = r["legs_shock"]
        print(f"   {label:<22s}"
              + "".join(f"{L[f'spread_{k}'] - b[f'spread_{k}']:15.2f}"
                        for k, _ in SOVEREIGN_LEGS)
              + f"{L['spread'] - b['spread']:10.2f}")

    # E2 -- which leg of the bond price moves
    if decompose:
        print("\n  E2  BOND-PRICE DECOMPOSITION of the SAME shock, impact leg (log %)",
              flush=True)
        names = [k for k, _ in BOND_CHANNELS]
        print("   scenario              " + "".join(f"{k[:12]:>14s}" for k in names),
              flush=True)
        for (label, _, _), d in zip(scenarios, decs):
            print(f"   {label:20s}" + "".join(f"{d[k][0]:14.4f}" for k in names),
                  flush=True)
        print("   -> the CLAIM is that the risk-premium leg carries the compression on"
              " the never-fired\n      path. If the liquidity leg carries it instead,"
              " the mechanism is direct relief,\n      not credibility.", flush=True)

    # THE CERTAINTY CURVE: the ergodic point AGAINST the announced probability. The
    # overlay figure shows responses over time at each phi; this one shows what the
    # announcement itself buys, which is the object of the experiment.
    if len(rest) > 1:
        b = rest[0]
        curve = dict(
            mu=[r["mu"] for r in rest],
            cred=[r["spread_bp"] for r in rest],
            sov=[r["sov_bp"] - b["sov_bp"] for r in rest],
            Q=[100 * (r["Q_bD"] / b["Q_bD"] - 1) for r in rest],
            Y=[100 * (r["Y"] / b["Y"] - 1) for r in rest],
            I=[100 * (r["I"] / b["I"] - 1) for r in rest])
        print(f"  figure -> "
              f"{plot_certainty_curve([a for a, _, _ in specs], curve, solved_ok)}",
              flush=True)
    print(f"  figure -> {OUTDIR}/ltro_activation.png", flush=True)
    cal["phi_ltro"] = phi0
    return scenarios, rest, decs


def main():
    cal = get_calibration()
    # match main.py: the net-worth floor is what keeps the deep default corners feasible,
    # so a standalone run without it is solving a different model from the pipeline
    cal["nw_floor_frac"] = 0.15
    ss = solve_steady_state(cal, verbose=False)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    run(cal, ss, sproc)


if __name__ == "__main__":
    main()
