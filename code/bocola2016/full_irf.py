# BOCOLA PASS-THROUGH IRF: A ONE-TIME PRICED SOVEREIGN-RISK SHOCK.
# The risk factor s jumps at t=0 (raising the priced default probability) then
# decays via its AR(1); the economy is traced forward on the NO-DEFAULT path
# (d=0 throughout -- risk priced, never realized, exactly Bocola's experiment).
# Deterministic (eps_z=eps_g=0); the two-branch expectations still price the
# default risk at every date via p^d(s_t). Net worth deteriorates as bond
# prices fall, the constraint tightens, and capital/output respond dynamically
# -- the response the static policy slice cannot show.
import numpy as np

from period_map import statics
from expectations_full import default_prob
from time_iteration_full import point_block_full


def s_path(cal, s_impact, T):
    # s JUMPS TO s_impact AT t=0, THEN REVERTS VIA ITS AR(1) (rho_s).
    s = np.empty(T)
    s[0] = s_impact
    for t in range(1, T):
        s[t] = (1 - cal["rho_s"]) * cal["s_star"] + cal["rho_s"] * s[t - 1]
    return s


def simulate_det(sol, cal, s_series, T):
    # FORWARD PATH ON THE d=0 RULES GIVEN A DETERMINISTIC s SERIES.
    grid, coef, quad = sol["grid"], sol["coef"], sol["quad"]
    K, B, P = cal["K_bg"], cal["B_bg"], cal["P_bg"]
    rec = {k: np.empty(T) for k in
           ("Y", "I", "C", "Q_b", "n", "mu", "lev", "l", "rk_spread", "s", "pd")}
    for t in range(T):
        st6 = np.array([K, B, P, cal["gamma"], cal["g_star"], s_series[t]])
        C = grid.eval(coef["C"][0], st6)[0]
        R = grid.eval(coef["R"][0], st6)[0]
        QB = grid.eval(coef["QB"][0], st6)[0]
        res, mu, av, sm = point_block_full([C, R, QB], st6, 0, grid, coef, cal, quad)
        rec["Y"][t] = sm["Y"]; rec["I"][t] = sm["I"]; rec["C"][t] = C
        rec["Q_b"][t] = QB; rec["n"][t] = sm["N"]; rec["mu"][t] = mu
        rec["lev"][t] = sm["assets"] / sm["N"]; rec["l"][t] = sm["l"]
        rec["rk_spread"][t] = sm["Zrent"] - (R - 1.0)   # rough capital-deposit gap
        rec["s"][t] = s_series[t]; rec["pd"][t] = default_prob(s_series[t])
        K, B, P = sm["Kp"], sm["Bp"], sm["Pp"]
    return rec


def pass_through_irf(sol, cal, pd_target=0.025, T=40):
    # IRF = SHOCKED PATH (s raised to hit pd_target) MINUS BASELINE (s=s*).
    s_impact = np.log(pd_target / (1.0 - pd_target))     # logit
    base = simulate_det(sol, cal, np.full(T, cal["s_star"]), T)
    shock = simulate_det(sol, cal, s_path(cal, s_impact, T), T)
    return base, shock, s_impact


def report_irf(base, shock, cal, s_impact):
    # PRINT THE HEADLINE PASS-THROUGH NUMBERS (% DEV, IMPACT AND PEAK).
    def dev(v):
        # PERCENT DEVIATION OF THE SHOCKED PATH FROM BASELINE.
        return 100.0 * (shock[v] / base[v] - 1.0)
    print(f"  s_impact={s_impact:+.2f}  pd_impact={shock['pd'][0]:.4f} "
          f"(from {base['pd'][0]:.4f})")
    print(f"  {'quarter':>7} {'Q_b%':>8} {'n%':>8} {'I%':>8} {'Y%':>8} "
          f"{'C%':>8} {'mu':>8}")
    for t in [0, 1, 2, 4, 8, 16, 30]:
        print(f"  {t:>7} {dev('Q_b')[t]:8.2f} {dev('n')[t]:8.2f} {dev('I')[t]:8.2f} "
              f"{dev('Y')[t]:8.3f} {dev('C')[t]:8.3f} {shock['mu'][t]:8.5f}")
    print(f"\n  impact: Q_b {dev('Q_b')[0]:+.1f}%  n {dev('n')[0]:+.1f}%  "
          f"I {dev('I')[0]:+.2f}%  Y {dev('Y')[0]:+.3f}%  C {dev('C')[0]:+.3f}%")
    print(f"  trough: I {dev('I').min():+.2f}%  Y {dev('Y').min():+.3f}%  "
          f"n {dev('n').min():+.1f}%  mu peak {shock['mu'].max():.5f}")
    print("  Bocola Fig.5 (pd 2.5%): I ~ -2%, Y ~ -0.3%, C ~ +0.3% on impact")
