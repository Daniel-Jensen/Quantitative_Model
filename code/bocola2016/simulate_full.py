# FULL-MODEL STOCHASTIC SIMULATOR: EVOLVES BOTH REGIMES WITH DEFAULT DRAWS.
# Draws (eps_z, eps_g, eps_s) and d' ~ Bernoulli(p^d(s)) each period, switching
# the active regime's rules. Used to trace the ergodic set of the endogenous
# states (for the PCA-rotated grid) and for IRFs with realized default.
import numpy as np

from period_map import statics
from expectations_full import next_exog_full, default_prob


def simulate_full(sol, cal, T=8000, burn=1500, seed=0):
    # FORWARD PATH OF THE 6 STATES UNDER DRAWN SHOCKS AND DEFAULT REALIZATIONS.
    grid, coef = sol["grid"], sol["coef"]
    rng = np.random.default_rng(seed)
    K, B, P = cal["K_bg"], cal["B_bg"], cal["P_bg"]
    Dz, g, s = cal["gamma"], cal["g_star"], cal["s_star"]
    d_cur = 0
    path = np.empty((T, 6))
    for t in range(T):
        st6 = grid.clip(np.array([K, B, P, Dz, g, s]))[0]
        C = grid.eval(coef["C"][d_cur], st6)[0]
        R = grid.eval(coef["R"][d_cur], st6)[0]
        QB = grid.eval(coef["QB"][d_cur], st6)[0]
        sm = statics(st6[0], st6[1], st6[2], st6[3], st6[4], C, R, QB, cal,
                     d=float(d_cur))
        ez, eg, es = rng.standard_normal(3)
        Dz, g, s = next_exog_full(st6[3], st6[4], st6[5],
                                  np.array([ez]), np.array([eg]),
                                  np.array([es]), cal)
        Dz, g, s = Dz[0], g[0], s[0]
        d_cur = int(rng.random() < default_prob(st6[5]))
        K, B, P = sm["Kp"], sm["Bp"], sm["Pp"]
        path[t] = [K, B, P, Dz, g, s]
    return path[burn:]
