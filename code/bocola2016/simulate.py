# ERGODIC SIMULATION AND BOX ADAPTATION FOR THE RESTRICTED MODEL.
# The solved rules are evaluated forward under drawn shocks to trace the
# ergodic set; the state box is then tightened to that cloud so the collocation
# grid spends its resolution where the model actually lives (and drops the
# deep-infeasible corners a BGP-centered box overcovers).
import numpy as np

from period_map import statics
from expectations import next_exog


def simulate(grid, coef, cal, T=20000, burn=2000, seed=0):
    # FORWARD-SIMULATE THE STATE PATH FROM THE BGP UNDER DRAWN SHOCKS.
    rng = np.random.default_rng(seed)
    s = np.array([cal["K_bg"], cal["B_bg"], cal["P_bg"],
                  cal["gamma"], cal["g_star"]], dtype=float)
    path = np.empty((T, 5))
    for t in range(T):
        sc = grid.clip(s)[0]
        C = grid.eval(coef["C"], sc)[0]
        R = grid.eval(coef["R"], sc)[0]
        QB = grid.eval(coef["QB"], sc)[0]
        st = statics(sc[0], sc[1], sc[2], sc[3], sc[4], C, R, QB, cal, d=0.0)
        ez, eg = rng.standard_normal(), rng.standard_normal()
        Dz_n, g_n = next_exog(sc[3], sc[4], np.array([ez]), np.array([eg]), cal)
        s = np.array([st["Kp"], st["Bp"], st["Pp"], Dz_n[0], g_n[0]])
        path[t] = s
    return path[burn:]


def ergodic_box(path, margin=0.10, qlo=0.001, qhi=0.999):
    # BOX COVERING THE ERGODIC CLOUD WITH A RELATIVE MARGIN ON EACH DIMENSION.
    lo = np.quantile(path, qlo, axis=0)
    hi = np.quantile(path, qhi, axis=0)
    span = hi - lo
    return lo - margin * span, hi + margin * span


def adapt_box(solve_fn, grid_fn, cal, lo0, hi0, rounds=4, verbose=True):
    # SOLVE -> SIMULATE -> RE-BOX UNTIL THE ERGODIC CLOUD SITS INSIDE THE BOX.
    # solve_fn(grid) -> sol dict with 'coef'; grid_fn(lo,hi) -> SmolyakGrid.
    lo, hi = np.array(lo0, float), np.array(hi0, float)
    for r in range(1, rounds + 1):
        grid = grid_fn(lo, hi)
        sol = solve_fn(grid)
        path = simulate(grid, sol["coef"], cal)
        lo_e, hi_e = ergodic_box(path)
        # keep the box no wider than needed but always covering the cloud
        lo_new = np.minimum(lo, lo_e)
        hi_new = np.maximum(hi, hi_e)
        # then shrink toward the cloud (the corners we want to drop)
        lo_new = 0.5 * (lo_new + lo_e)
        hi_new = 0.5 * (hi_new + hi_e)
        interior = (path.min(0) > lo_new).all() and (path.max(0) < hi_new).all()
        if verbose:
            frac_bad = (sol["mu"] >= 0.10).mean()
            print(f"  box round {r}: conv={sol['conv']:.2e}  "
                  f"grid_bad(mu>=.1)={frac_bad:.0%}  interior={interior}")
        if np.allclose(lo_new, lo, rtol=0.02) and np.allclose(hi_new, hi, rtol=0.02):
            lo, hi = lo_new, hi_new
            grid = grid_fn(lo, hi)
            sol = solve_fn(grid)
            return grid, sol, (lo, hi), simulate(grid, sol["coef"], cal)
        lo, hi = lo_new, hi_new
    grid = grid_fn(lo, hi)
    sol = solve_fn(grid)
    return grid, sol, (lo, hi), simulate(grid, sol["coef"], cal)
