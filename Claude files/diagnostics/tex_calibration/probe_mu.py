# PROBE (2026-10-27): mu, z AND alpha FOR BOTH BANKS AT EVERY NODE OF A SAVED d0 SOLUTION; NODES WHERE mu NEARS f
# MAKE THE FRANCHISE FIXED POINT alpha = f/(f - mu) DIVERGE.   python3 -u probe_mu.py <ckpt>
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import numpy as np


def main():
    import global_projection.solver_recursive.point_map as pm
    from global_projection.solver_recursive.decision_rules import SOLVE
    from global_projection.solver_recursive.state_grid import s_process_params, STATE_NAMES
    sv = pickle.load(open(HERE / "runs" / sys.argv[1], "rb"))
    rules, cal, ss = sv["rules"], sv["cal"], sv["ss"]
    sp = s_process_params(cal)
    u = rules.grid.points_unit
    print(f"f_D {cal['f_D']}  mu_eps {cal['mu_eps']}  lambda {cal['lambda_K_D']:.4f}")
    print(" node  axis          mu_D     z_D     alpha_D  n_D/n_ss   mu_F    alpha_F  n_F/n_ss")
    for i, S in enumerate(rules.grid.points):
        x = np.array([rules.vals[k][0][i] for k in SOLVE])
        v = pm._read_state_and_regime(S, 0, x, rules, cal, ss)
        pm._production(v, cal); pm._government(v, cal, ss); pm._balance_sheets(v, rules, cal)
        pm._continuation(v, rules, cal, sp, rules.n_gh, True)
        pm._discount_kernels(v, rules, cal); pm._incentive_constraint(v, cal)
        zD = 1.0 - v.E_Om_D * (1.0 + v.rdep_D) * v.n_D / v.lev_D
        ax = ",".join(f"{a}{'+' if u[i, j] > 0 else '-'}" for j, a in enumerate(STATE_NAMES) if abs(u[i, j]) > 1e-9) or "centre"
        print(f"  {i:2d}  {ax:12s} {v.mu_D:7.4f} {zD:+8.4f} {v.alpha_D_new:8.3f} {v.n_D / ss['ss_bank_D']['n_ss']:8.3f}"
              f"  {v.mu_F:7.4f} {v.alpha_F_new:8.3f} {v.n_F / ss['ss_bank_F']['n_ss']:8.3f}")


if __name__ == "__main__":
    main()
