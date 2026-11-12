# T1.4 (2026-10-26): A CONTINUATION PATH FOR THE d0 STAGE (pi = 0, eps 1e-2, SS CONSISTENT WITH IT, 15% FLOOR)
# FROM THE PRE-2026-10-26 CALIBRATION TO APPENDIX B, GROUP BY GROUP, EACH RUNG SEEDED WITH THE LAST
# SOLUTION'S DEVIATIONS FROM ITS OWN SS (IN THE RULES' FIT SPACE), A FAILED RUNG HALVED DOWN TO 1/16.
#   python3 -u "Claude files/diagnostics/tex_calibration/diag_path.py" > runs/diag_path.log
# Writes runs/diag_path_<group>.pkl after each completed group (rules, cal, ss).
import os
import pickle
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)
sys.path.insert(0, str(HERE))

import numpy as np

from diag_homotopy import OLD, economy

# the order: the household block first (Frisch 0.5 restores the GHH composite that G shrinks), the
# hard moves (G, Bohn, home bias) after it, the wedge last
PATH = (("households", ("sigma_D", "sigma_F", "frisch_D", "frisch_F")),
        ("markup", ("epsilon_D", "epsilon_F")),
        ("size+eta", ("size_F", "epsilon_trade")),
        ("foreign share", ("foreign_share_D",)),
        # after 2026-10-26's first pass broke inside delta_b_F once home bias was in, the fiscal groups go
        # first and home bias after them
        ("delta_b_F", ("delta_b_F",)),
        ("bohn", ("bohn_annual",)),
        ("G", ("G_D", "G_F")),
        # home bias last: from the G endpoint even 1/8 of its move fails (2026-10-27)
        ("banker beta", ("beta_inter_D", "beta_inter_F")),
        ("wedge", ("credit_spread_target_D", "credit_spread_target_F")),
        ("home bias", ("omega_home_D",)))
# the funded floor of the run (2026-10-26 decision: 50%); d0 does not reach it, the default regime does
FLOOR = 0.5


def relative_seed(prev, cal_p, ss_p, cal_k, ss_k):
    # THE NEW ECONOMY'S RULES: ITS SS PLUS THE LAST SOLUTION'S DEVIATION FROM ITS SS, NODE BY NODE.
    from global_projection.solver_recursive.decision_rules import RuleSet, STORE_RULES, to_fit, from_fit
    from global_projection.solver_recursive import recursive_experiment as RE
    ss_old = RuleSet.from_ss(prev.grid, ss_p, cal_p, n_regimes=2)
    new = RuleSet.from_ss(RE.normal_box(ss_k, cal_k, mu=1), ss_k, cal_k, n_regimes=2)
    for k in STORE_RULES:
        for d in range(new.n_regimes):
            dev = to_fit(k, prev.vals[k][d]) - to_fit(k, ss_old.vals[k][d])
            new.set_values(k, d, from_fit(k, to_fit(k, new.vals[k][d]) + dev))
    new.n_gh = prev.n_gh
    return new


# a rung whose Newton stalls at max|F| <= LOOSE (above the solver's 1e-8 stall acceptance) is kept and
# flagged: at appendix home bias the d0 system floors at 6.03e-8 from two different rungs
LOOSE = float(os.environ.get("PATH_LOOSE", "1e-7"))


def main():
    # THE PATH, OPTIONALLY RESUMED FROM A CHECKPOINT: python3 diag_path.py [checkpoint.pkl first_group]
    from calibration.global_projection import get_calibration
    from global_projection.solver_recursive.decision_rules import RuleSet
    from global_projection.solver_recursive import recursive_experiment as RE
    base = get_calibration()
    base.update(nw_floor_frac=FLOOR, mu_eps=1e-2, book_eps=1e-4, floor_eps=1e-3)
    new = {k: float(base[k]) for k in OLD}
    names = [g[0] for g in PATH]
    if len(sys.argv) > 2:
        sv = pickle.load(open(HERE / "runs" / sys.argv[1], "rb"))
        rules, cal, ss, cur = sv["rules"], sv["cal"], sv["ss"], sv["vals"]
        first = names.index(sys.argv[2])
        print(f"[resume] from {sys.argv[1]} at group {sys.argv[2]}", flush=True)
    else:
        cur, first = dict(OLD), 0
        cal, ss, sp = economy(base, cur)
        rules = RuleSet.from_ss(RE.normal_box(ss, cal, mu=1), ss, cal, n_regimes=2)
        rules.n_gh = RE.N_GH
        ok, it, w = RE._stage_wc(rules, cal, ss, sp, (0,), True, "d0 old calibration", False)
        print(f"[start] old calibration: ok {ok}  its {it}  max|F| {w:.2e}", flush=True)
        if not ok:
            return
    loose = []
    for name, keys in PATH[first:]:
        t0 = (cur[keys[0]] - OLD[keys[0]]) / (new[keys[0]] - OLD[keys[0]]) if new[keys[0]] != OLD[keys[0]] else 1.0
        t, h = float(np.clip(t0, 0.0, 1.0)), 1.0
        while t < 1.0 - 1e-12:
            h = min(h, 1.0 - t)  # a halved step must stay short of the end, or it re-tries the end
            nxt = min(1.0, t + h)
            vals = dict(cur, **{k: OLD[k] + nxt * (new[k] - OLD[k]) for k in keys})
            cal_k, ss_k, sp_k = economy(base, vals)
            trial = relative_seed(rules, cal, ss, cal_k, ss_k)
            ok, it, w = RE._stage(trial, cal_k, ss_k, sp_k, (0,), True, f"{name} {nxt:.4f}", False,
                                  warm=0, maxit=30)
            flag = "" if ok else ("  LOOSE (kept)" if w <= LOOSE else "")
            print(f"[{name}] fraction {nxt:.4f}: ok {ok}  its {it}  max|F| {w:.2e}{flag}", flush=True)
            if ok or w <= LOOSE:
                if not ok:
                    loose.append((name, nxt, w))
                rules, cal, ss, t, cur = trial, cal_k, ss_k, nxt, vals
                h = min(1.0, 2 * h)
                with open(HERE / "runs" / "diag_path_ckpt.pkl", "wb") as fh:
                    pickle.dump(dict(rules=rules, cal=cal, ss=ss, vals=cur), fh)
                continue
            h *= 0.5
            if h < 1 / 16:
                print(f"[{name}] BREAKS at fraction {t:.4f} of its move", flush=True)
                return
        print(f"[{name}] done", flush=True)
        with open(HERE / "runs" / f"diag_path_{name.replace(' ', '_').replace('+', '_')}.pkl", "wb") as fh:
            pickle.dump(dict(rules=rules, cal=cal, ss=ss, vals=cur), fh)
    print(f"[end] the d0 stage reached appendix B along the path; loose rungs {loose}", flush=True)


if __name__ == "__main__":
    main()
