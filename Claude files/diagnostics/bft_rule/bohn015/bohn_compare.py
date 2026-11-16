# BOHN 0.15 AGAINST 0.25 ON THE SAVED, CONVERGED NO-POLICY RULES OF THE 2026-11-06 CONTINUATION (12% BOX,
# RECAPITALISATION FLOOR 0.5). NO GLOBAL SOLVE: FITTED READS, EXACT ONE-QUARTER AND JOINT TWO-QUARTER CLEARING, THE
# EXACT REST POINT, EXACT DYNAMIC IRFS, THE SURPRISE 200 BP CAP, AND OFF-GRID ACCURACY ON A SIMULATED ERGODIC SET.
#   python3 -u bohn_compare.py > <out>/bohn_compare.log
import contextlib, copy, io, json, pickle, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import root, brentq

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from calibration.global_projection import get_calibration
from global_projection.steady_state import solve_steady_state
from global_projection.solver_recursive.state_grid import s_process_params, default_prob, IS, IBDD, IP_D, IK_D, IV
from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
from global_projection.solver_recursive.decision_rules import SOLVE7
from global_projection.solver_recursive import point_map as PM
from global_projection.solver_recursive import recursive_experiment as RE
from global_projection.solver_recursive import tpi_experiment as TE
from global_projection.solver_recursive import accuracy as AC

RUNS = ROOT / "Claude files" / "diagnostics" / "tpi_root_cause" / "runs"
OUT = ROOT.parent / "QUANTITATIVE_MODEL" / "results" / "GLOBAL_2026-11-12_pr38_bft" / "bohn015"
OUT.mkdir(parents=True, exist_ok=True)
PR38 = ROOT.parent / "QUANTITATIVE_MODEL" / "results" / "GLOBAL_2026-11-12_pr38" / "data"
CONT_KEYS = ("Q_bD", "alpha_D", "C_D", "N_D", "Kp_D", "p", "Q_bF", "alpha_F", "C_F", "N_F", "Kp_F")
SHOCKS = (0.015, 0.06)
t0 = time.time()


def log(*a):
    print(f"[{time.time() - t0:6.0f}s]", *a, flush=True)


class Econ:
    # ONE ECONOMY: THE CALIBRATION AT A BOHN COEFFICIENT, ITS STEADY STATE AND RISK PROCESS, AND ITS SAVED RULES
    def __init__(self, bohn, rules_file):
        base = get_calibration()
        base.update(nw_floor_frac=0.5, jacobian="complex", pd_box_top=0.12, tpi_on=False)
        self.cal = dict(base, bohn_annual=bohn)
        self.ss = solve_steady_state(self.cal, verbose=False)
        self.sp = s_process_params(self.cal)
        calibrate_household_anchors(self.cal, self.ss, self.sp)
        self.rules = pickle.load(open(rules_file, "rb"))
        for tag in ("_rest_point_exact", "_rest_point"):
            if hasattr(self.rules, tag):
                delattr(self.rules, tag)
        self.ngh = self.rules.n_gh or RE.N_GH
        self.bohn = bohn


class Proxy:
    def __init__(self, base):
        self._b = base

    def eval(self, name, d, x):
        return self._b.eval(name, d, x)

    def __getattr__(self, k):
        return getattr(self._b, k)


class ContProxy(Proxy):
    # continuation values replaced by supplied values at the stored next states
    def __init__(self, base, over):
        super().__init__(base)
        self._c = over

    def eval_cont(self, j, Sn):
        r = self._b.eval_cont(j, Sn)
        if j in self._c:
            r = dict(r)
            for k, val in self._c[j].items():
                r[k] = np.array(val, dtype=float)
        return r


def run_map(E, S, d, x, cont=None):
    cont = E.rules if cont is None else cont
    cal, ss = E.cal, E.ss
    v = PM._read_state_and_regime(S, d, x, cont, cal, ss)
    PM._production(v, cal); PM._government(v, cal, ss); PM._balance_sheets(v, cont, cal)
    PM._continuation(v, cont, cal, E.sp, E.ngh, False); PM._discount_kernels(v, cont, cal)
    PM._incentive_constraint(v, cal); PM._bond_demands(v, cal); PM._tpi_rule(v, cal, ss)
    PM._wages_and_dividends(v, cal); PM._households(v, cal, ss)
    v.res = PM._residual_vector(v, cal, ss)
    v.x = np.array(x, dtype=float)
    return v


def x_fit(E, S, d=0):
    return np.array([float(E.rules.eval_read(k, d, np.atleast_2d(S))[0]) for k in SOLVE7])


def fitted(E, S):
    v = run_map(E, S, 0, x_fit(E, S))
    v.err = float(np.max(np.abs(v.res)))
    return v


def clear(E, S, d=0, cont=None, x0=None):
    cont = E.rules if cont is None else cont
    x0 = x_fit(E, S, d) if x0 is None else x0
    f = lambda x: PM.point_residuals(S, d, x, cont, E.cal, E.ss, E.sp, n_gh=E.ngh)[0]
    sol = root(f, x0, method="hybr", tol=1e-13)
    v = run_map(E, S, d, sol.x, cont)
    v.err = float(np.max(np.abs(v.res)))
    return v


def two_quarter(E, S, x_today):
    # TODAY AND EVERY NEXT-QUARTER QUADRATURE STATE (BOTH REGIMES) CLEARED JOINTLY; THE QUARTER AFTER FROM THE RULES
    eps, _ = PM.gh_nodes(E.ngh)
    nreg, ncell, nx = E.rules.n_regimes, len(eps), len(x_today)

    def states(x):
        v = PM._read_state_and_regime(S, 0, x, E.rules, E.cal, E.ss)
        PM._production(v, E.cal); PM._government(v, E.cal, E.ss); PM._balance_sheets(v, E.rules, E.cal)
        return PM._next_states(v, E.sp, eps)[0]

    Sn0 = states(x_today)
    z0 = [x_today] + [clear(E, Sn0[i], d=j).x for j in range(nreg) for i in range(ncell)]
    z0 = np.concatenate(z0)

    def F(z, want=False):
        x, xn = z[:nx], z[nx:].reshape(nreg, ncell, nx)
        Sn = states(x)
        over, res_n = {j: {k: [] for k in CONT_KEYS} for j in range(nreg)}, []
        for j in range(nreg):
            for i in range(ncell):
                r, o = PM.point_residuals(Sn[i], j, xn[j, i], E.rules, E.cal, E.ss, E.sp, n_gh=E.ngh)
                res_n.append(r)
                vals = dict(Q_bD=o["Q_bD"], alpha_D=o["alpha_D"], C_D=o["C_D"], N_D=xn[j, i][0], Kp_D=xn[j, i][2],
                            p=xn[j, i][6], Q_bF=o["Q_bF"], alpha_F=o["alpha_F"], C_F=o["C_F"], N_F=xn[j, i][1],
                            Kp_F=xn[j, i][3])
                for k in CONT_KEYS:
                    over[j][k].append(vals[k])
        cp = ContProxy(E.rules, over)
        if want:
            return run_map(E, S, 0, x, cp), over
        return np.concatenate([PM.point_residuals(S, 0, x, cp, E.cal, E.ss, E.sp, n_gh=E.ngh)[0]] + res_n)

    sol = root(F, z0, method="hybr", tol=1e-12, options=dict(maxfev=200000))
    v, over = F(sol.x, want=True)
    v.err = float(np.max(np.abs(F(sol.x))))
    v.over = over
    return v


def summ(E, v):
    cal = E.cal
    f, lK, db = cal["f_D"], cal["lambda_K_D"], cal["delta_b_D"]
    psi = 1 - v.E_Om_D * (1 + v.rdep_D) * v.n_D / v.lev_D
    Y_ss = E.ss["ss_firm_D"]["Y_ss"]
    return dict(
        Q_bD=float(v.Q_bD), sov_bp=float(4e4 * (db * (1 / v.Q_bD - 1) - cal["delta_b_F"] * (1 / v.Q_bF - 1))),
        ls_bp=float(4e4 * lK * v.mu_D / v.E_Om_D), mu=float(v.mu_D), psi=float(psi),
        n=float(v.n_D), lev=float(v.lev_D), EOm=float(v.E_Om_D), R=float(1 + v.rdep_D), rdep_bp=float(4e4 * v.rdep_D),
        alpha=float(v.alpha_D_new), Y=float(v.Y_D), I=float(v.I_D), C=float(v.C_D), QK=float(v.Q_D),
        Kp=float(v.Kp_D), K=float(v.K_D),
        B=float(v.B_D), Bp=float(v.Bp_D), debt_gdp=float(v.Bp_D / (4 * Y_ss)), new=float(v.new_D),
        Tax=float(v.Tax_D), b_DD=float(v.b_D_D_new), b_lag=float(v.b_D_D_lag), P=float(v.P_D),
        exposure=float(v.Q_bD * v.b_D_D_new / v.assets_D), assets=float(v.assets_D),
        lev_ratio=float(v.assets_D / v.n_D),
        n_sov=float((1 - f) * v.Xi_D * v.b_D_D_lag),
        n_cap=float((1 - f) * (v.mpk_D + (1 - cal["delta_D"]) * v.Q_D) * v.K_D),
        n_dep=float(-(1 - f) * v.P_D), n_ent=float(cal["omega_ent_D"] * v.assets_D),
        lev_K=float(lK * (v.Q_D * v.Kp_D + v.L_wc_D)), lev_b=float(cal["lambda_bD_D"] * v.Q_bD * v.b_D_D_new),
        Epay=float(v.E_payD), nfa=float(v.nfa_dep_D / Y_ss), err=float(v.err))


def shocked(S, pd):
    S = S.copy()
    S[IS] = RE.s_from_pd(pd)
    return S


def psi_of(E, S_rest, pd, mode, x0=None):
    S = shocked(S_rest, pd)
    if mode == "fitted":
        v = fitted(E, S)
    else:
        v = clear(E, S, x0=x0)
        if mode == "2q":
            v = two_quarter(E, S, v.x)
    return summ(E, v)["psi"]


def threshold(E, S_rest, mode, lo=0.002, hi=0.06):
    g = lambda pd: psi_of(E, S_rest, pd, mode)
    glo, ghi = g(lo), g(hi)
    if glo > 0:
        return dict(pd=None, note=f"binds already at {lo}", psi_lo=glo)
    if ghi < 0:
        return dict(pd=None, note=f"slack at {hi}", psi_hi=ghi)
    return dict(pd=float(brentq(g, lo, hi, xtol=2e-5)), psi_lo=glo, psi_hi=ghi)


def cont_accuracy(E, v):
    # FITTED CONTINUATION AGAINST EXACT CLEARING AT EVERY NEXT-QUARTER QUADRATURE STATE
    eps, _ = PM.gh_nodes(E.ngh)
    Sn, _ = PM._next_states(v, E.sp, eps)
    out = dict(outside_now=float(np.max(E.rules.grid.outside(np.atleast_2d(v.S)))),
               outside_next=float(np.max(E.rules.grid.outside(Sn))))
    for j in range(E.rules.n_regimes):
        fit = E.rules.eval_cont(j, Sn)
        rel = {k: [] for k in ("Q_bD", "alpha_D", "C_D", "Kp_D")}
        for i in range(Sn.shape[0]):
            w = clear(E, Sn[i], d=j)
            ex = dict(Q_bD=w.Q_bD, alpha_D=w.alpha_D_new, C_D=w.C_D, Kp_D=w.Kp_D)
            for k in rel:
                rel[k].append(abs(float(ex[k]) / float(np.real(fit[k][i])) - 1))
        out[f"reg{j}"] = {k: float(max(val)) for k, val in rel.items()}
    return out


def quiet(fn, *a, **k):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        r = fn(*a, **k)
    return r, buf.getvalue()


if __name__ == "__main__":
    E = {0.25: Econ(0.25, RUNS / "rules_12pct_bohn0.250.pkl"), 0.15: Econ(0.15, RUNS / "rules_12pct_bohn0.150.pkl")}
    res = {}
    texts = {}
    # the economy check: the 0.25 rung against the reproduced PR #38 solution
    rules_pr38 = pickle.load(open(PR38 / "rules_risk.pkl", "rb"))
    for tag in ("_rest_point_exact", "_rest_point"):
        if hasattr(rules_pr38, tag):
            delattr(rules_pr38, tag)
    for b, e in E.items():
        r = res[str(b)] = {}
        cal, ss, sp, rules = e.cal, e.ss, e.sp, e.rules
        r["gamma_tau"] = ss["gs_D"]["gamma_tau"]
        r["debt_root"] = 1 - cal["delta_b_D"] + (cal["delta_b_D"] - ss["gs_D"]["gamma_tau"]) / ss["gs_D"]["Q_B_ss"]
        # 1. rest points: fitted iteration and exact fixed point
        S_fit = RE.stochastic_rest_point(rules, cal, ss, sp, verbose=False, exact=False)
        for tag in ("_rest_point_exact", "_rest_point"):
            if hasattr(rules, tag):
                delattr(rules, tag)
        S_rest = RE.stochastic_rest_point(rules, cal, ss, sp, verbose=False, exact=True)
        e.S_rest = S_rest
        S_ss = RE.ss_state(ss, cal, sp)
        r["rest_state_vs_ss_pct"] = {nm: float(100 * (S_rest[i] / S_ss[i] - 1)) for nm, i in
                                     (("K_D", IK_D), ("P_D", IP_D), ("b_DD", IBDD))}
        r["rest_state"] = S_rest.tolist()
        r["rest_fitted_state"] = S_fit.tolist()
        r["rest"] = dict(fitted_at_fitted=summ(e, fitted(e, S_fit)), exact=summ(e, clear(e, S_rest)))
        r["rest"]["two_quarter"] = summ(e, two_quarter(e, S_rest, clear(e, S_rest).x))
        r["rest"]["pd"] = float(default_prob(S_rest[IS]))
        log(b, "rest", r["rest"]["exact"]["mu"], r["rest"]["exact"]["n"])
        # 2. thresholds and a scan in the default probability
        r["threshold"] = {m: threshold(e, S_rest, m) for m in ("fitted", "exact", "2q")}
        log(b, "threshold", {m: t.get("pd") for m, t in r["threshold"].items()})
        scan = {}
        for pd in (0.005, 0.008, 0.010, 0.012, 0.015, 0.02, 0.03, 0.04, 0.06, 0.09, 0.12):
            scan[str(pd)] = dict(fitted=summ(e, fitted(e, shocked(S_rest, pd))),
                                 exact=summ(e, clear(e, shocked(S_rest, pd))))
        r["scan"] = scan
        # 3. impact responses: fitted, exact, two-quarter, each against the rest point by the same method
        base = dict(fitted=summ(e, fitted(e, S_rest)), exact=r["rest"]["exact"], two_quarter=r["rest"]["two_quarter"])
        r["impact"] = {}
        for pd in SHOCKS:
            S = shocked(S_rest, pd)
            vx = clear(e, S)
            lv = dict(fitted=summ(e, fitted(e, S)), exact=summ(e, vx), two_quarter=summ(e, two_quarter(e, S, vx.x)))
            r["impact"][str(pd)] = dict(levels=lv, base=base, cont_accuracy=cont_accuracy(e, vx))
            log(b, pd, {m: (round(l["ls_bp"], 1), round(100 * (l["Y"] / base[m]["Y"] - 1), 4)) for m, l in lv.items()})
        # exact dynamic IRFs
        r["irf"] = {}
        for pd in SHOCKS:
            p, txt = quiet(RE.dynamic_irf, rules, cal, ss, sp, pd_shock=pd, T=40, rest_verbose=False, exact=True)
            texts[f"irf_{b}_{pd}"] = txt
            r["irf"][str(pd)] = {k: p[k].tolist() for k in ("Y", "I", "C", "spread", "sov_bp", "dQ_bD", "n", "mu", "K",
                                                             "rdep", "at_wall")}
        log(b, "irfs done")
        # 4. the surprise 200 bp cap
        r["cap"] = {}
        for pd in SHOCKS:
            try:
                out, txt = quiet(TE.intervention, cal, ss, sp, rules, pd_shock=pd, mode="cap", cap_bp=200.0)
                texts[f"cap_{b}_{pd}"] = txt
                r["cap"][str(pd)] = dict(purchase_pct=100 * out["purchase"], relief=out["relief"],
                                         Y_off=out["irf_off"]["Y"].tolist(), Y_on=out["irf_on"]["Y"].tolist(),
                                         I_off=out["irf_off"]["I"].tolist(), I_on=out["irf_on"]["I"].tolist(),
                                         n_off=out["irf_off"]["n"].tolist(), n_on=out["irf_on"]["n"].tolist(),
                                         sp_off=out["irf_off"]["spread"].tolist(), sp_on=out["irf_on"]["spread"].tolist(),
                                         sov_off=out["irf_off"]["sov_bp"].tolist(), sov_on=out["irf_on"]["sov_bp"].tolist(),
                                         C_off=out["irf_off"]["C"].tolist(), C_on=out["irf_on"]["C"].tolist())
            except Exception as ex:
                r["cap"][str(pd)] = dict(error=repr(ex))
            log(b, "cap", pd, r["cap"][str(pd)].get("purchase_pct", r["cap"][str(pd)].get("error")))
        # 5. off-grid accuracy on a simulated ergodic set (fitted rules, as the solver's own report)
        states, off = AC.simulate(rules, cal, ss, sp, S_rest, T=1200, burn=200, seed=0)
        err = AC.euler_errors(rules, cal, ss, sp, states)
        ok = np.isfinite(err).all(axis=1)
        l10 = np.log10(np.maximum(err[ok], 1e-16))
        pdv = default_prob(states[:, IS])
        r["accuracy"] = dict(n=int(ok.sum()), l10_mean=float(l10.mean()), l10_max=float(l10.max()),
                             l10_mean_by_eq=l10.mean(axis=0).tolist(), l10_max_by_eq=l10.max(axis=0).tolist(),
                             off_frac=(off > 0).mean(axis=0).tolist(), pd_p95=float(np.percentile(pdv, 95)),
                             pd_max=float(pdv.max()), share_binding=None)
        # where the simulated economy binds: exact clearing on every 10th visited state
        bind = [summ(e, clear(e, S))["mu"] > 0 for S in states[::10]]
        r["accuracy"]["share_binding_exact"] = float(np.mean(bind))
        log(b, "accuracy", r["accuracy"]["l10_mean"], r["accuracy"]["l10_max"])
    # the economy check: PR #38's reproduced rules at 0.25
    e25 = E[0.25]
    S38 = RE.stochastic_rest_point(rules_pr38, e25.cal, e25.ss, e25.sp, verbose=False, exact=True)
    e38 = copy.copy(e25); e38.rules = rules_pr38; e38.ngh = rules_pr38.n_gh or RE.N_GH
    pr = dict(rest_exact=summ(e38, clear(e38, S38)), impact={str(pd): summ(e38, clear(e38, shocked(S38, pd))) for pd in SHOCKS})
    res["pr38_check_at_0.25"] = pr
    # why the threshold moves: each economy's rules at the other's rest-point balance sheets (exact), at 1.5% and 3%
    cross = {}
    for pd in (0.015, 0.03):
        for rb in (0.25, 0.15):
            for sb in (0.25, 0.15):
                cross[f"{pd}_rules{rb}_state{sb}"] = summ(E[rb], clear(E[rb], shocked(E[sb].S_rest, pd)))
    res["crossover"] = cross
    json.dump(res, open(OUT / "bohn_compare.json", "w"), indent=1,
              default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o))
    with open(OUT / "bohn_compare_irf_text.log", "w") as fh:
        for k, t in texts.items():
            fh.write(f"===== {k}\n{t}\n")
    log("saved", OUT)
