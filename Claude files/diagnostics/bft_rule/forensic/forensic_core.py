# Shared loaders and one-quarter clearing for the forensic decomposition of the k=2 BFT solution.
# No global solve: every object here is a one-quarter clearing on frozen, already-solved rules.
import sys, pickle, itertools
from pathlib import Path
import numpy as np
from scipy.optimize import root

ROOT = Path("/Users/Adam/Documents/uni/phd/research/QUANTITATIVE_MODEL_pr38")
sys.path.insert(0, str(ROOT))
from results_io import load
from global_projection.solver_recursive.state_grid import s_process_params, default_prob, IS, IM, IO, IBDD, IBDF, IP_D
from global_projection.solver_recursive.decision_rules import SOLVE7
from global_projection.solver_recursive import point_map as PM
from global_projection.solver_recursive import recursive_experiment as RE

MAIN = ROOT.parent / "QUANTITATIVE_MODEL"
DATA = MAIN / "results" / "GLOBAL_2026-11-12_pr38" / "data"
OUT = MAIN / "results" / "GLOBAL_2026-11-12_pr38_bft"

cal0 = load(DATA / "calibration"); cal0["tpi_on"] = False
ss = load(DATA / "steady_state")
sproc = s_process_params(cal0)
summ = load(OUT / "summary")
rules_off = pickle.load(open(DATA / "rules_risk.pkl", "rb"))
rules_pol = pickle.load(open(OUT / "rules_bft_k2.pkl", "rb"))
cal_pol = dict(cal0, tpi_on=True, tpi_rule="bft", tpi_gz=False, tpi_bft_er_ref=summ["er_ref_bp"] / 4e4,
               tpi_bft_k=2.0)
B = cal0["B_gov_D_ss"]
NGH = rules_off.n_gh or RE.N_GH


class Proxy:
    # the frozen rules, with today's C, alpha and r_wc (read from the fitted rules by the map) overridable
    def __init__(self, base, over=None):
        self._b, self._o = base, dict(over or {})

    def eval(self, name, d, x):
        if name in self._o:
            return np.array([self._o[name]])
        return self._b.eval(name, d, x)

    def __getattr__(self, k):
        return getattr(self._b, k)


class ContProxy(Proxy):
    # continuation values replaced by supplied (exactly cleared) values at the stored next states
    def __init__(self, base, cont_over):
        super().__init__(base)
        self._c = cont_over

    def eval_cont(self, j, Sn):
        r = self._b.eval_cont(j, Sn)
        if j in self._c:
            r = dict(r)
            for k, val in self._c[j].items():
                r[k] = np.array(val, dtype=float)
        return r


def run_map(S, d, x, cont, cal):
    v = PM._read_state_and_regime(S, d, x, cont, cal, ss)
    PM._production(v, cal); PM._government(v, cal, ss); PM._balance_sheets(v, cont, cal)
    PM._continuation(v, cont, cal, sproc, NGH, False); PM._discount_kernels(v, cont, cal)
    v._C_cur = float(np.real(cont.eval("C_D", v.d, v.Sm)[0]))
    PM._incentive_constraint(v, cal); PM._bond_demands(v, cal); PM._tpi_rule(v, cal, ss)
    PM._wages_and_dividends(v, cal); PM._households(v, cal, ss)
    v.res = PM._residual_vector(v, cal, ss)
    v.x = np.array(x, dtype=float)
    return v


def x_fit(rules, S, d=0):
    return np.array([float(rules.eval_read(k, d, np.atleast_2d(S))[0]) for k in SOLVE7])


def clear(rules, cal, S, share=None, d=0, cont=None, x0=None):
    # exact one-quarter clearing at S; share fixes today's purchase (TPI rule switched off for this quarter only)
    cont = rules if cont is None else cont
    c = cal if share is None else dict(cal, tpi_on=False)
    x0 = x_fit(rules, S, d) if x0 is None else x0
    if share is None:
        f = lambda x: PM.point_residuals(S, d, x, cont, c, ss, sproc, n_gh=NGH)[0]
        sol = root(f, x0, method="hybr", tol=1e-13); x = sol.x
    else:
        f = lambda z: PM.point_residuals(S, d, np.append(z, share), cont, c, ss, sproc, n_gh=NGH)[0][:-1]
        sol = root(f, x0[:-1], method="hybr", tol=1e-13); x = np.append(sol.x, share)
    v = run_map(S, d, x, cont, c)
    v.err = float(np.max(np.abs(v.res if share is None else v.res[:-1])))
    v.cal = c
    return v


def fitted(rules, cal, S, d=0):
    v = run_map(S, d, x_fit(rules, S, d), rules, cal)
    v.err = float(np.max(np.abs(v.res))); v.cal = cal
    return v


def consistent(rules, cal, S, share=None, tol=1e-11, maxit=80):
    # clearing in which today's C, alpha and r_wc used inside the map equal the cleared ones (fixed point)
    over, x0 = {}, None
    for it in range(maxit):
        v = clear(rules, cal, S, share, cont=Proxy(rules, over), x0=x0)
        new = dict(C_D=v.C_D, C_F=v.C_F, alpha_D=v.alpha_D_new, alpha_F=v.alpha_F_new, r_wc_D=v.r_wc_D,
                   r_wc_F=v.r_wc_F)
        gap = max(abs(new[k] - over.get(k, np.inf)) for k in new) if over else np.inf
        over, x0 = new, v.x
        if gap < tol:
            break
    v.cons_gap, v.cons_it = gap, it + 1
    return v


def rest_point(rules, cal):
    return RE.stochastic_rest_point(rules, cal, ss, sproc, verbose=False)


def shocked(S, pd):
    S = S.copy()
    if pd is not None:
        S[IS] = RE.s_from_pd(pd)
    return S


# ---------- the objects of the decomposition ----------
def legs(v):
    cal = v.cal
    f, om, db = cal["f_D"], cal["omega_ent_D"], cal["delta_b_D"]
    lK, lb, lbF, lcb = cal["lambda_K_D"], cal["lambda_bD_D"], cal["lambda_bF_D"], cal.get("lambda_cb_D", 0.0)
    payF_now = cal["delta_b_F"] + (1 - cal["delta_b_F"]) * v.Q_bF
    L = dict(
        # net worth n = (1-f)(X - P) + omega*assets
        n=v.n_D,
        n_cap=(1 - f) * (v.mpk_D + (1 - cal["delta_D"]) * v.Q_D) * v.K_D,
        n_sov=(1 - f) * v.Xi_D * v.b_D_D_lag,
        n_fb=(1 - f) * v.p * payF_now * v.b_F_D_lag,
        n_O=(1 - f) * v.O_lag, n_P=-(1 - f) * v.P_D, n_ent=om * v.assets_D, n_recap=v.recap_D,
        # divertible base
        lev=v.lev_D, lev_K=lK * (v.Q_D * v.Kp_D + v.L_wc_D), lev_b=lb * v.Q_bD * v.b_D_D_new,
        lev_bF=lbF * v.p * v.Q_bF * v.b_F_D_new, lev_cb=lcb * v.Z_cb,
        EOm=v.E_Om_D, R=1 + v.rdep_D, mu=v.mu_D, psi=1 - v.E_Om_D * (1 + v.rdep_D) * v.n_D / v.lev_D,
        ls_bp=4e4 * cal["lambda_K_D"] * v.mu_D / v.E_Om_D,
        Q=v.Q_bD, b_DD=v.b_D_D_new, b_DF=v.b_D_F_new, Bp=v.Bp_D, Mnew=v.M_cb_new, m=v.m_cb, Z=v.Z_cb,
        Epay=v.E_payD, EOmpay=v.E_Om_payD, foc=v.foc_D, assets=v.assets_D, dep=v.dep_D,
        Y=v.Y_D, I=v.I_D, C=v.C_D, C_F=v.C_F, slack_fit=v.slack_D, alpha_cur=v.alpha_D_cur,
        alpha_new=v.alpha_D_new, Kp=v.Kp_D, QK=v.Q_D, mpk=v.mpk_D, rdep_bp=4e4 * v.rdep_D, err=v.err,
        sov_bp=4e4 * (db * (1 / v.Q_bD - 1) - cal["delta_b_F"] * (1 / v.Q_bF - 1)), er_bp=4e4 * v.er_D,
        b_lag=v.b_D_D_lag, surv=v.surv, Xi=v.Xi_D, Pi_cb=v.Pi_cb, Op=v.Op_cb, O_lag=v.O_lag, M_lag=v.M_lag,
    )
    # bond-price factors: Q = [Epay/R] * [EOmpay/(EOm*Epay)] * [EOm*R/(EOm*R + lb*mu)]
    L["Q_rn"] = v.E_payD / (1 + v.rdep_D)
    L["cov"] = v.E_Om_payD / (v.E_Om_D * v.E_payD)
    L["gk"] = v.E_Om_D * (1 + v.rdep_D) / (v.E_Om_D * (1 + v.rdep_D) + lb * v.mu_D)
    L["Q_id"] = L["Q_rn"] * L["cov"] * L["gk"]
    # Epay = delta*E[hc] + (1-delta)*E[hc*Q'], by regime
    w = v.wgt
    hc = [cal["recovery_rate_D"] if r else 1.0 for r in v.reg]
    L["Epay_cpn"] = PM._expect(w, [db * hc[j] * np.ones_like(v.R[j]["Q_bD"]) for j in range(v.nreg)])
    L["Epay_Q_nd"] = PM._expect([w[0], 0 * w[1]], [(1 - db) * hc[0] * v.R[0]["Q_bD"], 0 * v.R[1]["Q_bD"]])
    L["Epay_Q_d"] = PM._expect([0 * w[0], w[1]], [0 * v.R[0]["Q_bD"], (1 - db) * hc[1] * v.R[1]["Q_bD"]])
    # SDF and franchise pieces of Omega_j = Lam_j (f + (1-f) alpha'_j)
    xc = PM._ghh(v_cur_C(v), v.N_D, cal["chi_D"], cal["frisch_D"], v.X_FLOOR_D, v.X_EPS)
    L["Lam"] = [cal["beta_inter_D"] * (xc / v.XN_D[j]) ** cal["sigma_D"] for j in range(v.nreg)]
    L["alp"] = [v.R[j]["alpha_D"] for j in range(v.nreg)]
    L["w"] = w
    return L


def v_cur_C(v):
    return v._C_cur if hasattr(v, "_C_cur") else None


def EOm_mix(La, Lb, f):
    # E[Lam_a (f + (1-f) alpha_b)]
    return PM._expect(La["w"], [La["Lam"][j] * (f + (1 - f) * Lb["alp"][j]) for j in range(len(La["Lam"]))])


def shapley(fun, n):
    # exact Shapley allocation of fun(all on) - fun(all off) over n factors
    from math import factorial
    vals = {}
    for mask in range(1 << n):
        vals[mask] = fun([(mask >> i) & 1 for i in range(n)])
    phi = np.zeros(n)
    for i in range(n):
        for mask in range(1 << n):
            if mask & (1 << i):
                continue
            k = bin(mask).count("1")
            wgt = factorial(k) * factorial(n - k - 1) / factorial(n)
            phi[i] += wgt * (vals[mask | (1 << i)] - vals[mask])
    return phi, vals[(1 << n) - 1] - vals[0]
