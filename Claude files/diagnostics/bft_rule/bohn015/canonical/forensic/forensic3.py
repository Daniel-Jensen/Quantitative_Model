# CANONICAL (Bohn 0.15) copy of bft_rule/forensic/forensic3.py, with two-quarter clearing of the controlled
# counterfactuals (fixed purchase today, the rule or no policy from next quarter) added.
# Forensic decomposition, part 3: the price of capital. Exact capital-pricing identity, the source of the Q_K rise
# (continuation groups switched one at a time, with today re-cleared), the revaluation identity, the derivative of mu
# with respect to capital, announced-rule paths against the delayed-start counterfactual, and one- vs two-quarter
# exact clearing. Saved rules only; no global solve.
import json, time
from forensic_core import *

OUTF = OUT / "forensic" / "forensic3_out.json"
t0 = time.time()
S_off = rest_point(rules_off, cal0)
S_pol = rest_point(rules_pol, cal_pol)
f, dK, lK, lb, om = cal0["f_D"], cal0["delta_D"], cal0["lambda_K_D"], cal0["lambda_bD_D"], cal0["omega_ent_D"]
IKP = 2
CONT_KEYS = ("Q_bD", "alpha_D", "C_D", "N_D", "Kp_D", "p", "Q_bF", "alpha_F", "C_F", "N_F", "Kp_F")
GROUPS = dict(Kcap=("Kp_D", "N_D"), Bond=("Q_bD",), Fran=("alpha_D",),
              Hh=("C_D", "C_F", "N_F", "p", "Q_bF", "Kp_F", "alpha_F"))
GNAMES = list(GROUPS) + ["Today", "Buy"]


class Hybrid(Proxy):
    # continuation from the no-policy rules, with the groups in `on` taken from the policy rules; today's fitted
    # objects (C, alpha, r_wc read by the map) from the policy rules when "Today" is on
    def __init__(self, on):
        super().__init__(rules_off)
        self.on = set(on)

    def eval(self, name, d, x):
        return (rules_pol if "Today" in self.on else rules_off).eval(name, d, x)

    def eval_cont(self, j, Sn):
        r = dict(rules_off.eval_cont(j, Sn))
        rp = rules_pol.eval_cont(j, Sn)
        for g, keys in GROUPS.items():
            if g in self.on:
                for k in keys:
                    r[k] = rp[k]
        return r


def kprice(v):
    # Q_K = E[Om payK'] / (E[Om](1+r) + lambda_K mu), payK' = mpk' + (1-delta) Q_K'
    cal = v.cal
    Zn = (1 - sproc["rho_z"]) * sproc["z_star"] + sproc["rho_z"] * v.Z_D
    mpk, QKn = [], []
    for j in range(v.nreg):
        m_, q_, _ = PM._cont_capital(np.maximum(v.R[j]["N_D"], 0.05), v.Kp_D, v.R[j]["Kp_D"], Zn, cal, "D")
        mpk.append(m_); QKn.append(q_)
    pay = [mpk[j] + (1 - dK) * QKn[j] for j in range(v.nreg)]
    E = lambda xs: float(np.real(PM._expect(v.wgt, xs)))
    EOm, R = float(np.real(v.E_Om_D)), 1 + float(np.real(v.rdep_D))
    Epay, EOmpay = E(pay), E([v.Om_D[j] * pay[j] for j in range(v.nreg)])
    out = dict(QK=float(np.real(v.Q_D)), Epay=Epay, Empk=E(mpk), EQKn=E(QKn), EQKn_nd=float(np.mean(QKn[0])),
               R=R, EOm=EOm, mu=float(np.real(v.mu_D)), cov=EOmpay / (EOm * Epay),
               gk=EOm * R / (EOm * R + lK * float(np.real(v.mu_D))), Kp=float(np.real(v.Kp_D)),
               I=float(np.real(v.I_D)), Kpp_nd=float(np.mean(np.real(v.R[0]["Kp_D"]))))
    out["QK_id"] = Epay / R * out["cov"] * out["gk"]
    return out


def kdecomp(a, b):
    # exact log split of Q_K(b)/Q_K(a)
    lQ = np.log(b["QK"] / a["QK"])
    dEp = b["Epay"] - a["Epay"]
    lEp = np.log(b["Epay"] / a["Epay"])
    p = dict(dividend=lEp * (b["Empk"] - a["Empk"]) / dEp, cont_price=lEp * (1 - dK) * (b["EQKn"] - a["EQKn"]) / dEp,
             deposit_rate=-np.log(b["R"] / a["R"]), covariance=np.log(b["cov"] / a["cov"]),
             gk_premium=np.log(b["gk"] / a["gk"]))
    pct = 100 * (b["QK"] / a["QK"] - 1)
    return dict(QK_pct=pct, pieces={k: pct * v / lQ for k, v in p.items()}, err=sum(p.values()) - lQ,
                EQKn_pct=100 * (b["EQKn"] / a["EQKn"] - 1), Empk_pct=100 * (b["Empk"] / a["Empk"] - 1),
                Kpp_nd_pct=100 * (b["Kpp_nd"] / a["Kpp_nd"] - 1))


def summ(v):
    L = legs(v)
    k = kprice(v)
    return dict(QK=k["QK"], mu=float(L["mu"]), psi=float(L["psi"]), ls_bp=float(L["ls_bp"]), n=float(L["n"]),
                Q=float(L["Q"]), I=k["I"], Kp=k["Kp"], rdep_bp=float(L["rdep_bp"]), EOm=float(L["EOm"]),
                m_pct=100 * float(L["m"]) / B, Y=float(L["Y"]), err=float(v.err))


def loop_gain(v):
    # analytic within-quarter gain of the capital demand price in its own level, holding K' and every other price
    L = legs(v)
    QK, K, Kp, n, lev, mu = (float(np.real(z)) for z in (v.Q_D, v.K_D, v.Kp_D, v.n_D, v.lev_D, v.mu_D))
    EOm, R = float(np.real(v.E_Om_D)), 1 + float(np.real(v.rdep_D))
    s_n = ((1 - f) * (1 - dK) * QK * K + om * QK * Kp) / n
    s_lev = lK * QK * Kp / lev
    dmu = -(1 - mu) * (s_n - s_lev) if mu > 0 else 0.0
    # sovereign leg: mu lowers the GK premium in Q_bD, which revalues the legacy bond book
    s_b = ((1 - f) * float(np.real(v.surv)) * (1 - cal0["delta_b_D"]) * float(np.real(v.Q_bD)) * float(np.real(v.b_D_D_lag))
           + om * float(np.real(v.Q_bD)) * float(np.real(v.b_D_D_new))) / n
    s_blev = lb * float(np.real(v.Q_bD)) * float(np.real(v.b_D_D_new)) / lev
    dlnQb_dmu = -lb / (EOm * R + lb * mu)
    return dict(s_n=s_n, s_lev=s_lev, dmu_dlnQK=dmu, gain_K=-lK * dmu / (EOm * R + lK * mu),
                s_b=s_b, s_blev=s_blev, dlnQb_dmu=dlnQb_dmu,
                dmu_dlnQb=-(1 - mu) * (s_b - s_blev) if mu > 0 else 0.0,
                dlnQK_dlnKp=cal0["ksi_D"] * QK * Kp / float(np.real(v.I_D)), lev_capital=QK * K / n)


def clear_fixK(rules, cal, S, Kp, x0, share=None):
    c = cal if share is None else dict(cal, tpi_on=False)
    keep = [i for i in range(len(x0)) if i != IKP and not (share is not None and i == len(x0) - 1)]

    def full(z):
        x = x0.copy(); x[keep] = z; x[IKP] = Kp
        if share is not None:
            x[-1] = share
        return x

    def fres(z):
        r = PM.point_residuals(S, 0, full(z), rules, c, ss, sproc, n_gh=NGH)[0]
        return np.delete(r, [0] + ([len(r) - 1] if share is not None else []))
    sol = root(fres, x0[keep], method="hybr", tol=1e-13)
    v = run_map(S, 0, full(sol.x), rules, c)
    v.err = float(np.max(np.abs(fres(sol.x)))); v.cal = c
    return v


def mu_derivative(rules, cal, S, v0, steps=(-1e-3, -3e-4, -1e-4, -3e-5, 3e-5, 1e-4, 3e-4, 1e-3)):
    # mu and Q_K as functions of K' (fractional deviation), every other condition of the quarter re-cleared
    rows = []
    for h in steps:
        w = clear_fixK(rules, cal, S, v0.Kp_D * (1 + h), v0.x)
        rows.append(dict(h=h, mu=float(w.mu_D), psi=float(1 - w.E_Om_D * (1 + w.rdep_D) * w.n_D / w.lev_D),
                         QK=float(w.Q_D), n=float(w.n_D), rdep_bp=float(4e4 * w.rdep_D), EOm=float(w.E_Om_D),
                         capeul=float(w.res[0]), err=w.err))
    return rows


# ---------- two-quarter joint clearing: today and next quarter at every (node, regime) cell solved together ----------
def two_quarter(rules, cal, S, x_today, share=None):
    # share fixes today's purchase (today's rule switched off); next quarter runs under cal as solved
    ct = cal if share is None else dict(cal, tpi_on=False)
    if share is not None:
        x_today = x_today.copy(); x_today[-1] = share
    eps, _ = PM.gh_nodes(NGH)
    nreg, ncell = rules.n_regimes, len(eps)
    nx = len(x_today)

    def states(x):
        v = PM._read_state_and_regime(S, 0, x, rules, ct, ss)
        PM._production(v, ct); PM._government(v, ct, ss); PM._balance_sheets(v, rules, ct)
        return PM._next_states(v, sproc, eps)[0]

    Sn0 = states(x_today)
    z0 = [x_today]
    for j in range(nreg):
        for i in range(ncell):
            z0.append(clear(rules, cal, Sn0[i], d=j).x)
    z0 = np.concatenate(z0)

    nt = nx if share is None else nx - 1
    if share is not None:
        z0 = z0[np.r_[0:nx - 1, nx:len(z0)]]

    def unpack(z):
        x = z[:nt] if share is None else np.append(z[:nt], share)
        return x, z[nt:].reshape(nreg, ncell, nx)

    def F(z, want=False):
        x, xn = unpack(z)
        Sn = states(x)
        over, res_n = {j: {k: [] for k in CONT_KEYS} for j in range(nreg)}, []
        for j in range(nreg):
            for i in range(ncell):
                r, o = PM.point_residuals(Sn[i], j, xn[j, i], rules, cal, ss, sproc, n_gh=NGH)
                res_n.append(r)
                vals = dict(Q_bD=o["Q_bD"], alpha_D=o["alpha_D"], C_D=o["C_D"], N_D=xn[j, i][0], Kp_D=xn[j, i][2],
                            p=xn[j, i][6], Q_bF=o["Q_bF"], alpha_F=o["alpha_F"], C_F=o["C_F"], N_F=xn[j, i][1],
                            Kp_F=xn[j, i][3])
                for k in CONT_KEYS:
                    over[j][k].append(vals[k])
        cp = ContProxy(rules, over)
        r0 = PM.point_residuals(S, 0, x, cp, ct, ss, sproc, n_gh=NGH)[0]
        if share is not None:
            r0 = r0[:-1]
        if want:
            v = run_map(S, 0, x, cp, ct); v.cal = ct
            return v
        return np.concatenate([r0] + res_n)

    sol = root(F, z0, method="hybr", tol=1e-12, options=dict(maxfev=200000))
    v = F(sol.x, want=True)
    v.err = float(np.max(np.abs(F(sol.x))))
    return v


res = {}
for pd in ((0.02, 0.015) if __name__ == "__main__" else ()):
    key = str(pd)
    So, Sp = shocked(S_off, pd), shocked(S_pol, pd)
    vA = clear(rules_off, cal0, So); vD = clear(rules_pol, cal_pol, Sp); vDc = clear(rules_pol, cal_pol, So)
    share = vDc.m_cb / B
    vB = clear(rules_off, cal0, So, share=share); vC = clear(rules_pol, cal_pol, So, share=0.0)
    r = {}
    # 1. exact capital pricing
    K = {nm: kprice(v) for nm, v in dict(A=vA, B=vB, C=vC, D_common=vDc, D_own=vD).items()}
    r["kprice"] = K
    r["kprice_id_err"] = {nm: k["QK_id"] / k["QK"] - 1 for nm, k in K.items()}
    r["kdec"] = {f"{b}_vs_A": kdecomp(K["A"], K[b]) for b in ("B", "C", "D_common", "D_own")}
    # 2. source of the Q_K rise: continuation groups and today's purchase switched to policy, today re-cleared exactly
    eq, imp = {}, {}
    for mask in range(1 << len(GNAMES)):
        on = [g for i, g in enumerate(GNAMES) if mask >> i & 1]
        sh = share if "Buy" in on else 0.0
        hc = Hybrid(on)
        w = clear(rules_off, cal0, So, share=sh, cont=hc, x0=(vDc.x if len(on) > 3 else vA.x).copy())
        eq[mask] = summ(w)
        # the impulse: today's allocation held at A's (purchase set), today's capital demand price re-evaluated
        xa = vA.x.copy(); xa[-1] = sh
        wi = run_map(So, 0, xa, hc, dict(cal0, tpi_on=False)); wi.cal = dict(cal0, tpi_on=False); wi.err = 0.0
        k = kprice(wi)
        imp[mask] = dict(Qd=k["QK_id"], mu=k["mu"], EOm=k["EOm"],
                         Epay=k["Epay"], R=k["R"], cov=k["cov"], gk=k["gk"],
                         psi=float(1 - wi.E_Om_D * (1 + wi.rdep_D) * wi.n_D / wi.lev_D))
    full = (1 << len(GNAMES)) - 1
    r["hybrid_check"] = dict(A_QK=eq[0]["QK"] - K["A"]["QK"], D_QK=eq[full]["QK"] - K["D_common"]["QK"],
                             D_mu=eq[full]["mu"] - float(vDc.mu_D), maxerr=max(e["err"] for e in eq.values()))
    from math import factorial
    def shap(vals):
        n = len(GNAMES); phi = np.zeros(n)
        for i in range(n):
            for m in range(1 << n):
                if m >> i & 1:
                    continue
                c = bin(m).count("1")
                phi[i] += factorial(c) * factorial(n - c - 1) / factorial(n) * (vals[m | 1 << i] - vals[m])
        return dict(zip(GNAMES, phi.tolist()))
    r["src_eq"] = {q: shap({m: (100 * np.log(e["QK"]) if q == "QK_pct" else e[q]) for m, e in eq.items()})
                   for q in ("QK_pct", "mu", "ls_bp", "psi")}
    r["src_eq"]["n_pct"] = shap({m: 100 * np.log(e["n"]) for m, e in eq.items()})
    r["src_imp"] = {"Qd_pct": shap({m: 100 * np.log(e["Qd"]) for m, e in imp.items()}),
                    "psi": shap({m: e["psi"] for m, e in imp.items()})}
    r["solo_eq"] = {g: {q: eq[1 << i][q] - eq[0][q] for q in ("QK", "mu", "ls_bp", "n")} for i, g in enumerate(GNAMES)}
    r["imp_total"] = dict(Qd_pct=100 * np.log(imp[full]["Qd"] / imp[0]["Qd"]), psi=imp[full]["psi"] - imp[0]["psi"],
                          eq_QK_pct=100 * np.log(eq[full]["QK"] / eq[0]["QK"]), eq_psi=eq[full]["psi"] - eq[0]["psi"])
    # 3. revaluation identity and loop gains
    r["loop"] = {nm: loop_gain(v) for nm, v in dict(A=vA, D_common=vDc, D_own=vD, C=vC).items()}
    nA, nD = legs(vA), legs(vDc)
    dQK = float(vDc.Q_D - vA.Q_D)
    r["reval"] = dict(dlnQK_pct=100 * dQK / float(vA.Q_D), legacy_capital=float((1 - f) * (1 - dK) * vA.K_D * dQK),
                      ratio=float((1 - f) * (1 - dK) * vA.Q_D * vA.K_D / vA.n_D),
                      legacy_pct=float(100 * (1 - f) * (1 - dK) * vA.K_D * dQK / vA.n_D),
                      mpk_pct=float(100 * (1 - f) * (vDc.mpk_D - vA.mpk_D) * vA.K_D / vA.n_D),
                      entry_QK_pct=float(100 * om * (vDc.Q_D * vDc.Kp_D - vA.Q_D * vA.Kp_D) / vA.n_D),
                      K_same=float(vDc.K_D - vA.K_D), n_cap_reported=float(100 * (nD["n_cap"] - nA["n_cap"]) / nA["n"]))
    # 4. derivative of mu in K'
    r["deriv"] = {nm: mu_derivative(rl, cl, S, v) for nm, rl, cl, S, v in
                  (("A", rules_off, cal0, So, vA), ("D_common", rules_pol, cal_pol, So, vDc),
                   ("D_own", rules_pol, cal_pol, Sp, vD))}
    # 5. one-quarter vs two-quarter exact clearing
    tq = {}
    for nm, rl, cl, S, v, shr in (("A", rules_off, cal0, So, vA, None), ("D_own", rules_pol, cal_pol, Sp, vD, None),
                                  ("D_common", rules_pol, cal_pol, So, vDc, None),
                                  ("B", rules_off, cal0, So, vB, share), ("C", rules_pol, cal_pol, So, vC, 0.0)):
        t1 = time.time()
        try:
            w = two_quarter(rl, cl, S, v.x, share=shr)
            nl1, nl2 = legs(v), legs(w)
            tq_n = dict(n_cap=(float(nl1["n_cap"]), float(nl2["n_cap"])), n_sov=(float(nl1["n_sov"]), float(nl2["n_sov"])),
                        Xi=(float(nl1["Xi"]), float(nl2["Xi"])), err2=w.err)
            tq[nm] = dict(one=summ(v), two=summ(w), sec=time.time() - t1, kp_one=kprice(v), kp_two=kprice(w), nw=tq_n)
        except Exception as e:
            tq[nm] = dict(one=summ(v), two=str(e))
        print(" two-quarter", nm, f"{time.time() - t1:.0f}s", flush=True)
    r["two_q"] = tq
    res[key] = r
    print(key, f"{time.time() - t0:.0f}s", flush=True)
    OUTF.parent.mkdir(exist_ok=True)
    json.dump(res, open(OUTF, "w"), default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o), indent=1)
print("saved")
