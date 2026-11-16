# Forensic decomposition, part 2: the capital-price channel, fixed-investment counterfactuals, uniqueness of the
# one-quarter equilibrium, fitted vs cleared vs fully consistent reads, continuation accuracy and box checks.
import json, time
from forensic_core import *

OUTF = OUT / "forensic" / "forensic2_out.json"
t0 = time.time()
S_off = rest_point(rules_off, cal0)
S_pol = rest_point(rules_pol, cal_pol)
f, lK = cal0["f_D"], cal0["lambda_K_D"]
IKP = 2  # Kp_D in the unknown vector; the D capital Euler is residual 0


def clear_fixK(rules, cal, S, Kp, share=None, x0=None):
    # clear every condition except the D capital Euler, with D investment (hence the capital price) held at Kp
    c = cal if share is None else dict(cal, tpi_on=False)
    x0 = x_fit(rules, S) if x0 is None else x0.copy()
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
    v.err = float(np.max(np.abs(fres(sol.x)))); v.cal = c; v.capeul = float(v.res[0])
    return v


def capsplit(L, v):
    return dict(QK=float(v.Q_D), mpk=float(v.mpk_D), Kp=float(v.Kp_D), I=float(v.I_D),
                n_QK=float((1 - f) * (1 - cal0["delta_D"]) * v.Q_D * v.K_D), n_mpk=float((1 - f) * v.mpk_D * v.K_D),
                K=float(v.K_D), n=float(v.n_D))


def short(L):
    return {k: float(L[k]) for k in ("mu", "psi", "ls_bp", "Q", "sov_bp", "n", "EOm", "rdep_bp", "Y", "I", "C", "C_F",
                                     "m", "err")}


def cont_check(rules, cal, v):
    # exact clearing at each next-state quadrature node, in both regimes, against the fitted continuation
    eps, _ = PM.gh_nodes(NGH)
    Sn, _ = PM._next_states(v, sproc, eps)
    out = dict(outside_now=float(np.max(rules.grid.outside(np.atleast_2d(v.S)))),
               outside_next=rules.grid.outside(Sn).max(axis=0).tolist())
    over, errs = {}, {}
    for j in range(v.nreg):
        if not np.any(v.wgt[j] > 0):
            continue
        vals = {k: [] for k in ("Q_bD", "alpha_D", "C_D", "N_D", "Kp_D", "p", "Q_bF", "alpha_F", "C_F", "N_F", "Kp_F")}
        ok = []
        for i in range(Sn.shape[0]):
            try:
                w = clear(rules, cal, Sn[i], d=j)
                ok.append(w.err)
                vals["Q_bD"].append(w.Q_bD); vals["alpha_D"].append(w.alpha_D_new); vals["C_D"].append(w.C_D)
                vals["N_D"].append(w.N_D); vals["Kp_D"].append(w.Kp_D); vals["p"].append(w.p)
                vals["Q_bF"].append(w.Q_bF); vals["alpha_F"].append(w.alpha_F_new); vals["C_F"].append(w.C_F)
                vals["N_F"].append(w.N_F); vals["Kp_F"].append(w.Kp_F)
            except Exception as e:
                ok.append(np.inf)
                for k in vals:
                    vals[k].append(np.nan)
        fit = v.R[j]
        errs[j] = dict(clear_err=ok, **{k: (np.nanmax(np.abs(np.array(vals[k]) / np.real(fit[k]) - 1)))
                                       for k in ("Q_bD", "alpha_D", "C_D", "N_D", "Kp_D")})
        errs[j]["mu_next_exact"] = None
        if np.all(np.isfinite(ok)) and max(ok) < 1e-8:
            over[j] = vals
    out["errs"] = errs
    # the continuation replaced by its exactly cleared values, today's allocation held, then re-cleared
    cp = ContProxy(rules, over)
    w0 = run_map(v.S, 0, v.x, cp, v.cal)
    psi0 = 1 - w0.E_Om_D * (1 + w0.rdep_D) * w0.n_D / w0.lev_D
    Qimp = w0.E_Om_payD / (w0.E_Om_D * (1 + w0.rdep_D) + cal0["lambda_bD_D"] * w0.mu_D)
    out["subst_fixed_alloc"] = dict(psi=float(psi0), mu=float(max(psi0, 0)), ls_bp=float(4e4 * lK * max(psi0, 0) / w0.E_Om_D),
                                    Q_implied=float(Qimp), max_resid=float(np.max(np.abs(w0.res))), regimes=list(over))
    try:
        w1 = clear(rules, cal, v.S, cont=cp, x0=v.x)
        L1 = legs(w1)
        out["subst_recleared"] = short(L1)
    except Exception as e:
        out["subst_recleared"] = str(e)
    return out


res = {}
for pd in (0.015, 0.06):
    key = str(pd)
    So, Sp = shocked(S_off, pd), shocked(S_pol, pd)
    vA = clear(rules_off, cal0, So); vD = clear(rules_pol, cal_pol, Sp); vDc = clear(rules_pol, cal_pol, So)
    share = vDc.m_cb / B
    vB = clear(rules_off, cal0, So, share=share); vC = clear(rules_pol, cal_pol, So, share=0.0)
    r = {}
    # 1. capital price and returns in each scenario
    r["capital"] = {nm: capsplit(legs(v), v) for nm, v in dict(A=vA, B=vB, C=vC, D_common=vDc, D_own=vD).items()}
    # 2. fixed-investment counterfactuals at the common state: D capital (so the capital price) held at A's
    KpA = vA.Kp_D
    fx = dict(B=clear_fixK(rules_off, cal0, So, KpA, share=share, x0=vB.x),
              C=clear_fixK(rules_pol, cal_pol, So, KpA, share=0.0, x0=vC.x),
              D=clear_fixK(rules_pol, cal_pol, So, KpA, x0=vDc.x),
              A=clear_fixK(rules_off, cal0, So, KpA, x0=vA.x))
    r["fixK"] = {nm: dict(**short(legs(v)), capeul=v.capeul, QK=float(v.Q_D)) for nm, v in fx.items()}
    # 3. uniqueness: the D capital Euler as a function of D investment, everything else cleared
    scan = {}
    for nm, rules, cal, S, sh in (("A", rules_off, cal0, So, None), ("D_common", rules_pol, cal_pol, So, None),
                                  ("D_own", rules_pol, cal_pol, Sp, None), ("C", rules_pol, cal_pol, So, 0.0)):
        base = vA.Kp_D
        rows, x0 = [], (vA.x if nm == "A" else (vC.x if nm == "C" else (vD.x if nm == "D_own" else vDc.x)))
        for frac in np.linspace(-0.004, 0.004, 17):
            try:
                w = clear_fixK(rules, cal, S, base * (1 + frac), share=sh, x0=x0)
                rows.append((float(frac), w.capeul, float(w.mu_D), w.err))
                x0 = w.x
            except Exception as e:
                rows.append((float(frac), np.nan, np.nan, np.inf))
        scan[nm] = rows
    r["scan"] = scan
    # alternative starting points for the policy clearing: from the no-policy allocation
    alt = clear(rules_pol, cal_pol, So, x0=np.append(vA.x[:-1], vDc.x[-1]))
    r["alt_root_D_common"] = dict(mu=float(alt.mu_D), Kp=float(alt.Kp_D), err=alt.err, Kp_main=float(vDc.Kp_D))
    alt2 = clear(rules_off, cal0, So, x0=np.append(vDc.x[:-1], 0.0))
    r["alt_root_A"] = dict(mu=float(alt2.mu_D), Kp=float(alt2.Kp_D), err=alt2.err, Kp_main=float(vA.Kp_D))
    # 4. fitted, exact and fully consistent reads
    reads = {}
    for nm, rules, cal, S in (("A_own", rules_off, cal0, So), ("D_own", rules_pol, cal_pol, Sp),
                              ("D_common", rules_pol, cal_pol, So)):
        vf, vx, vc = fitted(rules, cal, S), clear(rules, cal, S), consistent(rules, cal, S)
        reads[nm] = dict(fitted=short(legs(vf)), exact=short(legs(vx)), consistent=dict(**short(legs(vc)),
                         gap=vc.cons_gap, it=vc.cons_it),
                         fitted_today=dict(C_D=float(vx._C_cur), alpha=float(vx.alpha_D_cur), r_wc=float(4e4 * vx.r_wc_D_cur)),
                         cleared_today=dict(C_D=float(vx.C_D), alpha=float(vx.alpha_D_new), r_wc=float(4e4 * vx.r_wc_D)))
    r["reads"] = reads
    # 5. continuation accuracy and box checks
    r["cont"] = {nm: cont_check(rules, cal, v) for nm, rules, cal, v in
                 (("A_own", rules_off, cal0, vA), ("D_own", rules_pol, cal_pol, vD), ("D_common", rules_pol, cal_pol, vDc))}
    res[key] = r
    print(key, f"{time.time() - t0:.0f}s", flush=True)

json.dump(res, open(OUTF, "w"), default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o), indent=1)
print("saved")
