# Forensic decomposition of the CANONICAL k=2 BFT solution (Bohn 0.15; copy of bft_rule/forensic/forensic.py). One-quarter clearings on frozen rules only.
import json, time
from forensic_core import *

OUTF = OUT / "forensic" / "forensic_out.json"
t0 = time.time()
S_off = rest_point(rules_off, cal0)
S_pol = rest_point(rules_pol, cal_pol)
lb, lK, f = cal0["lambda_bD_D"], cal0["lambda_K_D"], cal0["f_D"]
res = {}

FAC = ["n_sovgain", "n_sovqty", "n_cap", "n_other", "Lam", "alpha", "R", "lev_swap", "lev_bprice", "lev_bqty",
       "lev_K", "lev_bF"]


def factors(a, d):
    # additive pieces of the change from a to d in n and in the divertible base
    dn = dict(n_sovgain=(1 - f) * (d["Xi"] - a["Xi"]) * a["b_lag"],
              n_sovqty=(1 - f) * d["Xi"] * (d["b_lag"] - a["b_lag"]),
              n_cap=d["n_cap"] - a["n_cap"])
    dn["n_other"] = d["n"] - a["n"] - sum(dn.values())
    dl = dict(lev_swap=-lb * d["Q"] * (d["Mnew"] - a["Mnew"]),
              lev_bprice=lb * (d["Q"] - a["Q"]) * a["b_DD"],
              lev_K=d["lev_K"] - a["lev_K"], lev_bF=d["lev_bF"] - a["lev_bF"] + d["lev_cb"] - a["lev_cb"])
    dl["lev_bqty"] = d["lev"] - a["lev"] - sum(dl.values())
    return dn, dl


def g_builder(a, d):
    dn, dl = factors(a, d)

    def state(z):
        zz = dict(zip(FAC, z))
        n = a["n"] + sum(zz[k] * dn[k] for k in dn)
        lev = a["lev"] + sum(zz[k] * dl[k] for k in dl)
        EOm = EOm_mix(d if zz["Lam"] else a, d if zz["alpha"] else a, f)
        R = d["R"] if zz["R"] else a["R"]
        psi = 1 - EOm * R * n / lev
        mu = max(psi, 0.0)
        return psi, mu, 4e4 * lK * mu / EOm
    return state, dn, dl


def decompose(a, d, tag):
    st, dn, dl = g_builder(a, d)
    chk_a, chk_d = st([0] * len(FAC)), st([1] * len(FAC))
    out = dict(tag=tag, mu=(a["mu"], d["mu"]), psi=(a["psi"], d["psi"]), ls=(a["ls_bp"], d["ls_bp"]),
               rebuild_err=(chk_a[2] - a["ls_bp"], chk_d[2] - d["ls_bp"], chk_d[0] - d["psi"]))
    for i, nm in enumerate(("psi", "mu", "ls")):
        phi, tot = shapley(lambda z: st(z)[i], len(FAC))
        out["sh_" + nm] = dict(zip(FAC, phi)); out["tot_" + nm] = tot
    # one-at-a-time bounds on the lending spread: each factor switched alone from a, and removed alone from d
    out["solo_from_a"] = {k: st([int(j == i) for j in range(len(FAC))])[2] - a["ls_bp"] for i, k in enumerate(FAC)}
    out["solo_from_d"] = {k: d["ls_bp"] - st([int(j != i) for j in range(len(FAC))])[2] for i, k in enumerate(FAC)}
    # bond price: exact log split
    dlQ = np.log(d["Q"] / a["Q"])
    dEp = d["Epay"] - a["Epay"]
    lEp = np.log(d["Epay"] / a["Epay"])
    pieces = dict(default_warranted=lEp * (d["Epay_cpn"] - a["Epay_cpn"]) / dEp if dEp else 0.0,
                  cont_nodefault=lEp * (d["Epay_Q_nd"] - a["Epay_Q_nd"]) / dEp if dEp else 0.0,
                  cont_default=lEp * (d["Epay_Q_d"] - a["Epay_Q_d"]) / dEp if dEp else 0.0,
                  deposit_rate=-np.log(d["R"] / a["R"]), covariance=np.log(d["cov"] / a["cov"]),
                  gk_premium=np.log(d["gk"] / a["gk"]))
    pct = 100 * (d["Q"] / a["Q"] - 1)
    out["Q_pct"] = pct
    out["Q_pieces_pct"] = {k: pct * v / dlQ for k, v in pieces.items()}
    out["Q_logsum_err"] = sum(pieces.values()) - dlQ
    # net worth: exact additive split, % of a's n; the sovereign capital gain allocated to the price pieces
    out["n_pct"] = 100 * (d["n"] / a["n"] - 1)
    nparts = dict(sov_capital_gain=(1 - f) * (d["Xi"] - a["Xi"]) * a["b_lag"],
                  sov_legacy_qty=(1 - f) * d["Xi"] * (d["b_lag"] - a["b_lag"]),
                  capital_return=d["n_cap"] - a["n_cap"], f_bonds=d["n_fb"] - a["n_fb"],
                  cb_claim_repaid=d["n_O"] - a["n_O"], deposits_owed=d["n_P"] - a["n_P"],
                  entry_transfer=d["n_ent"] - a["n_ent"], recap=d["n_recap"] - a["n_recap"])
    out["n_parts_pct"] = {k: 100 * v / a["n"] for k, v in nparts.items()}
    out["n_sum_err"] = sum(nparts.values()) - (d["n"] - a["n"])
    gain = nparts["sov_capital_gain"]
    dXi_from_Q = (1 - f) * d["surv"] * (1 - cal0["delta_b_D"]) * (d["Q"] - a["Q"]) * a["b_lag"]
    out["n_sovgain_by_Qpiece_pct"] = {k: 100 * dXi_from_Q * v / pct / a["n"] for k, v in out["Q_pieces_pct"].items()} if pct else {}
    out["levels"] = {k: (a[k], d[k]) for k in ("Q", "n", "lev", "EOm", "R", "mu", "psi", "ls_bp", "sov_bp", "rdep_bp",
                                              "b_DD", "b_DF", "Bp", "Mnew", "m", "Z", "Y", "I", "C", "C_F", "Kp", "QK",
                                              "alpha_cur", "alpha_new", "slack_fit", "err", "Epay", "cov", "gk", "er_bp")}
    out["dn"] = {k: 100 * v / a["n"] for k, v in dn.items()}
    out["dl"] = {k: 100 * v / a["lev"] for k, v in dl.items()}
    return out


def brief(L):
    return {k: float(L[k]) for k in ("mu", "psi", "ls_bp", "Q", "sov_bp", "n", "lev", "EOm", "rdep_bp", "Y", "I", "C",
                                     "C_F", "m", "Mnew", "b_DD", "b_DF", "Bp", "err", "slack_fit", "alpha_cur",
                                     "alpha_new", "er_bp")}


OUTF.parent.mkdir(exist_ok=True)
for pd in (None, 0.02, 0.015):
    key = "rest" if pd is None else str(pd)
    So, Sp = shocked(S_off, pd), shocked(S_pol, pd)
    v_A = clear(rules_off, cal0, So)
    v_D = clear(rules_pol, cal_pol, Sp)
    v_Dc = clear(rules_pol, cal_pol, So)            # policy at the no-policy initial state
    v_Ac = clear(rules_off, cal0, Sp)               # no policy at the policy initial state
    share = v_Dc.m_cb / B
    v_B = clear(rules_off, cal0, So, share=share)   # the same purchase as a surprise (no-policy expectations)
    v_C = clear(rules_pol, cal_pol, So, share=0.0)  # announced rule, no purchase this quarter
    v_Dchk = clear(rules_pol, cal_pol, So, share=share)
    A, D, Dc, Ac, Bv, C = (legs(v) for v in (v_A, v_D, v_Dc, v_Ac, v_B, v_C))
    r = dict(scen={k: brief(L) for k, L in dict(A_own=A, D_own=D, D_at_Soff=Dc, A_at_Spol=Ac, B_surprise=Bv,
                                                  C_announce_nobuy=C).items()},
             Dchk_gap=float(abs(legs(v_Dchk)["mu"] - Dc["mu"]) + abs(legs(v_Dchk)["Q"] - Dc["Q"])))
    r["dec_own"] = decompose(A, D, "own rest points: D_own vs A_own")
    r["dec_common"] = decompose(A, Dc, "common state S_off: policy vs no policy")
    r["dec_state"] = decompose(Dc, D, "policy rules: S_pol vs S_off (announcement-shifted initial state)")
    r["dec_purchase_surprise"] = decompose(A, Bv, "surprise purchase (B) vs A")
    r["dec_anticip"] = decompose(A, C, "announced, no purchase today (C) vs A")
    r["dec_purchase_given_rule"] = decompose(C, Dc, "purchase today under the rule (D) vs C")
    # divertibility counterfactual at D's allocation: claim counted at lambda_bD
    for nm, L in (("D_own", D), ("D_at_Soff", Dc)):
        lev_cf = L["lev"] + lb * L["Z"]
        psi_cf = 1 - L["EOm"] * L["R"] * L["n"] / lev_cf
        r["divertible_claim_" + nm] = dict(psi=psi_cf, mu=max(psi_cf, 0), ls_bp=4e4 * lK * max(psi_cf, 0) / L["EOm"],
                                           Z_over_lev=L["Z"] / L["lev"])
    # same, re-cleared for this quarter with lambda_cb = lambda_bD (continuation rules unchanged)
    vcf = clear(rules_pol, dict(cal_pol, lambda_cb_D=lb), Sp)
    r["divertible_claim_recleared_D_own"] = brief(legs(vcf))
    # accounting of the swap
    for nm, L, v in (("D_own", D, v_D), ("D_at_Soff", Dc, v_Dc)):
        r["acct_" + nm] = dict(
            swap_invariance=float(L["Q"] * L["b_DD"] + L["Z"] - L["Q"] * (L["Bp"] - L["b_DF"])),
            claim_vs_book=float(L["Z"] - L["Q"] * L["Mnew"]), next_claim=float(v.Op_cb - (1 + v.rdep_D) * v.Z_cb),
            clearing=float(v.Bp_D - v.b_D_D_new - v.b_D_F_new - v.M_cb_new),
            balance=float(v.assets_D - v.n_D - v.dep_D), lambda_cb=cal_pol.get("lambda_cb_D", 0.0),
            ic_closed_form=float(L["psi"]), ic_slack_fitted_alpha=float(L["slack_fit"]),
            ic_slack_own_alpha=float(L["alpha_new"] * L["n"] - L["lev"]))
    # expected Eurosystem P&L next quarter on the book carried out of D_own (by regime, per unit of SS D output)
    v = v_D
    db, rec = cal0["delta_b_D"], cal0["recovery_rate_D"]
    pnl = [((rec if v.reg[j] else 1.0) * (db + (1 - db) * v.R[j]["Q_bD"]) * v.M_cb_new - v.Op_cb) for j in range(v.nreg)]
    Yss = ss["ss_firm_D"]["Y_ss"]
    r["cb_pnl_next"] = dict(expected=float(PM._expect(v.wgt, pnl)) / Yss,
                            nodefault=float(np.mean(pnl[0])) / Yss, default=float(np.mean(pnl[1])) / Yss,
                            pd=float(default_prob(v.s)))
    res[key] = r
    print(f"{key} done {time.time() - t0:.0f}s", flush=True)

json.dump(res, open(OUTF, "w"), default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o), indent=1)
print("saved", OUTF)
