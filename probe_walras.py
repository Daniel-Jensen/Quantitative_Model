"""
Probe Walras leak on TFP shock IRF (walras-plumbing branch).

Strategy: evaluate each accounting identity at IRF time points
{0, 1, 5, 10, 50} and report which is shock-correlated (non-zero
when shock is on, zero otherwise).

If `goods_mkt_D` leaks but every targeted residual is zero, then by
Walras's law one of the *implicit* accounting identities (HH budget,
gov budget, bank wealth flow, trade-balance / NFA, capital producer
flow) must be misaligned with how the explicit blocks are wired.

We compute each identity from FIRST PRINCIPLES using IRF-pathed
variables, NOT from the block outputs — so we see whether the block
formulas agree with the underlying accounting they claim to enforce.
"""

import numpy as np
import sequence_jacobian as sj
from sequence_jacobian import simple, combine
import sys
sys.setrecursionlimit(5000)

import pickle
import pathlib

# Try to reuse a cached SS+IRF from the notebook to avoid re-solving.
CACHE = pathlib.Path('/tmp/walras_probe_cache.pkl')


def lvl(ss, irf, var, t):
    """Level at time t: ss + irf shift (if var has an IRF), else ss."""
    base = float(ss[var]) if var in ss.toplevel or var in ss.internals else 0.0
    if var in irf:
        return base + irf[var][t]
    return base


def report_identity(name, vals_by_t):
    print(f"\n{name}")
    print("  " + " ".join(f"t={t:>3}: {v:+.3e}" for t, v in vals_by_t.items()))


def check_identities(ss, irf, ts=(0, 1, 5, 10, 50)):
    """For each identity, compute residual at given time points."""
    results = {}

    # Convenience: get level of var at time t. var(-1) means lag.
    def L(var, t, lag=0):
        if t - lag < 0:
            return float(ss[var])
        return lvl(ss, irf, var, t - lag)

    for country in ['D', 'F']:
        c = country

        # ---- 1. Capital accumulation identity ----
        #   K = (1-δ)·K(-1) + Φ(I/K(-1))·K(-1)
        # where Φ(x) = γ0·x^(1-ξ) + γ1
        # This is enforced by capital_res_D (target).
        cap_law = {}
        for t in ts:
            K   = L(f'K_{c}', t)
            Kl  = L(f'K_{c}', t, lag=1)
            I   = L(f'I_{c}', t)
            d   = float(ss[f'delta_{c}'])
            ksi = float(ss[f'ksi_{c}'])
            g0  = float(ss[f'gamma0_{c}'])
            g1  = float(ss[f'gamma1_{c}'])
            iota = I / Kl
            phi_K = g0 * iota ** (1 - ksi) + g1
            res = K - (1 - d) * Kl - phi_K * Kl
            cap_law[t] = res
        results[f'cap_law_{c}'] = cap_law

        # ---- 2. cap_profit consistency ----
        #   cap_profit = Q·(K - (1-δ)·K(-1)) - I
        cap_prof = {}
        for t in ts:
            Q  = L(f'Q_{c}', t)
            K  = L(f'K_{c}', t)
            Kl = L(f'K_{c}', t, lag=1)
            I  = L(f'I_{c}', t)
            d  = float(ss[f'delta_{c}'])
            cp_check = Q * (K - (1 - d) * Kl) - I
            cp_block = L(f'cap_profit_{c}', t)
            cap_prof[t] = cp_block - cp_check
        results[f'cap_profit_def_{c}'] = cap_prof

        # ---- 3. Goods market residual (block formula) ----
        gm = {}
        for t in ts:
            Y  = L(f'Y_{c}', t)
            C  = L(f'C_{c}', t)
            I  = L(f'I_{c}', t)
            G  = L(f'G_{c}', t)
            Phi = L(f'Phi_{c}', t)
            T  = L(f'T_{c}', t)
            cp = L(f'cap_profit_{c}', t)
            NX = L(f'NX_{c}', t)
            P  = L(f'P_CES_{c}', t)
            gm_with_cp = Y - (P * C + I + G + Phi + T + cp) - NX
            gm_no_cp   = Y - (P * C + I + G + Phi + T) - NX
            gm[t] = (gm_with_cp, gm_no_cp)
        results[f'goods_mkt_{c}'] = gm

        # ---- 4. Government budget identity ----
        #   coupon + G = P·TAX + net_issuance
        gov = {}
        for t in ts:
            bg  = L(f'b_gov_{c}', t)
            bgl = L(f'b_gov_{c}', t, lag=1)
            G   = L(f'G_{c}', t)
            TAX = L(f'TAX_{c}', t)
            P   = L(f'P_CES_{c}', t)
            qb  = L(f'q_b_{c}', t)
            dr  = L(f'def_rate_{c}', t)
            rr  = float(ss[f'recovery_rate_{c}'])
            zw  = float(ss[f'zeta_writeoff_{c}'])
            db  = float(ss[f'delta_b_{c}'])
            h   = 1 - rr
            sc  = 1 - zw * dr * h
            coupon = db * (1 - dr * h) * bgl
            net_iss = qb * (bg - sc * (1 - db) * bgl)
            res = coupon + G - P * TAX - net_iss
            gov[t] = res
        results[f'gov_budget_{c}'] = gov

        # ---- 5. Bank wealth flow identity ----
        #   n = (1-f)·[(1+rn)·n(-1) + cap_profit] + m - Phi - T
        bank = {}
        for t in ts:
            n   = L(f'n_inter_{c}', t)
            nl  = L(f'n_inter_{c}', t, lag=1)
            rn  = L(f'rn_{c}', t)
            f   = float(ss[f'f_{c}'])
            m   = L(f'm_{c}', t)
            Phi = L(f'Phi_{c}', t)
            T   = L(f'T_{c}', t)
            cp  = L(f'cap_profit_{c}', t)
            gross = (1 + rn) * nl + cp
            res = (1 - f) * gross + m - Phi - T - n
            bank[t] = res
        results[f'bank_flow_{c}'] = bank

        # ---- 6. Bank balance sheet identity ----
        #   Q·K + q_b_D·b_D + q_b_F·b_F = θ·n
        bs = {}
        for t in ts:
            Q   = L(f'Q_{c}', t)
            K   = L(f'K_{c}', t)
            qbD = L('q_b_D', t)
            qbF = L('q_b_F', t)
            if c == 'D':
                bD = L('b_D_D', t)
                bF = L('b_F_D', t)
            else:
                bD = L('b_D_F', t)
                bF = L('b_F_F', t)
            th = L(f'theta_{c}', t)
            n  = L(f'n_inter_{c}', t)
            bs[t] = Q * K + qbD * bD + qbF * bF - th * n
        results[f'bank_bs_{c}'] = bs

        # ---- 7. Banker dividend identity ----
        #   div = f·gross - m
        divres = {}
        for t in ts:
            div = L(f'div_{c}', t)
            rn  = L(f'rn_{c}', t)
            nl  = L(f'n_inter_{c}', t, lag=1)
            f   = float(ss[f'f_{c}'])
            m   = L(f'm_{c}', t)
            cp  = L(f'cap_profit_{c}', t)
            net_div = f * ((1 + rn) * nl + cp) - m
            divres[t] = div - net_div
        results[f'div_def_{c}'] = divres

        # ---- 8. HH aggregate budget identity (in home goods units) ----
        # Per-agent: dep + c = Rgross·dep(-1) + z
        # Aggregate: DEP + C = Rgross·DEP(-1) + Z_agg
        # where Z_agg = sum_i z_i = sum_i (y_pre_i - t_paid_i)
        # In aggregate: Z_agg = (w·N + div)/P_CES - TAX
        # Multiply by P_CES:
        #   P·DEP + P·C = (1+rdep)·P(-1)·DEP(-1) + w·N + div - P·TAX
        # Using P·DEP = D_supply:
        #   D_supply + P·C = (1+rdep)·D_supply(-1) + w·N + div - P·TAX
        hh = {}
        for t in ts:
            DEP = L(f'DEP_{c}', t)
            DEPl = L(f'DEP_{c}', t, lag=1)
            C   = L(f'C_{c}', t)
            TAX = L(f'TAX_{c}', t)
            w   = L(f'w_{c}', t)
            N   = L(f'N_{c}', t)
            div = L(f'div_{c}', t)
            P   = L(f'P_CES_{c}', t)
            Pl  = L(f'P_CES_{c}', t, lag=1)
            rdep = L(f'rdep_{c}', t)
            Rgross = (1 + rdep) * Pl / P
            # Per-agent identity → bundle units:
            lhs = DEP + C
            rhs = Rgross * DEPl + (w * N + div) / P - TAX
            hh[t] = lhs - rhs
        results[f'hh_budget_{c}'] = hh

        # ---- 9. Cobb-Douglas production identity ----
        prod = {}
        for t in ts:
            Y = L(f'Y_{c}', t)
            Z = L(f'Z_{c}', t)
            K = L(f'K_{c}', t)
            N = L(f'N_{c}', t)
            a = float(ss[f'alpha_{c}'])
            prod[t] = Y - Z * K ** a * N ** (1 - a)
        results[f'prod_{c}'] = prod

    # ---- 10. Trade balance (closure) ----
    #   NX_D + p·NX_F = 0
    tb_close = {}
    for t in (0, 1, 5, 10, 50):
        NXD = L('NX_D', t)
        NXF = L('NX_F', t)
        p   = L('p', t)
        tb_close[t] = NXD + p * NXF
    results['trade_close'] = tb_close

    # ---- 11. Trade balance vs imports ----
    #   NX_D = IM_F - p·IM_D ;  NX_F = IM_D - IM_F / p
    tb_D = {}
    tb_F = {}
    for t in (0, 1, 5, 10, 50):
        IMD = L('IM_D', t)
        IMF = L('IM_F', t)
        p   = L('p', t)
        tb_D[t] = L('NX_D', t) - (IMF - p * IMD)
        tb_F[t] = L('NX_F', t) - (IMD - IMF / p)
    results['NX_D_def'] = tb_D
    results['NX_F_def'] = tb_F

    # ---- 12. External account ----
    ca = {}
    for t in (0, 1, 5, 10, 50):
        NXD = L('NX_D', t)
        qbD = L('q_b_D', t)
        qbF = L('q_b_F', t)
        qbDl = L('q_b_D', t, lag=1)
        qbFl = L('q_b_F', t, lag=1)
        bFD = L('b_F_D', t)
        bDF = L('b_D_F', t)
        bFDl = L('b_F_D', t, lag=1)
        bDFl = L('b_D_F', t, lag=1)
        rbD = L('rb_actual_D', t)
        rbF = L('rb_actual_F', t)
        receipts = (1 + rbF) * qbFl * bFDl
        payments = (1 + rbD) * qbDl * bDFl
        nfa = qbF * bFD - qbD * bDF
        ca[t] = NXD + receipts - payments - nfa
    results['ca_res_D'] = ca

    return results


def main():
    print("Loading notebook model + cached SS/IRF...")
    if CACHE.exists():
        with open(CACHE, 'rb') as f:
            data = pickle.load(f)
        ss = data['ss']
        irf_Z = data['irf_Z']
        irf_def = data['irf_def']
        print("Loaded from cache.")
    else:
        print("Cache not found. Run notebook first to populate cache, "
              "or use cache_walras_probe.py to dump SS + IRF.")
        return

    print("\n" + "=" * 80)
    print("PROBE: TFP shock on Z_D (1% AR(0.8))")
    print("=" * 80)

    res = check_identities(ss, irf_Z)
    for name, vals in res.items():
        if name.startswith('goods_mkt_'):
            print(f"\n{name}  (with cap_profit, without cap_profit)")
            for t, (a, b) in vals.items():
                tag_a = "*** LEAKS" if abs(a) > 1e-10 else "ok"
                tag_b = "*** LEAKS" if abs(b) > 1e-10 else "ok"
                print(f"  t={t:>3}: with_cp={a:+.3e} [{tag_a}]   "
                      f"no_cp={b:+.3e} [{tag_b}]")
        else:
            print(f"\n{name}")
            for t, v in vals.items():
                tag = "*** LEAKS" if abs(v) > 1e-10 else "ok"
                print(f"  t={t:>3}: {v:+.3e}  [{tag}]")

    print("\n" + "=" * 80)
    print("PROBE: Default shock on shock_def_D (1% AR(0.8))")
    print("=" * 80)

    res = check_identities(ss, irf_def)
    for name, vals in res.items():
        if name.startswith('goods_mkt_'):
            print(f"\n{name}  (with cap_profit, without cap_profit)")
            for t, (a, b) in vals.items():
                tag_a = "*** LEAKS" if abs(a) > 1e-10 else "ok"
                tag_b = "*** LEAKS" if abs(b) > 1e-10 else "ok"
                print(f"  t={t:>3}: with_cp={a:+.3e} [{tag_a}]   "
                      f"no_cp={b:+.3e} [{tag_b}]")
        else:
            print(f"\n{name}")
            for t, v in vals.items():
                tag = "*** LEAKS" if abs(v) > 1e-10 else "ok"
                print(f"  t={t:>3}: {v:+.3e}  [{tag}]")


if __name__ == '__main__':
    main()
