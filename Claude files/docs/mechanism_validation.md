# HHBANK v12 — Mechanism Validation

**Status:** complete. Method: IRF decompositions on the audited re-run (baseline + W-1/W-2/W-3-fixed variants), not channel-kill re-solves; the decompositions identify channel sizes directly.

Experiment of record: `shock_def_D` = +1pp quarterly default probability, ρ = 0.8 (the notebook's own experiment). Baseline calibration (`writeoff_enabled = 0`).

---

## 1. Sovereign–bank doom loop

**Intended:** default risk ↑ → bond prices ↓ → bank net worth ↓ → credit/leverage ↓ → output ↓ → fiscal position ↓ → default risk ↑.

**Implemented links:** def_rate(+1) → FOC spread (ψ_spread = 0.778) → q_b_D ↓ (−2.7% on impact) → rb_actual_D < 0 → rn → n_inter (P2) → θ via IC → K via balance sheet; fiscal: q_b ↓ → issuance revenue ↓ → b_gov → def_rate (def_scale).

**Verdict: NOT OPERATIVE — three of five links break or flip sign.**

1. **Bank net worth RISES** on the default shock: n_inter_D(0) = +4.6%. Decomposition of rn_D at t=0 (audit run):
   - bond capital loss contribution: −0.25% (φ_bD = 0.25, drb_D = −2.5%)
   - capital-return contribution: +6.76% — dominated by **rdep_D falling −1.49pp at t=0**, leveraged ×3.6
   - net: rn_D(0) = +5.2%
   The deposit rate is a *contemporaneous* unknown paid on the *predetermined* deposit stock (9× quarterly GDP): the equilibrium flight-to-safety fall in rdep transfers ≈ 1.49pp × 9 ≈ 0.134 (= the entire +0.139 net-worth rise) from households to banks at impact. Deposits are implicitly state-contingent claims (audit.md T-2). Git history confirms the authors fought this ("Fix perverse GDP increase under default shock: fix rdep exogenously", then reverted).
2. **Output RISES**: Y_D(0) = +0.012% (perverse; survives the Walras fixes at +0.009%).
3. **Fiscal feedback is anti-doom-loop**: b_gov_D *falls* (−0.18% at t0, −0.29% peak) after the shock (windfall-boosted dividends raise taxable income; lower rdep lowers debt service via the GE response), so endogenous def_rate contribution is *negative* from t=1 onward: def_rate(1) = 0.8% (pure AR shock) − ~0.02% (endogenous).
4. Bond-loss link is intrinsically small: bonds = 10% of bank assets, duration ~2.5y (δ_b = 0.1), and `writeoff_enabled = 0` ⇒ zero cash-flow losses; only mark-to-market q movements bite.
5. Spread link works (below).

**Missing links / fixes:** make deposits non-contingent (pay rdep_{t−1} on Dep_{t−1} — re-date `deposit_return_*` and `bank_return_*`); enable write-offs or lengthen duration to give default actual teeth; reconsider whether capital-return windfall (Q_D rises on impact) is intended.

## 2. Endogenous spread amplification

**Implemented:** spread_t+1 required = xs_ss + ψ_spread·def_rate(t+1), ψ_spread = λ·ψ_λB/(βΩ) = 0.778; def_rate endogenous in b_gov(−1) with slope ≈ 0.11 per unit quarterly-debt ratio.

**Verdict: PRESENT but essentially exogenous.** spread_rb(0) = +26bp/q (≈ +1.04pp annualized) per 1pp default shock — entirely the direct ψ_spread·def term. The *endogenous* (debt-feedback) component is negative (see above). Amplification loop spread→debt→spread has the wrong sign in equilibrium. The diagnostic cell TPI-11 itself documents that def_rate is invariant to who holds the bonds.

## 3. Intermediary net-worth channel

**Verdict: PRESENT for TFP shocks, INVERTED for default shocks.** On the TFP shock the channel is hyperactive: n_inter_D(0) = −20.6% on a +1% Z_D shock (θ jumps +22%) — implausibly large, driven by the same contemporaneous-rdep mechanism (rdep_D(0) = +9.4pp quarterly). Any business-cycle moment discipline would reject this volatility. On default shocks the channel runs backwards (windfall). The channel's *transmission* to real activity (n → θ → K) is intact mechanically (financial_solved IC works), but its *input* is dominated by the deposit-repricing transfer.

## 4. External adjustment channel

**Verdict: QUALITATIVELY PRESENT, QUANTITATIVELY UNRELIABLE pre-fix.** NX_D, p, IM respond with correct signs; but goods_mkt_F leaks up to 2% of GDP (TFP shock), and post-fix p and NX peaks change by 25–100% (walras_forensics.md §7). CA = ΔNFA fails pre-fix (W-1); holds to 4e−7 post-fix. All open-economy quantities must be quoted from the fixed model only.

## 5. Portfolio reallocation channel

**Verdict: PRESENT but calibrated to near-zero.** Cross-border holdings move ≤ 0.29% of their SS level at peak (b_F_D: −0.0013 on SS 0.461) on the default shock. ψ_bF_D = ψ_bD_F = 0.5 was chosen by an IRF-half-life grid search (model_notes), not data; it pins levels tightly. The "contagion through cross-border portfolios" story is thus assumed away in calibration. n_inter_F still falls (−1.5% baseline; −1.1% fixed) — but via q_b_F valuation and ToT, not reallocation.

---

## Summary table

| Mechanism | Implemented? | Sign as intended? | Quantitatively material? |
|---|---|---|---|
| Doom loop (full circle) | partially | **No** (3 links flip) | No |
| Spread amplification | yes | direct leg yes; feedback leg **No** | direct leg yes |
| Net-worth channel | yes | **inverted for default shocks** | yes (overactive for TFP) |
| External adjustment | yes | yes | unreliable pre-fix |
| Portfolio reallocation | yes | yes | **negligible by calibration** |

## Highest-value model improvements (mechanism-side)

1. Re-date deposit rate to be predetermined (kills the windfall; likely restores doom-loop sign).
2. Apply W-1/W-2/W-3 fixes (restores accounting; halves spurious spillovers).
3. Set `writeoff_enabled = 1` with ζ > 0, δ_b ≈ 0.04 (real losses + realistic duration → real doom loop).
4. Re-derive λ_gk multi-asset (audit.md C-1) so Δ's are meaningful and ≤ 1.
5. Calibrate ψ portfolio penalties to observed cross-border holding volatilities, not IRF half-lives.
