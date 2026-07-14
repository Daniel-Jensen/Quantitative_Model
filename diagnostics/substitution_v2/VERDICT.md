# VERDICT v2 — deleveraging vs substitution at psi_lambda_B = 2.8 (macro-pru-fix)

**(a) Yes — `Y_D` falls at the data-disciplined dial** (impact −0.0283%SS, trough
−0.0323%SS at t=26), matching the sweep's interpolated −0.033%SS; the peak spread is
147.2bp (≈150bp target), so the impulse transmits and the sign is not a config
artifact. **(b) Yes — the model is deleveraging-dominant by the decisive classifier:
the capital LEVEL `K_D` is negative throughout (impact −0.0106%SS, trough −0.0657%SS
at t=36) — it never rises — alongside net worth −1.26%SS and total assets (θ·n)
−0.34%SS.** Substitution is present but dominated, and — importantly — it points
**toward sovereign bonds, not toward capital**.

Branch `macro-pru-fix` `c6f5707`; `EL_price=0.102491`, `psi_lambda_B=2.8`,
`psi_spread=0.7259`. SS market clearing ≈0 (`goods_mkt_D=-3.9e-7`). Evidence:
`run_log.md`, `ss_values.json`, figures `v2_01/02/03`.

---

## Step 2 — regime: deleveraging-dominant, with a *quantity* tilt into sovereigns

| balance-sheet item | impact | extremum (t) | direction |
|---|---:|---:|---|
| capital LEVEL `K_D` | −0.011%SS | **−0.066%SS (t=36)** | **falls, never rises** |
| net worth `n_inter_D` | −1.26%SS | −1.26%SS (t=0) | falls (worst on impact) |
| total assets `θ·n` | −0.34%SS | −0.34%SS (t=0) | falls |
| Greek-sov **quantity** `b_D_D` | +0.75%SS | **+1.35%SS (t=4)** | **rises** |
| bond **market-value** share | −0.12pp | −0.12pp | falls (MTM) |

The bank **accumulates more Greek sovereign bonds in quantity** (`b_D_D`↑) while its
**capital level and net worth shrink** — the classic contraction-with-sovereign-
absorption pattern. But because `q_b_D` crashes −2.88%, the bond book's *market
value* share actually edges **down**: this is deleveraging in market value, with a
quantity tilt toward sovereigns (banks absorbing the issuance foreigners shed), not
a market-value renationalisation. The decisive number — the sign of the capital
level — is unambiguously negative.

## Step 3a — ΔY is investment/deleveraging, NOT a terms-of-trade NX effect

Goods-market decomposition (impact, identity residual 5e−17):

| contribution | impact | note |
|---|---:|---|
| investment `I_D` | **−1.145e-3** | dominant drag (deleveraging) |
| bundle consumption `P·C_D` | −0.325e-3 | falls |
| net exports `NX_D` | **+1.188e-3** | **positive — cushions, does not drive** |
| government `G_D` | 0 | fixed |
| **net `Y_D`** | **−0.283e-3** | |

The negative aggregate is driven by **investment and consumption falling**; net
exports move the *other* way (+1.19e-3) and nearly offset the investment drag on
impact. So the small headline ΔY (−0.028%) **masks a larger real-economy/credit
contraction cushioned by rising net exports** — the two-country terms-of-trade
channel works against the result, it does not manufacture it.

## Step 3b — substitution push vs deleveraging pull (EL_price on/off)

| capital response | impact | extremum |
|---|---:|---:|
| net `K_D` (EL on) | −0.011%SS | −0.066%SS |
| substitution-leaning (EL off) | −0.009%SS | −0.047%SS |
| deleveraging pull (on−off) | −0.002%SS | **−0.020%SS** |

`EL_price` **deepens** the capital contraction by ~30–40% (−0.047% → −0.066%). And
even with `EL_price` OFF, capital still **falls** (−0.047%): the pre-existing agency
(`psi_spread`) + IC-collateral + MTM channels already deleverage. **A substitution
push that lifts the capital level never emerges in any counterfactual.** So: the
substitution push into capital is ≈0/negative; the deleveraging pull dominates at
every horizon; `EL_price` is a genuine (not double-counted) additional deleveraging
force on top of the existing ones.

## Step 4 — validation moment holds at the calibrated dial

Net-worth-to-spread pass-through at `psi_lambda_B=2.8`: **−0.856%SS per 100bp**
(peak spread 147.2bp, peak Δn_inter −1.26%SS) — matches the sweep's dial-invariant
flag (≈−0.85). This is the number to validate against bank-equity/sovereign event
studies; it is stable at the calibrated dial.

## Step 6 — consistency & the renationalisation nuance

Capital/credit **contracting** (`K_D`↓, `I_D`↓); Greek-sov **quantity** `b_D_D`↑;
bond **market-value** share ↓. The model does **not** produce the doubly-counterfactual
`bonds↓ & capital↑ & credit↑` (checked: NO). One nuance worth flagging against the
prior that "the model won't reproduce renationalisation": it **does** reproduce the
renationalisation *direction in quantity* (`b_D_D`↑ as capital/credit fall) — but
mechanically, via government financing + foreign flight clearing, **not** via a
risk-shifting motive, and **not** in market value. So the sign pattern is consistent
with a Greek-crisis contraction; it just arrives through the financing/clearing
channel rather than the (absent) political motive.

## Step 5 — not triggered

Capital levels do not rise, so there is no channel inconsistency to localise and no
`recommended_fix.md` is warranted.

---

**Bottom line.** At the data-disciplined dial, deleveraging cleanly dominates
substitution (capital level falls, never rises; net worth and total assets fall),
the output decline is a genuine investment/credit contraction (cushioned, not caused,
by net exports), `EL_price` strengthens rather than distorts this, and the doom-loop
pass-through validation moment is stable at −0.86%/100bp. The only refinement to the
prior: the bank tilts into sovereigns in *quantity* (renationalisation direction) via
the financing/clearing channel, while deleveraging in market value.
