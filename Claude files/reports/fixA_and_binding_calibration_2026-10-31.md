# Fix A (Gertler–Karadi reserves) and the search for a binding calibration (B), 2026-10-31

Global solver only, branch `calibration`. Follows `tpi_mechanism_and_cross_border_2026-10-30.md` (§6 options A–E).

## 1. Fix A: the Eurosystem claim is no longer divertable

**Change.**
- When the Eurosystem buys D bonds, the D bank receives a safe claim on it (`Z_cb`, the reserves the purchase pays with).
- That claim used to count against the D bank's incentive constraint at the bond's divertable share (`lambda_bD_D`), so a
  purchase freed no room under the constraint at given prices.
- It now counts at `lambda_cb_D` (the divertable share of the Eurosystem claim), set to 0 as in Gertler–Karadi (2011),
  where central-bank intermediation is not subject to the bankers' diversion problem.
- Files: `calibration/global_projection.py` (the new key) and `point_map._incentive_constraint` (one term).
- `lambda_cb_D` = `lambda_bD_D` restores the old design; `tpi_purchase_tests.py` uses it for the comparison.

**What it changes.**
- **Baseline: nothing.** With the TPI off, the economy never holds a Eurosystem claim, and the rest point, impulse
  responses, decompositions and calibration table are unchanged to the printed digits.
  - The solved rules do change at the grid points that place a claim on the bank's balance sheet (the edges of the
    Eurosystem-book dimension): the old solution's residual there is 3.8e-2 under fix A.
  - The re-solve converges in 4 Newton steps (2.2e-10). The rest point and the p^d-2% impact read identically.
- **Purchases (one-off, rules held fixed):** each euro bought now frees `lambda_bD_D` euros of the constraint directly,
  so purchases bite wherever the constraint binds. The ceiling is unchanged: the most purchases can remove is the
  liquidity premium.

A one-off purchase on the fix-A solution (`runs/b_sweep_5bp.pkl`), rules held fixed, with the other states at the rest
point. Each cell gives the D bond price, the lending spread and D output:

| p^d | old design: buy 5% | fix A: buy 5% | either design: buy 15% |
|---|---|---|---|
| 2% (slack) | +0.43%, 0 bp, +0.001% | +0.43%, 0 bp, +0.001% | +2.01%, 0 bp, +0.005% |
| 5% | +0.60%, +12 bp, −0.010% | +0.79%, −118 bp, +0.059% | +3.17%, −173 bp (all of it), +0.078% |
| 12% | +1.26%, −167 bp, +0.078% | +1.38%, −239 bp, +0.116% | +5.58%, −693 bp (all of it), +0.327% |

- **Fix A where the constraint binds:** a 5% purchase now cuts the lending spread by 118 bp at p^d 5% and 239 bp at
  12%. Under the old design, on the same rules, the 5% purchase even raised the spread at 5% (+12 bp), because the
  claim it paid with crowded the constraint.
- **Where the constraint is slack** (2%), neither design has a real effect.
- **At 15%** both designs coincide: the constraint is gone (μ_D = 0) either way.
- These numbers come from the fix-A rules. They differ slightly from the 2026-10-30 report, which held the old rules
  fixed and patched the constraint.
- **`run.py global` under fix A (2026-11-01):**
  - the impact table and the bond split are identical;
  - the rest point and calibration table move by at most 3e-9 relative, and the risk IRF by at most 6e-8 (solver
    tolerance).
  - `tpi_purchase_tests.log`, `tpi_mechanism_check.log` (purchase columns only; the premium split is unchanged) and
    `model_side.log` (unchanged) were rerun on the new rules.

## 2. B: is there a pre-crisis lending-spread target at which the paper goes through?

**What "the paper goes through" means here.** SPEC.md's TPI story needs a constrained marginal holder: "the wedge
exists because the marginal holder is constrained; TPI relieves the constraint". So the sweep tests three things, by
exact reads on each re-solved model:
1. the D bank's incentive constraint binds at the stochastic rest point (μ_D, the multiplier on that constraint, > 0);
2. the acceptance signs at the p^d 2% headline shock: Q_bD (D bond price), n_D and n_F (D and F bank net worth), Y_D
   (D output) and C_D (D consumption) fall, and the lending spread rises;
3. a purchase of 5% of the D debt at that shock lowers the lending spread by at least 1 bp and raises output (fix A).

**The experiment.**
- `credit_spread_target_D` (the deterministic-SS D lending spread, the instrument appendix B sets to 5 bp/yr) is raised
  in rungs from 5 bp.
- Each rung is a warm, hard, fine-grid re-solve, with the other instruments (leverage, exposures, the F spread) and f_D
  (the banker exit rate) held.
- Driver: `Claude files/diagnostics/tex_calibration/b_spread_sweep.py`. Log: `runs/b_spread_sweep.log` and `.json`;
  solutions in `runs/b_sweep_<bp>bp.pkl`.

| SS spread instrument (bp/yr) | slack ratio at rest | E[Ω] at rest | IC binds from p^d | sovereign spread at rest | lev_D | 2% shock: lending spread | Y_D | C_D | n_F | buy 5% at 2%: lending spread |
|---|---|---|---|---|---|---|---|---|---|---|
| 5 (appendix B) | 1.125 | 1.117 | 3.93% | 62 bp | 7.98 | 0 | −0.031% | +0.120% | +0.026% | 0 |
| 10 | 1.118 | 1.122 | 3.72% | 66 bp | 7.96 | 0 | −0.032% | +0.122% | +0.031% | 0 |
| 20 | 1.110 | 1.140 | 3.49% | 74 bp | 7.93 | 0 | −0.032% | +0.122% | +0.036% | 0 |
| 40 | 1.108 | 1.194 | 3.44% | 93 bp | 7.89 | 0 | −0.031% | +0.118% | +0.040% | 0 |
| 60 | 1.114 | 1.265 | 3.66% | 112 bp | 7.86 | 0 | −0.028% | +0.111% | +0.039% | 0 |

- The slack ratio is E[Ω](1 + r)n / divertable assets; 1 means binding.
- E[Ω] is the D bank's expected value of a marginal euro of equity.
- The 80 bp rung stalled (Newton steps cut to 1/32 at max|F| 3.6e-5), and the sweep was stopped there.
- On the coarse grid (2026-10-28) the rest-point slack kept rising up to a 160 bp instrument, and the SS became
  infeasible near 320 bp.

**Result: no SS spread target makes the D bank bind at the rest point, or at the 2% shock.**
- The slack ratio bottoms out at about 1.108 near 40 bp, then rises again.
- The binding threshold never falls below p^d 3.4%, against the 2% headline shock.
- None of the three conditions holds at any rung.

**Why.** Raising the SS spread makes banking more profitable, above all in the future crisis states where the
constraint binds and banks earn large liquidity premia.
- That raises the franchise value today (E[Ω] 1.117 → 1.265 over 5 → 60 bp). A more valuable bank has more to lose by
  diverting, so the constraint relaxes.
- In the deterministic SS the two effects are consistent with binding, since that is how the instrument is defined.
- At the stochastic rest point, the anticipation of profitable crisis states outweighs the tighter SS from about 40 bp
  on.
- This is the same anticipation effect that leaves the appendix-B economy 12% slack at rest.

## 3. Two wrong signs that the calibration does not fix

At the 2% shock, D consumption (+0.12%) and F bank net worth (+0.03%) rise at every rung. Exact reads on the 5 bp
solution at larger shocks:

| p^d | μ_D | lending spread | Y_D | C_D | n_F | D deposit rate |
|---|---|---|---|---|---|---|
| 2% | 0 | 0 | −0.031% | +0.120% | +0.026% | −13 bp |
| 4% | 0.003 | +14 bp | −0.049% | +0.175% | −0.006% | −24 bp |
| 5% | 0.041 | +173 bp | −0.135% | +0.087% | +0.013% | −84 bp |
| 8% | 0.114 | +472 bp | −0.294% | −0.069% | +0.009% | −195 bp |
| 12% | 0.173 | +693 bp | −0.409% | −0.170% | −0.039% | −276 bp |

- **D consumption** falls only from about 8%.
- **F bank net worth** stays within ±0.04%, economically zero, at every level.
- **The common factor is the union deposit rate,** which falls with D risk. That lowers F banks' funding cost and
  supports D consumption.
- This needs its own diagnosis (a consumption decomposition via `measure.py`) and is independent of B.

## 4. What this leaves for the paper

The D constraint binds only in crisis states (p^d above about 3.4–3.9%), whatever the spread target. Three routes
remain, all author decisions:
1. **Crisis-sized headline shocks.** No recalibration. The 2010–11 Greek average p^d is 9.1%/qtr, inside the binding
   region. There the Bocola pass-through operates (at 8%: lending spread +472 bp, Y_D −0.29%), purchases bite under
   fix A (§1), and C_D falls. n_F stays near zero.
2. **Another calibration lever** to move the threshold below 2%. f_D is the candidate. On the earlier narrow box it moved
   the threshold little (1.27–1.36% over f_D 0.06–0.10); untested on the final box.
3. **Accept occasional binding as the result:** the TPI repairs transmission only in crises, which matches its ECB
   mandate ("unwarranted, disorderly" spread moves).

SPEC.md's always-binding framing (written for the SSJ model) does not carry over to the global model at the appendix-B
targets.

## 5. Does fix A make the TPI work?

Partly. It repairs the transmission channel but not the cap.
- **Works:** wherever the D constraint binds (p^d above about 3.9%), purchases now relieve it.
  - At p^d 5%, a 5% purchase moves the lending spread −118 bp instead of +12 bp, and output +0.059% instead of −0.010%.
  - At 12%: −239 bp instead of −167 bp.
- **Does not work:**
  - **The 200 bp cap is still out of reach.** Purchases remove only the liquidity premium. At p^d 12% the cap needs a
    D bond price about 30% higher, and a 5% purchase adds about 1.4%.
  - **No effect at the 2% headline shock.** The constraint is slack there, and §2 shows no spread calibration changes
    that.
- **Not yet tested:** the full TPI equilibrium. `RUN_TPI` is still False, and the cap walk has not been re-run under
  fix A. The announcement effect of a credible floor would add to the one-off numbers, but it cannot remove the
  expected default loss.
