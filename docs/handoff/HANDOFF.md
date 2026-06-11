---

## Handoff: 2026-06-11T08:44:22Z (auto-saved before compaction)

### Compaction Metadata
- Trigger: auto
- Custom instructions: (none)
- Transcript: /Users/Adam/.claude/projects/-Users-Adam-Documents-uni-phd-research-QUANTITATIVE-MODEL/1074cac6-aaa6-48c0-b68f-807e7b331a27.jsonl
- CWD: /Users/Adam/Documents/uni/phd/research/QUANTITATIVE_MODEL

### Last User Message (transcript tail)
(unavailable)

### Last Assistant Message (transcript tail)
Note for DJ review: PR #26 against `main` will show the full divergence since `audit` branched from current main — the fix commits are clean, but the docs commit references `bank-cal` calibration porting as a *follow-up* (not done here, per your scope).

### Git Snapshot
- Branch: audit
- Status:
 M fig_tpi_bond_price.png
 M fig_tpi_effectiveness.png
 M fig_tpi_spread_mitigation.png
 M fig_tpi_utility.png
 M fig_tpi_welfare_bar.png
 M fig_tpi_welfare_macro.png
 M fig_tpi_welfare_spread.png
 D model_notes
 D model_v11.ipynb
?? "OLD models /model_v11.ipynb"
?? audit_artifacts/bankcal_stability_log.txt
?? audit_artifacts/fix_log.txt
?? audit_artifacts/irfs.npz
?? audit_artifacts/irfs_fixed.npz
?? audit_artifacts/irfs_stable.npz
?? audit_artifacts/nb_cells.py
?? audit_artifacts/philamb_log.txt
?? audit_artifacts/philamb_log2.txt
?? audit_artifacts/run1_log.txt
?? audit_artifacts/run2_postfix_log.txt
?? audit_artifacts/run3_allfix_log.txt
?? audit_artifacts/run4_kt_log.txt
?? audit_artifacts/tpi_log.txt
?? audit_artifacts/tpi_log2.txt
?? docs/handoff/
?? docs/model_notes
- Recent commits:
52f17d5 Add forensic audit reports and reproduction scripts
4c810e1 Rework W-1 fix: keep Y = F(K_t), pay new-capital product to capital producer
396cbd9 Fix accounting and timing bugs found in forensic audit (W-1, W-2, W-3, T-2, A-2, TPI-1)
7e94175 Merge pull request #25 from Daniel-Jensen/Sunday-fixes
754ea86 Add TPI policy analysis: closed-loop spread rule, welfare figures, and equation fixes

### Model Summary

- Two-country HANK model (D=Greece, F=Germany) with GK bank intermediation, sovereign default risk, and TPI policy extension.
- Six structural bugs fixed on `audit` branch (now `AB-audit`): W-1 capital income timing, W-2 F-bank ToT revaluation, W-3 F-bank FOC p-conversion, T-2 predetermined deposit rate, A-2 smart-steady m formula, TPI-1 CB budget remittance.
- Post-fix Walras: goods_mkt_F ≤8e−10, ca_res_D ≤5.8e−8 on 1% TFP shock (were 2.0e−2 / 1.5e−4 pre-fix).
- Default shock signs correct post-T-2-fix: n_inter_D[0]=−3.5%, Y_D[0]=−2.5e−4 (were positive/perverse).
- TPI post-fix: ΔW_D=+1.88%, ΔW_F=−1.90% at γ=10 (zero-sum burden shift; spread not compressed). All pre-fix TPI figures stale.
- phi_lamb=0.15 required for stability at current amplification (psi_lambda_B=3, def_scale=0.25); Bohn=0.60/yr vs literature 0.10–0.15/yr — tension documented, re-mapping required.
- Open author decisions: C-1 (Delta_cross=1.45>1; hardcode 0.2/0.4 per bank-cal), S-1 (writeoff_enabled=0; no realized losses).
- Repo restructured: equations → `code/`, figures → `plots/`, canonical docs (STATE/SPEC/PROCESS/HANDOFF) added to `docs/`.
- PR #26 open (AB-audit → main); bank-cal branch must not be merged (port values only).

### Handoff Context (paste into next session)

Active branch: `AB-audit`. Active notebook: `code/model_v12.ipynb`.

Key calibration (cell `f042652d` or equivalent):
- `phi_lamb_D = phi_lamb_F = 0.15` — stable at current amplification; needs re-mapping after calibration port
- `def_scale_D = 0.25`, `psi_lambda_B = 3.0`, `delta_b = 0.10`, `f = 0.12`, `theta = 4.0`
- `writeoff_enabled = 0`, `Delta_cross = 1.45` (both open; see C-1, S-1 in STATE.md)

Next task: port bank-cal calibration values onto AB-audit (see `docs/bank_cal_review.md` §Recommendation):
- `delta_b_D=0.036`, `delta_b_F=0.038`, `f=0.03`, `Delta_D=0.2`/`Delta_F=0.4` (hardcoded, resolves C-1)
- `recovery_rate=0.40`, `zeta_writeoff=1.0`, EBA bilateral shares, decide S-1
- After porting: re-map (phi_lamb, def_scale) bifurcation using `audit_artifacts/philamb_test.py` as template

Run environment: `/opt/anaconda3/envs/ssj/bin/python`. Each Jacobian solve ~3 min. Regression test: `audit_artifacts/run_audit.py`.

---
