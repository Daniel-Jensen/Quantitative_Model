---
---


## Handoff: 2026-06-11T15:17:13Z (auto-saved before compaction)


### Compaction Metadata
- Trigger: auto
- Custom instructions: (none)
- Transcript: /Users/Adam/.claude/projects/-Users-Adam-Documents-uni-phd-research-QUANTITATIVE-MODEL/1074cac6-aaa6-48c0-b68f-807e7b331a27.jsonl
- CWD: /Users/Adam/Documents/uni/phd/research/QUANTITATIVE_MODEL


### Last User Message (transcript tail)
(unavailable)


### Last Assistant Message (transcript tail)
(unavailable)


### Git Snapshot
- Branch: endo-def
- Recent commits:
  - 1da134e Endogenous default: use actual Y in Debt/GDP denominator
  - edf1115 Consolidate to single HANDOFF.md; fix stale branch/PR references
  - d5f70e8 Fill HANDOFF.md model summary and handoff context placeholders
  - a0ddc18 Restructure repo layout and add canonical post-audit documentation
  - 52f17d5 Add forensic audit reports and reproduction scripts


### Model Summary

- Two-country HANK model (D=Greece, F=Germany), SSJ framework, GHH preferences, GK financial intermediaries, Hatchondo-Martinez perpetuity bonds (decay δ_b=0.10)
- All 6 audit fixes applied on `AB-audit` branch (PR #27); `endo-def` branches from `AB-audit`
- Key fixes: predetermined deposit rate (T-2 timing), Walras residuals both satisfied post-fix, `financial_solved_D/F` inner `.solved()` blocks handle `{nu_K, nu_b, eta, theta}` — these must NOT appear in outer `unknowns_tp`
- `endo-def` feature: `def_rate_D = shock_def_D + def_scale_D * (f(b_gov_D(-1)/Y_D) - f(ss_ratio))` — actual `Y_D` in denominator replaces constant `Y_ss_D`, so output contractions raise Debt/GDP endogenously; SS unchanged (Y_D=Y_ss_D=1.0 at SS)
- Repo structure: model code in `code/`, figures in `plots/`, docs in `docs/`; canonical equations in `code/equations_D.py`, `code/equations_F.py`, `code/equations_global.py`
- SS Walras residuals: `ca_res_D` and `goods_mkt_F` both hold to solver tolerance (verified post-fix)
- IRF comparison script at `audit_artifacts/endo_def_comparison.py` — compares endo-def (Y_D) vs AB-audit (Y_ss) Jacobians
- `unknowns_tp` (23 vars): `['K_D','n_inter_D','div_D','I_D','Q_D','b_gov_D','N_D','b_F_D','w_D','rdep_D','K_F','n_inter_F','div_F','I_F','Q_F','b_gov_F','N_F','b_D_F','w_F','rdep_F','p','q_b_D','q_b_F']`
- `targets_tp` (23): `['deposit_mkt_D','K_res_D','n_inter_val_D','div_res_D','capital_res_D','q_res_D','b_gov_res_D','b_F_D_res','labor_mkt_res_D','w_res_D','deposit_mkt_F','K_res_F','n_inter_val_F','div_res_F','capital_res_F','q_res_F','b_gov_res_F','b_D_F_res','labor_mkt_res_F','w_res_F','goods_mkt_D','rb_D_res','rb_F_res']`
- Savefig paths in `code/model_v12.ipynb` still use bare filenames (save to `code/` not `plots/`) — needs fix


### Handoff Context (paste into next session)

**Branch**: `endo-def` (contains all AB-audit fixes + endogenous default feature)

**Immediate pending tasks:**
1. **IRF comparison**: `audit_artifacts/endo_def_comparison.py` was running when session ended. Check `audit_artifacts/endo_def_log.txt` for results. If failed, rerun with `/opt/anaconda3/envs/ssj/bin/python audit_artifacts/endo_def_comparison.py > audit_artifacts/endo_def_log.txt 2>&1`. The script compares `def_rate_D = f(b_gov(-1)/Y_D)` vs `f(b_gov(-1)/Y_ss)` IRFs under `shock_def_D` and `Z_D` shocks.
2. **Fix savefig paths**: In `code/model_v12.ipynb`, all `fig*.savefig('fig_tpi_*.png', ...)` calls need `../plots/` prefix, since the notebook is now in `code/` not the repo root.
3. **Port bank-cal calibration values**: `bank-cal` branch has updated `phi_lamb` and `def_scale` parameters. Do not merge; port values manually into `code/model_v12.ipynb`.

**Key invariants to preserve:**
- `theta_D`/`theta_F` must NOT appear in `unknowns_tp` — handled inside `financial_solved_D/F` inner blocks
- SS rest point of `def_rate_D/F = 0` is maintained because `debt_ratio - ss_ratio = 0` at SS
- Python env: `/opt/anaconda3/envs/ssj/bin/python` (base env has broken liblapack symlink)
- `code` is a reserved Python module name — use `sys.path.insert(0, str(CODE_DIR))` + bare imports, not `from code.equations_D import`

---
