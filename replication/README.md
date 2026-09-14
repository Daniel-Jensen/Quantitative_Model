# replication/

Third-party replication packages, kept separate from this project's own code
(`code/`) and write-up (`docs/`). Nothing here is authored by this project.

## `bocola2016/`

Luigi Bocola, "The Pass-Through of Sovereign Risk", *Journal of Political
Economy* 124(4), 2016 — the author's official MATLAB replication package,
downloaded unmodified. It is the reference implementation for this project's
default mechanism (see the root `CLAUDE.md`).

Layout:

| Path | Contents |
|------|----------|
| `Model/` | Model solution, IRFs, LTRO experiment, risk-premium decomposition |
| `Estimation_Step1/`, `Estimation_Step2/` | Two-step Bayesian estimation, particle filter, Smolyak solution files |
| `Cross_section/` | Compustat/Fama-French cross-sectional evidence |
| `Figures and Tables/` | Scripts reproducing each published figure and table |
| `Data.xls` | Source data |
| `Readme.txt` | The author's own readme — start there |

**Not in git:** `Matfiles/` (165 MB of solved-model output, one file above
GitHub's 100 MB per-file limit) is gitignored. The `.m` scripts in `Model/`
and `Estimation_Step*/` regenerate it. The smaller per-subdirectory
`Matfiles/` folders are tracked.
