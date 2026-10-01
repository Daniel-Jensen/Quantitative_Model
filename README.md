Code for first project 

Link to overleaf: 


https://www.overleaf.com/project/698b4f88aeef1d0e1d08cc0c

## Running the models

Two solvers of the same model, one entry point, [run.py](run.py):

```bash
python3 run.py                    # MODEL set in run.py's configuration block
python3 run.py ssj                # linear sequence-space model   (linear_ssj/)
python3 run.py global [--quick]   # nonlinear global projection   (global_projection/)
python3 run.py both   [--quick]   # both, then the side-by-side comparison
python3 run.py both --plots-only  # redraw every figure from the saved data
python3 run.py compare            # re-render the comparison from the saved data
```

Each run computes first, saves everything to `results/<SSJ|GLOBAL>/data/`, then draws the
figures from that saved data into `results/<SSJ|GLOBAL>/figures/`; the console log is
`results/<SSJ|GLOBAL>/run.log` and the comparison is `results/COMPARISON/`. The SSJ model
needs its own interpreter: set `SSJ_PYTHON` at the top of `run.py` or in the environment.

Documentation, notes, reports and tests written by Claude are kept apart in
[`Claude files/`](<Claude files>), at the path each would otherwise have had (start with
[REFACTOR_ARCHITECTURE.md](<Claude files/REFACTOR_ARCHITECTURE.md>) for the repository map).

## Notebook hygiene setup

`.gitattributes` declares `nbstripout` as the clean filter and `nbdime` as the
diff/merge driver for `*.ipynb` files. The filters only take effect after each
collaborator runs the setup commands once per local clone:

```bash
pip install nbstripout nbdime
nbstripout --install        # registers clean filter in .git/config
nbdime config-git --enable  # registers diff/merge drivers in .git/config
```

After setup, committed notebooks will have outputs, execution counts, and
volatile metadata stripped automatically, and `git diff` / `git merge` on
notebooks will use nbdime's cell-aware rendering.
