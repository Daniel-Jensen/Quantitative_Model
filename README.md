Code for first project 

Link to overleaf: 


https://www.overleaf.com/project/698b4f88aeef1d0e1d08cc0c


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
