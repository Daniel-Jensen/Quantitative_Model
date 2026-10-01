"""The two solvers must stay isolated: linear_ssj/ and global_projection/ never import
each other, and the shared calibration/ package imports neither (so either interpreter
can load it). Pure AST, no model solve; runs under any python3:

    python3 -m pytest "Claude files/test_isolation.py"   (or run it as a script)
"""
import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]   # repo root
FORBIDDEN = {"linear_ssj": {"global_projection"},
             "global_projection": {"linear_ssj"},
             "calibration": {"linear_ssj", "global_projection"}}


def imported_tops(path):
    for n in ast.walk(ast.parse(path.read_text())):
        if isinstance(n, ast.Import):
            yield from (a.name.split(".")[0] for a in n.names)
        elif isinstance(n, ast.ImportFrom) and n.level == 0 and n.module:
            yield n.module.split(".")[0]


def test_solvers_do_not_import_each_other():
    hits = [f"{f.relative_to(ROOT)} imports {top}"
            for pkg, banned in FORBIDDEN.items()
            for f in sorted((ROOT / pkg).rglob("*.py"))
            for top in imported_tops(f) if top in banned]
    assert not hits, "cross-solver imports:\n" + "\n".join(hits)


if __name__ == "__main__":
    test_solvers_do_not_import_each_other()
    print("test_isolation: PASSED")
