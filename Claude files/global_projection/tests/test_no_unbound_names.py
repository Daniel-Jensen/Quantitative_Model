# STATIC GATE: NO FUNCTION MAY USE A NAME THAT NOTHING BINDS.
# This is here because that bug shipped and cost a four-hour run. An earlier multi-part
# edit failed its last assertion, so NONE of its parts were written; the retry re-applied
# only the part that had failed, leaving `solved_ok` referenced in ltro_experiment.run
# and defined nowhere. Every import succeeded, every test passed, and the pipeline ran
# for an hour before reaching the line and raising NameError. Nothing else in the suite
# looks at a function that is never called by a test.
#
# The check must handle four things a naive version gets wrong -- each produced a false
# positive on this package, and a checker that cries wolf is worse than none:
#   `import x as y` binds y not x; tuple unpacking binds every target; a NESTED function
#   sees its enclosing function's locals; and a LAMBDA is its own scope.
import ast, builtins, os, pathlib, sys

ROOT = pathlib.Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PKGS = ("solver_recursive", "blocks", "config", "reporting")
BUILTINS = set(dir(builtins))


def targets(node, out):
    if isinstance(node, ast.Name): out.add(node.id)
    elif isinstance(node, (ast.Tuple, ast.List)):
        for e in node.elts: targets(e, out)
    elif isinstance(node, ast.Starred): targets(node.value, out)

def arg_names(a):
    out = {x.arg for x in list(a.args) + list(a.kwonlyargs) + list(a.posonlyargs)}
    for v in (a.vararg, a.kwarg):
        if v: out.add(v.arg)
    return out

def free_in(node, bound):
    u, stack = set(), [node]
    while stack:
        x = stack.pop()
        if isinstance(x, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue
        if isinstance(x, ast.Lambda):
            u |= free_in(x.body, bound | arg_names(x.args))
            continue
        if isinstance(x, ast.Name) and isinstance(x.ctx, ast.Load): u.add(x.id)
        stack.extend(ast.iter_child_nodes(x))
    return u - bound

def bound_here(fn):
    b = arg_names(fn.args)
    stack = list(fn.body)
    while stack:
        x = stack.pop()
        if isinstance(x, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            b.add(x.name); continue
        if isinstance(x, ast.Lambda): continue
        if isinstance(x, ast.Name) and isinstance(x.ctx, ast.Store): b.add(x.id)
        elif isinstance(x, (ast.Import, ast.ImportFrom)):
            b |= {n.asname or n.name.split('.')[0] for n in x.names}
        elif isinstance(x, ast.ExceptHandler) and x.name: b.add(x.name)
        elif isinstance(x, ast.comprehension): targets(x.target, b)
        elif isinstance(x, (ast.Global, ast.Nonlocal)): b |= set(x.names)
        stack.extend(ast.iter_child_nodes(x))
    return b

def used_here(fn):
    u = set()
    for st in fn.body: u |= free_in(st, set())
    return u

def walk_fns(body, enclosing, path, bad):
    for x in body:
        if isinstance(x, (ast.FunctionDef, ast.AsyncFunctionDef)):
            scope = enclosing | bound_here(x)
            free = used_here(x) - scope - BUILTINS
            if free: bad.append(f"  {path}::{x.name}  UNBOUND {sorted(free)}")
            walk_fns(x.body, scope, path, bad)
        elif isinstance(x, ast.ClassDef):
            walk_fns(x.body, enclosing, path, bad)



def _scan(paths):
    bad = []
    for path in paths:
        tree = ast.parse(pathlib.Path(path).read_text())
        mod = set()
        for x in ast.walk(tree):
            if isinstance(x, (ast.Import, ast.ImportFrom)):
                mod |= {n.asname or n.name.split(".")[0] for n in x.names}
            elif isinstance(x, ast.Assign):
                for tg in x.targets: targets(tg, mod)
            elif isinstance(x, (ast.AnnAssign, ast.AugAssign)): targets(x.target, mod)
            elif isinstance(x, (ast.FunctionDef, ast.ClassDef)): mod.add(x.name)
            elif isinstance(x, ast.For): targets(x.target, mod)
            elif isinstance(x, ast.withitem) and x.optional_vars:
                targets(x.optional_vars, mod)
        walk_fns(tree.body, mod, os.path.relpath(path, ROOT), bad)
    return bad


def test_no_unbound_names():
    # EVERY MODULE ON THE EXECUTION PATH, not just the ones a test imports.
    files = [str(ROOT / "main.py")]
    for pkg in PKGS:
        files += [str(p) for p in sorted((ROOT / pkg).glob("*.py"))]
    bad = _scan(files)
    assert not bad, "unbound names:\n" + "\n".join(bad)
    print(f"  {len(files)} modules scanned, no unbound names: PASSED")


def test_the_check_actually_catches_one():
    # A GATE THAT CANNOT FAIL IS NOT A GATE. Plant the exact defect that shipped -- a
    # name used in a function body and bound nowhere -- and require the scan to find it.
    import tempfile
    src = ("def run():\n"
           "    acc = []\n"
           "    for x in (1, 2):\n"
           "        acc.append(x)\n"
           "    return [v for v in never_bound]\n")
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
        fh.write(src); tmp = fh.name
    try:
        bad = _scan([tmp])
        assert bad and "never_bound" in bad[0], f"the check missed a planted defect: {bad}"
    finally:
        os.unlink(tmp)
    print("  planted defect detected: PASSED")


if __name__ == "__main__":
    test_no_unbound_names()
    test_the_check_actually_catches_one()
    print("test_no_unbound_names: ALL PASSED")
