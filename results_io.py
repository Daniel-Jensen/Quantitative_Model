"""Pipeline results on disk as plain data: arrays -> <stem>.npz, everything else -> <stem>.json.

Shared by run.py (global model, in-process) and linear_ssj/main.py (SSJ model, its own
interpreter), so it needs only the standard library and numpy and imports neither solver.

    save(stem, obj)   nested dicts / lists / tuples of arrays and scalars
    load(stem)        the same structure back; arrays and floats round-trip bit-identically

Dicts whose keys are all strings stay plain JSON objects (readable as-is); any other
dict keeps its key types (e.g. float TPI gammas) as a list of [key, value] pairs.
"""
import json
import shutil
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def model_dirs(base, fresh=False):
    """<base>/data and <base>/figures (created) and the <base>/run.log path.

    fresh=True empties data/ and figures/ first, so a run's outputs are one consistent
    set: a quick preview must not leave a previous full run's policy figures beside it.
    """
    base = Path(base)
    for sub in ("data", "figures"):
        d = base / sub
        if fresh and d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True, exist_ok=True)
    return base / "data", base / "figures", base / "run.log"


def _encode(obj, arrays, path):
    if isinstance(obj, np.ndarray) and obj.dtype != object:
        arrays[path] = obj
        return {"__array__": path}
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    if isinstance(obj, (list, tuple, np.ndarray)):
        items = [_encode(v, arrays, f"{path}/{i}") for i, v in enumerate(obj)]
        return items if isinstance(obj, list) else {"__tuple__": items}
    if hasattr(obj, "items"):   # dict, and dict-likes such as SSJ's ImpulseDict
        pairs = list(obj.items())
        if all(isinstance(k, str) for k, _ in pairs):
            return {k: _encode(v, arrays, f"{path}/{k}") for k, v in pairs}
        return {"__pairs__": [[_encode(k, arrays, path), _encode(v, arrays, f"{path}/{k}")]
                              for k, v in pairs]}
    return {"__repr__": repr(obj)}   # not restorable; recorded so the file still says what it was


def _decode(obj, arrays):
    if isinstance(obj, list):
        return [_decode(v, arrays) for v in obj]
    if not isinstance(obj, dict):
        return obj
    if "__array__" in obj:
        return arrays[obj["__array__"]]
    if "__tuple__" in obj:
        return tuple(_decode(v, arrays) for v in obj["__tuple__"])
    if "__pairs__" in obj:
        return {_decode(k, arrays): _decode(v, arrays) for k, v in obj["__pairs__"]}
    return {k: _decode(v, arrays) for k, v in obj.items()}


def save(stem, obj):
    """Write obj as <stem>.json (+ <stem>.npz when it holds arrays); return the json path."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    arrays = {}
    tree = _encode(obj, arrays, "")
    stem.with_suffix(".json").write_text(json.dumps(tree, indent=1))
    npz = stem.with_suffix(".npz")
    if arrays:
        np.savez(npz, **arrays)
    elif npz.exists():
        npz.unlink()   # a stale array file from an earlier run must not pair with this json
    return stem.with_suffix(".json")


def load(stem):
    """Read what save() wrote."""
    stem = Path(stem)
    npz = stem.with_suffix(".npz")
    arrays = dict(np.load(npz, allow_pickle=False)) if npz.exists() else {}
    return _decode(json.loads(stem.with_suffix(".json").read_text()), arrays)


if __name__ == "__main__":
    # round-trip self-check on every container shape the pipelines save
    import tempfile
    x = {"irf": {"Y": np.linspace(0, 1, 7), "C": np.arange(3.0)},
         "gammas": [0, 2.0, 5], "by_gamma": {0: {"q": np.ones(2)}, 2.0: {"q": np.zeros(2)}},
         "scen": [("label", {"Y_D": np.array([0.1, np.nan])}, (0.1, 0.2, 0.3, 1.0))],
         "scalar": np.float64(1 / 3), "flag": True, "none": None, "nested": {"a": {"b": [1, "x"]}}}
    with tempfile.TemporaryDirectory() as d:
        save(Path(d) / "t", x)
        y = load(Path(d) / "t")
    assert np.array_equal(y["irf"]["Y"], x["irf"]["Y"]) and y["irf"]["C"].dtype == np.float64
    assert set(y["by_gamma"]) == {0, 2.0} and isinstance(list(y["by_gamma"])[1], float)
    assert y["scen"][0][0] == "label" and isinstance(y["scen"][0], tuple)
    assert np.array_equal(y["scen"][0][1]["Y_D"], x["scen"][0][1]["Y_D"], equal_nan=True)
    assert y["scen"][0][2] == (0.1, 0.2, 0.3, 1.0) and y["scalar"] == 1 / 3
    assert y["gammas"] == [0, 2.0, 5] and y["nested"] == {"a": {"b": [1, "x"]}} and y["none"] is None
    print("results_io: round-trip OK")
