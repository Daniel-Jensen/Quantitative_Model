# MEASURE THE GLOBAL MODEL'S STATE BOXES FROM A SOLVED RULE SET AND WRITE THEIR SPEC (2026-10-21).
# The normal box covers normal times: the stochastic rest point, the headline IRF paths and one-step
# successors from interior points. The post-default box covers the state a default produces from
# any of those points and the straight segment from there back to the rest point (the recovery).
# Each box puts the D-bank block (capital, deposit obligations, the D bank's D-bond book, the
# cross-border deposit position) on the principal axes of its cloud, in net-worth units
# (state_grid.align_weights), with a margin; the other stocks get fractional bands from the cloud,
# never narrower than the old V5 box. The blend hands the continuation from the post-default rules
# to the normal ones along the line from the post-default cluster to the SS, inside the normal box.
# usage: python3 "Claude files/diagnostics/state_box/measure_box.py" RULES.pkl [OUT.json] [--hard-ss]
import json
import os
import pickle
import sys
from datetime import date

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

import numpy as np

from calibration.global_projection import get_calibration
from global_projection.steady_state import solve_steady_state
from global_projection.solver_recursive.state_grid import (s_process_params, align_weights, _box_centre,
                                                           ALIGN_BLOCK, IK_F, IP_F, IBDF, IBFD, IZ, IS,
                                                           IM, IO)
from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
from global_projection.solver_recursive.decision_rules import SOLVE
import global_projection.solver_recursive.point_map as PM
import global_projection.solver_recursive.recursive_experiment as RE

MARGIN = 0.15  # of each axis' range, either side
MIN_HALF = 0.01  # the narrowest half-width, in net-worth units
# the old V5 bands, the floor for the other stocks' bands
V5 = dict(K_F=0.02, P_F=0.0226, b_DF=0.12, b_FD=0.12, Z_D=0.03)
# the blend: the post-default rules hold up to 90% of the normal box's reach along the recovery
# line and the normal rules from 40% of it
BLEND_POST, BLEND_NORMAL = 0.9, 0.4


def _successors(rules, cal, ss, sp, S, d, eps):
    # THE NEXT STATES AT THE GIVEN s' INNOVATIONS FROM THE PERIOD MAP AT S IN REGIME d (FITTED RULES).
    x = np.array([float(rules.eval(k, d, np.atleast_2d(S))[0]) for k in SOLVE])
    v = PM._read_state_and_regime(S, d, x, rules, cal, ss)
    PM._production(v, cal)
    PM._government(v, cal, ss)
    PM._balance_sheets(v, rules, cal)
    return PM._next_states(v, sp, eps)[0]


def _clouds(rules, cal, ss, sp, seed=11):
    # THE NORMAL CLOUD AND THE POST-DEFAULT CLOUD, IN NATURAL UNITS.
    g = rules.grid
    eps, _ = PM.gh_nodes(5)
    S0 = RE.stochastic_rest_point(rules, cal, ss, sp, verbose=False)
    irf = []
    for pd0 in (0.01, 0.02):
        S, s0 = S0.copy(), RE.s_from_pd(pd0)
        for t in range(40):
            S[IS] = sp["s_star"] + (s0 - sp["s_star"]) * sp["rho_s"] ** t
            irf.append(S.copy())
            S = RE.advance(RE.read_at(rules, cal, ss, sp, S.copy()), S, sp)
    irf = np.array(irf)
    irf_in = irf[np.all(np.abs(g.to_unit(irf)) <= 1.0, axis=1)]
    U = np.random.default_rng(seed).uniform(-0.8, 0.8, size=(500, g.d))
    U[:, [IM, IO]] = 0.0
    interior = g.from_unit(U)
    starts = np.vstack([interior, irf_in])
    succ = np.vstack([_successors(rules, cal, ss, sp, S, 0, eps) for S in starts])
    normal = np.vstack([S0[None], irf, succ])
    post = []
    for S in np.vstack([g.points[np.abs(g.points[:, IM]) < 1e-12], irf_in, interior[:200]]):
        Sp = _successors(rules, cal, ss, sp, S, 1, np.zeros(1))[0]
        if np.all(np.isfinite(Sp)):
            post.append(Sp)
    post = np.array(post)
    # the recovery runs from each post-default state back towards the rest point
    seg = [Sp + t * (S0 - Sp) for Sp in post for t in (0.25, 0.5, 0.75, 1.0)]
    return normal, np.vstack([post, seg]), post, S0


def _box_spec(cloud, ss, cal, c):
    # PRINCIPAL AXES AND BOUNDS FOR THE D-BANK BLOCK, AND FRACTIONAL BANDS FOR THE OTHER STOCKS.
    blk = list(ALIGN_BLOCK)
    cen = np.array([c["K_D"], c["P_D"], c["b_DD"], 0.0])
    Wd = (cloud[:, blk] - cen) * align_weights(ss, cal)
    _, sv, Vt = np.linalg.svd(Wd - Wd.mean(axis=0), full_matrices=False)
    Z = Wd @ Vt.T
    lo, hi = Z.min(axis=0), Z.max(axis=0)
    w = hi - lo
    lo, hi = lo - MARGIN * w, hi + MARGIN * w
    mid = 0.5 * (lo + hi)
    half = np.maximum(0.5 * (hi - lo), MIN_HALF)
    lo, hi = mid - half, mid + half
    bands = {}
    for i, key in ((IK_F, "K_F"), (IP_F, "P_F"), (IBDF, "b_DF"), (IBFD, "b_FD"), (IZ, "Z_D")):
        ref = {"K_F": c["K_F"], "P_F": c["P_F"], "b_DF": c["b_DF"], "b_FD": c["b_FD"], "Z_D": c["Z_D"]}[key]
        f = cloud[:, i] / ref - 1.0
        flo, fhi = f.min(), f.max()
        wf = fhi - flo
        bands[key] = [min(flo - MARGIN * wf, -V5[key]), max(fhi + MARGIN * wf, V5[key])]
    # no box reaches short positions in the D bank's F book: the no-short-sale pair holds at every
    # node, so a state below zero is an interpolation artefact (read through the box edge instead)
    bands["b_FD"][0] = max(bands["b_FD"][0], -1.0)
    bands.update(M_cb=0.32, O_cb=0.25)
    return dict(axes=Vt.tolist(), lo=lo.tolist(), hi=hi.tolist(), bands=bands,
                axis_sd=(sv / np.sqrt(len(cloud))).tolist())


def _blend(normal_spec, post_states, ss, cal, c):
    # THE RECOVERY LINE (UNIT VECTOR FROM THE POST-DEFAULT CLUSTER TO THE SS) AND THE HANDOVER POINTS.
    blk = list(ALIGN_BLOCK)
    cen = np.array([c["K_D"], c["P_D"], c["b_DD"], 0.0])
    q2 = ((post_states[:, blk] - cen) * align_weights(ss, cal)).mean(axis=0)
    u = -q2 / np.linalg.norm(q2)
    # how far the normal box reaches from its centre towards the post-default cluster
    A = np.asarray(normal_spec["axes"])
    lo, hi = np.asarray(normal_spec["lo"]), np.asarray(normal_spec["hi"])
    a = A @ (-u)
    reach = np.min([(hi[i] if a[i] > 0 else lo[i]) / a[i] for i in range(len(a)) if abs(a[i]) > 1e-12])
    return dict(direction=u.tolist(), post_until=-BLEND_POST * reach, normal_from=-BLEND_NORMAL * reach,
                cluster=float(-np.linalg.norm(q2)), reach=float(reach))


def main():
    # MEASURE AND WRITE THE SPEC.
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    path = args[0]
    out = args[1] if len(args) > 1 else os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))),
        "global_projection", "solver_recursive", "box_specs.json")
    d_ = pickle.load(open(path, "rb"))
    rules = d_["rules"] if isinstance(d_, dict) else d_
    cal = get_calibration()
    cal["nw_floor_frac"] = 0.15
    if isinstance(d_, dict) and "mu_eps" in d_:
        cal["mu_eps"] = float(d_["mu_eps"])
    if "--hard-ss" in sys.argv:  # the rule set was solved on the SS for the hard max
        e = cal["mu_eps"]; cal["mu_eps"] = 0.0
        ss = solve_steady_state(cal, verbose=False)
        cal["mu_eps"] = e
    else:
        ss = solve_steady_state(cal, verbose=False)
    sp = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sp)
    c = _box_centre(ss, cal)
    normal, postc, post_states, S0 = _clouds(rules, cal, ss, sp)
    spec_n = _box_spec(normal, ss, cal, c)
    spec_p = _box_spec(postc, ss, cal, c)
    blend = _blend(spec_n, post_states, ss, cal, c)
    doc = dict(normal=spec_n, post=spec_p, blend=blend,
               measured_from=os.path.basename(path), date=str(date.today()),
               note="D-bank block (K_D, P_D, b_DD, V_dep) in net-worth units (state_grid.align_weights); "
                    "other bands are fractions of the SS stock; regenerate with this script")
    with open(out, "w") as fh:
        json.dump(doc, fh, indent=1)
    names = ("K_D", "P_D", "b_DD", "V_dep")
    for lab, spec in (("normal", spec_n), ("post-default", spec_p)):
        print(f"{lab} box:")
        for k in range(4):
            print(f"   axis {k} (sd {spec['axis_sd'][k]:.3f}): " + " ".join(
                f"{n} {spec['axes'][k][j]:+.3f}" for j, n in enumerate(names))
                  + f"   [{spec['lo'][k]:+.3f}, {spec['hi'][k]:+.3f}]")
        print("   bands: " + "  ".join(f"{k} {v}" if not isinstance(v, list) else f"{k} [{v[0]:+.3f},{v[1]:+.3f}]"
                                       for k, v in spec["bands"].items()))
    print(f"blend: recovery line {np.round(blend['direction'], 3).tolist()}, cluster at {blend['cluster']:+.3f}, "
          f"normal box reach {blend['reach']:.3f}: post rules up to {blend['post_until']:+.3f}, "
          f"normal from {blend['normal_from']:+.3f} (net-worth units)")
    print(f"written: {out}")


if __name__ == "__main__":
    main()
