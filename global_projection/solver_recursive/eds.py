# EPSILON-DISTINGUISHABLE-SET (EDS) GRID FOR THE RECURSIVE SOLVER (Maliar-Maliar).
# The Smolyak BOX collocates ~40-50% infeasible corners at mu=2 (the near-unit-root
# capital forces a wide box; the post-default d=1 regime has no equilibrium at its
# cross-corners). Those corners cannot be fitted, tightened, or masked away -- they
# must NOT be collocated. EDS does exactly that: SIMULATE the model to get the
# ergodic cloud, keep a well-spaced subset (farthest-point sampling), and fit the
# rules by COMPLETE-POLYNOMIAL LEAST SQUARES on that cloud. The solver then only ever
# visits states the model actually reaches -- no infeasible corners, no global-fit
# poisoning, no Gibbs (low complete degree). EDSGrid is a drop-in for SmolyakGrid
# (same points/n/d/fit/eval/basis/clip interface), so RuleSet / point_map / the
# time-iteration driver are reused UNCHANGED.
from itertools import combinations_with_replacement

import numpy as np

from solver_recursive.state_grid import (chebyshev_basis_1d, build_state_box,
                                          s_process_params, default_prob,
                                          IK_D, IK_F, IP_D, IP_F, IBDD, IBDF, IBFD, IV,
                                          IS, IZ, STATE_NAMES)
from solver_recursive.decision_rules import RuleSet, STORE_RULES
from solver_recursive.recursive_main import (time_iteration, ss_state, ss_x,
                                             calibrate_household_anchors)
from solver_recursive.point_map import point_residuals, SOLVE7

# The state indices come from state_grid, not local literals. They were hard-wired to
# IS, IZ = 5, 6 here -- a third private copy of the convention, which the 7 -> 9 state
# change would have silently turned into b_DF and V_dep.


class EDSGrid:
    # ARBITRARY COLLOCATION POINTS + COMPLETE-DEGREE CHEBYSHEV BASIS, RIDGE-LS FIT.
    def __init__(self, points, degree=2, lo=None, hi=None, ridge=1e-3, active=None):
        self.points = np.atleast_2d(np.asarray(points, dtype=float))
        self.n, self.d = self.points.shape
        self.lo = self.points.min(0) if lo is None else np.asarray(lo, float)
        self.hi = self.points.max(0) if hi is None else np.asarray(hi, float)
        self.hi = np.where(self.hi > self.lo + 1e-9, self.hi, self.lo + 1e-6)
        self.degree = int(degree)
        # ACTIVE dimensions = those that actually vary on the cloud. Some endogenous
        # stocks barely move on the ergodic path (near-zero range); their basis columns
        # are degenerate and inject a spurious slope. Keep only monomials that use active
        # dims -> the rules are a polynomial in the varying states, flat in the rest.
        self.active = (np.asarray(active, dtype=bool) if active is not None
                       else (self.hi - self.lo) > 1e-5 * (np.abs(self.lo) + 1.0))
        self._exps = [e for e in self._monomials(self.d, self.degree)
                      if all(e[j] == 0 or self.active[j] for j in range(self.d))]
        self.n_basis = len(self._exps)
        self._maxe = max((max(e) for e in self._exps), default=0)
        Phi = self._basis(self.points)                        # (n, n_basis)
        self._Phi = Phi
        A = Phi.T @ Phi
        ds = float(np.mean(np.diag(A))) + 1e-12
        td = np.array([sum(e) for e in self._exps])           # total degree per basis fn
        # ridge on ALL non-constant terms (degree>=1), scaled to the MEAN diagonal ds:
        # a well-resolved dimension (large own diagonal, e.g. s) is barely touched, but
        # a near-constant/degenerate dimension (small own diagonal, e.g. the F-side or Z
        # states that hardly move in the D-risk sim) is strongly damped -> coefficients
        # cannot explode. Only the constant term (degree 0) is unpenalised (level stays
        # unbiased). Floor for invertibility.
        A[np.diag_indices_from(A)] += ds * (ridge * (td >= 1) + 1e-8)
        self._A_lu = np.linalg.cholesky(A)
        self._PhiT = Phi.T

    @staticmethod
    def _monomials(d, degree):
        # ALL EXPONENT VECTORS OF TOTAL DEGREE <= degree (the complete polynomial).
        exps = [(0,) * d]
        for total in range(1, degree + 1):
            for combo in combinations_with_replacement(range(d), total):
                e = [0] * d
                for j in combo:
                    e[j] += 1
                exps.append(tuple(e))
        return exps

    def _to_unit(self, x):
        return 2.0 * (np.atleast_2d(x) - self.lo) / (self.hi - self.lo) - 1.0

    def _basis(self, x):
        u = np.clip(self._to_unit(x), -1.5, 1.5)              # mild extrapolation guard
        T = [chebyshev_basis_1d(u[:, j], self._maxe) for j in range(self.d)]
        B = np.ones((u.shape[0], self.n_basis))
        for b, e in enumerate(self._exps):
            for j in range(self.d):
                if e[j]:
                    B[:, b] *= T[j][:, e[j]]
        return B

    def basis(self, x):
        return self._basis(x)

    def fit(self, values):
        # RIDGE LEAST-SQUARES COEFFICIENTS (over-determined: n points > n_basis).
        rhs = self._PhiT @ np.asarray(values, dtype=float)
        y = np.linalg.solve(self._A_lu, rhs)
        return np.linalg.solve(self._A_lu.T, y)

    def eval(self, coeffs, x):
        return self._basis(x) @ coeffs

    def clip(self, x):
        return np.clip(np.atleast_2d(x), self.lo, self.hi)


def _next_state(S, x, out, s_next):
    # NEXT-PERIOD STATE FROM THE PERIOD-MAP OUTPUTS (Kp from SOLVE7, the rest from out).
    Sn = np.empty(len(STATE_NAMES))
    Sn[IK_D], Sn[IK_F] = x[2], x[3]
    Sn[IP_D], Sn[IP_F] = out["Pp_D"], out["Pp_F"]
    Sn[IBDD], Sn[IBDF] = out["b_D_D_new"], out["b_D_F_new"]
    Sn[IBFD] = out["b_F_D_new"]
    Sn[IV] = out["Vp_dep"]
    Sn[IS], Sn[IZ] = s_next, S[IZ]
    return Sn


def simulate_cloud(rules, cal, ss, sproc, d, T=3000, seed=0, no_default=False,
                   burn=200):
    # FORWARD-SIMULATE THE MODEL ON THE FITTED RULES (read the policy, step the state).
    rng = np.random.default_rng(seed)
    S = ss_state(ss, cal, sproc).copy()
    out_c = np.empty((T, len(STATE_NAMES)))
    for t in range(T):
        Sm = np.atleast_2d(S)
        x = np.array([float(rules.eval(k, d, Sm)[0]) for k in SOLVE7])
        try:
            _, o = point_residuals(S, d, x, rules, cal, ss, sproc, n_gh=5,
                                   no_default=no_default)
        except (ValueError, RuntimeError, FloatingPointError):
            S = ss_state(ss, cal, sproc).copy()               # reset on a bad step
            out_c[t] = S
            continue
        out_c[t] = S
        s_next = ((1.0 - sproc["rho_s"]) * sproc["s_star"] + sproc["rho_s"] * S[IS]
                  + sproc["sigma_s"] * rng.standard_normal())
        s_next = float(np.clip(s_next, -11.0, -2.5))
        S = _next_state(S, x, o, s_next)
        # CLIP to the seed box: the Chebyshev rules explode if extrapolated, so the sim
        # must stay where the seed is reliable. The cloud is then the in-box ergodic
        # path -- the feasible interior, never the infeasible corners.
        S = rules.grid.clip(S).ravel()
        if not np.all(np.isfinite(S)):
            S = ss_state(ss, cal, sproc).copy()
    return out_c[burn:]


def eds_select(cloud, n_target, seed=0):
    # FARTHEST-POINT SAMPLING: a well-spaced (epsilon-distinguishable) subset that
    # covers the cloud uniformly, on unit-normalised coordinates.
    cloud = np.asarray(cloud, dtype=float)
    lo, hi = cloud.min(0), cloud.max(0)
    U = (cloud - lo) / (hi - lo + 1e-12)
    rng = np.random.default_rng(seed)
    idx = [int(rng.integers(len(U)))]
    d2 = np.sum((U - U[idx[0]]) ** 2, axis=1)
    for _ in range(min(n_target, len(U)) - 1):
        j = int(np.argmax(d2))
        idx.append(j)
        d2 = np.minimum(d2, np.sum((U - U[j]) ** 2, axis=1))
    return cloud[idx]


def build_cloud(seed_rules, cal, ss, sproc, T=3000):
    # The d=0 ERGODIC cloud (clipped to the seed box) serves BOTH regimes: the d=1
    # haircut enters the period map's surv factor and its low-B continuation clips to
    # the grid bound -- the same approximation the Smolyak solve already makes. This
    # avoids simulating d=1 at post-haircut B, which is far outside the seed box (the
    # seed rules explode there). Z_D jittered (degenerate in the risk sim).
    cloud = simulate_cloud(seed_rules, cal, ss, sproc, d=0, T=T, no_default=False)
    cloud = cloud[np.all(np.isfinite(cloud), axis=1)]
    cloud[:, IZ] = cal["Z_ss_D"] * (1.0 + np.random.default_rng(7).uniform(
        -0.03, 0.03, len(cloud)))
    return cloud


def build_eds_grid(seed_rules, cal, ss, sproc, n_points=140, degree=2,
                   T=3000, verbose=False, cloud=None, ridge=1e-3):
    # BUILD (or reuse) THE CLOUD -> EDS subset (SS anchored) -> EDSGrid over both regimes.
    if cloud is None:
        cloud = build_cloud(seed_rules, cal, ss, sproc, T=T)
    pts = eds_select(cloud, n_points)
    pts = np.vstack([ss_state(ss, cal, sproc), pts])          # anchor the SS
    lo = np.minimum(pts.min(0), cloud.min(0))
    hi = np.maximum(pts.max(0), cloud.max(0))
    # ACTIVE dims from the CLOUD variance (not the grid extent -- the SS anchor would
    # otherwise make a degenerate dim look active). Rules are a polynomial in these.
    active = cloud.std(0) > 1e-4 * (np.abs(cloud.mean(0)) + 1.0)
    g = EDSGrid(pts, degree=degree, lo=lo, hi=hi, ridge=ridge, active=active)
    if verbose:
        print(f"  cloud={len(cloud)} pts -> EDS grid={g.n} pts, degree {degree} "
              f"({g.n_basis} basis fns), fit cond={np.linalg.cond(g._Phi.T @ g._Phi):.1e}",
              flush=True)
    return g


def solve_eds(cal, ss, sproc, n_points=140, degree=2, seed_mu=1, verbose=True,
              max_it=60, damp=0.2):
    # FULL EDS PIPELINE: mu=1 Smolyak seed (for the CLOUD only) -> EDS grid -> solve
    # on it from the FEASIBLE from_ss cold start (NOT the seed extrapolation, which
    # blows up outside the seed's box). Mirrors solve_recursive: d0 no-default, then
    # the d1 recovery homotopy, then the joint two-regime solve.
    from solver_recursive.recursive_experiment import solve_recursive
    if verbose:
        print("  [eds] solving mu=1 Smolyak seed (for the ergodic cloud) ...", flush=True)
    seed_rules = solve_recursive(cal, ss, sproc, mu=seed_mu, verbose=False, s_refine=0)
    grid = build_eds_grid(seed_rules, cal, ss, sproc, n_points=n_points,
                          degree=degree, verbose=verbose)
    rules = RuleSet.from_ss(grid, ss, cal)                    # feasible everywhere
    if verbose:
        print("  [eds] d0 no-default solve ...", flush=True)
    _, _, w0, _ = time_iteration(rules, cal, ss, sproc, regimes=(0,), no_default=True,
                              damp=0.25, tol=1e-6, max_it=max_it, n_gh=5, verbose=verbose)
    for k in STORE_RULES:                                     # warm d1 from converged d0
        rules.set_values(k, 1, rules.vals[k][0].copy())
    rec_t = cal["recovery_rate_D"]
    w1 = np.nan
    for rec in (0.85, 0.70, 0.55, rec_t):                    # d1 recovery homotopy
        cal["recovery_rate_D"] = rec
        _, _, w1, _ = time_iteration(rules, cal, ss, sproc, regimes=(1,), no_default=False,
                                  damp=0.2, tol=1e-6, max_it=40, n_gh=5, verbose=False)
        if verbose:
            print(f"  [eds] d1 homotopy recovery={rec:.2f}: worst={w1:.2e}", flush=True)
    cal["recovery_rate_D"] = rec_t
    if verbose:
        print("  [eds] joint two-regime solve ...", flush=True)
    ok, it, wj, _ = time_iteration(rules, cal, ss, sproc, regimes=(0, 1), no_default=False,
                                damp=damp, tol=1e-6, max_it=max_it, n_gh=5, verbose=verbose)
    if verbose:
        print(f"  [eds] DONE: d0 worst={w0:.2e}  d1 worst={w1:.2e}  joint worst={wj:.2e}",
              flush=True)
    return rules
