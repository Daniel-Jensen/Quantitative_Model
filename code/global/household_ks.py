# KRUSELL-SMITH HOUSEHOLD LAYER: HA POLICIES OVER THE AGGREGATE STATE GRID.
# The heterogeneous-agent block (household.py / distribution.py) is kept as the
# single source of household optimality; this layer makes it RECURSIVE in the
# aggregate state. Approximation (agreed Tier-3): each country's distribution
# is carried through its MEAN W_X only (1-moment KS); the cross-sectional
# SHAPE is anchored at the steady-state distribution, shifted proportionally
# to match W. Households forecast with the equilibrium rules themselves
# (internally consistent KS -- no separate perceived law of motion).
#
# ITERATION SCHEME (the outer KS loop around the price-rule Newton):
#   given a RuleSet (prices/valuations + current C rules):
#   1. For each grid point (d, S) and each country X:
#        income objects (w, Div, Tax) from the period map at that point;
#        deposit return locked one period: R_X = 1 + rdep_X(d, S).
#   2. EGM backward step on (a, e) x grid: continuation marginal utility is
#        E[u'(c(a', e'; S'))] with S' from the same two-branch quadrature as
#        expectations.py (eps_s nodes x d'), evaluating the CURRENT c-rules.
#        Iterate to policy convergence (the household problem is stationary
#        given the rules). GHH utility exactly as household.py.
#   3. Aggregate with the shifted-SS-shape distribution at mean W_X:
#        C_X(d, S), A'_X(d, S)  ->  W transition W_X' = R_X(d,S) * A'_X(d,S).
#   4. rules.set_values("C_X", d, C_X) for both regimes; return the W' maps
#        for the period map's state transition.
#   Outer loop: KS step (1-4) alternates with the price-rule Newton until the
#   C rules and price rules are jointly fixed. The AGGREGATE-EULER closure is
#   the degenerate first iterate of this loop (skip 2-3, close C with the
#   rep-agent Euler) -- the validation ladder inside the same architecture.
#
# REUSE: the EGM backward step and distribution forward step call the SAME
# kernels as the transition solver (fast_kernels.py numba + numpy fallback);
# this module only re-indexes them from time-paths to grid-points.
import numpy as np

from decision_rules import RuleSet   # noqa: F401  (API typing/reference)


def shifted_ss_distribution(D_ss, a_grid, W_target, W_ss):
    # 1-MOMENT KS DISTRIBUTION: THE SS SHAPE SHIFTED TO MATCH MEAN WEALTH.
    # Proportional shift on the asset dimension (mass stays on the grid via
    # the same lottery reassignment distribution.py uses).
    raise NotImplementedError("assembled with the KS iteration (step 3)")


def egm_on_grid(rules, cal, ss, quad, country="D"):
    # POLICY FUNCTIONS c(a, e; d, S) BY EGM WITH QUADRATURE CONTINUATIONS.
    raise NotImplementedError("assembled with the KS iteration (step 2)")


def ks_update(rules, cal, ss, quad, mode="ks"):
    # ONE OUTER KS STEP: HOUSEHOLD POLICIES -> AGGREGATES -> UPDATED C RULES.
    # mode="aggregate_euler" is the degenerate first iterate (no EGM): close
    # C_X on the grid with the households' aggregate deposit Euler -- used to
    # initialize and as the Tier-1 validation ladder.
    raise NotImplementedError("assembled after recursive_residual.py lands")
