# DECISION-RULE LAYER: PER-REGIME CHEBYSHEV COEFFICIENTS ON THE STATE GRID.
# Every equilibrium object is approximated as a rule x(d, S) with SEPARATE
# coefficient sets for the default regime d_D in {0,1}. Two kinds:
#   SOLVE7   -- the pointwise Newton unknowns (transition.py's stacked system):
#               N, Kp, rdep per country + the terms of trade p.
#   DERIVED4 -- the banker valuations alpha, Q_b, READ OFF the Euler recursions
#               given the frozen continuation (point_map computes them).
# Both are stored and interpolated (continuation + simulation need them). The
# two consumption rules (HH_RULES) are updated by the household layer, not the
# price solve. Storage layout is the one place the stacking convention lives.
import numpy as np

SOLVE7 = ("N_D", "N_F", "Kp_D", "Kp_F", "rdep_D", "rdep_F", "p")
# banker valuations + household aggregates, all READ OFF the recursions/closure
DERIVED = ("alpha_D", "alpha_F", "Q_bD", "Q_bF", "C_D", "C_F", "A_D", "A_F")
STORE_RULES = SOLVE7 + DERIVED            # interpolated for continuation/sim
ALL_RULES = STORE_RULES
# back-compat alias (older imports)
DERIVED4 = DERIVED


class RuleSet:
    # COEFFICIENTS AND POINT VALUES FOR EVERY RULE IN BOTH REGIMES.

    def __init__(self, grid):
        # EMPTY CONTAINER BOUND TO ONE GRID (values (n,) PER RULE PER REGIME).
        self.grid = grid
        self.vals = {k: [np.empty(grid.n), np.empty(grid.n)] for k in ALL_RULES}
        self.coef = {k: [None, None] for k in ALL_RULES}

    def set_values(self, name, d, values, weights=None, ridge=0.0):
        # SET POINT VALUES FOR ONE RULE IN ONE REGIME AND REFIT ITS COEFFICIENTS.
        # weights/ridge (optional) switch to the masked ridge-LS fit; the default
        # (both absent) keeps the exact square collocation solve.
        self.vals[name][d] = np.asarray(values, dtype=float).copy()
        if weights is None and ridge == 0.0:
            self.coef[name][d] = self.grid.fit(self.vals[name][d])
        else:
            self.coef[name][d] = self.grid.fit_weighted(self.vals[name][d],
                                                        weights, ridge)

    def eval(self, name, d, x):
        # EVALUATE ONE RULE IN REGIME d AT NATURAL-COORDINATE POINTS x.
        return self.grid.eval(self.coef[name][d], x)

    def eval_all(self, d, x):
        # EVALUATE EVERY RULE IN REGIME d AT POINTS x -> DICT OF ARRAYS.
        B = self.grid.basis(x)
        return {k: B @ self.coef[k][d] for k in ALL_RULES}

    def copy(self):
        # DEEP-ENOUGH COPY FOR A FROZEN CONTINUATION (values + coefficients).
        rs = RuleSet(self.grid)
        for k in ALL_RULES:
            for d in (0, 1):
                rs.vals[k][d] = self.vals[k][d].copy()
                rs.coef[k][d] = (None if self.coef[k][d] is None
                                 else self.coef[k][d].copy())
        return rs

    @classmethod
    def from_ss(cls, grid, ss, cal):
        # STEADY-STATE COLD START IN BOTH REGIMES. Constant SS everywhere EXCEPT
        # Kp_D/Kp_F, which track the state's own K so the Jermann inversion is
        # feasible at every grid point (SS-constant Kp asks for an infeasible
        # I/K at high-K corners); N is normalised to 1 at the SS.
        rs = cls(grid)
        bk_D, bk_F = ss["ss_bank_D"], ss["ss_bank_F"]
        const = dict(N_D=1.0, N_F=1.0,
                     rdep_D=cal["r_dep_D_target"], rdep_F=cal["r_dep_F_target"],
                     p=ss["p_ss"],
                     alpha_D=bk_D["alpha_ss"], alpha_F=bk_F["alpha_ss"],
                     Q_bD=ss["Q_bD_ss"], Q_bF=ss["Q_bF_ss"],
                     C_D=ss["C_D_ss"], C_F=ss["C_F_ss"],
                     A_D=ss["A_D_ss"], A_F=ss["A_F_ss"])
        for k, v in const.items():
            for d in (0, 1):
                rs.set_values(k, d, np.full(grid.n, v))
        for d in (0, 1):                        # Kp tracks the K state (feasible)
            rs.set_values("Kp_D", d, grid.points[:, 0].copy())
            rs.set_values("Kp_F", d, grid.points[:, 1].copy())
        return rs
