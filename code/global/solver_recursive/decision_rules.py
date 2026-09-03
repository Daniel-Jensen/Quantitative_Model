# DECISION-RULE LAYER: PER-REGIME CHEBYSHEV COEFFICIENTS ON THE STATE GRID.
# Every equilibrium object is approximated as a rule x(j, S) with SEPARATE
# coefficient sets for the compound regime j (see regime_table). Two kinds:
#   SOLVE   -- the pointwise Newton unknowns: N, Kp, rdep and household saving A per
#              country, the terms of trade p, BOTH sovereign prices Q_bD/Q_bF and the
#              cross-border holdings b_DF/b_FD. (SOLVE7 is a back-compat alias; the
#              name is historical, the tuple is 13 long.)
#   DERIVED -- objects READ OFF the Euler recursions given the frozen continuation
#              (point_map computes them): the banker valuations alpha, the
#              consumptions, r_wc.
# Both are stored and interpolated (continuation + simulation need them). Storage
# layout is the one place the stacking convention lives -- every caller takes the
# order from SOLVE rather than repeating it, because a duplicated literal is exactly
# how the 7 -> 11 change silently broke recursive_main._sweep's x_ss.
import numpy as np

# BOTH SOVEREIGN MARKETS CLEAR. Each bond used to be FORCE-FED to the banks:
# b_D_D = (1-shareF)*B' at a FIXED SS share, with the price then read off the
# D bank's Euler given that imposed quantity. Neither intermediary had a demand
# schedule and no market cleared. Now both banks' D-bond FOCs are residuals, b_DF is
# an unknown, b_DD = B' - b_DF clears the market by construction, and Q_bD is the
# price that does it. Q_bD stays a STORED rule (the continuation needs Q_bD'), it is
# just no longer read off a recursion.
# A_D IS SOLVED TOO -- THE UNION DEPOSIT MARKET. Under the old NATIONAL clearing each
# household was force-fed its own bank's funding need (A_D = dep_D/P_CES_D) and the F
# household's Euler was computed and DROPPED, so consumption was the bookkeeping
# residual of the bank balance sheet: C_D = W/P_CES + inc - A_D with the two gross legs
# ~8 and C_D ~0.79, income contributing 2% of the movement against 42% from each gross
# leg. A_D is now the D household's CHOICE (euler_D), euler_F is restored as a residual,
# and union clearing is an EXPLICIT residual with A_F a genuine unknown too. A_F was
# briefly left as the residual OF the clearing identity -- algebraically the same system,
# but the Newton then had A_F absorbing 100% of any union funding swing (a 2% capital
# move is ~0.18 of funding, ~23% of a household's consumption), C_F slammed into its
# clip plateau and hybr made ZERO progress at 6/19 points. Both savings solved, clearing
# scaled by SS deposits, is the conditioned form of the identical equilibrium.
# THE CB BACKSTOP ADDS NO UNKNOWN. The LTRO facility is a fixed envelope, fully drawn
# whenever it is offered (weakly optimal: it is lent at the deposit rate and relaxes the
# constraint), so there is no quantity to solve for and no complementarity.
SOLVE = ("N_D", "N_F", "Kp_D", "Kp_F", "rdep_D", "rdep_F", "p",
         "Q_bD", "b_DF", "Q_bF", "b_FD", "A_D", "A_F")
SOLVE7 = SOLVE                              # back-compat alias (older imports)
# banker valuations + household aggregates, all READ OFF the recursions/closure
DERIVED = ("alpha_D", "alpha_F", "C_D", "C_F",
           "r_wc_D", "r_wc_F")
STORE_RULES = SOLVE7 + DERIVED            # interpolated for continuation/sim
ALL_RULES = STORE_RULES
# back-compat alias (older imports)
DERIVED4 = DERIVED

# RULES INTERPOLATED IN LOGS (Bocola parameterises every policy as exp(ss + gamma)).
# Fitting log x rather than x makes the interpolant positive by construction, so the
# hard clips that used to protect positivity are unnecessary and the fit stops having
# to represent a plateau. Values are STORED in levels throughout (warm starts, damping
# and every caller are unchanged); only the Chebyshev fit and evaluation go through
# the transform. rdep is a rate that may legitimately go negative, so it is carried as
# the GROSS rate 1 + r -- Bocola's R -- which is positive.
LOG_RULES = frozenset({"N_D", "N_F", "Kp_D", "Kp_F", "p",
                       "alpha_D", "alpha_F", "Q_bD", "Q_bF", "b_DF", "b_FD",
                       "C_D", "C_F", "A_D", "A_F"})
GROSS_RULES = frozenset({"rdep_D", "rdep_F", "r_wc_D", "r_wc_F"})
_FIT_FLOOR = 1e-12


def to_fit(name, v):
    # LEVELS -> THE QUANTITY ACTUALLY FITTED BY THE CHEBYSHEV COLLOCATION.
    if name in GROSS_RULES:
        return np.log(np.maximum(1.0 + np.asarray(v, dtype=float), _FIT_FLOOR))
    if name in LOG_RULES:
        return np.log(np.maximum(np.asarray(v, dtype=float), _FIT_FLOOR))
    return np.asarray(v, dtype=float)


# THE REGIME TABLE: the compound index j -> (default d', CB-active m'). The index used
# to BE the default indicator; it now carries the CB regime too, because the TPI
# backstop is a second discrete state the continuation has to be conditioned on.
#   n = 2  the pre-TPI model: j IS d, no central bank anywhere.
#   n = 3  adds a CB-active regime in the NO-DEFAULT states only. Right for an
#          instrument that is conditional on the sovereign -- a yield peg, or an LTRO
#          under the collateral-ineligibility rule the ECB applied to Greek paper in
#          2012 and 2015 -- and 25% cheaper than n = 4.
#   n = 4  makes the CB regime ORTHOGONAL to default: the facility is available in the
#          default state too. Right for an instrument that supports BANKS rather than
#          the sovereign, and NECESSARY for the risk-premium channel: the default branch
#          carries little probability mass but the largest payoff deviation, so it
#          dominates cov(Omega, payD), which is the term a credible backstop compresses.
# n = 2 is bit-for-bit the old layout, which is what makes phi = 0 nest exactly.
_REG_TABLE = {2: ((0, 0), (1, 0)),
              3: ((0, 0), (0, 1), (1, 0)),
              4: ((0, 0), (0, 1), (1, 0), (1, 1))}


def regime_table(n_regimes):
    # (d, m) PAIRS FOR EVERY REGIME INDEX, single-sourced.
    return _REG_TABLE[int(n_regimes)]


def from_fit(name, y):
    # FITTED QUANTITY -> LEVELS (the inverse of to_fit).
    # The exponent is clamped only to keep a diverging transient iterate FINITE
    # (exp overflows to inf at ~709, and inf propagates into every expectation);
    # +-50 spans 1e-22 .. 5e21, so it is unreachable by any admissible value.
    if name in GROSS_RULES:
        return np.exp(np.clip(y, -50.0, 50.0)) - 1.0
    if name in LOG_RULES:
        return np.exp(np.clip(y, -50.0, 50.0))
    return y


class RuleSet:
    # COEFFICIENTS AND POINT VALUES FOR EVERY RULE IN EVERY REGIME.

    def __init__(self, grid, n_regimes=2):
        # EMPTY CONTAINER BOUND TO ONE GRID (values (n,) PER RULE PER REGIME).
        self.grid = grid
        self.n_regimes = int(n_regimes)
        self.reg = regime_table(self.n_regimes)
        self.vals = {k: [np.empty(grid.n) for _ in self.reg] for k in ALL_RULES}
        self.coef = {k: [None for _ in self.reg] for k in ALL_RULES}
        # quadrature order the rules were SOLVED under; time_iteration stamps it and
        # every reader (IRFs, decompositions, accuracy) uses it, so a solve at n_gh=5
        # is never read back at n_gh=7 and charged the difference as approximation error
        self.n_gh = None

    def set_values(self, name, d, values, weights=None, ridge=0.0):
        # SET POINT VALUES FOR ONE RULE IN ONE REGIME AND REFIT ITS COEFFICIENTS.
        # Values are stored in LEVELS; the fit is on to_fit(name, .) so the log rules
        # are collocated in logs. weights/ridge (optional) switch to the masked
        # ridge-LS fit; the default (both absent) keeps the exact square solve.
        self.vals[name][d] = np.asarray(values, dtype=float).copy()
        y = to_fit(name, self.vals[name][d])
        if weights is None and ridge == 0.0:
            self.coef[name][d] = self.grid.fit(y)
        else:
            self.coef[name][d] = self.grid.fit_weighted(y, weights, ridge)

    def eval(self, name, d, x):
        # EVALUATE ONE RULE IN REGIME d AT NATURAL-COORDINATE POINTS x.
        return from_fit(name, self.grid.eval(self.coef[name][d], x))

    def eval_all(self, d, x):
        # EVALUATE EVERY RULE IN REGIME d AT POINTS x -> DICT OF ARRAYS.
        B = self.grid.basis(x)
        return {k: from_fit(k, B @ self.coef[k][d]) for k in ALL_RULES}

    def copy(self):
        # DEEP-ENOUGH COPY FOR A FROZEN CONTINUATION (values + coefficients).
        rs = RuleSet(self.grid, self.n_regimes)
        rs.n_gh = self.n_gh
        for k in ALL_RULES:
            for d in range(self.n_regimes):
                rs.vals[k][d] = self.vals[k][d].copy()
                rs.coef[k][d] = (None if self.coef[k][d] is None
                                 else self.coef[k][d].copy())
        return rs

    @classmethod
    def from_ss(cls, grid, ss, cal, n_regimes=2):
        # STEADY-STATE COLD START IN EVERY REGIME. Constant SS everywhere EXCEPT
        # Kp_D/Kp_F, which track the state's own K so the Jermann inversion is
        # feasible at every grid point (SS-constant Kp asks for an infeasible
        # I/K at high-K corners); N is normalised to 1 at the SS.
        rs = cls(grid, n_regimes)
        bk_D, bk_F = ss["ss_bank_D"], ss["ss_bank_F"]
        const = dict(N_D=1.0, N_F=1.0,
                     rdep_D=cal["r_dep_D_target"], rdep_F=cal["r_dep_F_target"],
                     p=ss["p_ss"],
                     alpha_D=bk_D["alpha_ss"], alpha_F=bk_F["alpha_ss"],
                     Q_bD=ss["Q_bD_ss"], Q_bF=ss["Q_bF_ss"],
                     b_DF=cal["b_D_F_ss"], b_FD=cal["b_F_D_ss"],
                     C_D=ss["C_D_ss"], C_F=ss["C_F_ss"],
                     A_D=ss["A_D_ss"], A_F=ss["A_F_ss"],
                     # r_wc = rdep + lambda*mu/Omega, constant at the SS
                     r_wc_D=cal["r_dep_D_target"] + cal["credit_spread_target_D"],
                     r_wc_F=cal["r_dep_F_target"] + cal["credit_spread_target_F"])
        for k, v in const.items():
            for d in range(rs.n_regimes):
                rs.set_values(k, d, np.full(grid.n, v))
        for d in range(rs.n_regimes):           # Kp tracks the K state (feasible)
            rs.set_values("Kp_D", d, grid.points[:, 0].copy())
            rs.set_values("Kp_F", d, grid.points[:, 1].copy())
        return rs
