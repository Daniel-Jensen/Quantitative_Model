# DECISION RULES: PER-REGIME CHEBYSHEV FITS ON THE STATE GRID.
import numpy as np

# the Newton unknowns at each grid point; the TPI variable x_cb is last
SOLVE = ("N_D", "N_F", "Kp_D", "Kp_F", "rdep_D", "rdep_F", "p",
         "Q_bD", "b_DF", "Q_bF", "b_FD", "A_D", "A_F", "x_cb")
SOLVE7 = SOLVE  # back-compat alias
# objects read off the recursions, stored for the continuation
DERIVED = ("alpha_D", "alpha_F", "C_D", "C_F",
           "r_wc_D", "r_wc_F")
STORE_RULES = SOLVE7 + DERIVED
ALL_RULES = STORE_RULES
# back-compat alias
DERIVED4 = DERIVED

# rules fitted in logs, or as gross rates, so the fit stays positive; the cross-border books
# b_DF, b_FD are fitted in levels (2026-10-20): each has a no-short-sale corner, where a log
# unknown's effect vanishes (the D bank's default-regime book: Jacobian column ~1e-10)
LOG_RULES = frozenset({"N_D", "N_F", "Kp_D", "Kp_F", "p",
                       "alpha_D", "alpha_F", "Q_bD", "Q_bF",
                       "C_D", "C_F", "A_D", "A_F"})
GROSS_RULES = frozenset({"rdep_D", "rdep_F", "r_wc_D", "r_wc_F"})
_FIT_FLOOR = 1e-12


def to_fit(name, v):
    # LEVELS -> THE QUANTITY THE CHEBYSHEV FIT USES.
    if name in GROSS_RULES:
        return np.log(np.maximum(1.0 + np.asarray(v, dtype=float), _FIT_FLOOR))
    if name in LOG_RULES:
        return np.log(np.maximum(np.asarray(v, dtype=float), _FIT_FLOOR))
    return np.asarray(v, dtype=float)


# the regime index is the default indicator: 0 = no default, 1 = default
REGIMES = (0, 1)


def regime_table(n_regimes):
    # THE DEFAULT INDICATOR OF EACH REGIME.
    assert int(n_regimes) == len(REGIMES), "the model has exactly two regimes: d = 0, 1"
    return REGIMES


def from_fit(name, y):
    # FITTED QUANTITY -> LEVELS; THE CLIP ONLY KEEPS A DIVERGING ITERATE FINITE.
    if name in GROSS_RULES:
        return np.exp(np.clip(y, -50.0, 50.0)) - 1.0
    if name in LOG_RULES:
        return np.exp(np.clip(y, -50.0, 50.0))
    return y


class RuleSet:
    # COEFFICIENTS AND POINT VALUES FOR EVERY RULE IN EVERY REGIME.
    # class defaults, so rule sets pickled before the post-default box (or the travelling box
    # spec, recursive_experiment._build_rules) still load
    post = home = blend = box = None

    def __init__(self, grid, n_regimes=2):
        # AN EMPTY RULE SET ON ONE GRID.
        self.grid = grid
        self.n_regimes = int(n_regimes)
        self.reg = regime_table(self.n_regimes)
        self.vals = {k: [np.empty(grid.n) for _ in self.reg] for k in ALL_RULES}
        self.coef = {k: [None for _ in self.reg] for k in ALL_RULES}
        # the quadrature order the rules were solved with, reused by every reader
        self.n_gh = None
        # the post-default rule set on its own box, its link back here, and the blend between
        # the two (attach_post); all None for a single box
        self.post, self.home, self.blend = None, None, None

    def set_values(self, name, d, values, weights=None, ridge=0.0):
        # SET ONE RULE'S VALUES IN ONE REGIME AND REFIT IT.
        self.vals[name][d] = np.asarray(values, dtype=float).copy()
        y = to_fit(name, self.vals[name][d])
        if weights is None and ridge == 0.0:
            self.coef[name][d] = self.grid.fit(y)
        else:
            self.coef[name][d] = self.grid.fit_weighted(y, weights, ridge)

    def eval(self, name, d, x):
        # EVALUATE ONE RULE IN REGIME d AT POINTS x.
        return from_fit(name, self.grid.eval(self.coef[name][d], x))

    def eval_all(self, d, x):
        # EVALUATE EVERY RULE IN REGIME d AT POINTS x.
        B = self.grid.basis(x)
        return {k: from_fit(k, B @ self.coef[k][d]) for k in ALL_RULES}

    def attach_post(self, post, blend):
        # ATTACH THE POST-DEFAULT RULE SET AND THE BLEND THAT HANDS THE CONTINUATION BETWEEN THEM.
        self.post, self.blend = post, blend
        post.home, post.blend = self, blend

    def blend_weight(self, x):
        # THE POST-DEFAULT RULES' WEIGHT AT POINTS x: 1 TOWARDS THE POST-DEFAULT CLUSTER, 0 NEAR
        # THE SS, A C2 QUINTIC STEP BETWEEN, ALONG THE RECOVERY LINE IN NET-WORTH UNITS.
        b = self.blend
        p = ((np.atleast_2d(x)[:, b["block"]] - b["centre"]) * b["weights"]) @ b["direction"]
        t = np.clip((p - b["post_until"]) / (b["normal_from"] - b["post_until"]), 0.0, 1.0)
        return 1.0 - t * t * t * (10.0 - 15.0 * t + 6.0 * t * t)

    def eval_cont(self, d, x):
        # THE CONTINUATION RULES IN REGIME d AT NEXT STATES x, EACH SET READ THROUGH ITS OWN BOX
        # EDGE; ONE FUNCTION WHICHEVER SET'S NODE ASKS.
        base = self if self.home is None else self.home
        r = base.eval_all(d, base.grid.clip(x))
        if base.post is None:
            return r
        w = base.blend_weight(x)
        if not np.any(w > 0.0):
            return r
        rp = base.post.eval_all(d, base.post.grid.clip(x))
        return {k: (1.0 - w) * r[k] + w * rp[k] for k in r}

    def eval_read(self, name, d, x):
        # ONE RULE AT STATES x FOR THE READERS (IRFs, REST POINT), BLENDED AS THE CONTINUATION IS.
        base = self if self.home is None else self.home
        v = base.eval(name, d, x)
        if base.post is None:
            return v
        w = base.blend_weight(x)
        return v if not np.any(w > 0.0) else (1.0 - w) * v + w * base.post.eval(name, d, x)

    def copy(self):
        # A COPY FOR A FROZEN CONTINUATION (THE POST-DEFAULT SET IS COPIED AND RE-LINKED).
        rs = RuleSet(self.grid, self.n_regimes)
        rs.n_gh, rs.box = self.n_gh, self.box
        for k in ALL_RULES:
            for d in range(self.n_regimes):
                rs.vals[k][d] = self.vals[k][d].copy()
                rs.coef[k][d] = (None if self.coef[k][d] is None
                                 else self.coef[k][d].copy())
        if self.post is not None:
            rs.attach_post(self.post.copy(), self.blend)
        return rs

    @classmethod
    def from_ss(cls, grid, ss, cal, n_regimes=2):
        # STEADY-STATE COLD START IN EVERY REGIME.
        rs = cls(grid, n_regimes)
        bk_D, bk_F = ss["ss_bank_D"], ss["ss_bank_F"]
        const = dict(N_D=1.0, N_F=1.0,
                     rdep_D=cal["r_dep_D_target"], rdep_F=cal["r_dep_F_target"],
                     p=ss["p_ss"],
                     alpha_D=bk_D["alpha_ss"], alpha_F=bk_F["alpha_ss"],
                     Q_bD=ss["Q_bD_ss"], Q_bF=ss["Q_bF_ss"],
                     b_DF=cal["b_D_F_ss"], b_FD=cal["b_F_D_ss"],
                     C_D=ss["C_D_ss"], C_F=ss["C_F_ss"],
                     A_D=ss["A_D_ss"], A_F=ss["A_F_ss"], x_cb=0.0,
                     # r_wc = rdep + the credit spread at the SS
                     r_wc_D=cal["r_dep_D_target"] + cal["credit_spread_target_D"],
                     r_wc_F=cal["r_dep_F_target"] + cal["credit_spread_target_F"])
        for k, v in const.items():
            for d in range(rs.n_regimes):
                rs.set_values(k, d, np.full(grid.n, v))
        for d in range(rs.n_regimes):  # Kp tracks the K state so investment is feasible everywhere
            rs.set_values("Kp_D", d, grid.points[:, 0].copy())
            rs.set_values("Kp_F", d, grid.points[:, 1].copy())
        return rs
