# ALL MODEL PARAMETERS FOR THE TWO-COUNTRY HANK-GK MONETARY UNION.
# D-bonds are D-good claims (priced with rdep_D), F-bonds F-good claims
# (rdep_F); cross-border legs convert via p (D-goods per F-good).
# Bocola (2016) Tables 1-2 anchor every parameter with a direct counterpart.
# The deliberate divergences are flagged "vs Bocola" at the parameter itself —
# each one is load-bearing, so read the note before importing his value.


def get_calibration():
    # BUILD THE PARAMETER DICT CONSUMED BY EVERY BLOCK.
    cal = dict(
        # Household preferences. Bocola §II.A.1 uses log utility (NOT Epstein-Zin)
        # => sigma = 1. His nu = 0.5 (Table 1) is the INVERSE Frisch chosen for a
        # Frisch elasticity of 2; `frisch` here is the elasticity itself.
        # GHH (not his separable form) is his own §V.C open-economy fix.
        sigma_D=1.0,   sigma_F=1.0,
        frisch_D=2.0,  frisch_F=2.0,
        chi_D=0.5417,  chi_F=0.5417,   # warm start; SS solve overwrites to pin N_ss=1

        # Idiosyncratic income (Rouwenhorst)
        n_e_D=2,       n_e_F=2,
        rho_e_D=0.9,   sigma_e_D=0.2,
        rho_e_F=0.9,   sigma_e_F=0.2,

        # Asset grids
        a_min_D=0.0, a_max_D=87.2, n_a_D=250, a_curve_D=2.0,
        a_min_F=0.0, a_max_F=87.2, n_a_F=250, a_curve_F=2.0,

        # Firms (Cobb-Douglas, flexible prices)
        epsilon_D=6.0, epsilon_F=6.0,   # demand elasticity -> mc = (eps-1)/eps
        Z_ss_D=0.45,   Z_ss_F=0.45,     # warm start; SS solve overwrites to pin Y_ss=1

        # Capital (Jermann 1998 adjustment cost). ksi = elasticity of Tobin's q wrt I/K.
        #
        # 0.50, inside Bocola's Table 2 posterior [0.324, 0.525] (mean 0.426) and
        # replacing the 0.15 that the 2026-08-28 audit priced.
        #
        # WHY THE 0.15 WENT. It was chosen to kill a positive output hump from q2, but
        # the audit shows Bocola's OWN closed-economy IRF crosses zero at q12 and settles
        # at +0.07%, and his net worth crosses at q8 and settles at +2.9% -- the pattern
        # 0.15 was removing is in the reference model. Meanwhile ksi scales the capital
        # -price leg of bank-equity destruction one-for-one, because dlogQ_K = ksi *
        # dlog(I/K) in both codes: at his shock that leg is -2.87pp of his -6.74pp total
        # (43%), against -0.51pp of our -4.17pp at ksi = 0.15 (12%). Bank net worth is
        # what drives mu, the credit spread and the working-capital wedge, so the
        # departure was paying for a persistence pattern it did not need with a third of
        # the amplification.
        alpha_D=0.30,  alpha_F=0.30,    # capital share (Bocola Table 1)
        delta_D=0.025, delta_F=0.025,
        ksi_D=0.50,    ksi_F=0.50,

        # Financial intermediary. f = exit/payout share, so the Omega kernel puts
        # weight (1-f) on the franchise value.
        #
        # 0.08, NOT Bocola's psi = 0.9646 survival (f = 0.0354). This is the only bank
        # parameter still departing from him: leverage (5.0, Table 1), the sovereign
        # exposure (7.6% of assets, Table B1) and the recovery (0.45, Greek PSI) are all
        # his, and the spread target has its own note above. f is the least identified
        # number in the block -- a banker SURVIVAL rate is not observable -- which is why
        # it is the one that moves.
        #
        # WHY IT IS STILL 0.08 AT THE 100 bp TARGET. Re-measured on the coarse grid at
        # spread = 100 bp/yr, reading against the model's own rest point:
        #
        #     f     lambda  alpha_ss  mu_rest  Y fitted  Y exact  |gap|   d_r_wc
        #   0.0354  0.3065   1.5325   0.0112   -0.048%  -0.062%  0.014   +23 bp/yr
        #   0.0800  0.2363   1.1817   0.0098   -0.102%  -0.123%  0.021   +54 bp/yr  <- ships
        #   0.1200  0.2228   1.1142   0.0100   -0.137%  -0.149%  0.013   +68 bp/yr
        #
        # Two things this shows that the old (8 bp) measurements could not. First, at a
        # 100 bp target Bocola's OWN f = 0.0354 is the WORST of the three on fidelity as
        # well as on transmission: alpha_ss = f/(f - mu_ss), so with mu_ss = 0.0123 his f
        # implies alpha_ss = 1.53 against his own 1.026, while f = 0.08 gives 1.18 and
        # f = 0.12 gives 1.11. Raising f moves alpha and lambda TOWARD his values here.
        # Second, the identification gap is now small at every f -- the kink, not f, was
        # what made it large.
        #
        # 0.12 is better on every number in that table and essentially reaches his
        # open-economy -0.157%. It is NOT shipped because the cost lands on the one thing
        # f actually means: 88%/quarter survival is an average banker horizon of 8.3
        # quarters, ~2 years, against Bocola's 28 and Gertler-Kiyotaki's ~36. 0.08 is
        # 12.5 quarters, ~3 years -- already short, and the smaller departure. The
        # f = 0.12 row is recorded so the trade-off is visible rather than assumed.
        #
        # MECHANISM (unchanged, and it is why f is the right dial for TRANSMISSION even
        # though it is the wrong one for the kink): the risk shock raises leverage ~7%
        # and the franchise value alpha' by about as much, and Omega = beta[f + (1-f)a']
        # passes (1-f)alpha/[f+(1-f)alpha] of the latter through -- 0.961 at f = 0.04,
        # 0.921 at f = 0.08. That coefficient is the only thing standing between
        # d log(leverage) and d log E[Om] in dmu = d log(lev) - d log E[Om].
        # NB (2026-08-29): at mu_ss = 0.001 the model's STOCHASTIC rest point has
        # mu = EXACTLY 0 -- the economy sits ON the KKT kink, where a polynomial
        # interpolant of a C0 multiplier is least reliable (see CLAUDE.md). Raising f
        # helps the pass-through partly BECAUSE it lifts the rest point off the kink.
        # A calibration with a materially binding ergodic constraint would make the
        # level of the output response identified rather than bracketed.
        f_D=0.08,               f_F=0.08,
        # R^bg = 1.003 quarterly (Bocola Table 1 sample-average risk-free rate).
        r_dep_D_target=0.003,   r_dep_F_target=0.003,
        # beta*R = 1 at the SS under log utility => beta_inter = 1/R^bg. This MUST
        # move with r_dep_target: leaving it at 0.99 collapses the alpha fixed
        # point to a near-tangency and the SS solve fails.
        beta_inter_D=0.997,     beta_inter_F=0.997,
        # lambda and omega_ent are SOLVED (calibrate_bank_targets) to hit these.
        # BOCOLA'S OWN BANK CALIBRATION, decoded from his solved parameter vector
        # (Model/Matfiles/model_solution_mean.mat + model_param.m): lev = 5,
        # lambda = 0.20513 => alpha_ss = lambda*lev = 1.0256, psi = 0.96 => mu_ss =
        # 0.00100, excess return = (lambda/alp)(1/beta)(mu/(1-mu)) = 2.0bp/qtr = 8bp/yr,
        # omega_ent = 0.00745. This calibration reproduces all four to 3-4 digits.
        #
        # THE SPREAD TARGET IS 100 bp/yr, NOT BOCOLA'S 8 (2026-08-29). The leverage (5.0,
        # Table 1), the sovereign exposure (7.6% of assets, Table B1) and the recovery
        # (0.45, Greek PSI) are data-anchored and stay his. The spread does not.
        #
        # WHY. At 8 bp the model's STOCHASTIC REST POINT has mu = EXACTLY 0: the economy
        # sits ON the KKT kink. mu = max(.,0) is C0, so a Chebyshev interpolant returns
        # mu > 0 in a neighbourhood where the truth is 0, and reading the fitted rules
        # against clearing the period map exactly at the same state then disagree by MORE
        # than the response being measured. THE LEVEL OF THE OUTPUT RESPONSE IS NOT
        # IDENTIFIED THERE. Raising the target lifts the rest point off the kink; nothing
        # else does, and in particular f does NOT -- calibrate_bank_targets forces
        # alpha_ss = lambda*theta exactly, so the SS is marginally binding whatever f is,
        # and the binding margin in leverage units is ~theta*mu_ss with mu_ss set by the
        # SPREAD. Measured (coarse grid, rest point, p^d = 1.98% shock):
        #
        #   spread    mu_rest  spread_rest   Y fitted  Y exact   |gap|
        #      8 bp   0.00000       0.0 bp    -0.094%  -0.007%   0.087   <- on the kink
        #     25 bp   0.00000       0.0 bp    -0.097%  -0.061%   0.036
        #    100 bp   0.00983      79.0 bp    -0.102%  -0.123%   0.021   <- BINDS, shipped
        #    250 bp   0.02860     228.7 bp    -0.096%  -0.118%   0.022   <- no further gain
        #    400 bp   infeasible at f = 0.08 (franchise fold)
        #
        # 100 bp is where the constraint starts binding ergodically and the identification
        # gap collapses 4x; 250 buys nothing more and pushes alpha_ss to 1.61.
        #
        # IT IS ALSO A DATA-ANCHORED NUMBER, not a fitted one: it is Gertler-Kiyotaki's
        # OWN target ("an average credit spread of 100 basis points per year and an
        # economy-wide leverage ratio of 4"), and euro-area periphery bank lending spreads
        # ran 100-300 bp over the risk-free rate in 2011-12. The 200bp/leverage-4 pair
        # this file used to carry was labelled "(GK11)" but is 2x their own target.
        #
        # WHAT IT COSTS. Bocola ESTIMATES mu^bg = 0.00087 (Step-1 posterior mean, from
        # his constructed multiplier series), which is 8 bp/yr; 100 bp is 12x that. The
        # defence is that his estimate belongs to his CLOSED benchmark, where the
        # constraint does not have to transmit anything -- computed exactly from his own
        # solved coefficients, mu = 0 along his entire benchmark IRF except one quarter
        # and binds on 1.2% of his ergodic set, and his output response comes from the
        # labour-supply wealth effect instead. This model uses his SS V.C transmission
        # (GHH + working capital), where the wedge lambda*mu/E[Om] is the ONLY channel
        # into output, so a constraint that is slack in the ergodic region transmits
        # nothing. The two cannot both be imported.
        #
        # THE BENCHMARK IS NOT -1.05/-1.44/-1.53. That is a cumulated quarterly GROWTH
        # gap x400 over an 8-quarter estimated shock sequence; its output LEVEL
        # equivalent is -0.26/-0.36/-0.38%. The like-for-like single-shock IRF targets
        # at p^d = 1.98%/qtr are -0.157% (his SS V.C open economy) and -0.222% (his
        # closed benchmark) -- see reporting/prints.py.
        leverage_target_D=5.0,            leverage_target_F=5.0,
        credit_spread_target_D=0.0025,    credit_spread_target_F=0.0025,    # 100 bp/yr
        # warm starts; overwritten by calibrate_bank_targets in steady_state.py
        lambda_K_D=0.22,        lambda_K_F=0.22,
        lambda_bD_D=0.22,       lambda_bD_F=0.22,
        lambda_bF_D=0.22,       lambda_bF_F=0.22,
        omega_ent_D=0.002,      omega_ent_F=0.002,

        # SGU debt-elastic premium on the cross-border (net external) position, Bocola's
        # ONLY foreign friction (residual_model_open.m: R = 1/beta + 0.01*B_for/gdp). Enters
        # the deposit-UIP keyed to the P_D-P_F wealth imbalance; 0 at the symmetric SS, so
        # it is undistorting there and only induces stationarity off-SS (fixes the F-side
        # wealth quasi-unit-root). kappa_nfa=0 nests the frictionless UIP exactly. Verified
        # (scratchpad): kappa=0.01 => external-position AR(1) root ~0.92, half-life ~2y.
        # ON at Bocola's own 0.01 since 2026-08-25. At kappa = 0 the model has no force
        # returning wealth to the symmetric SS, so the sovereign-risk IRF never comes
        # back: 24 quarters after a shock that has 88% decayed, Y_D was still -0.373% and
        # FALLING, with K_D -1.26%, P_D -4.66%, B_D +8.5% and bank net worth +25.7% all
        # permanently displaced. That is the documented F-side quasi-unit-root, not a
        # persistence result, and reading a trough off a non-stationary path is
        # meaningless. This is the standard Schmitt-Grohe-Uribe (2003) stationarity
        # induction and it is Bocola's own foreign friction, so turning it on costs no
        # fidelity; it is undistorting at the symmetric SS (P_D = P_F => premium 0).
        kappa_nfa=0.01,
        # Cross-border bond portfolio adjustment costs (legacy transition solver only)
        # BOTH CROSS-BORDER PORTFOLIO ADJUSTMENT COSTS ARE LOAD-BEARING NOW. They used to
        # be dead (transition-solver only). With each sovereign split solved from the two
        # banks' FOCs, psi is the ONLY thing giving the foreign demand schedule a slope:
        # lambda_b*Q*b is ~7% of divertable assets, so both schedules are otherwise
        # near-flat and the split is numerically indeterminate.
        # D bond, measured: at 0.05 a 0.5% price wedge supported a 54% position swing and
        # the rest point walked b_DF from 0.196 to 0.360. At 2.0 the split holds to 4
        # decimals and every collocation point clears at 1e-14.
        # F bond, measured on the 21-point 10-state grid: 0.05 -> 10/21 points fail at
        # |F| = 8.7e-02; 0.5 -> 1/21 at 3.9e-03; 2.0 -> 0/21 at 6.2e-13. Same value for
        # both legs, since the two markets are structurally identical.
        psi_bF_D=2.0,           psi_bD_F=2.0,
        b_F_D_ss=0.196,         b_D_F_ss=0.196,   # ~20% of each supply (contagion leg)
        excess_return_F_D_ss=0.0,                 # overwritten after the SS solve
        excess_return_D_F_ss=0.0,                 # overwritten after the SS solve

        # Government bonds. delta_b = 0.056 = Bocola's pi (Table 1, fraction of the
        # HM/CE perpetuity maturing each quarter). Duration ~1/pi (long duration is
        # what makes priced risk generate large MTM losses; at delta_b = 0.25 the
        # repricing shrinks ~6x and the risk channel turns EXPANSIONARY).
        delta_b_D=0.056,        delta_b_F=0.056,
        # B_gov set so the D-bank's holdings of D-sovereign are 7.6% of its total
        # assets -- Bocola's exp^bg (Table 1; Table B1 gives 160/2093 = 0.076).
        # The earlier 3.722 ("93% of GDP") MISREAD Bocola's "93% of bank EQUITY"
        # holdings figure as a debt/GDP ratio, giving a 3x-too-high exposure and a
        # huge default fiscal windfall (55% haircut on ~93% of GDP) that made
        # default EXPANSIONARY. Banks hold ~all modeled debt, so B_gov ~= 0.076 x
        # bank assets ~= 0.97 (about 24% of annual GDP).
        B_gov_D_ss=0.98,        B_gov_F_ss=0.98,

        # Default risk (Bocola 2016): the PRICED default probability pi_t is an
        # exogenous input path to the solver (his s-shock), built per experiment
        # in main.py. Only D is risky; the feared event is a pure haircut.
        recovery_rate_D=0.45,   # 55% haircut (Greek PSI 2012; Bocola D = 0.55)
        # s-process innovation sd. Bocola Table 2 posterior mean (his param(23)); with
        # rho_s = 0.95 the unconditional sd of s is 2.02, and the box covers +-2.16 of
        # those (his own coverage). The pre-2026-08-15 value 1.5075 was NOT an estimate:
        # it was backed out of "a single +2sd innovation must lift p^d from 0.1% to 2%",
        # which no collocation box can hold -- the box clip then cut the EFFECTIVE
        # persistence of s from 0.95 to 0.80 and removed ~62% of the long bond's
        # repricing.
        # 0.4455 = 0.63/sqrt(2), NOT the reported posterior mean, and the difference is a
        # QUADRATURE CONVENTION, not a recalibration. point_map.gh_nodes uses numpy's
        # hermegauss -- the PROBABILISTS' rule, whose nodes are already in sd units, so
        # s' = mean + sigma*node makes sigma a genuine sd. Bocola's GaussHermite.m returns
        # the PHYSICISTS' rule and applies the SAME sigma*node map, which silently
        # rescales his innovation by 1/sqrt(2): his SOLVED model behaves as if sigma =
        # 0.4455. Feeding his reported 0.63 into a correct probabilists' rule reproduces
        # his parameter but not his model. Measured cost of the mismatch: unconditional
        # sd of s 2.02 vs 1.43, ergodic E[p^d] 0.66% vs 0.27%, D-bond marked 0.9057 vs
        # 0.9178, and a steady-state credit spread 32.8 bp vs 25.5 bp. The comment here
        # reached this conclusion on 2026-08-15 and the value below was never changed.
        sigma_s=0.4455,

        # THE CENTRAL-BANK BACKSTOP: A STOCHASTIC LTRO (Bocola's own instrument).
        # With per-period probability phi_ltro the central bank offers collateralised
        # credit of size ltro_D/ltro_F. It is lent at the DEPOSIT RATE, so it changes the
        # COMPOSITION of the bank's funding -- divertable deposits for non-divertable
        # central-bank credit -- at an unchanged rate: every budget identity in the model
        # is untouched and the whole effect is in the incentive constraint,
        #     mu_ratio = N'/(lambda*A')  ->  (N' + m)/(lambda*(A' - m)).
        # See point_map.py and docs/ltro_backstop_plan.md.
        #
        # phi_ltro is a per-experiment SCALAR, not a state: a phi dimension would centre
        # its box at 0.5 and the steady state would stop being a collocation node. One
        # solve per activation, each EXACT at its own phi. 0.0 is the nesting value and
        # every non-backstop run sits there, so phi = 0 reproduces the no-backstop model
        # EXACTLY rather than to tolerance (test_recursive_nesting N3).
        phi_ltro=0.0,
        #
        # THE ENVELOPE IS SIZED TO RELIEVE THE CONSTRAINT, NOT TO UNBIND IT, and that is
        # the whole calibration decision. Backing E[Om]R = 1.16714 out of mu_ss = 0.01231
        # and solving for the facility that drives mu to zero:
        #
        #     m (share of quarterly GDP)     mu at the SS     mu in the crisis state
        #        0.0%                          0.01231              0.02339
        #        1.0%                          0.00617              0.01725
        #        2.0%                          0.00003              0.01110
        #        3.4%                          0                    0.00248
        #       40.0%  (Bocola's own)          0                    0
        #
        # 2.0% of quarterly GDP already unbinds the constraint at the steady state and
        # 3.4% unbinds it in the crisis. Bocola's 40% is roughly TWELVE times the size
        # that fully neutralises the crisis, and at his size mu = 0 with enormous margin
        # in every relieved state -- so the entire m = 1 coefficient set sits ON the KKT
        # kink, where mu = max(.,0) is C0, a Chebyshev interpolant returns mu > 0 where
        # the truth is 0, and the fitted-versus-exact read disagrees by more than the
        # response being measured. That is the identification pathology the 100 bp spread
        # target was adopted to escape; re-entering it through the facility size would
        # give it back.
        # Run 40% as a documented upper-bound variant WITH the caveat, not as a headline.
        #
        # 2.0% SHIPS, and it is chosen against the SOLVED rules, not the deterministic
        # algebra. On the no-facility four-regime solve the multiplier is 0.01979 at the
        # grid centre and 0.05637 at the crisis corner (p^d = 4.82%/qtr) -- both well
        # above the deterministic-SS 0.01231, because those are risk-priced rules. The
        # facility multiplies E[Om]R*n/lev by (1 + m/n)/(1 - lambda*m/lev), so holding
        # n, lev and E[Om] fixed:
        #
        #     m       mu at the centre        mu at the crisis corner
        #    0.6%     0.01614  (-18%)         0.05285  (-6%)
        #    1.2%     0.01248  (-37%)         0.04933  (-13%)
        #    2.0%     0.00760  (-62%)         0.04463  (-21%)
        #    4.0%     0        (ON THE KINK)  0.03285  (-42%)
        #
        # THE TABLE ABOVE IS THE WRONG TEST, and 2.0% was shipped on it and had to be
        # withdrawn. It reads mu at the GRID CENTRE, but the object that has to stay off
        # the kink is the STOCHASTIC REST POINT, where mu is 0.00983 -- half the centre's
        # 0.01979. Measured on the solved four-regime rules at 2.0%, mu at the rest point
        # ran 0.00983 (phi=0) -> 0.00627 (phi=0.5) -> 0.00000 (phi=1): the facility put
        # the ergodic point ON mu = 0 at full credibility, and the identification went
        # with it -- the fitted-vs-exact output bracket at phi = 1 was
        # [-0.1524%, +0.0148%], wider than the response and straddling zero.
        #
        # SIZE AGAINST THE REST POINT INSTEAD. The measured slope is
        # d(mu_rest)/d(phi*m) = -0.4915 per unit of facility, so
        #     mu_rest(phi=1) = 0.00983 - 0.4915*m,
        # and keeping a comfortable margin above the kink at FULL credibility -- the
        # worst case, since phi and m only ever enter as a product here -- gives
        #     m = 1.0%  ->  mu_rest(phi=1) = 0.0049, half the no-backstop value.
        # That is the largest envelope for which every reported phi stays off the kink.
        # A bigger facility is not "more policy", it is less identification.
        ltro_D=0.010,           ltro_F=0.0,
        # STATE-CONTINGENT ACTIVATION (None = the constant-phi design, unchanged).
        # With a constant phi the facility is offered in EVERY state, and measured, that
        # is where it does most of its work: at phi = 0.5 the multiplier falls 37% at the
        # ergodic rest point against 8.9% at the headline shock. A facility that bites
        # hardest in normal times is a permanent liquidity subsidy, not a backstop -- and
        # it MOVES THE STEADY STATE, which is why every activation rests somewhere
        # different and the cross-phi IRFs are not directly comparable.
        # Setting ltro_s_thr makes the offer probability logistic in the exogenous risk
        # factor:  phi(s) = phi_ltro / (1 + exp(-(s - ltro_s_thr)/ltro_s_width)).
        # The threshold is s at p^d = 1%/qtr. Measured on that profile:
        #     rest point     p^d = 0.10%/qtr  ->  phi = 0.010 x phi_ltro  (off)
        #     headline shock p^d = 1.98%/qtr  ->  phi = 0.800 x phi_ltro  (on)
        #     crisis corner  p^d = 4.82%/qtr  ->  phi = 0.962 x phi_ltro  (on)
        # That is a backstop: the rest point stays put, every activation shares one
        # steady state, and the whole effect lands in the states where the spread
        # actually is elevated.
        # ON BY DEFAULT. The constant-phi design distorts the steady state, and it does
        # so in the wrong direction: the facility is offered in EVERY state, so a fixed
        # envelope relieves a large FRACTION of a small multiplier in normal times and a
        # small fraction of a large one in a crisis. Measured at phi = 0.5:
        #     ergodic rest point (p^d = 0.10%/qtr)   mu 0.00983 -> 0.00619   -37.0%
        #     headline shock     (p^d = 1.98%/qtr)   mu 0.02365 -> 0.02154    -8.9%
        # Four times more relief where it is least needed. That is a permanent liquidity
        # subsidy rather than a backstop, and it is why every activation used to rest
        # somewhere different -- which in turn made the cross-phi IRFs incomparable,
        # because each was differenced against its own moved rest point.
        # ltro_s_thr = log(0.01/0.99) puts the trigger at p^d = 1%/qtr:
        #     rest point      p^d = 0.10%  ->  phi = 0.010 x phi_ltro   (off)
        #     headline shock  p^d = 1.98%  ->  phi = 0.800 x phi_ltro   (on)
        #     crisis corner   p^d = 4.82%  ->  phi = 0.962 x phi_ltro   (on)
        # so the steady state is left where it was and the whole effect lands in the
        # states where the spread actually is elevated. Set None for the constant-phi
        # design, which nests exactly and is what the earlier results were run under.
        ltro_s_thr=-4.59512,    ltro_s_width=0.5,
        # ltro_F = 0 targets the facility on the country in crisis. The actual 3-year
        # LTROs were euro-area-wide, which is ltro_F = ltro_D; that is a different
        # experiment (it also relieves the F bank, and the union deposit market carries
        # the difference), not a robustness check on this one.

        # Working capital (Neumeyer-Perri): firms pre-finance zeta x wage bill at
        # r_wc = rdep(-1) + lambda*mu/Omega. The only spread->output channel;
        # zeta = 0 nests it off exactly.
        zeta_wc_D=1.0,          zeta_wc_F=1.0,

        # Fiscal. Bocola's rule tau(S) = tau* exp{z} + gamma_tau*B with gamma_tau =
        # 1.0 (Table 1): taxes respond one-for-one to the debt LEVEL, so a default
        # (lower B) lowers taxes -- the fiscal-relief leg of the default event. The
        # bank-loss leg dominates (default contractionary) once exposure is 7.6%.
        # BOHN RULE STRENGTH, stated as the DEBT ROOT it delivers rather than as a
        # coefficient. government.govt_steady_state inverts
        #     root = (1-delta_b) + (delta_b - gamma_tau)/Q_B_ss
        # for gamma_tau, so the calibrated object is the persistence of the debt stock.
        # WHY THIS IS NOT Bocola's gamma_tau = 1 (Table 1), in either reading:
        #   - as a LEVEL coefficient, gamma_tau = 1 makes a 55% haircut a 53%-of-GDP
        #     tax windfall and leaves the sovereign MORE indebted after default;
        #   - as an ELASTICITY (Tax = Tax_ss*(B/B_ss)^phi) it does not stabilise debt at
        #     all here, because G_D = 0 leaves Tax_ss = 0.00297 against B_ss = 0.98, so
        #     dTax/dB = 0.003 and the root is 1.0002. Bocola's gamma_tau = 1 acts on a
        #     tax that is a real share of GDP; this one on 0.3% of it.
        #   - raising the ELASTICITY instead (phi = 15 for root 0.955) swings taxes
        #     0.29x-3.17x across the +-8% B band and to 6e-6 at the default node -- a
        #     convexity the mu=1 Chebyshev basis cannot carry. Measured: the SS spread
        #     went 129bp -> 250bp and the C_D response to risk flipped POSITIVE.
        # The linear level rule is exactly representable in the basis, and at recovery
        # 0.45 its relief leg is ~2% of GDP. At the old root 0.9929 (half-life 97
        # quarters) B_D marched past +20% and hit the collocation wall at q7 -- the
        # figure's "trough then flat recovery" WAS that wall. root = 0.93 (half-life 9.5
        # quarters) holds the peak debt deviation near +6%, inside the b_band = 0.12 the
        # risk grid now carries; root = 0.95 peaks at +8.3% and still escaped a 0.08
        # band. Peak deviation scales as impulse/(1-root), so the root and the band are
        # calibrated as a PAIR against the IRF's own law of motion.
        # The steady state is unchanged for any root (at B = anchor the rule returns
        # Tax_ss identically), so no bank or household recalibration follows.
        debt_root_D=0.93,       debt_root_F=0.93,
        G_D=0.0,                G_F=0.0,

        # COUNTRY SIZE. F is EIGHT TIMES D (2026-08-28): the audit found that a
        # symmetric union makes D half the union, so D's own sovereign shock moves the
        # union real deposit rate 45 bp/yr and cancels 78% of the credit-spread rise
        # before it reaches any firm's wage bill. Bocola's SV.C open economy has no such
        # feedback -- his R = 1/beta + 0.01*(B_for/gdp) is a WORLD rate. Every country
        # variable stays PER CAPITA of its own country and is unchanged at the steady
        # state; size_F enters ONLY where D and F quantities are aggregated (goods
        # market, union deposit clearing, the two sovereign markets, the union wealth
        # identity). p_ss = 1 and every SS ratio is preserved -- see steady_state.py.
        size_D=1.0,             size_F=8.0,

        # Trade / CES basket. HOME BIAS MUST SCALE WITH SIZE or trade cannot balance:
        # at p = 1 balanced trade needs size_D*(1-omega_D)*C = size_F*(1-omega_F)*C, so
        # (1-omega_F) = (1-omega_D)*size_D/size_F. omega_home_F is DERIVED from
        # omega_home_D and the sizes in get_calibration below rather than set here, so
        # the two can never fall out of step. The small country is the open one: D
        # imports 15% of its basket, F imports 1.875% of its.
        omega_home_D=0.85,      epsilon_trade=0.5,

        # Solver settings
        T=300,                 # risk-shock horizon (T=100 truncates, T=500 identical)
        tol_hh=1e-12,
        tol_dist=1e-12,
        tol_mkt=1e-12,         # SS stage-1 hybr xtol
        tol_transition=1e-10,  # 7T acceptance; do NOT tighten (hybr plateaus ~5e-11)
        n_jobs=0,              # FD-Jacobian workers; 0 -> os.cpu_count()
        use_numba=True,        # JIT EGM/distribution kernels; numpy fallback otherwise
    )
    # DERIVED: the size-consistent foreign-goods weight (see omega_home_D above).
    # size_F = size_D reproduces the symmetric calibration exactly.
    cal["omega_home_F"] = 1.0 - (1.0 - cal["omega_home_D"]) * cal["size_D"] / cal["size_F"]
    cal["omega_home"] = cal["omega_home_D"]      # legacy key: D's weight
    return cal
