# BOCOLA (2016) CALIBRATION: TABLE 1/2 TARGETS -> BACKED-OUT STRUCTURAL PARAMS.
# The paper (footnote 14) reparameterizes [beta, lambda, omega, delta, chi,
# iota, tau_star, a1, a2] by balanced-growth-path targets. This module holds
# the fixed/estimated values and solves the BGP relations for the structural
# parameters, then reports every implied great ratio so the mapping is
# checkable against the paper rather than asserted.
import numpy as np


def base_params():
    # FIXED (TABLE 1) AND ESTIMATED-POSTERIOR-MEAN (TABLE 2) INPUTS.
    return dict(
        # Table 1 (fixed / external)
        exp_bg=0.076,        # Q_B*B / total bank assets
        pi_mat=0.056,        # fraction of bonds maturing each quarter
        g_star=0.198,        # spending-output ratio on BGP (exp{g*})
        rho_g=0.894, sig_g=0.012,
        iy_bg=0.213,         # investment-output ratio
        l_bg=0.318,          # hours worked
        R_bg=1.003,          # quarterly gross risk-free/deposit rate
        alpha=0.30,          # capital share
        nu=0.5,              # inverse Frisch elasticity
        lev_bg=5.0,          # bank leverage
        gamma_tau=1.0,       # fiscal rule coefficient on debt
        adj_bg=0.0,          # capital adjustment cost on BGP (Q_K normalized 1)
        qb_bg=1.0,           # bond price normalized to 1 on BGP
        # Table 2 (estimated posterior means)
        mu_bg=0.001,         # Lagrange multiplier on BGP  (mu^bg x100 = .10)
        psi=0.96,            # banker survival probability
        xi=0.42,             # elasticity of Tobin's q wrt I/K
        gamma=0.43 / 400.0,  # deterministic tech growth per quarter
        rho_z=0.51, sig_z=0.0059,
        s_star=-7.06, rho_s=0.95, sig_s=0.63,
        # default
        haircut=0.55,        # D, Greek PSI 2012
    )


def solve_bgp(p=None, delta_mode="standard"):
    # SOLVE THE BALANCED-GROWTH-PATH SYSTEM FOR THE STRUCTURAL PARAMETERS.
    # delta_mode="standard" (PRIMARY, user-chosen 2026-07-23): delta fixed at
    #   a conventional 0.025 -> sensible K/Y=2.66yr and a well-conditioned
    #   global solve; i/y is then an implied CHECK (0.277 vs 0.213 target).
    # delta_mode="euler" (ROBUSTNESS VARIANT): delta backed out jointly from
    #   the capital Euler + Jermann replacement + i/y target (the literal
    #   footnote-14 reading). Hits every stated target but implies a
    #   near-frictionless K/Y=10yr whose global solve is markedly harder to
    #   converge (see docs/bocola2016_replication.md).
    if p is None:
        p = base_params()
    G = np.exp(p["gamma"])                      # gross tech growth factor
    beta = G / p["R_bg"]                        # household Euler on BGP
    # marginal value of wealth av solves av = [(1-psi)+psi*av]/(1-mu) on BGP
    av = (1.0 - p["psi"]) / (1.0 - p["mu_bg"] - p["psi"])
    lam = av / p["lev_bg"]                       # binding leverage = av/lambda
    Lam_hat = av * (1.0 - p["mu_bg"]) / p["R_bg"]   # E[Lambda_hat'] on BGP
    # capital credit spread from the capital Euler on BGP
    RK_minus_R = lam * p["mu_bg"] / Lam_hat
    R_K = p["R_bg"] + RK_minus_R                 # gross capital return on BGP

    alpha = p["alpha"]
    if delta_mode == "euler":
        # R_K = (1-delta) + alpha*Y/K ; I/K = G-1+delta (Jermann, Phi(x*)=x*);
        # i/y = (I/K)/(Y/K). Two eqs in {delta, Y/K}: eliminate Y/K.
        # alpha*Y/K = R_K-(1-delta) = delta + (R_K-1)
        #   Y/K = (delta + R_K-1)/alpha
        # I/K = G-1+delta ;  i/y = (G-1+delta) / ((delta+R_K-1)/alpha)
        #   = alpha*(G-1+delta)/(delta+R_K-1) = iy_bg  ->  linear in delta
        a = alpha - p["iy_bg"]
        b = alpha * (G - 1.0) - p["iy_bg"] * (R_K - 1.0)
        delta = -b / a
    else:
        delta = 0.025
    YK = (delta + (R_K - 1.0)) / alpha           # Y/K on BGP
    IK = G - 1.0 + delta                          # I/K on BGP (Jermann repl.)
    iy = IK / YK                                  # investment-output (check)
    KY = 1.0 / YK

    # Jermann Phi(x)=a1 x^(1-xi)+a2 normalized: Phi'(x*)=1 (Q_K=1) and
    # Phi(x*)=x* (adj_bg=0). x*=I/K.
    xstar = IK
    a1 = xstar ** p["xi"] / (1.0 - p["xi"])       # from a1(1-xi)x*^(-xi)=1
    a2 = xstar - a1 * xstar ** (1.0 - p["xi"])    # from Phi(x*)=x*

    # labor FOC: chi*l^nu * c = (1-alpha)*Y/l  (GHH-free separable log utility)
    #   -> chi = (1-alpha)*(Y/l)/(c * l^nu). Use BGP shares per unit output:
    #   with Y=1 units, c/Y = 1 - i/y - g, w = (1-alpha)Y/l.
    cy = 1.0 - iy - p["g_star"]
    chi = (1.0 - alpha) / (cy * p["l_bg"] ** (1.0 + p["nu"]))

    # bond coupon iota from q_b^bg=1: bond Euler on BGP with no default risk,
    # R_B = [pi + (1-pi)(iota+Q_B)]/Q_B = R_B_required. On BGP bonds are held
    # by banks at the same shadow value as capital: R_B-R = lam*mu/Lam_hat too
    # (single lambda). So R_B = R_K. With Q_B=1: pi+(1-pi)(iota+1) = R_B.
    R_B = R_K
    iota = (R_B - p["pi_mat"]) / (1.0 - p["pi_mat"]) - 1.0

    # bond/capital quantities from exp_bg = Q_B*B/(Q_K*K+Q_B*B), Q_K=Q_B=1
    #   B/K = exp/(1-exp)
    BK = p["exp_bg"] / (1.0 - p["exp_bg"])

    # --- BGP ALLOCATION LEVELS (units: Y_bg = 1 via a TFP normalization) ---
    # Production Y = Zbar * K^alpha (e^z l)^(1-alpha); Zbar picked so Y_bg = 1
    # at K_bg = KY, l_bg, Dz = gamma. This only sets units for the detrended
    # levels; ratios above are Zbar-invariant.
    Zbar = 1.0 / (KY ** alpha * (G * p["l_bg"]) ** (1.0 - alpha))
    Y_bg = 1.0
    K_bg = KY                       # K/Y = KY with Y = 1
    C_bg = cy                       # c/y
    I_bg = iy                       # i/y
    B_bg = BK * K_bg                # B/K * K
    # bank block on BGP (Q_K = Q_B = 1, d = 0)
    Zrent = alpha * Y_bg / K_bg     # profit per unit capital Z = alpha Y/K
    X_bg = ((1.0 - delta) + Zrent) * K_bg \
        + (p["pi_mat"] + (1.0 - p["pi_mat"]) * (iota + 1.0)) * B_bg
    N_bg = (K_bg + B_bg) / p["lev_bg"]     # binding leverage (K+B)/N = lev
    dep_bg = K_bg + B_bg - N_bg            # deposits b' = assets - net worth
    P_bg = p["R_bg"] * dep_bg              # gross payment owed next period
    omega = (N_bg - p["psi"] * (X_bg - P_bg)) / X_bg   # entrant endowment share
    # fiscal intercept tau* from the government budget on BGP
    #   B*G = B[1+(1-pi)iota] + g*Y - tau ; tau = tau* e^gamma + gamma_tau*B
    tau_bg = B_bg * (1.0 + (1.0 - p["pi_mat"]) * iota - G) + p["g_star"] * Y_bg
    tau_star = np.exp(-p["gamma"]) * (tau_bg - p["gamma_tau"] * B_bg)

    out = dict(**p, G=G, beta=beta, av_bg=av, lam=lam, Lam_hat_bg=Lam_hat,
               R_K_bg=R_K, RK_minus_R=RK_minus_R, delta=delta, YK=YK, KY=KY,
               IK=IK, iy_implied=iy, cy=cy, a1=a1, a2=a2, xstar=xstar,
               chi=chi, R_B_bg=R_B, iota=iota, BK=BK, delta_mode=delta_mode,
               Zbar=Zbar, Y_bg=Y_bg, K_bg=K_bg, C_bg=C_bg, I_bg=I_bg,
               B_bg=B_bg, Zrent_bg=Zrent, X_bg=X_bg, N_bg=N_bg, dep_bg=dep_bg,
               P_bg=P_bg, omega=omega, tau_bg=tau_bg, tau_star=tau_star)
    return out


def bgp_report(bgp):
    # HUMAN-READABLE BGP TABLE FOR CHECKING AGAINST THE PAPER.
    lines = [
        "BALANCED-GROWTH-PATH  (Bocola 2016 calibration)",
        f"  delta_mode                 = {bgp['delta_mode']}",
        f"  G (tech growth)            = {bgp['G']:.6f}   (gamma={bgp['gamma']:.6e})",
        f"  beta                       = {bgp['beta']:.6f}",
        f"  av (marg. value wealth)    = {bgp['av_bg']:.6f}",
        f"  lambda (divertable share)  = {bgp['lam']:.6f}",
        f"  E[Lambda_hat']             = {bgp['Lam_hat_bg']:.6f}",
        f"  R_K - R (credit spread)    = {bgp['RK_minus_R']:.6e}  "
        f"({bgp['RK_minus_R']*400*100:.1f} bps ann)",
        f"  delta (depreciation)       = {bgp['delta']:.6f}  "
        f"({bgp['delta']*4*100:.2f}%/yr)",
        f"  K/Y  (capital-output)      = {bgp['KY']:.3f} quarters "
        f"({bgp['KY']/4:.2f} yrs)",
        f"  I/K                        = {bgp['IK']:.6f}",
        f"  i/y implied                = {bgp['iy_implied']:.4f}  "
        f"(target {bgp['iy_bg']:.4f})",
        f"  c/y                        = {bgp['cy']:.4f}",
        f"  a1, a2 (Jermann)           = {bgp['a1']:.4f}, {bgp['a2']:.4f}",
        f"  chi (labor weight)         = {bgp['chi']:.4f}",
        f"  iota (bond coupon)         = {bgp['iota']:.6f}",
        f"  B/K (bond/capital)         = {bgp['BK']:.4f}",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    for mode in ("euler", "standard"):
        print(bgp_report(solve_bgp(delta_mode=mode)))
        print()
