

var
    c        // consumption
    N        // labor
    w        // wage
    Y        // output
    K        // end-of-period capital stock
    I        // investment
    Q        // Tobins Q (price of installed capital)
    iota     // investment rate I/K(-1)
    mpk      // marginal product of capital
    rk       // ex-post return on capital
    r        // ex-ante real deposit/bond rate (set at t, paid at t+1)
    rn       // ex-post return on bank net worth
    nu       // GK auxiliary: marginal value of assets
    eta      // GK auxiliary: marginal value of net worth
    theta    // bank leverage Q*K / nB
    nB       // bank net worth
    m        // transfer to newly born bankers
    B        // government bonds outstanding
    T        // lump-sum tax
    G        // government spending
    Z        // TFP
    spread   // credit spread rk - r (reporting only)
    lev      // leverage = Q*K/nB (reporting only, equals theta)
;

varexo
    eps_G    // government spending innovation
    eps_Z    // TFP innovation (set to zero in baseline experiment)
;

parameters
    betta        // HH discount factor
    sigma        // CRRA coefficient (= 1/EIS in notebook)
    frisch       // Frisch elasticity
    vphi         // labor disutility
    alpha        // capital share
    delta        // depreciation
    ksi          // Jermann adjustment cost curvature
    gamma0       // Jermann adjustment cost level
    gamma1       // Jermann adjustment cost shifter
    f            // banker exit probability  (1-f survive each period)
    lambda_gk    // fraction of assets divertible
    betta_b      // banker discount factor
    omega_m      // banker entry transfer / Y ratio
    rho_G        // persistence of G shock
    rho_Z        // persistence of TFP
    G_ss         // SS gov spending
    B_ss         // SS gov debt
    phi_b        // tax response to debt deviations
    T_ss         // SS taxes
    nB_ss        // SS bank net worth (used as initval helper)
;

/* ---------- Calibration (standard quarterly) ---------------------------- */

betta     = 0.99;          // -> r_ss ~ 1.01% per quarter
sigma     = 1.00;          // log utility (notebook used eis=0.5 i.e. sigma=2; we pick the standard value)
frisch    = 0.5;           // notebook value
alpha     = 0.33;          // capital share (notebook 0.35, we use canonical 0.33)
delta     = 0.025;         // 10% annual depreciation
ksi       = 0.5;           // notebook value (Jermann)
f         = 0.028;         // banker exit (survival 0.972 - GK 2011)
lambda_gk = 0.381;         // GK 2011 baseline
betta_b   = 0.99;          // banker discount = HH discount in baseline
omega_m   = 0.002;         // small entry transfer (calibrated to keep nB>0)
rho_G     = 0.90;          // notebook value
rho_Z     = 0.95;          // standard

// Jermann adjustment-cost coefficients (so that in SS Q=1 and Phi(delta)=delta)
gamma0 = delta^ksi / (1 - ksi);
gamma1 = -delta * ksi / (1 - ksi);

// Government targets (pin G_ss & B_ss; T from gov BC in SS)
G_ss = 0.20;               // ~20% of Y (Y normalized to ~1)
B_ss = 0.40;               // ~40% of Y
phi_b = 0.10;              // strong enough for non-explosive debt
T_ss  = G_ss + (1/betta - 1) * B_ss;   // T = G + r*B in SS

// Filled by initval/ss block solver; harmless placeholder
nB_ss = 0.15;

// Labor disutility - calibrated by hand so that Y_ss = 1 (notebook normalization).
// With the rest of the calibration this implies N_ss ~ 0.34.
vphi = 28.3882;

/* ---------- Model ------------------------------------------------------- */

model;

/* 1. Household Euler (single liquid asset, return r set at t for t+1) */
c^(-sigma) = betta * (1 + r) * c(+1)^(-sigma);

/* 2. Labor supply */
vphi * N^(1/frisch) = w * c^(-sigma);

/* 3. Production */
Y = Z * K(-1)^alpha * N^(1-alpha);

/* 4. Marginal products */
w   = (1 - alpha) * Y / N;
mpk = alpha * Y / K(-1);

/* 5. Investment rate */
iota = I / K(-1);

/* 6. Capital LoM with Jermann adjustment cost
      K = (1 - delta) K_{-1} + Phi(iota) K_{-1},
      Phi(iota) = gamma0 * iota^(1-ksi) + gamma1                          */
K = (1 - delta) * K(-1) + (gamma0 * iota^(1 - ksi) + gamma1) * K(-1);

/* 7. Tobin's Q from firm's FOC:  Q * Phi'(iota) = 1                       */
Q = 1 / ( gamma0 * (1 - ksi) * iota^(-ksi) );

/* 8. Ex-post return on capital (notebook's spec: ignores adj-cost wedge,
      consistent with banks earning the gross capital return)              */
1 + rk = ( mpk + (1 - delta) * Q ) / Q(-1);

/* ---- Gertler-Karadi banking block (matches notebook timing) ------------ */

/* 9. nu : marginal value of assets
      nu_t = beta_b * (f + (1-f)*lambda_gk*theta_{t+1}) * (rk_{t+1}-r_t)   */
nu = betta_b * ( f + (1 - f) * lambda_gk * theta(+1) ) * ( rk(+1) - r );

/* 10. eta : marginal value of net worth                                   */
eta = betta_b * ( f + (1 - f) * lambda_gk * theta(+1) ) * ( 1 + r );

/* 11. Incentive constraint -> leverage                                    */
theta = eta / (lambda_gk - nu);

/* 12. Bank balance sheet : value of capital = leverage x equity           */
Q * K = theta * nB;

/* 13. Net worth law of motion                                             */
nB = (1 - f) * (1 + rn) * nB(-1) + m;

/* 14. Return on bank equity (notebook: rn = theta_{-1}*(rk - r_{-1}) + r_{-1}) */
rn = theta(-1) * ( rk - r(-1) ) + r(-1);

/* 15. Banker entry (notebook cell 20: m = omega * Y)                      */
m = omega_m * Y;

/* ---- Government -------------------------------------------------------- */

/* 16. Debt LoM */
B = (1 + r(-1)) * B(-1) + G - T;

/* 17. Tax rule (debt-stabilizing feedback) */
T = T_ss + phi_b * ( B(-1) - B_ss );

/* 18. G follows AR(1) */
G = (1 - rho_G) * G_ss + rho_G * G(-1) + eps_G;

/* 19. TFP (kept exogenous; baseline shock is on G) */
log(Z) = rho_Z * log(Z(-1)) + eps_Z;

/* 20. Goods market clearing (closes the system by Walras' law) */
Y = c + I + G;

/* Reporting */
spread = rk - r(-1);
lev    = Q * K / nB;

end;

/* ---------- Steady state ------------------------------------------------ */
/*
   We let Dynare's nonlinear solver find the steady state from initial
   guesses below.  The block is mutually consistent:

     r        = 1/betta - 1
     Q        = 1, iota = delta, K from mpk = alpha*Y/K
     rk       = mpk - delta
     theta    = eta/(lambda_gk - nu),
                with nu  = betta_b*(f+(1-f)*lambda_gk*theta)*(rk-r)
                     eta = betta_b*(f+(1-f)*lambda_gk*theta)*(1+r)
                so   theta = (1+r)/(rk - r + lambda_gk/[betta_b*(f+(1-f)*lambda_gk*theta)] ) ...
     nB       = K/theta  (since Q=1)
     m        = nB * (1 - (1-f)*(1+rn))   from nB LoM
     omega_m  is fixed; m must equal omega_m*Y -> this pins down Y given
              the rest.  In the baseline calibration omega_m=0.002 and the
              numbers work out to Y close to 1.
*/

// Initial values are the analytical/numerical SS computed offline,
// so Dynare's nonlinear solver converges in a handful of iterations.
initval;
    Z      = 1.0;
    r      = 0.010101;
    Q      = 1.0;
    iota   = 0.025;
    Y      = 1.0000;
    N      = 0.3437;
    K      = 8.7446;
    I      = 0.2186;
    c      = 0.5814;
    w      = 1.9495;
    mpk    = 0.037738;
    rk     = 0.012738;

    theta  = 6.5134;
    nu     = 0.006369;
    eta    = 2.440129;
    nB     = 1.3426;
    rn     = 0.027274;
    m      = 0.002000;

    G      = G_ss;        // 0.20
    B      = B_ss;        // 0.40
    T      = T_ss;        // 0.2040

    spread = 0.002637;
    lev    = 6.5134;
end;

steady;
check;

/* ---------- Shock and IRFs ---------------------------------------------- */

shocks;
    var eps_G; stderr 0.01;     // 1% of SS Y, like notebook's dG[0]=1
    var eps_Z; stderr 0.01;
end;

stoch_simul(order=1, irf=100) Y c I N r rk spread Q K nB theta G B T;
