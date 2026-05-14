import numpy as np
import matplotlib.pyplot as plt
import sequence_jacobian as sj
from sequence_jacobian import simple, solved, combine, create_model
from sequence_jacobian import grids, hetblocks


calibration_start = {

    #==> Household parameters (single deposit asset) efe
    'frisch': 1.0,    # Frisch elasticity of labour supply (1/varphi)
    'eis': 0.5,       # Intertemporal elasticity of substitution (1/sigma)
    'habit': 0.00,    # Habit formation parameter (h)
    'C_lag': 0.0,     # Previous Consumption

    #==> Deposit rate 
    'rdep': 0.01,

    #==> Bond rate
    'rb': 0.01,

    #==> Government bonds
    'B_supply': 0.6,

    #==> Policy parameters
    'T': 0.15,         # Linear tax on labour income
    'phi_T': 0.06,     # Fiscal adjustment rule

    ##==> Aggregate targets
    'Y':  1.00,
    'N':  1.00,
    'w':  0.70,

    ##==> Financial intermediary
    'f': 0.06,
    'lambda_gk': 0.116,
    'ksi': 0.5,                  # Capital adjustment cost
    'n_inter': 0.75,             # Financial intermediary net worth
    'theta':   4,                # Capital leverage (k_inter / n_inter)

    #==> Production
    'alpha': 0.35,               # Capital share in production
    'delta': 0.0125,             # Depreciation rate (quarterly)

    #==> Deposit grid (single asset)
    'nZ': 2,                       # Matrix loads with 19 Number of income states
    'nDep': 49,                    # Number of deposit grid points
    'Depmax': 100,                 # Maximum deposits

    #==> Not used unless nZ != 19
    'rho_z': 0.9,
    'sigma_z': 0.5,

    #==> Tobin's q (= 1 in steady state with standard adjustment cost calibration)
    'Q': 1.0,

    #==> Government bonds outstanding (= B_supply at SS) and fiscal-rule anchors
    'b_gov':   0.6,    # = B_supply
    'b_gov_ss': 0.6,   # SS anchor used in tax rule
    'T_ss':    0.15,   # SS tax rate anchor used in tax rule
}

calibration_hh = {
    **calibration_start,

    #==> Calibrated so deposit market clears
    'beta': 0.9696492489377851,

    #==> Dividend income (bank equity return to HH) in steady state
    'div': 0.11,
}



# 1. Initialize Household
def hh_init(dep_grid, z, rdep, eis):
    coh = (1 + rdep) * dep_grid[np.newaxis, :] + z[:, np.newaxis]
    Vdep = (1 + rdep) * coh ** (-1 / eis)
    return Vdep

# 2. Backward Step (HH Block)
@sj.het(exogenous='Pi', policy='dep', backward='Vdep', backward_init=hh_init)
def hh(Vdep_p, dep_grid, z, rdep, beta, eis):
    uc_nextgrid = beta * Vdep_p
    c_nextgrid = uc_nextgrid ** (-eis)
    coh = (1 + rdep) * dep_grid[np.newaxis, :] + z[:, np.newaxis]
    
    # Policy iteration
    dep = sj.interpolate.interpolate_y(c_nextgrid + dep_grid, coh, dep_grid)
    sj.misc.setmin(dep, dep_grid[0])
    
    c = coh - dep
    uce = c ** (-1 / eis)        # marginal utility; UCE aggregate needed by sdf block
    Vdep = (1 + rdep) * uce
    return Vdep, dep, c, uce

# 3. Grids and Income
def make_grids(rho_z, sigma_z, nZ, Depmax, nDep):
    e_grid, _, Pi = sj.grids.markov_rouwenhorst(rho_z, sigma_z, nZ)
    dep_grid = sj.grids.asset_grid(0, Depmax, nDep)
    return e_grid, Pi, dep_grid

def income(e_grid, w):
    # Standard HANK income: idiosyncratic productivity × real wage.
    # Z reaches households through w (labor block: w=(1-α)Y/N, N depends on Z),
    # so Z→labor→w→income→C is the transmission path without Z here directly.
    z = w * e_grid
    return z

hh_extended = hh.add_hetinputs([make_grids, income])

# 4. Steady State & Market Blocks
@simple
def smart_steady(theta, Y, T, n_inter, rdep, alpha, delta, f, N, B_supply, rb):
    K        = theta * n_inter
    phi_b    = B_supply / n_inter
    rk       = alpha * Y / K - delta
    rn       = theta * (rk - rdep) + phi_b * (rb - rdep) + rdep
    m        = n_inter * (1 - (1 - f) * (1 + rn))
    k_inter  = K
    I        = K * delta
    D_supply = (theta - 1) * n_inter + B_supply
    G        = T - rb * B_supply
    Z        = Y / ((K ** alpha) * (N ** (1 - alpha))) # Fixed: Changed Ze to Z
    rdep_ante  = rdep
    return K, rk, rn, m, k_inter, I, D_supply, G, Z, rdep_ante

@simple
def market_clearing(Y, C, I, G, DEP, n_inter, theta, B_supply):
    # 'DEP' comes from the hh_extended block aggregate
    goods_mkt   = Y - C - I - G
    deposit_mkt = DEP - ((theta - 1) * n_inter + B_supply)
    return goods_mkt, deposit_mkt

@simple
def steady_auxilliary(theta, rk, rdep, delta, alpha, Y, K, N, lambda_gk, beta, ksi, rn):
    iota   = delta
    mpk    = alpha * (Y / K)
    w      = (1 - alpha) * Y / N
    Omega  = theta * lambda_gk / (beta * (1 + rn))
    nu     = beta * Omega * (rk - rdep)
    eta    = beta * Omega * (1 + rdep)
    gamma0 = delta ** ksi / (1 - ksi)
    gamma1 = -delta * ksi / (1 - ksi)
    return iota, mpk, w, Omega, nu, eta, gamma0, gamma1

@simple
def banker_div(rn, n_inter):
    div = rn * n_inter
    return div

@simple
def sdf(beta, UCE):
    SDF = beta * UCE(+1) / UCE
    return SDF

# 5. Create Model
ha = sj.create_model([hh_extended, smart_steady, market_clearing, steady_auxilliary, banker_div, sdf], name="Simple HA Model")

# 6. Solve Steady State
# Ensure calibration_start contains all required inputs (theta, Y, n_inter, etc.)
unknowns_ss = {'beta': (0.95, 0.98)} 
targets_ss  = ['deposit_mkt']

ss = ha.solve_steady_state(
    calibration_start, 
    unknowns_ss, 
    targets_ss, 
    solver='bisect' # Use bisect for 1D stable solving
)

print(f"Deposit market clearing: {ss['beta']}")

D = ss.internals['hh']['D'].sum(axis=0)
dep_grid = ss.internals['hh']['dep_grid']
plt.plot(dep_grid, D.cumsum())
plt.ylim([0.2, 1])
plt.xlim([0, 100])
plt.xlabel('Assets')
plt.ylabel('Cumulative distribution')
plt.show()


cali = ss
calibration = dict(ss)



# ── Off-Steady-State Equations ────────────────────────────────────────────────
@simple
def capital_adj(Y, K, Q, I, alpha, delta, gamma0, gamma1, ksi):
    iota        = I / K(-1)
    mpk         = alpha * Y / K(-1)
    rk          = (mpk + (1 - delta) * Q) / Q(-1) - 1
    q_res       = Q - 1 / (gamma0 * (1 - ksi) * iota ** (-ksi))
    capital_res = K - (1 - delta) * K(-1) - (gamma0 * iota ** (1 - ksi) + gamma1) * K(-1)
    return iota, mpk, rk, q_res, capital_res


@simple
def intermediation_IC(nu, eta, lambda_gk):
    theta = eta / (lambda_gk - nu)
    return theta


@simple
def bank_return(theta, rk, rdep, b_gov, n_inter, rb):
    # phi_b is predetermined (lagged balance-sheet ratio)
    phi_b_lag = b_gov(-1) / n_inter(-1)
    rn = theta(-1) * (rk - rdep) + phi_b_lag * (rb - rdep) + rdep
    return rn


@simple
def intermediation_P1(rk, rdep, nu, lambda_gk, eta, theta, SDF, f):
    # Bellman equations: bankers fund at deposit rate rd
    Omega_p1 = f + (1 - f) * lambda_gk * theta(+1)
    nu_res   = nu  - SDF * Omega_p1 * (rk(+1) - rdep(+1))
    eta_res  = eta - SDF * Omega_p1 * (1 + rdep(+1))
    return nu_res, eta_res

# K = theta * n_inter is the bank balance-sheet identity.
# Keeping it separate (not inside intermediation_P1) lets K be an outer unknown
# so capital_adj can use K without creating a cycle with intermediation_ne_solved.
@simple
def k_balance_sheet(theta, n_inter, K):
    K_res = K - theta * n_inter
    return K_res


@simple
def intermediation_P2(rn, n_inter, m, f):
    n_inter_val = (1 - f) * (1 + rn) * n_inter(-1) + m - n_inter
    return n_inter_val


@simple
def intermediation_P3(theta, n_inter, b_gov):
    D_supply = (theta - 1) * n_inter + b_gov
    return D_supply


@simple
def government(T, rb, b_gov):
    G = T - (1 + rb) * b_gov(-1) + b_gov
    return G


@simple
def mon_pol(rdep_ante):
    rdep = rdep_ante(-1)
    return rdep


# fiscal is split into two blocks to break the government <-> fiscal cycle.
# T only depends on lagged b_gov so it can be computed before government.
# The budget residual needs G (from government) so it comes after.
@simple
def tax_rule(b_gov, T_ss, b_gov_ss, phi_T):
    T = T_ss + phi_T * (b_gov(-1) - b_gov_ss)
    return T

@simple
def budget_residual(b_gov, G, T, rb):
    b_gov_res = (1 + rb) * b_gov(-1) + G - T - b_gov
    return b_gov_res


@simple
def rk_from_production(Y, K, alpha, delta):
    # rk at t uses capital chosen at t-1 (standard timing)
    rk = alpha * Y / K(-1) - delta
    return rk


@simple
def habit_law(C, C_lag):
    habit_res = C_lag - C(-1)
    return habit_res


@simple
def labor(Y, K, Z, alpha):
    N = (Y / Z / K(-1) ** alpha) ** (1 / (1 - alpha))
    w = (1 - alpha) * (Y / N)
    return N, w



print("Got here!1")



# Merge all financial-intermediary blocks into one solved block.
#
# Why: three separate blocks created a ring in the outer DAG —
#   bank_networth (needs rk) → capital_adj (needs K) → intermediation_ne (needs n_inter) → bank_networth
#
# By combining them the dependencies are internal to the solver:
#   - intermediation_IC uses (nu, eta) → theta
#   - bank_return uses theta(-1), n_inter(-1), rk (external) → rn
#   - intermediation_P2 uses rn → n_inter_val  [target for n_inter unknown]
#   - intermediation_P1 uses (nu, eta, theta, rk, SDF) → nu_res, eta_res  [targets for nu, eta]
#   - intermediation_P3 uses theta, n_inter → D_supply
# No cycle exists within this DAG.
financial_solved = combine([
    intermediation_IC,
    bank_return,
    intermediation_P2,
    intermediation_P1,
    intermediation_P3,
]).solved(
    unknowns={
        'n_inter': (1e-6, 10 * float(cali['n_inter'])),
        'nu':  float(cali['nu']),
        'eta': float(cali['eta']),
    },
    targets=['n_inter_val', 'nu_res', 'eta_res'],
    solver='broyden_custom'
)

print("financial_solved inputs: ", financial_solved.inputs)
print("financial_solved outputs:", financial_solved.outputs)



# ── 1. Finalize the Model Assembly ───────────────────────────────────────────

# Assemble all blocks into the full dynamic model
# We include the 'solved' blocks which handle internal feedback loops (like bank net worth)
ha_full = sj.create_model([
    hh_extended,
    financial_solved,
    k_balance_sheet,
    capital_adj,
    sdf,
    mon_pol,
    tax_rule,
    budget_residual,
    labor,
    habit_law,
    banker_div,
    market_clearing
], name="Full HANK Model")

print("Full Model Inputs:", ha_full.inputs)
print("Full Model Outputs:", ha_full.outputs)

# ── 2. Define Transition Path Unknowns and Targets ───────────────────────────

# For the transition path, we need to specify which variables adjust to clear markets.
# Typical choices: rdep_ante (monetary policy) and Y (output/goods market).
# SSJ's H_Z only includes targets with a computational path from Z (partial Jacobian).
# q_res has no Z-path (I is at SS, Q and K held fixed as outer unknowns in H_Z).
# b_gov_res has no Z-path (G is exogenous, b_gov is an outer unknown held fixed in H_Z).
# Only goods_mkt, deposit_mkt, K_res are reachable from Z; use a 3×3 system.
# Q and b_gov stay at their SS values for the transition path (standard first-pass simplification).
unknowns_tp = ['rdep_ante', 'Y', 'K']
targets_tp  = ['goods_mkt', 'deposit_mkt', 'K_res']

# ── 3. Jacobian Calculation ──────────────────────────────────────────────────

T = 300 
exogenous = ['Z'] # We want to see how the model reacts to TFP shocks

print(f"Computing Jacobian for horizon T={T}...")
G = ha_full.solve_jacobian(ss, unknowns_tp, targets_tp, exogenous, T=T)

# ── 4. Impulse Response Functions (Example: TFP Shock) ───────────────────────

# Define a persistent TFP shock
rho_Z = 0.9
dZ = 0.01 * rho_Z ** np.arange(T)

# Compute linear IRFs using the Jacobian
irf = {}
for var in ['Y', 'C', 'K', 'rdep']:
    irf[var] = G[var]['Z'] @ dZ

# ── 5. Visualization of Results ──────────────────────────────────────────────

plt.figure(figsize=(10, 6))

plt.subplot(2, 2, 1)
plt.plot(irf['Y'][:50])
plt.title('Output (Y)')
plt.grid(True)

plt.subplot(2, 2, 2)
plt.plot(irf['C'][:50])
plt.title('Consumption (C)')
plt.grid(True)

plt.subplot(2, 2, 3)
plt.plot(irf['K'][:50])
plt.title('Capital (K)')
plt.grid(True)

plt.subplot(2, 2, 4)
plt.plot(irf['rdep'][:50])
plt.title('Deposit Rate (rdep)')
plt.grid(True)

plt.tight_layout()
plt.show()

# ── 6. Final Calibration Export ──────────────────────────────────────────────

# Storing the final steady state for further use
final_calibration = ss
print("Steady State beta:", ss['beta'])
print("Steady State Output:", ss['Y'])

