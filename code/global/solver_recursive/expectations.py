# TWO-BRANCH QUADRATURE EXPECTATIONS FOR THE RECURSIVE SOLVER.
# The banker/household expectations that enter bank_backward's FOCs are
# computed as GENUINE integrals over next period's state: Gauss-Hermite nodes
# for the s innovation x the default realization d' in {0,1} weighted by
# p^d(s_t). Next-period objects come from the fitted decision rules evaluated
# at the implied next state -- the bad state is the SAME rules at a reachable
# point on the grid. Z and G are deterministic
# for the D-risk experiment (agreed Tier cuts), so eps_s is the only
# continuous innovation: 7 GH nodes x 2 branches = 14 continuation points.
import numpy as np

from solver_recursive.state_grid import default_prob

# state indices, matching state_grid.STATE_NAMES (7-state)
IK_D, IK_F, IP_D, IP_F, IB_D, IS, IZ = range(7)


def gh_nodes(n=7):
    # GAUSS-HERMITE NODES/WEIGHTS FOR ONE STANDARD-NORMAL INNOVATION.
    x, w = np.polynomial.hermite_e.hermegauss(n)
    return x, w / w.sum()


def next_states(S_end, s_now, sproc, eps):
    # NEXT-PERIOD STATE POINTS GIVEN END-OF-PERIOD STOCKS AND s DRAWS.
    # S_end = [K_D', K_F', P_D', P_F', B_D'] (period-map outputs, d'-independent
    # except B_D' -- the haircut applies on ARRIVAL in d'=1); eps (m,) draws.
    m = eps.size
    s_next = (1.0 - sproc["rho_s"]) * sproc["s_star"] \
        + sproc["rho_s"] * s_now + sproc["sigma_s"] * eps
    out = np.empty((m, 6))
    out[:, :5] = S_end
    out[:, IS] = s_next
    return out


def continuation(rules, S_next, d_next, cal, haircut_states=None):
    # RULE-IMPLIED CONTINUATION OBJECTS IN REGIME d_next AT THE NEXT STATES.
    # Returns what bank_backward's two-branch FOCs consume per node: alpha',
    # C', Q_bD', Q_bF', p', rdep' -- evaluated with the d_next coefficient set.
    # haircut_states: in d'=1 the D-bond stock arrives haircut (recovery), so
    # the evaluation point differs; the caller passes the adjusted states.
    S_eval = rules.grid.clip(haircut_states if haircut_states is not None
                             else S_next)
    r = rules.eval_all(d_next, S_eval)
    return dict(alpha_D=r["alpha_D"], alpha_F=r["alpha_F"],
                C_D=np.maximum(r["C_D"], 1e-8),
                C_F=np.maximum(r["C_F"], 1e-8),
                Q_bD=np.maximum(r["Q_bD"], 1e-4),
                Q_bF=np.maximum(r["Q_bF"], 1e-4),
                p=np.maximum(r["p"], 1e-3),
                rdep_D=r["rdep_D"], rdep_F=r["rdep_F"])


def two_branch_weights(s_now, n_gh=7):
    # QUADRATURE WEIGHTS FOR THE (eps_s x d') MIXTURE AT CURRENT s.
    eps, w = gh_nodes(n_gh)
    pd = float(default_prob(s_now))
    return eps, w * (1.0 - pd), w * pd, pd     # (nodes, w_nd, w_def, p^d)


def bank_expectations(point, S_end, rules, cal, sproc, ss, n_gh=7):
    # THE EXPECTATION TERMS OF bank_backward'S TWO-BRANCH FOCs, AS INTEGRALS.
    # Produces, for ONE grid point in ONE current regime, every E[.] the
    # occasionally-binding FOC block needs:
    #   E_Lamhat_X            = E[Lambda_hat_X']                     (X = D,F)
    #   E_Lamhat_rk_X         = E[Lambda_hat_X' * rk_X']
    #   E_Lamhat_paydD        = E[Lambda_hat' * D-bond payoff']  (haircut in d'=1)
    #   E_Lamhat_paybF        = E[Lambda_hat' * F-bond payoff']
    #   (cross-border legs with p'/p conversions, same structure)
    # Lambda_hat_X' = beta_X*(C_X/C_X')*[(1-f_X)+... ] uses the SAME kernel
    # convention as bank.py (Omega = beta*[f + (1-f)*alpha']); rk' is built
    # from the same capital-return formula as capital.py on the continuation
    # rules. IMPLEMENTATION LANDS WITH recursive_residual.py: the formulas are
    # assembled there against the untouched bank.py conventions and verified
    # by the nesting test (pi = 0 must reproduce the deterministic path).
    raise NotImplementedError(
        "bank_expectations is assembled alongside recursive_residual.py; "
        "see the module header for the integral list and conventions.")
