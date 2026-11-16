# THE HEADLINE SHOCKS IN UNITS OF THE RISK PROCESS, WITH THE CODE'S OWN TRANSFORMATION (LOGISTIC p^d = 1/(1+e^-s),
# s' = (1-rho) s* + rho s + sigma eps, s* = logit(pd_rest)); AND THE PREREQUISITES OF THE BOHN 0.15 k = 2 SOLVE.
import json, pickle, sys
from pathlib import Path
import numpy as np
from scipy.stats import norm
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(ROOT))
from calibration.global_projection import get_calibration
from global_projection.solver_recursive.state_grid import s_process_params, default_prob, IS, IM, _s_bounds
from global_projection.solver_recursive import recursive_experiment as RE

cal = dict(get_calibration(), pd_box_top=0.12)
sp = s_process_params(cal)
rho, sig, s0 = sp["rho_s"], sp["sigma_s"], sp["s_star"]
sd_u = sig / np.sqrt(1 - rho ** 2)
lo, hi = _s_bounds(cal, None, None, None)
mean_pd = quad(lambda s: default_prob(s) * norm.pdf(s, s0, sd_u), s0 - 12 * sd_u, s0 + 12 * sd_u)[0]
out = dict(rho_s=rho, sigma_s=sig, s_star=s0, pd_rest=float(default_prob(s0)), sd_uncond=sd_u, mean_pd=mean_pd,
           box_pd=[float(default_prob(lo)), float(default_prob(hi))], shocks={})
for pd in (0.015, 0.02, 0.06):
    s = RE.s_from_pd(pd)
    d = s - s0
    tail = float(norm.sf(d / sd_u))
    # a level shock held to the impact quarter, then decaying at rho: p^d after h quarters
    path = [float(default_prob(s0 + rho ** h * d)) for h in (1, 2, 4, 8)]
    out["shocks"][str(pd)] = dict(s=s, ds=d, uncond_sd=d / sd_u, innov_sd=d / sig, tail_prob=tail,
                                  return_years=1 / tail / 4, pd_annual=1 - (1 - pd) ** 4, pd_after_h_1_2_4_8=path)
print(json.dumps(out, indent=1))
