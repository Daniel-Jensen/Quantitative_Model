"""Run the pipeline through build_and_solve and save IRFs for comparison.

Usage:
    /opt/anaconda3/envs/ssj/bin/python code/dump_irfs.py OUT.npz
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from calibration import get_calibration
from steady_state import solve_steady_state
from ic_delta_calibration import calibrate_ic_delta
from depreciation_calibration import calibrate_depreciation
from full_model import build_and_solve

KEYS = ('Y_D', 'C_D', 'I_D', 'n_inter_D', 'K_D', 'b_gov_D', 'w_D', 'N_D',
        'p', 'q_b_D', 'spread_rb', 'Y_F', 'C_F', 'I_F', 'n_inter_F')


def main(out_path):
    r = calibrate_depreciation(calibrate_ic_delta(
        solve_steady_state(get_calibration())))
    m = build_and_solve(r)
    payload = {}
    for tag in ('irfs_def_D', 'irfs_Z_D'):
        for k in KEYS:
            if k in m[tag]:
                payload[f'{tag}__{k}'] = np.asarray(m[tag][k])
    # Store the SS level for EVERY dumped series. A missing ss__ entry tempts a
    # consumer into a divisor of 1.0, which silently reports a level deviation
    # as a percentage -- exactly the mislabelling CLAUDE.md records (n_inter_ss
    # = 2.138 and K_ss = 10.8 are not ~1, and a past bug mislabelled those by
    # 2.1x and 10x). I_D in particular has no SS entry in ss_final and must be
    # reconstructed as delta*K.
    ss = m['ss_final']
    for k in KEYS:
        if k in ('I_D', 'I_F'):
            suf = k[-1]
            payload[f'ss__{k}'] = np.asarray(
                float(ss[f'delta_{suf}']) * float(ss[f'K_{suf}']))
        else:
            payload[f'ss__{k}'] = np.asarray(float(ss[k]))
    missing = [k for k in KEYS if f'ss__{k}' not in payload]
    assert not missing, f'no SS level stored for: {missing}'
    np.savez(out_path, **payload)
    print(f'wrote {out_path} with {len(payload)} arrays '
          f'({len(KEYS)} series, all with SS levels)')


if __name__ == '__main__':
    main(sys.argv[1])
