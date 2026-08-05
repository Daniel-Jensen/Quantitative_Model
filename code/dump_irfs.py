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
    for k in ('Y_D', 'C_D', 'n_inter_D', 'K_D'):
        payload[f'ss__{k}'] = np.asarray(float(m['ss_final'][k]))
    np.savez(out_path, **payload)
    print(f'wrote {out_path} with {len(payload)} arrays')


if __name__ == '__main__':
    main(sys.argv[1])
