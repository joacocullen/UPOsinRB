"""
Newton-GMRes-Hookstep convergence of an invariant solution via spooky.

Reads two YAML config files (--solver, --newton) and loads the initial
condition from the directory given by newton.yaml's `input_dir` field.
"""

import sys
import yaml
import argparse
import numpy as np
from types import SimpleNamespace

import spooky as sp
from spooky.solvers import SPECTER
from spooky.methods import DynSys


def load_config(path):
    with open(path) as f:
        return SimpleNamespace(**yaml.safe_load(f))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--solver', required=True, help='Path to solver.yaml')
    parser.add_argument('--newton', required=True, help='Path to newton.yaml')
    args = parser.parse_args()

    pm_solver = load_config(args.solver)
    pm_newton = load_config(args.newton)

    print(f"Ra = {pm_solver.ra:.2e}  Pr = {pm_solver.pr}  "
          f"Nx={pm_solver.Nx} Nz={pm_solver.Nz}")

    grid = sp.Grid2D_semi(Lx=pm_solver.Lx, Lz=pm_solver.Lz,
                          Nx=pm_solver.Nx, Nz=pm_solver.Nz, dt=pm_solver.dt)
    solver = SPECTER(grid,
                     nprocs=pm_solver.nprocs,
                     ra=pm_solver.ra,
                     pr=pm_solver.pr,
                     gamma=pm_solver.gamma,
                     solver='BOUSS',
                     ftypes=['vx', 'vz', 'th'],
                     precision='double')

    newt = DynSys(pm_newton, solver)

    if pm_newton.restart_iN == 0:
        fields = solver.load_fields(pm_newton.input_dir, pm_newton.start_idx)
        T  = pm_newton.T
        sx = pm_newton.sx
    else:
        import os
        restart_path = os.path.join(pm_newton.output_dir,
                                    f'iN{pm_newton.restart_iN:02}')
        fields = solver.load_fields(restart_path, 0)
        T, sx = newt.get_restart_values(pm_newton.restart_iN)

    U = newt.flatten(fields)
    if sx is not None:
        X = np.append(U, [T, sx])
    elif T is not None:
        X = np.append(U, T)
    else:
        X = U  # steady state with fixed Tconst

    print('Running Newton-Krylov-Hookstep solver...')
    newt.run_newton(X)
    print('Done.')


if __name__ == '__main__':
    main()
