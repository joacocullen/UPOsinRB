"""
Compute Floquet multipliers of a converged periodic orbit via spooky's
Arnoldi iteration (DynSys.floq_exp).

Usage
-----
    python floquet_spectrum.py \\
        --solver ../convergence/params/PO1_Ra1.25e6/solver.yaml \\
        --newton ../convergence/params/PO1_Ra1.25e6/newton.yaml \\
        --state-dir ../data/invariant_solutions/PO1 \\
        --n-eigs 100 \\
        --output PO1_Ra1.25e6_floquet.npz
"""

import sys
import yaml
import argparse
import numpy as np
from types import SimpleNamespace
from pathlib import Path

import spooky as sp
from spooky.solvers import SPECTER
from spooky.methods import DynSys

sys.path.insert(0, str(Path(__file__).parent.parent / 'convergence'))
from io_utils import make_solver


def load_config(path):
    with open(path) as f:
        return SimpleNamespace(**yaml.safe_load(f))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--solver',    required=True)
    parser.add_argument('--newton',    required=True)
    parser.add_argument('--state-dir', required=True,
                        help='Directory containing SPECTER .out files for the converged state')
    parser.add_argument('--state-idx', type=int, default=0)
    parser.add_argument('--n-eigs',   type=int, default=100)
    parser.add_argument('--tol',      type=float, default=1e-10)
    parser.add_argument('--output',   default='floquet.npz')
    args = parser.parse_args()

    pm_solver = load_config(args.solver)
    pm_newton = load_config(args.newton)

    solver = make_solver(pm_solver.ra, pm_solver.pr,
                         pm_solver.Nx, pm_solver.Nz, pm_solver.dt,
                         nprocs=pm_solver.nprocs)

    newt = DynSys(pm_newton, solver)

    fields = solver.load_fields(args.state_dir, args.state_idx)
    T  = pm_newton.T
    sx = pm_newton.sx

    U = newt.flatten(fields)
    if sx is not None:
        X = np.append(U, [T, sx])
    else:
        X = np.append(U, T)

    print(f'Computing {args.n_eigs} Floquet exponents (Ra={pm_solver.ra:.2e}, T={T})...')
    eigval_H, eigvec_H, Q = newt.floq_exp(X, args.n_eigs, tol=args.tol)

    multipliers = np.exp(eigval_H * T)
    idx = np.argsort(-np.abs(multipliers))
    multipliers = multipliers[idx]
    eigval_H    = eigval_H[idx]
    eigvec_H    = eigvec_H[:, idx]

    print('Top 5 Floquet multipliers:')
    for i, mu in enumerate(multipliers[:5]):
        print(f'  mu_{i+1} = {mu:.6f}   |mu| = {abs(mu):.6f}')

    np.savez(args.output,
             multipliers=multipliers,
             floquet_exponents=eigval_H,
             modes=eigvec_H,
             Ra=pm_solver.ra,
             T=T)
    print(f'Saved to {args.output}')


if __name__ == '__main__':
    main()
