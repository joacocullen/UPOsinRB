"""
Recurrence analysis on a DNS trajectory saved as SPECTER binary files.

Computes the recurrence function
    G(t, tau) = min_s  ||T_s X(t+tau) - X(t)|| / ||X(t)||
over a grid of (t, tau) pairs to identify near-recurrences as initial
guesses for the Newton solver.

The horizontal shift T_s is optimised using the phase of the x-direction
cross-correlation in Fourier space (exact for periodic boundaries).

"""

import argparse
import numpy as np

from io_utils import make_solver


def flatten(fields):
    return np.concatenate([f.ravel() for f in fields])


def optimal_shift(U_a, U_b, Nx, n_fields, Nz):
    """
    Find the integer x-shift s* (in grid points) that minimises
    ||roll(U_b, s*, axis=0) - U_a||.
    Uses cross-correlation via FFT; O(Nx log Nx) per call.
    """
    # Reshape to (n_fields*Nz, Nx) so we can FFT along x at once
    A = U_a.reshape(n_fields, Nx, Nz).sum(axis=(0, 2))  # project onto x
    B = U_b.reshape(n_fields, Nx, Nz).sum(axis=(0, 2))
    corr = np.real(np.fft.ifft(np.fft.fft(A).conj() * np.fft.fft(B)))
    s_star = int(np.argmax(corr))
    return s_star


def shift_fields(fields, s):
    """Roll all fields by s grid points in the x-direction."""
    return [np.roll(f, s, axis=0) for f in fields]


def recurrence_distance(fields_t, fields_t_tau):
    """||T_{s*} X(t+tau) - X(t)|| / ||X(t)|| minimised over integer shifts."""
    Nx, Nz = fields_t[0].shape
    n_fields = len(fields_t)
    U_t     = flatten(fields_t)
    U_t_tau = flatten(fields_t_tau)

    s_star = optimal_shift(U_t, U_t_tau, Nx, n_fields, Nz)
    U_shifted = flatten(shift_fields(fields_t_tau, s_star))

    norm_t = np.linalg.norm(U_t)
    g = np.linalg.norm(U_shifted - U_t) / norm_t if norm_t > 0 else 1.0
    return g, s_star


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--opath',     required=True,
                        help='Directory containing SPECTER output files')
    parser.add_argument('--start-idx', type=int, default=1)
    parser.add_argument('--end-idx',   type=int, required=True,
                        help='Last snapshot index to load')
    parser.add_argument('--Nx',   type=int, required=True)
    parser.add_argument('--Nz',   type=int, required=True)
    parser.add_argument('--ra',   type=float, required=True)
    parser.add_argument('--pr',   type=float, default=1.0)
    parser.add_argument('--dt',   type=float, required=True,
                        help='SPECTER time step (free-fall units)')
    parser.add_argument('--ostep', type=int, required=True,
                        help='Steps between saved snapshots (tstep in parameter.inp)')
    parser.add_argument('--tau-min',   type=float, default=2.0)
    parser.add_argument('--tau-max',   type=float, default=50.0)
    parser.add_argument('--threshold', type=float, default=0.05)
    args = parser.parse_args()

    solver = make_solver(args.ra, args.pr, args.Nx, args.Nz, args.dt)

    dt_save = args.dt * args.ostep  # physical time between snapshots
    tau_steps_min = max(1, round(args.tau_min / dt_save))
    tau_steps_max = round(args.tau_max / dt_save)

    indices = list(range(args.start_idx, args.end_idx + 1))
    print(f'Loading {len(indices)} snapshots from {args.opath}...')
    snapshots = [solver.load_fields(args.opath, idx) for idx in indices]
    print('Done.')

    hits = []
    n = len(snapshots)
    for i in range(n):
        for dtau in range(tau_steps_min, min(tau_steps_max + 1, n - i)):
            g, s = recurrence_distance(snapshots[i], snapshots[i + dtau])
            if g < args.threshold:
                t   = (args.start_idx + i) * dt_save
                tau = dtau * dt_save
                hits.append((t, tau, s, g))

    hits.sort(key=lambda x: x[3])
    print(f'\nFound {len(hits)} near-recurrences (G < {args.threshold}):')
    print(f'{"t":>10}  {"tau":>10}  {"s (px)":>8}  {"G":>10}')
    for t, tau, s, g in hits[:20]:
        print(f'{t:10.3f}  {tau:10.4f}  {s:8d}  {g:10.4e}')


if __name__ == '__main__':
    main()
