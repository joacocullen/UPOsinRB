import numpy as np
import spooky as sp
from spooky.solvers import SPECTER


def make_solver(ra, pr, Nx, Nz, dt, nprocs=1, Lx=2*np.pi, Lz=np.pi, gamma=1.):
    grid = sp.Grid2D_semi(Lx=Lx, Lz=Lz, Nx=Nx, Nz=Nz, dt=dt)
    return SPECTER(grid, nprocs=nprocs, ra=ra, pr=pr, gamma=gamma,
                   solver='BOUSS', ftypes=['vx', 'vz', 'th'], precision='double')
