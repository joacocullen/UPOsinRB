# Invariant solutions

One converged state per solution family at the representative Ra used in the paper.

| Solution | Ra        | Period            |
| ST       | 8×10^5    |   —               |
| PO1      | 1.25×10^6 | 11.33853009569369 |
| PO2      | 5×10^6    | 32.95424449318564 |
| PO3      | 9×10^6    |  2.012760757483008|

Each subdirectory (`ST/`, `PO1/`, `PO2/`, `PO3/`) contains SPECTER binary
files `vx.00000.out`, `vz.00000.out`, `th.00000.out` (float64, physical space,
shape [Nx, Nz]). The `vy.00000.out` and `pr.00000.out` files are unused by
the Newton solver.
