# UPOs in Rayleigh-Bénard convection

Code and data accompanying:

> **Characterizing periodic orbits in two-dimensional Rayleigh-Bénard flows**  
> J. Cullen, M. Y. Vinograd, P. Clark Di Leoni (2026)

We compute and characterize a steady state (ST) and three families of periodic
orbits (PO1, PO2, PO3) in 2D Rayleigh-Bénard convection at Pr = 1, using
Newton-GMRes-Hookstep continuation from the full PDE.

## Dependencies

- **SPECTER** — pseudospectral DNS solver. Must be compiled and available as `./BOUSS`.  
  Repository: https://github.com/PClarkDiLeoni/SPECTER
- **spooky** — Newton-GMRes-Hookstep library:

```bash
conda env create -f environment.yml
conda activate uposinrb
```

## Physical setup

| Parameter | Value |
|-----------|-------|
| Pr | 1 |
| Domain [Lx, Lz] | [2π, π] |
| Boundary conditions | Periodic in x, no-slip at z = 0, π |
| Ra range | 10⁵ – 2×10⁷ |
| Resolutions [Nx, Nz] | [256,103], [256,231], [512,512] |

Time in free-fall units t' = h/U, U = √(gαΔT h).

## Repository structure

```
convergence/
  io_utils.py           — solver factory
  find_recurrences.py   — recurrence analysis to seed the Newton solver
  converge_upo.py       — Newton-GMRes-Hookstep convergence script
  params/
    ST_Ra8e5/           — solver.yaml + newton.yaml
    PO1_Ra1.25e6/
    PO2_Ra5e6/
    PO3_Ra9e6/
floquet/
  floquet_spectrum.py   — Floquet multipliers via Arnoldi iteration
data/
  invariant_solutions/  — converged fields as SPECTER binaries (one dir per family)
```

## Running a convergence

The newton.yaml for each solution already points `input_dir` at the provided
converged state. From inside `convergence/`:

```bash
python converge_upo.py \
    --solver params/PO1_Ra1.25e6/solver.yaml \
    --newton params/PO1_Ra1.25e6/newton.yaml
```

Convergence reports are written to `output/` and `reports/` as configured in
`newton.yaml`. To restart from Newton iteration N, set `restart_iN: N`.

**Finding new initial guesses from a DNS trajectory:**

```bash
python find_recurrences.py \
    --opath /path/to/dns/output \
    --start-idx 1 --end-idx 500 \
    --Nx 256 --Nz 103 --ra 1.25e6 --pr 1 \
    --dt 5e-4 --ostep 200 \
    --tau-min 5 --tau-max 20 --threshold 0.05
```

Copy the snapshot at the best (t, τ) into a new `input_dir` and set `T` in
`newton.yaml`.

## Floquet analysis

```bash
cd floquet
python floquet_spectrum.py \
    --solver ../convergence/params/PO1_Ra1.25e6/solver.yaml \
    --newton ../convergence/params/PO1_Ra1.25e6/newton.yaml \
    --state-dir ../data/invariant_solutions/PO1 \
    --n-eigs 100 \
    --output PO1_Ra1.25e6_floquet.npz
```

Output `.npz` contains `multipliers`, `floquet_exponents`, `modes`, `Ra`, `T`.

## Converged solutions

See `data/invariant_solutions/README.md`.

## Citation

```bibtex
@article{cullen2026upos,
  title   = {Characterizing periodic orbits in two-dimensional {Rayleigh-B\'enard} flows},
  author  = {Cullen, Joaqu\'in and Vinograd, Melisa Y. and {Clark Di Leoni}, Patricio},
  journal = {},
  year    = {2026},
}
```
