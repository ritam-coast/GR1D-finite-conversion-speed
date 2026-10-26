# GR1D with TOV initial data, piecewise-polytrope EoS, and finite-speed phase-transition conversion

This branch builds on upstream GR1D (https://github.com/evanoconnor/GR1D, commit `20837489`, branch `main` of this repository).
It is the code used for the Letter "Quark deconfinement as a natural trigger of Type-I critical collapse:
black holes from neutron stars deemed stable" (Borthakur, Saha, Mallick).

## Additions

- **TOV initial data** (`src/TOVstar.F90`), read from a profile file, with optional radial velocity
  perturbations (uniform `U=V0` or homologous `U=V0 r/R`). The metric is re-solved from the constraints.
- **Piecewise-polytrope cold EoS + Gamma_th thermal part** (`eoskey=4` with a `pp_eos.in` file),
  P = P_cold(rho) + (Gamma_th - 1) rho [eps - eps_cold(rho)]. Up to 128 pieces.
- **Robust atmosphere**, primitive-recovery fallback at stellar surfaces, and a CFL guard.
- **Diagnostics**: central density and lapse, unbound mass, baryon mass, and quark- and mixed-phase masses
  (`out/tovdiag.dat`, `out/tovdiag2.dat`); a BH is declared when the central lapse falls below a threshold (`out/BH`).
- **Finite-speed conversion (optional)**: if `conv.in` is present in the run directory, the cold EoS in each
  zone is a blend of the hybrid table `pp_eos.in` and a hadronic table `pp_eos_had.in`,
  P_cold = (1 - f) P_had + f P_hyb (eps_cold likewise), with converted fraction f(r,t):
  - mode 1, uniform relaxation: f = 1 - exp(-t/tau), `conv.in` = `1 tau_ms width_km`;
  - mode 2, outward front: f = 1 behind r_f = v_f t, linear ramp of width w, `conv.in` = `2 v_f_in_c width_km`.

  - mode 3, Lagrangian relaxation: f is carried by the fluid in the advected scalar `ye` (conserved `D*ye`),
    starts at 0, and after every time step grows as df/dt = (1 - f)/tau wherever rho > rho_conv;
    `conv.in` = `3 tau_ms width_km` on the first line and `rho_conv` [g/cm^3] on the second line;
  - mode 4, Lagrangian front: as mode 3, but f = max(f, ramp) behind the outward front r_f = v_f t
    for matter with rho > rho_conv; `conv.in` = `4 v_f_in_c width_km` and `rho_conv` on the second line.

  Conversion is at fixed rest-mass density and specific internal energy, so the released energy becomes heat.
  In modes 3 and 4 conversion is irreversible and attached to the matter; matter compressed above rho_conv
  later in the evolution starts unconverted and converts at the same rate. Modes 1 and 2 prescribe f(r,t)
  directly (Eulerian, not advected) and are kept for comparison.
  Without `conv.in` the code behaves exactly as without this option.

## Input files in the run directory

- `parameters` (standard GR1D, `initial_data = "TOV"`, `eoskey = 4`)
- `tov_params.in`: profile file, kick amplitude, kick time, second kick amplitude, diagnostic cadence,
  BH lapse threshold, outer boundary (code units), perturbation type (0 uniform, 2 homologous)
- `pp_eos.in`: first line `n_pieces Gamma_th`, then `rho_lo K Gamma a` per piece (cgs, P/c^2 = K rho^Gamma)
- optional `conv.in` and `pp_eos_had.in` (same format as `pp_eos.in`)

## Build

```
cp make.inc.example make.inc   # or adapt your own make.inc (HDF5 Fortran, LAPACK)
make
```
