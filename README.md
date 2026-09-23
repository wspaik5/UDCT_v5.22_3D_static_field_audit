# UDCT v5.22 - 3D Static-Field Audit

Methods-correction and reproducibility note. It does not claim a Bullet Cluster fit, a new MOND mechanism, a relativistic lensing solution, or evidence that UDCT outperforms MOND or dark-matter models.

**Author:** Won Shik Paik  
**Paper:** *UDCT v5.22, Audited Three-Dimensional Static-Field Note*  
**Date:** 22 September 2026

## Version genealogy

UDCT v5.22 supplements and does not replace the locked v5.4 RAR, v5.5 BTFR, or v5.6 dwarf records. This repository contains only the v5.22 three-dimensional static-field audit.

## Scope and locked definitions

The frozen acceleration scale is

```text
xi0 = 1.082401e-10 m s^-2
```

The reported statistic is a sampled contrast, not a global argmax:

```text
C_samp = [Sigma_eff(-d, 0) + Sigma_eff(+d, 0)] / [2 Sigma_eff(0, 0)]
```

The scripts also report the full-map argmax separately.

The locked toy geometries are:

| Geometry | Total mass | Gas component | Each galaxy component | Galaxy coordinates |
|---|---:|---|---|---:|
| Baseline | 1e14 Msun | fraction 0.85, Plummer scale 300 kpc | fraction 0.075, Plummer scale 50 kpc | +/-500 kpc |
| Favorable | 3e14 Msun | fraction 0.85, Plummer scale 500 kpc | fraction 0.075, Plummer scale 70 kpc | +/-300 kpc |

No additional geometry is included.

## Field equations

QUMOND uses

```text
Laplacian(Phi) = div[nu(|grad(Phi_N)| / xi0) grad(Phi_N)]
nu(x) = sqrt(1 + 1/x)
Laplacian(Phi_N) = 4 pi G rho_b
```

AQUAL uses

```text
div[mu(|grad(Phi)| / xi0) grad(Phi)] = 4 pi G rho_b
mu(y) = 2y / (sqrt(1 + 4y^2) + 1)
```

The numerical constants, Plummer normalization, Picard update, and FFT conventions are preserved in the source files.

## Environment

From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

## Exact run commands

Recompute the Newtonian and QUMOND audit matrix:

```bash
python3 src/v522_audit_recompute.py \
  --out data/v522_qumond_newton_audit.csv \
  --json /tmp/v522_qumond_newton_audit.json
```

Recompute the baseline AQUAL audit:

```bash
python3 src/aqual_3d_compare.py \
  --N 96 \
  --box 3000 \
  --floor 1e-5 \
  --geometry baseline \
  --out data/v522_aqual_baseline_N96_B3000.json
```

Recompute the favorable AQUAL audit:

```bash
python3 src/aqual_3d_compare.py \
  --N 96 \
  --box 3000 \
  --floor 1e-5 \
  --geometry favorable \
  --out data/v522_aqual_favorable_N96_B3000.json
```

The AQUAL calculations are materially slower than the FFT audit and print Picard-iteration progress to standard output.

## Locked outputs

For the favorable toy at 100 kpc Gaussian smoothing:

| Operator | Locked C_samp |
|---|---:|
| Newtonian baryons | 1.328 |
| QUMOND | 1.210 |
| AQUAL | 1.160 |

At 150 kpc smoothing, all three sampled contrasts are below 1.

The machine-readable locked values are in `expected/table2_favorable.json`. The QUMOND audit matrix is in `expected/table4_qumond_matrix.json`.

## Grid and boundary-condition disclosure

- Newtonian and QUMOND headline audit: periodic FFT Poisson calculation, `N=128`, box size `3000 kpc`, `Delta x=23.4375 kpc` (reported as `23.44 kpc`).
- AQUAL headline audit: Picard iteration plus conjugate gradients, isolated monopole boundary condition, `N=96`, box size `3000 kpc`, `Delta x=31.25 kpc`.
- The grids are unequal. Their classifications may be compared, but the difference between QUMOND `1.210` and AQUAL `1.160` is not a precision ranking of the operators.
- This package does not claim continuum convergence.

## Interpretation lock

The Newtonian favorable contrast is larger than the QUMOND and AQUAL contrasts. The toy result is therefore dominated by compactness, projection, and smoothing; it is not evidence for a uniquely nonlinear modified-gravity displacement.

## Reproduction acceptance rule

A third party should compare the recomputed favorable 100 kpc `C_samp` values with `expected/table2_favorable.json`. If any reprinted `C_samp` differs by more than `0.02` from its locked value, the package has failed reproduction. `C_samp` and the global argmax must remain separate diagnostics.

## Related links

- v5.4 RAR: https://github.com/wspaik5/UDCT_v5.4_SPARC_RAR_Reproducibility
- v5.5 BTFR: https://github.com/wspaik5/UDCT_v5.5_SPARC_Q1_BTFR_Reproducibility

## Exclusions

This repository contains no SPARC rotmod files, real Bullet Cluster maps, or v2.x unification material.

## License

MIT License. See `LICENSE`.
