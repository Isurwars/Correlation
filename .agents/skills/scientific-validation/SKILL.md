---
name: scientific-validation
description: Enforce physical invariants, sum rules, literature benchmark calibration, and single-precision numerical stability for atomistic structural and dynamic calculators. Use when adding or calibrating calculators, order parameters, or phase markers.
---

# Scientific Validation & Physical Invariants

## 1. Mandatory Physical Sanity Checks
Verify every calculator against analytical limits and canonical synthetic structures in Google Test:
- **RDF $g(r)$:** $g(r) \to 1$ at large $r$; ideal crystal peaks match exact lattice shell distances and coordination numbers.
- **Hyperuniformity ($\sigma_N^2(R), \chi_H(R), S(Q)$):**
  - Ideal lattice: $\sigma_N^2(R) \propto R^\beta$ with $\beta \approx 2.0$, $\chi_H(R) \to 0$, Torquato-Stillinger $H = S(0^+)/S(q_{\text{peak}}) \le 10^{-3}$.
  - Uncorrelated Poisson / liquid: $\beta \approx 3.0$, $\chi_H(R_{\max}) \approx \rho k_B T \kappa_T > 0$.
- **Dynamics ($D_{\text{MSD}}, D_{\text{VACF}}, \tau, De$):**
  - Units must match across C++ core, CLI, Python bindings, and Slint tooltips (`Å²/fs` for $D$, `fs` for $\tau$, dimensionless $De = \tau / t_{\text{obs}}$).
  - Regime thresholds: $De < 0.1$ (Liquid-like), $0.1 \le De \le 10$ (Viscoelastic), $De > 10$ (Solid-like).

## 2. Literature Benchmark Calibration
Before setting or modifying phase-search reference markers (`PhaseSearchTypes.hpp`) or canonical preset bounds:
- Cross-check peer-reviewed literature (e.g., Errington & Debenedetti *Nature* 2001, Lechner & Dellago *JCP* 2008, Piaggi & Parrinello *JCP* 2017, Martelli et al. *PRB* 2017 / *PNAS* 2018, Torquato *Phys. Rep.* 2018).
- Document the exact citation and physical ranges in a comment above the marker array.
- Account for all-atom vs coarse-grained cells (e.g., `examples/a-Ice/ice.outmol` contains both O and H atoms).

## 3. Numerical Precision (`USE_SINGLE_PRECISION=ON`)
- `real_t` defaults to `float`. Use `double` internally for OLS regressions, log-log slopes, variance accumulators ($\langle N^2 \rangle - \langle N \rangle^2$), and trapezoidal integrals before casting outputs to `real_t`.
- Guard all divisions, `std::log`, and `std::sqrt` against zero/negative/non-finite inputs (`std::isfinite`).
