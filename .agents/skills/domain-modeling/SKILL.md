---
name: domain-modeling
description: Establish, refine, and enforce ubiquitous domain language and scientific terminology in docs/CONTEXT.md. Use when clarifying scientific concepts, adding new atomic properties, or reconciling naming ambiguities.
---

# Domain Modeling

Builds and maintains a shared ubiquitous language for `Correlation` to eliminate ambiguity in atomic structural analysis, condensed matter physics, and parallel computing.

## Phase 1: Identify Domain Concepts

When designing or refactoring algorithms, identify terms across the four primary domains:

| Domain | Core Concepts | Common Confusions / Anti-Patterns |
| :--- | :--- | :--- |
| **Atomic Trajectories** | Frame, Timestep, Simulation Cell, Periodic Boundary Conditions (PBC), Minimum Image Convention (MIC) | Confusing fractional coordinates with cartesian coordinates; ignoring triclinic box tilting. |
| **Pair Correlations** | Radial Distribution Function ($g(r)$), Pair Distribution Function ($G(r)$), Structure Factor ($S(q)$) | Conflating $g(r)$ (density-normalized) with $G(r)$ (reduced/neutron-weighted) or partial vs. total correlations. |
| **Angular & Topology** | Planar Angle Distribution (PAD), Bond Angle Distribution (BAD), Ring Statistics (Primitive vs. King's rings) | Ambiguous cut-off distance definitions for nearest neighbor graph construction. |
| **Dynamical Properties** | Mean Squared Displacement (MSD), Velocity Autocorrelation Function (VACF), Vibrational Density of States (VDOS) | Confusing instantaneous velocities with finite-difference displacements; window averaging artifacts. |

## Phase 2: Updating `docs/CONTEXT.md`

Whenever a new domain entity or mathematical convention is introduced, update or create `docs/CONTEXT.md`:

```markdown
# Domain Context & Ubiquitous Language

## [Entity / Metric Name]
- **Definition:** Exact physical/mathematical meaning.
- **Mathematical Formulation:** Formula (LaTeX style).
- **Coordinate Conventions:** Cartesian (Ångströms), fractional, reciprocal ($1/\text{Å}$).
- **Normalizations:** Volume density $\rho_0$, number of pairs, Faber-Ziman vs. Ashcroft-Langreth weights.
- **Code Representations:** Associated C++ structs, typedefs, and Slint properties.
```

## Phase 3: Enforce Ubiquitous Language in Code

1. Verify that C++ class, struct, and variable names reflect the exact physical definitions (e.g. `RadialDistributionFunction`, not `DistCalc`).
2. Ensure comments and Doxygen documentation use standard physical units (Å, ps, eV, THz, radians/degrees).
