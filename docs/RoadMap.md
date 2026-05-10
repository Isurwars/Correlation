# Correlation 3.0.0 Roadmap

This document outlines the strategic direction and planned features for the `3.0.0` release of `Correlation`. The goal of this release is to expand the software's reach by supporting more file formats, adding requested analysis methods, and improving the underlying architecture for better integration into modern scientific workflows.

## 1. Supported Software & File Formats
Expanding support to attract users from biomolecular, solid-state, and soft-matter communities.

- [x] **VASP (`POSCAR`, `CONTCAR`, `XDATCAR`)**: Essential for the inorganic materials, ceramics, and solid-state physics communities.
- [x] **GROMACS (`.gro`, `.xtc`, `.trr`)**: Opens `Correlation` to the soft-matter and biophysics communities.
- [x] **Protein Data Bank (`.pdb`)**: The universal standard for biological and organic structures.
- [ ] **Quantum ESPRESSO & CP2K**: Support or documentation for these popular *ab initio* molecular dynamics alternatives.

## 2. New Distribution Functions and Analyses
Enhancing the suite of structural and dynamic tools.

- [x] **Mean Squared Displacement (MSD)**: The standard way to compute diffusion coefficients ($D$) via the Einstein relation. Natural addition alongside the existing VACF.
- [x] **Common Neighbor Analysis (CNA) / Polyhedral Template Matching (PTM)**: Crucial for classifying atoms into localized crystallographic environments (FCC, BCC, HCP, Icosahedral) during freezing/melting or in metallic glasses.
- [x] **Steinhardt Bond-Orientational Parameters ($Q_4, Q_6, W_6$)**: The gold standard for distinguishing rotationally invariant local structures in liquids vs amorphous/crystalline solids.
- [ ] **Spatial Distribution Functions (SDF)**: Computing the 3D probability density of finding an atom around a central molecule (e.g., mapping hydration shells), moving beyond 1D $g(r)$.
- [x] **Hydrogen Bond Analysis**: Using geometric criteria (distance + angle cutoffs) to count hydrogen bonds, and dynamical H-bond autocorrelation functions to measure their lifetimes.

## 3. Code Architecture & "Quality of Life" Improvements
Modernizing the codebase and improving developer/user experience.

- [x] **Python Bindings (`pybind11` or `nanobind`)**: Wrap the C++ core into a Python module (e.g., `pip install correlation-analysis`) to allow researchers to bypass the UI and integrate the fast parallel engine directly into their Jupyter Notebooks and existing analysis pipelines (ASE, MDAnalysis).
- [x] **Plugin/Factory Pattern for Readers & Calculators**: Refactor the architecture so adding a new `Reader`, `Calculator` and `Writer` only requires dropping in a new class that automatically registers itself with the GUI, without modifying the core UI and `AppController`.
- [x] **Explicit SIMD Vectorization**: Explore using explicit SIMD (AVX2/AVX-512) for the inner loops of distance calculations to achieve even higher performance on top of the existing Intel TBB multithreading.
- [x] **Automated Release Pipelines**: Set up GitHub Actions to automatically compile binaries for Windows/Linux/macOS, generate installers, and update the AUR package upon a new GitHub Release.
|