# Correlation 3.0.0 Roadmap

This document outlines the strategic direction and planned features for the `3.0.0` milestone of `Correlation`. This release aims to transform the software into a state-of-the-art engine for large-scale structural analysis, focusing on high-performance computing, advanced 3D spatial analysis, and seamless integration into modern scientific ecosystems.

## 1. Supported Software & File Formats
Expanding support to attract users from biomolecular, solid-state, and soft-matter communities.

- [x] **VASP (`POSCAR`, `CONTCAR`, `XDATCAR`)**: Essential for the inorganic materials and solid-state physics communities.
- [x] **GROMACS (`.gro`, `.xtc`, `.trr`)**: Opens `Correlation` to the soft-matter and biophysics communities.
- [x] **Protein Data Bank (`.pdb`)**: The universal standard for biological and organic structures.
- [ ] **Quantum ESPRESSO & CP2K**: Support for popular *ab initio* molecular dynamics formats.
- [ ] **Extended XYZ**: Improved support for custom metadata and per-atom properties.

## 2. New Distribution Functions and Analyses
Enhancing the suite of structural and dynamic tools.

- [x] **Mean Squared Displacement (MSD)**: Standard diffusion coefficient ($D$) computation.
- [x] **Common Neighbor Analysis (CNA)**: Classification into localized crystallographic environments.
- [x] **Steinhardt Bond-Orientational Parameters ($Q_4, Q_6, W_6$)**: Gold standard for liquid vs. solid structure.
- [x] **Hydrogen Bond Analysis**: Geometric criteria and autocorrelation functions for H-bond lifetimes.
- [ ] **Spatial Distribution Functions (SDF)**: 3D probability density mapping for hydration shells and solvation structures.
- [ ] **Cluster Analysis**: Connectivity-based clustering algorithms to identify molecular fragments or crystalline grains.

## 3. High-Performance Computing (HPC) & Architecture
Modernizing the core engine for massive datasets.

- [x] **Explicit SIMD Vectorization**: AVX2/AVX-512 optimization for core calculation loops.
- [x] **Plugin/Factory Pattern**: Decoupled architecture for readers and calculators.
- [ ] **Spatial Partitioning (Cell-Lists)**: Implement $O(N)$ neighbor search to replace the current $O(N^2)$ distance tensor approach for large systems.
- [ ] **Task-Based Parallelism**: Advanced TBB task graph implementation for better multi-core scaling.
- [ ] **Memory-Mapped I/O**: Efficient loading of multi-gigabyte trajectories.

## 4. Ecosystem & User Experience
Improving integration and developer/user workflows.

- [x] **Python Bindings**: Direct C++ core integration for Jupyter and automated pipelines.
- [x] **Automated Release Pipelines**: CI/CD for cross-platform binaries and installers.
- [ ] **Refined Python API**: NumPy-native integration for seamless data exchange.
- [ ] **Interactive Plotting**: Enhanced GUI with interactive zoom/pan for distribution functions.
- [ ] **Analysis Comparison**: Built-in functionality to overlay and compare results from multiple runs.
- [ ] **WASM Web Interface**: Explore running the core engine in the browser for zero-install previews.