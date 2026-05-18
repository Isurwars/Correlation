# Correlation 3.0.0 Roadmap

This document outlines the strategic direction and planned features for the `3.0.0` milestone of `Correlation`. This release aims to transform the software into a state-of-the-art engine for large-scale structural analysis, focusing on high-performance computing, advanced 3D spatial analysis, and seamless integration into modern scientific ecosystems.

## 🚀 Best Options to Attack Next (Current Priorities)

Based on recent progress, the following features are the primary focus for the immediate next development sprints to achieve the 3.0.0 vision:

1. **Memory-Mapped I/O** (Core): Crucial for lazy loading massive trajectories (e.g., GROMACS, LAMMPS) without exceeding system memory limits.
2. **Analysis Comparison** (UI): Quality-of-life feature allowing users to overlay and compare plots from multiple frames or datasets directly in the GUI.
3. **Refined Python API** (Integration): NumPy-native integration to allow seamless, high-performance data exchange between the C++ engine and Python environments.
4. **Extended XYZ** (Integration): Support for the Extended XYZ file format to handle custom metadata and per-atom properties.

## 1. Core (HPC & Architecture)
Modernizing the core engine for massive datasets and pushing the limits of computational performance.

- [x] **Explicit SIMD Vectorization**: AVX2/AVX-512 optimization for core calculation loops.
- [x] **Plugin/Factory Pattern**: Decoupled architecture for readers and calculators.
- [x] **Spatial Partitioning (Cell-Lists)**: Implement $O(N)$ neighbor search to replace the current $O(N^2)$ distance tensor approach for large systems.
- [x] **Task-Based Parallelism**: Advanced TBB task graph implementation for better multi-core scaling.
- [x] **Performance Benchmarking Suite**: Google Benchmark framework integration with comprehensive benchmarks for calculators and analyzers to monitor and guide optimization.
- [ ] **Memory-Mapped I/O**: Efficient loading of multi-gigabyte trajectories.
- [ ] **GPU Acceleration**: Experimental CUDA/HIP support for heavy calculations like S(Q).

## 2. Analysis (Distributions & Tools)
Enhancing the suite of structural and dynamic analysis tools.

- [x] **Mean Squared Displacement (MSD)**: Standard diffusion coefficient ($D$) computation.
- [x] **Common Neighbor Analysis (CNA)**: Classification into localized crystallographic environments.
- [x] **Steinhardt Bond-Orientational Parameters ($Q_4, Q_6, W_6$)**: Gold standard for liquid vs. solid structure.
- [x] **Hydrogen Bond Analysis**: Geometric criteria and autocorrelation functions for H-bond lifetimes.
- [x] **Spatial Distribution Functions (SDF)**: 3D probability density mapping. The initial focus will be on exporting 3D density maps (e.g., `.cube` files) for external visualization.
- [x] **Cluster Analysis**: Connectivity-based clustering algorithms to identify molecular fragments or crystalline grains.

## 3. Integration (Interoperability & Python)
Expanding support to attract users from various computational chemistry and physics communities.

- [x] **VASP (`POSCAR`, `CONTCAR`, `XDATCAR`)**: Essential for the inorganic materials and solid-state physics communities.
- [x] **GROMACS (`.gro`, `.xtc`, `.trr`)**: Opens `Correlation` to the soft-matter and biophysics communities.
- [x] **Protein Data Bank (`.pdb`)**: The universal standard for biological and organic structures.
- [x] **Python Bindings**: Direct C++ core integration for Jupyter and automated pipelines.
- [x] **Quantum ESPRESSO**: High priority support for *ab initio* MD output.
- [x] **CP2K Support**: Support for this popular molecular dynamics format.
- [ ] **Extended XYZ**: Improved support for custom metadata and per-atom properties.
- [ ] **Refined Python API**: NumPy-native integration for seamless data exchange.

## 4. UI (User Experience & Workflows)
Improving developer and user workflows through the graphical interface.

- [x] **Automated Release Pipelines**: CI/CD for cross-platform binaries and installers.
- [x] **Static Plot Management**: In-memory SVG generation preventing disk-bound delays, with a native file-saving dialog built into the user interface.
- [ ] **Analysis Comparison**: Built-in functionality to overlay and compare results from multiple runs statically.
- [ ] **WASM Web Interface**: Explore running the core engine in the browser for zero-install previews.