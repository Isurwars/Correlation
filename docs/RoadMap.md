# Correlation 3.0.0 Roadmap

This document outlines the strategic direction and planned features for the `3.0.0` milestone of `Correlation`. This release aims to transform the software into a state-of-the-art engine for large-scale structural analysis, focusing on high-performance computing, advanced 3D spatial analysis, and seamless integration into modern scientific ecosystems.

> **Current version:** 2.9.0 — All originally planned 3.0.0 features have been implemented. The remaining work focuses on hardening, polish, and new stretch goals.

## 🚀 Stretch Goals & Next Sprints

All primary 3.0.0 features are complete. The following represent the next development focus:

1. **Analysis Comparison** (UI): Quality-of-life feature allowing users to overlay and compare plots from multiple frames or datasets directly in the GUI.
2. **GPU Acceleration Hardening**: Broader CUDA kernel coverage beyond S(Q), performance benchmarks on real hardware, and automated CI testing on GPU runners.
3. **WASM Web Interface Polish**: Expand the Emscripten bindings to support more calculators and add a standalone browser demo.
4. **Refined Python API — Documentation & Examples**: NumPy-native integration is implemented; publish Jupyter notebook examples and add type stubs.

## 1. Core (HPC & Architecture)
Modernizing the core engine for massive datasets and pushing the limits of computational performance.

- [x] **Explicit SIMD Vectorization**: AVX2/AVX-512 optimization for core calculation loops.
- [x] **Plugin/Factory Pattern**: Decoupled architecture for readers and calculators.
- [x] **Spatial Partitioning (Cell-Lists)**: Implement $O(N)$ neighbor search to replace the current $O(N^2)$ distance tensor approach for large systems.
- [x] **Task-Based Parallelism**: Advanced TBB task graph implementation for better multi-core scaling.
- [x] **Performance Benchmarking Suite**: Google Benchmark framework integration with comprehensive benchmarks for calculators and analyzers to monitor and guide optimization.
- [x] **Memory-Mapped I/O**: Efficient lazy loading of multi-gigabyte trajectories via `MappedFile`. All major readers (XYZ, LAMMPS, GROMACS, XDATCAR) use memory-mapped frames.
- [x] **GPU Acceleration**: CUDA-accelerated S(Q) via `GPUSQCalculator` (`-DBUILD_WITH_CUDA=ON`), with automatic CPU fallback when no compatible device is detected.

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
- [x] **Extended XYZ**: Full support for standard and Extended XYZ (`comment`-line key=value metadata, custom per-atom columns, lattice vectors, and multi-frame trajectories).
- [x] **Refined Python API**: NumPy-native zero-copy integration via `.positions` and `.velocities` buffer-protocol properties on `Cell`. Deprecation warnings guide users away from the old copy-based API.

## 4. UI (User Experience & Workflows)
Improving developer and user workflows through the graphical interface.

- [x] **Automated Release Pipelines**: CI/CD for cross-platform binaries and installers.
- [x] **Static Plot Management**: In-memory SVG generation preventing disk-bound delays, with a native file-saving dialog built into the user interface.
- [x] **WASM Web Interface**: Core engine runs in the browser via Emscripten (`-DBUILD_WASM=ON`). Exposes `Cell`, `Trajectory`, `DistributionFunctions`, RDF, and PAD to JavaScript via embind.
- [ ] **Analysis Comparison**: Built-in functionality to overlay and compare results from multiple runs statically.