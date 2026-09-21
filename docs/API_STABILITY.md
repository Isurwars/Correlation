# Correlation v4.0.0 API Stability & Compatibility Contract

*Last Updated: September 2026*  
*Version: 4.0.0*  
*Status: **FROZEN***

---

## 1. Executive Summary & Semantic Versioning Commitment

`Correlation` v4.0.0 formalizes a strict **Public API Freeze** governed by **Semantic Versioning 2.0.0**. Starting with v4.0.0, the public C++, Python, and WebAssembly interfaces are guaranteed backwards-compatible across all subsequent `v4.x.y` releases.

### 1.1 SemVer 2.0.0 Rules
- **Patch Releases (`v4.0.x`)**: Bug fixes, performance optimizations, documentation updates, and compiler warnings resolution. Zero API additions or behavioral alterations.
- **Minor Releases (`v4.x.0`)**: Backwards-compatible feature additions (new calculators, supplementary reader formats). Existing interfaces remain immutable.
- **Major Releases (`v5.0.0`)**: Breaking architectural changes (e.g. AoS to SoA memory restructuring for `Cell`). Prohibited within the 4.x lifecycle.

### 1.2 Deprecation & Migration Policy
1. No public class, struct, function, method, or binding will be deleted or altered incompatibly within the `4.x` lifecycle.
2. If an API is superseded, it will be marked `[[deprecated("...")]]` in C++ and `warnings.warn(..., DeprecationWarning)` in Python.
3. Deprecated APIs will continue functioning until the `5.0.0` major release.

---

## 2. C++23 Public API Surface Freeze

### 2.1 Core Domain Models
- **`correlation::core::Cell`**: Atom container, periodic boundary conditions, volume, minimum-image conventions, bulk loading (`addAtom`, `reserveAtoms`).
- **`correlation::core::Atom`**: Position, velocity (storing force vectors in MLIP trajectories), chemical symbol, mass, covalent radius.
- **`correlation::core::Trajectory`**: Frame indexing, lazy memory-mapped file offsets, time steps.
- **`correlation::core::MappedFile`**: POSIX `mmap` RAII abstraction for zero-copy file parsing.
- **`correlation::math::Vector3<real_t>` & `correlation::math::Matrix3<real_t>`**: Geometry primitives.
- **`real_t` Type Alias**: Configurable precision float alias (`float` in single-precision, `double` in standard).

### 2.2 Readers & Trajectory Parsers
- **`correlation::readers::BaseReader`**: Pure virtual interface (`getName`, `getExtensions`, `isTrajectory`, `readStructure`, `readTrajectory`).
- **`correlation::readers::ReaderFactory`**: Auto-registration mechanism and content-based format sniffing.
- **Supported Reader Suite**: CAR, CELL, CIF, ARC, LAMMPS dump, CASTEP MD, DMol3 OUTMOL, VASP POSCAR/XDATCAR, GROMACS, PDB, XYZ/EXTXYZ, Quantum ESPRESSO, CP2K, ORCA, GPAW, ABINIT, DFTB+, MACE, CHGNet, GAP, and NequIP.

### 2.3 Calculators & Analysis Suite
- **Radial & Angular Distributions**: `RDFCalculator`, `PADCalculator`, `CNCalculator`, `DistanceCalculator`.
- **Dynamics & Scattering**: `MSDCalculator`, `VACFCalculator`, `VDOSCalculator`, `XRDCalculator`, `StructureFactorCalculator`.
- **Structural Motifs & Order**: `SteinhardtCalculator`, `CNACalculator`, `VoronoiCalculator`, `HBondCalculator`, `SDFCalculator`, `ClusterCalculator`, `LocalEntropyCalculator`, `ChiralityCalculator`, `HyperuniformityCalculator`.
- **ML & Electronic Structure**: `MLIPCalculator`, `TorchGNNModel`, `PeriodicGraphBuilder`, `GraphDescriptors`, `TDOSCalculator`, `StructuralElectronicCorrelation`.
- **Unified Engine**: `DistributionFunctions` facade and `CalculatorFactory`.

### 2.4 C++ Quality Standards Gate
- Strict compilation: `-Wall -Wextra -Wpedantic -Werror` on GCC and Clang.
- Function Cognitive Complexity $\le 25$ (`readability-function-cognitive-complexity`).
- Zero `NOLINT` or inline suppression policy.
- Full Doxygen comments on all public interfaces.

---

## 3. Python Public API Surface Freeze

### 3.1 Module Namespace (`correlation`)
- `correlation.Cell`, `correlation.Trajectory`, `correlation.DistributionFunctions`, `correlation.Histogram`.
- `correlation.get_registered_calculators()`, `correlation.get_registered_readers()`, `correlation.get_registered_writers()`.

### 3.2 Zero-Copy Buffer Protocol
- **Read Access**: `Cell.positions` and `Cell.velocities` expose zero-copy strided NumPy views (`np.ndarray(shape=(N, 3), dtype=real_t)`).
- **Write Access**: `Cell.from_arrays(positions, symbols)` provides bulk SIMD loading directly from contiguous NumPy arrays, avoiding Python loop overhead.

### 3.3 Materials Science Ecosystem Bridges (`correlation.adapters`)
- **ASE**: `from_ase`, `to_ase`, `from_ase_trajectory`, `to_ase_trajectory`.
- **Pymatgen**: `from_pymatgen`, `to_pymatgen`.
- **PyG / Torch Geometric**: `to_torch_geometric`, `to_atom_graphs`.

---

## 4. WebAssembly & Web Worker API Freeze

### 4.1 WASM Module (`correlation_wasm.js` / `.wasm`)
- Embind exports: `createCorrelationModule()`, `readFromBuffer(strData, filename)`, `Cell`, `Trajectory`, `DistributionFunctions`.
- SIMD128 vectorization enabled by default.

### 4.2 Web Worker Message Protocol (`worker.js`)
- Commands: `LOAD_FILE` (`{ text, filename }`), `RUN_ANALYSIS` (`{ rMax, binWidth, calcPad }`).
- Responses: `READY`, `FILE_LOADED`, `STATUS`, `ANALYSIS_COMPLETE` (`{ results }`), `ERROR`.
