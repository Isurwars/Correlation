# ADR 0004: Pure Python Ecosystem Adapters (ASE & Pymatgen) over Pybind11

## Status
Accepted

## Context
In materials science and molecular dynamics, the Atomic Simulation Environment (ASE) and Pymatgen are de facto standard Python libraries for atomic data structures. Building C++ pybind11 modules with hard dependencies on ASE or Pymatgen introduces severe ABI coupling, requires Python packages to be available at C++ compile time, complicates cross-compilation (wheels for manylinux, macOS, Windows), and makes packaging fragile.

## Decision
1. **Separation of Layers**:
   - Keep the C++ pybind11 extension module (`_correlation`) minimal, fast, and purely dependent on NumPy buffer protocol and standard types.
   - Implement ecosystem conversion logic in pure Python (`python/correlation/adapters.py`).
2. **Lazy Dynamic Feature Detection**:
   - Detect ASE and Pymatgen at runtime when `correlation` is imported. If either is missing, core Correlation functionality remains 100% operational without error.
3. **Bidirectional API**:
   - Provide explicit functional APIs: `from_ase`, `to_ase`, `from_ase_trajectory`, `to_ase_trajectory`, `from_pymatgen`, `to_pymatgen`.
4. **Convenience Monkey-Patching**:
   - Safely monkey-patch `.to_correlation()` onto `ase.Atoms` and `pymatgen.core.Structure` when present, and inject `.to_ase()` and `.to_pymatgen()` onto `correlation.Cell`.

## Consequences
### Positive
- Zero compile-time or link-time dependency on ASE, Pymatgen, or third-party Python packages in C++ wheels.
- Fast wheel compilation and simplified CI build pipelines (`scikit-build-core`).
- Seamless developer user experience with standard duck-typing and method calls (`atoms.to_correlation()`).

### Negative / Trade-offs
- Coordinate conversions pass through NumPy arrays, introducing a microscopic copy overhead compared to theoretical in-memory pointer sharing, though completely negligible compared to subsequent neighbor searches and pair calculations ($< 0.1\%$ runtime).
