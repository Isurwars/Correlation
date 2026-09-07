[![Version](https://img.shields.io/badge/version-4.0.0-green)](https://img.shields.io/badge/version-4.0.0-green) [![License](https://img.shields.io/badge/license-AGPLv3-brightgreen)](https://img.shields.io/badge/license-AGPLv3-brightgreen) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)

# Correlation: Contributing Guidelines

We welcome contributions to `Correlation` and strive to maintain rigorous standards of scientific integrity, code reproducibility, and computational performance. As a research-grade atomic structural analysis suite, we enforce strict architectural and quality gates across all contributions.

---

## Code of Conduct

All contributors are expected to uphold our [**Code of Conduct**](CODE_OF_CONDUCT.md). Please report unacceptable behavior to the project maintainers.

---

## Development Standards & Quality Gates

### 1. C++23 Language Standards
- The core and calculator libraries target **C++23** (`-std=c++23`).
- Prefer modern constructs: `std::span`, `std::ranges`, concepts, `std::expected`, `constexpr`, and structured bindings.
- **Strict RAII**: Avoid raw allocations (`new`/`delete`). Use smart pointers (`std::unique_ptr`, `std::shared_ptr`) and RAII memory containers (`DeviceBuffer<T>`, `std::vector`).
- **Data-Oriented Design**: Keep hot loops contiguous in memory (`alignas(64)` where applicable for cache-line alignment) and free of branch mispredictions.

### 2. Cognitive Complexity Gate ($\le 25$)
- All C++ functions must maintain Cognitive Complexity $\le 25$ as enforced by `clang-tidy` (`readability-function-cognitive-complexity`).
- Decompose monolithic routines into clean, single-purpose helper functions and parameter structs.
- **Zero Suppression Policy**: Prohibit `NOLINT`, `NOLINTNEXTLINE`, or inline static analysis suppressions. Resolve root causes.

### 3. Formatting & Static Analysis
- **`clang-format`**: Run `clang-format -i <file>` on all modified C++ files before submitting.
- **`clang-tidy`**: All code must compile cleanly against `compile_commands.json` with zero warnings under `-Wall -Wextra -Wpedantic -Werror`.

### 4. Doxygen Documentation
- All public and protected classes, structs, methods, and functions in header files (`include/**`) must have complete Doxygen blocks.
- Document `@brief`, `@param[in,out]`, `@return`, and `@throws` explicitly.

### 5. Slint UI MVVM Architecture
- **Declarative UI**: Slint files (`ui/*.slint`) contain strictly UI layout, Material Design 3 design tokens, and property bindings. Scientific calculations must never be implemented in Slint scripts.
- **Thread Safety**: Long-running calculations must run asynchronously on worker threads. Never block the Slint event loop. UI updates must be dispatched via `slint::invoke_from_event_loop()`.

### 6. Python Bindings & Ecosystem Adapters
- The C++ extension (`_correlation`) is built with pybind11 and provides NumPy buffer protocol interop.
- Ecosystem adapters (`ase`, `pymatgen`) reside in pure Python (`python/correlation/adapters.py`) using dynamic feature detection. Do not add compile-time Python package dependencies to C++.

---

## Local Development Workflow

### Building the Project
```bash
# Clone repository with submodules
git clone https://github.com/Isurwars/Correlation.git
cd Correlation

# Configure desktop build (Release mode, GUI, Tests)
cmake -B build -S . \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_GUI=ON \
  -DBUILD_TESTING=ON \
  -DBUILD_PYTHON_BINDINGS=ON

# Compile with parallel workers
cmake --build build -j$(nproc)
```

### Running Tests
```bash
# Run C++ Google Test suite
ctest --test-dir build --output-on-failure

# Run Python unit and adapter tests
pytest tests/
```

### Building WebAssembly (Optional)
```bash
emcmake cmake -B build-wasm -S . \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_WASM=ON \
  -DBUILD_GUI=OFF \
  -DBUILD_TESTING=OFF \
  -DUSE_SINGLE_PRECISION=ON

cmake --build build-wasm --target correlation_wasm
```

---

## Git Workflow & Conventional Commits

We follow trunk-based development with Conventional Commits:
- `feat(<scope>): add new calculator or feature`
- `fix(<scope>): resolve calculation bug or UI issue`
- `refactor(<scope>): simplify complexity without changing behavior`
- `test(<scope>): add unit or regression test cases`
- `docs(<scope>): update documentation or tutorial notebooks`
- `ci(<scope>): update GitHub Actions workflows`

Valid scopes: `core`, `calculators`, `readers`, `writers`, `app`, `ui`, `bindings`, `wasm`, `packaging`.

---

## License & Attribution

By contributing to **Correlation**, you agree that your contributions will be licensed under the GNU Affero General Public License v3 (AGPLv3).
