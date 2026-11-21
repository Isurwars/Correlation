---
name: cpp23-scientific
description: Enforce modern C++23 standards, SOLID zero-cost abstractions, RAII, and cache-aligned OpenMP/TBB parallelism for scientific calculators. Use when writing, refactoring, or profiling C++23 compute code.
---

# C++23 Scientific Computing & SOLID Architecture

Unified guide for modern C++23 idioms, SOLID design, and multi-threaded parallelism in `Correlation`.

## 1. Modern C++23 Standards & RAII
- **Language:** Strictly C++23 (`std::span`, `std::ranges`, `std::expected`, `constexpr`, `std::jthread`, designated initializers).
- **Ownership:** Zero raw `new`/`delete`. Prefer value semantics, `std::unique_ptr`, and `std::shared_ptr<const T>` for immutable snapshots shared across threads.
- **Immutability & Safety:** Mark non-mutating queries `[[nodiscard]]` and `const noexcept` where applicable. Prefer `const auto` on `static_cast<T>(...)` (`modernize-use-auto`) and `.at()` for `std::array` indexing.
- **Precision Discipline:** Respect `real_t` (`USE_SINGLE_PRECISION=ON` → `float`). Perform sensitive accumulations (OLS regressions, variance sums, entropy integrals) in `double` internally, then cast to `real_t`.
- See [core-guidelines.md](references/core-guidelines.md) for detailed C++23 examples.

## 2. SOLID Adapted for High-Performance Compute
- **Zero-Cost Hot Loops (OCP/LSP):** Use C++20/23 concepts and templates—never virtual dispatch—inside per-atom or per-pair compute loops.
- **Boundary Interfaces (DIP/ISP):** Use pure abstract interfaces (`BaseCalculator`, `IRenderEngine`) and dependency injection (`std::function` providers) at architectural boundaries (`src/app/`, `src/render3d/`).
- **Single Responsibility (SRP):** Keep pure math stateless (`HyperuniformityCalculator::computeRealSpaceDiagnostics`), separate from UI formatting (`PlotController`) and persistence (`FileWriter`).
- **Cognitive Complexity:** Keep every function strictly `<= 25` (`readability-function-cognitive-complexity`), zero `NOLINT` suppressions.

## 3. Thread Safety & Cache Locality (TBB / OpenMP)
- **False Sharing Prevention:** Align shared/per-thread accumulator structs with `alignas(64)` (`DynamicProperties`, `HyperuniformityDiagnostics`).
- **Locking:** Always use `std::scoped_lock` / `std::unique_lock`; never hold a mutex across heavy computation or Slint UI event-loop dispatches.
- **TBB Parallelism:** Prefer `tbb::parallel_for` / `tbb::parallel_reduce` with `tbb::enumerable_thread_specific` for thread-local histograms, reducing deterministically after the parallel region.
- **Thresholding:** Gate parallel dispatch on workload size (e.g., `atoms.size() >= 2048`) to avoid task-scheduling overhead on tiny cells.
