---
trigger: always_on
---

# Rule: C++ Code Style, Memory Safety, and Modern Standards

*Activation Mode: Glob (`**/*.{hpp,cpp,cxx,h,cc,c}`)*

## 1. Static Analysis & Formatting (Clang Tooling)
- **Tooling Compliance:** All C++ code must strictly align with workspace `.clang-format` and `.clang-tidy` rules.
- **Verification Gate:** Run `clang-tidy` against `compile_commands.json` on all modified files prior to task completion. Resolve all warnings and linter findings before presenting code.
- **Strict Compilation:** Code must compile cleanly with `-Wall -Wextra -Wpedantic -Werror`.
- **Zero Narrowing:** Implicit narrowing conversions (e.g., `double` to `float`, `size_t` to `int`) are prohibited. Use `static_cast` or explicit conversions where appropriate.

## 2. Resource & Memory Safety (RAII & Smart Pointers)
- **RAII Mandate:** All resource lifetimes (memory, files, GPU buffers, threads) must be managed using RAII wrappers.
- **Smart Pointers:** 
  - Use `std::unique_ptr` by default for single ownership.
  - Use `std::shared_ptr` / `std::weak_ptr` only when true shared runtime ownership is required.
  - **Strictly Banned:** Raw `new` and `delete`. Use `std::make_unique` or `std::make_shared`.
- **Safe Containers:** C-style arrays (`type arr[N]`) and pointer arithmetic are prohibited. Use `std::array`, `std::vector`, or `std::span` (for non-owning contiguous memory views).
- **Header Guards:** Every header file must use `#pragma once` at the very top.

## 3. Modern C++ Architecture (C++20/C++23)
- **Const Correctness:** Apply `const` and `constexpr` aggressively. Mark variables, parameters, and member functions `const` unless mutation is required.
- **Error Handling & Nullability:**
  - Use `std::optional<T>` for optional values rather than sentinel values or nullptr.
  - Use `std::expected<T, E>` or `std::variant` for typed error reporting rather than raw integer return codes.
- **Move Semantics & Pass-by-Value:**
  - Pass heavy types by `const T&` for read-only access.
  - Pass by value and `std::move` when sinks take ownership.
- **Explicit Conversions:** Mark all single-argument constructors and conversion operators as `explicit` to prevent implicit type coercions.

## 4. Doxygen Documentation Standard
- **Scope:** Every `class`, `struct`, `enum`, `concept`, and public/protected function in header files (`.hpp`) must include Doxygen blocks.
- **Formatting:** Use block style `/** ... */`.
- **Tag Directives:** Explicitly annotate `@param[in]`, `@param[out]`, `@param[in,out]`, `@return`, and `@throws`.

### Example
```cpp
/**
 * @brief Computes the pair distribution function g(r) from simulation coordinates.
 * @param[in] coordinates Atomistic positions in Cartesian space.
 * @param[in] cutoff Radius limit for the correlation calculation.
 * @return Distribution profile containing radius bins and g(r) amplitudes.
 * @throws std::invalid_argument If cutoff is non-positive.
 */
[[nodiscard]] DistributionProfile computeCorrelation(
    std::span<const Vector3D> coordinates, 
    double cutoff
);