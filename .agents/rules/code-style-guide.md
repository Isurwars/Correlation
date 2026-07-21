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

## 2. Doxygen Documentation Standard
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
```

## Reference

- **See skill:** `cpp-coding-standards` for comprehensive C++ coding standards and guidelines (RAII, smart pointers, modern C++20/23 architecture, const correctness, move semantics, concurrency patterns).