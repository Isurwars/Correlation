# Rule: C++ Code Style, Memory Safety, and Modern Standards

*Activation Mode: Glob (`**/*.{hpp,cpp,cxx,h,cc,c}`)*

## 1. Static Analysis & Formatting (Clang Tooling)
- **Tooling Compliance:** All C++ code must strictly align with workspace `.clang-format` and `.clang-tidy` rules.
- **Verification Gate:** Run `clang-tidy` against `compile_commands.json` on all modified files prior to task completion. Resolve all warnings and linter findings before presenting code.
- **Strict Compilation:** Code must compile cleanly with `-Wall -Wextra -Wpedantic -Werror` and strict conversion-as-error gates:
  - **Clang:** `-Wfloat-conversion -Wimplicit-float-conversion -Wshorten-64-to-32 -Werror=float-conversion -Werror=implicit-float-conversion -Werror=shorten-64-to-32`
  - **GCC:** `-Wfloat-conversion -Werror=float-conversion`
  - **MSVC:** `/W4 /we4244 /we4267 /we4305`
- **Cognitive Complexity Gate:** Function cognitive complexity must be $\le 25$ (`readability-function-cognitive-complexity`). See [complexity-gate](file:///home/isurwars/Projects/Correlation/.agents/rules/complexity-gate.md).
- **Zero Narrowing & Precision Portability (`C4244`, `C4267`, `C4305`):**
  - Implicit narrowing conversions (`double` to `real_t`/`float`, `size_t` to `AtomID`/`uint32_t` or `int`, `int` to `char`) are hard compilation and `clang-tidy` errors.
  - Always cast floating-point literals assigned or passed to `real_t` with `static_cast<real_t>(...)` (and use `std::numbers::*_v<real_t>` for mathematical constants).
  - Always cast `size_t` indices to `correlation::core::AtomID` (`static_cast<AtomID>(...)`) when calling graph/atom index APIs.
  - Never pass raw `::tolower` or `std::tolower` to `std::ranges::transform`; always use `[](unsigned char chr) { return static_cast<char>(std::tolower(chr)); }`.
  - In `.cu` / GPU files, qualify `<cmath>` functions (`std::sqrt`, `std::cos`, `std::sin`, `std::atan2`) so `float` overloads are selected instead of C `double` functions.
- **Prohibit NOLINT Suppressions:** Never insert `NOLINT`, `NOLINTNEXTLINE`, `NOLINTBEGIN`, `NOLINTEND`, or any inline suppression comments. Fix the root cause in code design or types.

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

- **See rule:** [solid-principles](file:///home/isurwars/Projects/Correlation/.agents/rules/solid-principles.md) for C++23 SOLID architectural standards (concepts vs interfaces, Single Responsibility, Dependency Inversion).
- **See skill:** [solid-design](file:///home/isurwars/Projects/Correlation/.agents/skills/solid-design/SKILL.md) for practical C++23 SOLID patterns, anti-pattern catalogs, and refactoring guidelines.
- **See skill:** [cpp-standards](file:///home/isurwars/Projects/Correlation/.agents/skills/cpp-standards/SKILL.md) for comprehensive C++ coding standards and guidelines (RAII, smart pointers, modern C++20/23 architecture, const correctness, move semantics, concurrency patterns).