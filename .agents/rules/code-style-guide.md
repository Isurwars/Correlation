---
trigger: always_on
---

# Rule: C++ Code Style, Structure, and ECC Skill Standards
*Activation Mode: Glob (`src/**/*.{hpp,cpp,cxx,h,cc}`)*

## 1. Static Analysis & Formatting (Clang-Tidy & Clang-Format)
- **Strict Compliance:** All generated, modified, or refactored C++ code must strictly comply with the workspace `.clang-tidy` and `.clang-format` specifications.
- **Verification Gate:** Before declaring a coding task finished, you must run `clang-tidy` on the modified files using the project compile commands database. Fix any warnings or errors automatically before presenting the code to the user.
- **Compilation Hygiene:** Code modifications must compile cleanly under strict warning parameters (`-Wall -Wextra -Werror`). Avoid all implicit narrowings or type conversions that risk data loss.

## 2. Safety & Resource Management (ECC Standards)
- **RAII Mandate:** Enforce Resource Acquisition Is Initialization strictly for all resources, file descriptors, and lifetime management.
- **Smart Pointers:** 
  - Prefer `std::unique_ptr` by default for single ownership structures.
  - Use `std::shared_ptr` exclusively when shared runtime ownership is explicitly required.
  - **Strictly Banned:** Manual memory management via raw `new` or `delete` keywords.
- **Memory Boundaries:** Do not use raw C-style arrays. Utilize modern safe container wrappers such as `std::vector` or `std::array`.
- **Header Guards:** Always use `#pragma once` as the standard file header safety guard.

## 3. Modern C++ Architecture
- **Language Target:** Align implementations with C++20 or C++23 specifications.
- **Const Correctness:** Enforce aggressive `const` and `constexpr` optimization. Mark variables, parameters, and member functions as `const` unless state mutation is actively required by the logic.
- **Type Safety Primitives:** Use modern types for safe conditional management and robust error paths:
  - `std::optional` for handling optional/nullable values.
  - `std::variant` or `std::expected` for safe, typed error handling over raw status flags or error codes.

## 4. Doxygen Documentation Standard
Every class, struct, enum, and public/protected function must include valid Doxygen blocks directly above the declaration in the header (`.hpp`) file. 
- Use the `/** ... */` block style.
- Explicitly document all parameters (`@param[in,out]`), return values (`@return`), and potential exceptions (`@throws`).

### Example:
```cpp
/**
 * @brief Computes the structural correlation function.
 * @param[in] data Vector containing the raw atomistic simulation coordinates.
 * @param[in] cutoff Radius limit for the correlation calculations.
 * @return A structural distribution profile.
 */
DistributionProfile computeCorrelation(const std::vector<double>& data, double cutoff);