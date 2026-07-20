---
trigger: always_on
---

# Rule: C++ Code Style and Structure Guide
*Activation Mode: Glob (`src/**/*.{hpp,cpp,cxx,h,cc}`)*

## 1. Static Analysis & Formatting (Clang-Tidy & Clang-Format)
- **Strict Compliance:** All generated, modified, or refactored C++ code must strictly comply with the workspace `.clang-tidy` and `.clang-format` specifications.
- **Verification Gate:** Before declaring a coding task finished, you must run `clang-tidy` on the modified files using the project compile commands database. Fix any warnings or errors automatically before presenting the code to the user.

## 2. Doxygen Documentation Standard
Every class, struct, enum, and public/protected function must include valid Doxygen blocks directly above the declaration in the header (`.hpp`) file. 
- Use the `/** ... */` block style.
- Explicitly document all parameters (`@param[in,out]`), return values (`@return`), and potential exceptions (`@throws`).

### Example:
/**
 * @brief Computes the structural correlation function.
 * @param[in] data Vector containing the raw atomistic simulation coordinates.
 * @param[in] cutoff Radius limit for the correlation calculations.
 * @return A structural distribution profile.
 */
DistributionProfile computeCorrelation(const std::vector<double>& data, double cutoff);

## 3. Strict Signature Ordering (HPP vs CPP Alignment)
To maintain code scannability, the ordering of function parameters, local variable declarations, and initializer lists must follow a strict structural pattern between declaration and definition.

### Initializer Lists
In the `.cpp` constructors, initialize member variables in the exact order they are declared in the `.hpp` class definition to avoid compiler initialization order warnings.

### Parameter and Variable Ordering Template
When implementing or modifying functions, maintain the following logical grouping for input, processing, and output variables:

1. **Inputs First (Read-Only):** `const T&` or primitive inputs passed by value.
2. **In/Out Parameters:** References or pointers modified by the function.
3. **State/Iterators:** Local loop indices, flags, and structural tracking variables.
4. **Output/Return Cache:** The final variable holding the return value, declared right before the core logic block.