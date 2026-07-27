---
name: cpp-error-diagnosis
description: Structured compiler and runtime error parsing for C++20, OpenMP, linker failures, and sanitizer reports with root-cause remediation patterns.
---

# C++ Error Diagnosis Skill

This skill provides structured diagnostic protocols for interpreting and resolving C++ compilation, linkage, and runtime errors in the Correlation project.

## 1. Diagnostic Output Format

All error reports must follow the strict format:

```
[File:Line] → [Error Category] → [Root Cause] → [Fix Action]
```

**Example:**
```
[src/calculators/RDFCalculator.cpp:142] → Template Instantiation Error →
  Missing `#include <algorithm>` for `std::ranges::sort` →
  Add `#include <algorithm>` to file header
```

## 2. Common Error Patterns & Resolution

### Compiler Errors (Clang/GCC)

| Error Pattern | Root Cause | Fix |
| :--- | :--- | :--- |
| `no matching function for call` | Wrong argument types or missing overload | Check parameter types, add `static_cast` or correct overload |
| `use of undeclared identifier` | Missing `#include` or namespace qualification | Add include or `using` declaration |
| `no member named 'X' in 'Y'` | Accessing non-existent member or wrong struct type | Verify struct definition in header |
| `cannot convert 'X' to 'Y'` | Implicit narrowing or type mismatch | Use `static_cast<Y>(x)` or correct the source type |
| `constexpr variable must be initialized` | Missing initializer for `constexpr` | Provide compile-time constant value |
| `concept 'X' not satisfied` | Template argument doesn't meet concept constraints | Check concept requirements and add missing operations |

### Linker Errors

| Error Pattern | Root Cause | Fix |
| :--- | :--- | :--- |
| `undefined reference to 'X'` | Missing `target_link_libraries` or source file not compiled | Add library to CMake target or include `.cpp` in sources |
| `multiple definition of 'X'` | Symbol defined in header without `inline` | Add `inline` keyword or move definition to `.cpp` |
| `undefined symbol: vtable` | Missing virtual destructor implementation | Implement `virtual ~ClassName() = default;` in `.cpp` |

### OpenMP / Concurrency Errors

| Error Pattern | Root Cause | Fix |
| :--- | :--- | :--- |
| `data race on variable 'X'` (TSan) | Shared mutable state without synchronization | Use `std::atomic`, `std::mutex`, or OpenMP `reduction` clause |
| `false sharing` (performance) | Thread-local data on same cache line | Use `alignas(64)` on thread-local accumulators |
| `use-after-free` in parallel loop | Iterator invalidation or dangling reference | Copy data into thread-local storage before parallel region |

### Sanitizer Reports

| Sanitizer | Report Prefix | Action |
| :--- | :--- | :--- |
| ASan | `heap-use-after-free`, `stack-buffer-overflow` | Fix lifetime/bounds; replace raw pointers with smart pointers |
| LSan | `Detected memory leaks` | Ensure RAII wrappers for all allocations |
| UBSan | `undefined behavior` | Add bounds checks, fix signed overflow, validate pointer alignment |
| TSan | `Data race detected` | Add synchronization primitives |

## 3. Diagnostic Workflow

1. **Capture:** Run build command, capture full stderr output.
2. **Extract:** Parse first error only (cascading errors are usually symptoms).
3. **Classify:** Match error text against patterns in tables above.
4. **Locate:** Navigate to exact file and line from error message.
5. **Fix:** Apply the prescribed fix action.
6. **Verify:** Rebuild and confirm error is resolved before proceeding to next error.
