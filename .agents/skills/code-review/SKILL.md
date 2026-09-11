---
name: code-review
description: Comprehensive pre-commit and pre-merge code review checklist. Use when the user asks to "review code", "check changes", or before committing critical features.
---

# Code Review

Rigorous technical audit before staging, committing, or merging changes in `Correlation`. Ensures zero regressions in C++20 standards, thread safety, cognitive complexity, and build integrity.

## Phase 1: Context & Diff Inspection

1. **Inspect Working Tree & Diff:**
   ```bash
   git diff --stat
   git diff -U3
   ```
2. **Classify Affected Subsystems:**
   - **Calculators / Math:** (`src/calculators/`, `include/calculators/`)
   - **Core & Domain:** (`include/core/`, `src/core/`)
   - **UI & ViewModels:** (`ui/`, `src/app/`, `include/app/`)
   - **Build & Bindings:** (`CMakeLists.txt`, `cmake/`, `src/bindings/`)

## Phase 2: Quality Gates Audit

| Gate | Check | Command / Verification | Pass Criteria |
| :--- | :--- | :--- | :--- |
| **1. Compilation** | Strict compiler flags | `cmake --build build -j$(nproc)` | 0 warnings with `-Wall -Wextra -Wpedantic -Werror` |
| **2. Clang-Tidy** | Static analysis | `clang-tidy -p build <modified_files>` | 0 warnings. **Zero** `NOLINT` comments allowed. |
| **3. Complexity** | Cognitive complexity threshold | `clang-tidy -checks="readability-function-cognitive-complexity"` | Cognitive complexity $\le 25$ per function. |
| **4. Concurrency** | Thread safety & alignment | Code inspection | No shared mutable state without lock/reduction; `alignas(64)` on thread-local accumulators. |
| **5. RAII & Ownership** | Memory safety | Code inspection | No raw `new`/`delete`; proper use of `std::span`, `std::unique_ptr`, `std::shared_ptr`. |
| **6. Documentation** | Doxygen completeness | Check public headers | `@param`, `@return`, `@throws` present on all public API methods. |
| **7. Testing** | Unit & Functional tests | `ctest --test-dir build --output-on-failure` | 100% passing tests across all test suites. |

## Phase 3: Diagnostic Output

Present review findings using the strict diagnostic format:
```
[File:Line] → [Violation Category] → [Issue Description] → [Required Action]
```

### Review Decision Verdict
Conclude the review with one of three explicit states:
- **LGTM (Ready to Commit):** All gates pass, zero issues.
- **NEEDS WORK:** Minor warnings or styling/doc issues that must be remediated.
- **BLOCKED:** Critical defects (compilation error, thread race, memory leak, complexity $> 25$).
