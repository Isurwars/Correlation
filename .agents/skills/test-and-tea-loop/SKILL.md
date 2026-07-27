---
name: test-and-tea-loop
description: Self-correcting build-test loop that compiles, diagnoses failures, fixes source code, and re-verifies until all tests pass.
---

# Test-and-Tea Loop

Automated build → diagnose → fix → retest loop for the Correlation C++ project using CMake and Google Test.

## Phase 1: Build Target Selection

Select the appropriate build target based on the scope of changes:

| Changed Files | Build Target | Test Binary |
| :--- | :--- | :--- |
| `src/calculators/**`, `src/core/**`, `include/**` | `correlation_unit_tests` | `./build/tests/correlation_unit_tests` |
| `src/io/**`, `src/parsers/**`, `tests/functional/**` | `correlation_functional_tests` | `./build/tests/correlation_functional_tests` |
| `ui/**/*.slint`, `src/app/**` | `correlation_gui_tests` | `./build/tests/correlation_gui_tests` |
| `src/cli/**` | `correlation_cli_e2e_tests` | `./build/tests/correlation_cli_e2e_tests` |
| Multiple scopes or uncertain | Build all: `correlation_unit_tests correlation_functional_tests correlation_gui_tests` | Run each binary sequentially |

## Phase 2: The Self-Correction Compilation Loop

1. **Build:**
   ```bash
   cmake --build build --target <target> -j$(nproc)
   ```

2. **On Build Failure:**
   - Parse first compiler error only (cascading errors are symptoms).
   - Use `[File:Line] → [Error Type] → [Fix Action]` diagnostic format.
   - For **syntax/type errors**: navigate to exact file/line and fix.
   - For **linker errors**: inspect `CMakeLists.txt` target linkage.
   - For **Slint generation errors**: check `.slint` syntax and property types.

3. **Loop Limit:** Stop and alert the user after **3 consecutive failures** on the same root cause.

## Phase 3: Test Execution & Verification

Once the build compiles with 0 errors:

1. **Execute:** Run the compiled test binary.
2. **On Test Failure:**
   - Parse Google Test output for `FAILED` assertions.
   - Fix the **source code**, not the test (unless the test logic is fundamentally incorrect).
   - Re-run only the failing test suite, not the entire battery.
3. **Final Verification:**
   - Re-run the full test binary once all individual failures are resolved.
   - Confirm `[  PASSED  ]` with zero failures before declaring completion.

## Phase 4: Post-Verification Cleanup

1. Run `clang-format` on all modified `.cpp`/`.hpp` files if available.
2. Run `graphify update .` if source files were modified (per graphify-maintenance skill).
3. Present final success state with test count summary.
