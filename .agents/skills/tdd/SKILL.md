---
name: tdd
description: Test-driven development (Red-Green-Refactor) and self-correcting build-test loop. Use when implementing new functionality, writing tests before code, or iterating on test failures.
---

# Test-Driven Development (TDD) Loop

Automated build → diagnose → fix → retest loop for the Correlation C++ and Python suites using CMake, Google Test, and pytest.

## Phase 1: Build Target & Test Selection

Select the appropriate build target or test suite based on the scope of changes:

| Changed Files | Target / Suite | Execution Command |
| :--- | :--- | :--- |
| `src/calculators/**`, `src/core/**`, `src/readers/**`, `src/writers/**`, `src/analysis/**`, `src/math/**`, `include/**` | `correlation_unit_tests` | `ctest --test-dir build -R correlation_unit_tests --output-on-failure` |
| `tests/functional/**` | `correlation_functional_tests` | `ctest --test-dir build -R correlation_functional_tests --output-on-failure` |
| `ui/**/*.slint`, `src/app/**` | `correlation_gui_tests` | `ctest --test-dir build -R correlation_gui_tests --output-on-failure` |
| `src/bindings/**`, `python/**` | Python Bindings Suite | `.venv/bin/pytest tests/test_*.py` |
| Multiple scopes or uncertain | Build all C++ targets | `cmake --build build -j$(nproc) && ctest --test-dir build --output-on-failure` |

## Phase 2: The Self-Correction Compilation Loop

1. **Build:**
   ```bash
   cmake --build build --target <target> -j$(nproc)
   ```

2. **On Build Failure:**
   - Parse first compiler error only (cascading errors are symptoms).
   - Use `[File:Line] -> [Error Type] -> [Fix Action]` diagnostic format.
   - For **syntax/type errors**: navigate to exact file/line and fix.
   - For **linker errors**: inspect `CMakeLists.txt` target linkage.
   - For **Slint generation errors**: check `.slint` syntax and property types.

3. **Loop Limit:** Stop and alert the user after **3 consecutive failures** on the same root cause.

## Phase 3: Test Execution & Verification

Once the build compiles with 0 errors:

1. **Execute:** Run the compiled test binary or test runner.
2. **On Test Failure:**
   - Parse Google Test / pytest output for `FAILED` assertions.
   - Fix the **source code**, not the test (unless the test logic is fundamentally incorrect).
   - Re-run only the failing test suite, not the entire battery.
3. **Final Verification:**
   - Re-run the full target test once all individual failures are resolved.
   - Confirm `[  PASSED  ]` with zero failures before declaring completion.

## Phase 4: Post-Verification Cleanup

1. Run `clang-format -i` on all modified `.cpp`/`.hpp` files.
2. Run `graphify update .` if source files were added, renamed, or modified.
3. Present final success state with test count summary.
