# Rule: Testing Standards & Coverage

*Activation Mode: Glob (`**/tests/**`, `**/*Tests.cpp`, `**/*Test.cpp`)*

## 1. Test Targets (CMake)

| Target | Scope | Binary Path |
| :--- | :--- | :--- |
| `correlation_unit_tests` | Unit tests for calculators, parsers, utilities | `build/tests/correlation_unit_tests` |
| `correlation_functional_tests` | Integration tests for full analysis pipelines | `build/tests/correlation_functional_tests` |
| `correlation_gui_tests` | Slint UI property sync and controller tests | `build/tests/correlation_gui_tests` |
| `correlation_cli_e2e_tests` | End-to-end CLI argument and output validation | `build/tests/correlation_cli_e2e_tests` |

## 2. Naming Conventions

- **Test Files:** `<ComponentName>Tests.cpp` inside `tests/unit/`, `tests/functional/`, or `tests/e2e/`.
- **Test Suites:** PascalCase class name + `Tests` suffix (e.g., `RDFCalculatorTests`, `AppControllerTests`).
- **Test Cases:** Descriptive camelCase verbs (e.g., `ComputesCorrectBinCounts`, `HandlesEmptyTrajectory`, `RejectsNegativeCutoff`).

## 3. Framework & Assertions

- **Framework:** Google Test (`gtest`) exclusively. Do not introduce Catch2, Doctest, or other frameworks.
- **Assertions:** Prefer `EXPECT_*` over `ASSERT_*` unless early termination is required. Use `EXPECT_NEAR` for floating-point comparisons with explicit tolerance.
- **Fixtures:** Use `TEST_F` with `::testing::Test` base class for shared setup/teardown.

## 4. Verification Gate

- All modified source code must pass the **relevant** test target before task completion.
- New public functions in headers require at least one corresponding unit test.
- GUI property changes require `correlation_gui_tests` to pass.
- Never modify test assertions to make failing tests pass — fix the source code instead.
