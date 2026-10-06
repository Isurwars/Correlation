# Rule: Testing Standards & Coverage

*Activation Mode: Glob (`**/tests/**`, `**/*Tests.cpp`, `**/*Test.cpp`, `**/tests/test_*.py`)*

## 1. Test Targets & Suites

| Target / Suite | Scope | Execution Command |
| :--- | :--- | :--- |
| `correlation_unit_tests` | Unit tests for calculators, readers, writers, math, core | `ctest --test-dir build -R correlation_unit_tests --output-on-failure` |
| `correlation_functional_tests` | Integration tests for end-to-end analysis pipelines | `ctest --test-dir build -R correlation_functional_tests --output-on-failure` |
| `correlation_gui_tests` | Slint UI property sync and controller ViewModel tests | `ctest --test-dir build -R correlation_gui_tests --output-on-failure` |
| `Python Bindings Suite` | NumPy / PyBind11 bindings and adapters validation | `.venv/bin/pytest tests/test_*.py` |
| `Fuzz Suite` | LLVM libFuzzer targets for trajectory readers | Built under `tests/fuzz/` when fuzzing enabled |

## 2. Naming Conventions

- **Test Files:** `<ComponentName>Tests.cpp` inside `tests/unit/` or `tests/functional/`. Python tests named `test_<feature>.py`.
- **Test Suites:** PascalCase class/component name + `Tests` suffix (e.g., `RDFCalculatorTests`, `AppControllerTests`).
- **Test Cases:** Descriptive PascalCase verbs (e.g., `ComputesCorrectBinCounts`, `HandlesEmptyTrajectory`, `RejectsNegativeCutoff`).

## 3. Framework & Assertions

- **C++ Framework:** Google Test (`gtest`) exclusively for unit and integration testing. Do not introduce Catch2, Doctest, or other C++ test frameworks.
- **Python Framework:** `pytest` exclusively for Python bindings and adapter regression tests.
- **Assertions:** Prefer `EXPECT_*` over `ASSERT_*` unless early test termination is required. Use `EXPECT_NEAR` for floating-point comparisons with explicit epsilon tolerance.
- **Fixtures:** Use `TEST_F` with `::testing::Test` base class for shared setup/teardown.

## 4. Verification Gate

- All modified source code must pass the **relevant** test target before task completion.
- New public functions in headers require at least one corresponding unit test.
- GUI property and ViewModel modifications require `correlation_gui_tests` to pass.
- Python bindings edits (`src/bindings/`) require `pytest tests/test_*.py` to pass.
- Never modify test assertions to make failing tests pass — fix the source code instead.
