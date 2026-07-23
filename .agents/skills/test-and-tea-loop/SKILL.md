---
name: test-and-tea-loop
description: Self-correcting build-test loop that compiles, diagnoses failures, fixes source code, and re-verifies until all tests pass.
---

# Test-and-Tea Loop

Automated build → diagnose → fix → retest loop for C++ projects using CMake and Google Test.

## Phase 1: Environment & Build Detection
1. **Locate Build Directory**: Scan the workspace for `build/`, `out/`, or `target/`.
2. **Identify Test Runner**: Check the workspace files to see if the project uses Google Test (`gtest`) or Catch2. 
3. **Establish Build Commands**: Default to `cmake --build build/` unless a specific Makefile or Bazel configuration is detected.

## Phase 2: The Self-Correction Compilation Loop
Execute the compilation command in the local environment terminal. If the build fails:
1. **Parse Compiler Logs**: Capture stdout/stderr. Extract exact file paths, line numbers, and error codes (e.g., Clang/GCC errors, missing templates, unresolved linker symbols).
2. **Diagnose and Fix**: 
    - For **syntax/type errors**, navigate to the specific file/line and rewrite the broken code.
    - For **linker errors** (unresolved external symbols), inspect `CMakeLists.txt` or the corresponding header files to ensure targets are linked correctly.
3. **Recurse**: Re-run the build command. 
4. **Loop Limit**: Stop and alert the user if the build fails 3 consecutive times on the same root cause.

## Phase 3: Test Execution & Verification
Once the build compiles successfully with 0 warnings/errors:
1. **Execute Binaries**: Run the compiled test binary (e.g., `./build/bin/run_tests`).
2. **Evaluate Test Assertions**: If any tests fail, treat the test output logs as a failure condition. Modify the source code (not the tests, unless the test logic itself is fundamentally broken or outdated) to fix the logical bug.
3. **Final Polish**: Re-verify the build one last time, run `clang-format` if available, and present the final success state.
