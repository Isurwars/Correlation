---
name: sanitizer-validator
description: Instruments the C++ build with LLVM/GCC sanitizers to trap, diagnose, and fix memory leaks, undefined behavior, and race conditions at runtime.
triggers:
  - post_compilation
  - pre_test_execution
capabilities:
  - terminal_execution
  - filesystem_read_write
  - flag_injection
dependencies:
  - test-and-tea-loop
---

# Sanitizer Validator

This skill adds runtime memory and concurrency safety checks to your agentic workflow. It compiles code using LLVM/GCC Sanitizers, runs the test suite, parses runtime reports, and automatically patches bugs.

## Execution Pipeline

### Phase 1: Sanitizer Flag Injection
1. **Detect Compiler**: Verify that the workspace uses Clang or GCC.
2. **Determine Sanitizer Target**: Select the appropriate sanitizer based on the current testing context:
    - **AddressSanitizer (ASan) & LeakSanitizer (LSan)**: Default choice for memory safety (`-fsanitize=address,undefined`).
    - **ThreadSanitizer (TSan)**: Used for multi-threaded applications (`-fsanitize=thread`). Note: TSan cannot be combined with ASan.
3. **Inject Build Flags**: Modify `CMakeLists.txt` (or inject via environment variables `CXXFLAGS` and `LDFLAGS`) to include:
   ```flags
   -fsanitize=address,undefined -fno-omit-frame-pointer -g -O1
   ```

### Phase 2: Runtime Execution & Trap Capture
1. **Clean Rebuild**: Force a complete rebuild of the test suite with the injected instrumentation flags.
2. **Execute Instrumented Binaries**: Run the test suite under the sanitizer environment.
3. **Trap Runtime Interceptors**: Monitor stdout/stderr for explicit sanitizer error reports, including:
    - `ASan: heap-use-after-free` or `stack-buffer-overflow`
    - `LSan: Detected memory leaks`
    - `UBSan: undefined behavior` (e.g., null pointer dereferences, integer overflows)
    - `TSan: Data race detected`

### Phase 3: Root-Cause Remediation & Resolution Loop
If a sanitizer report is captured, treat it as a critical test failure:
1. **Parse the Stack Trace**: Extract the file path, line number, and function signature where the violation originated.
2. **Apply Architectural Fixes**:
    - **For ASan/LSan**: Replace raw pointers with smart pointers (`std::unique_ptr`, `std::shared_ptr`), fix lifecycle ownership, or correct array indexing bounds.
    - **For UBSan**: Add explicit boundary validation checks or use safe type casting.
    - **For TSan**: Introduce `std::mutex`, `std::lock_guard`, or convert shared variables to `std::atomic`.
3. **Re-Validate**: Re-run Phase 2. Continue this cycle until the entire test suite passes with zero sanitizer violations.
4. **Clean Exit**: Strip the sanitizer compilation flags before final artifact delivery to ensure production builds remain uninstrumented.
