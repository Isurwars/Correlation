---
name: sanitizers
description: Instrument C++ build with LLVM/GCC sanitizers (ASan, UBSan, TSan) to trap and diagnose memory leaks, undefined behavior, and data races at runtime. Use when testing memory or thread safety.
---

# Runtime Sanitizers (ASan / UBSan / TSan)

This skill directs runtime memory and concurrency safety checks using dedicated presets defined in `CMakePresets.json`. It builds instrumented binaries, runs the test suite, parses runtime reports, and guides root-cause remediation.

## Execution Pipeline

### Phase 1: Sanitizer Preset Execution
Use the repository's pre-configured CMake presets (`asan` and `tsan`):

1. **AddressSanitizer (ASan) & UndefinedBehaviorSanitizer (UBSan):**
   ```bash
   cmake --preset asan
   cmake --build --preset asan -j$(nproc)
   ctest --test-dir build-asan --output-on-failure
   ```
2. **ThreadSanitizer (TSan):**
   ```bash
   cmake --preset tsan
   cmake --build --preset tsan -j$(nproc)
   ctest --test-dir build-tsan --output-on-failure
   ```
*Note: Never modify `CMakeLists.txt` to inject sanitizer flags manually. Use presets.*

### Phase 2: Trap Capture & Diagnostic Parsing
Monitor test output for explicit sanitizer error reports:
- `ASan: heap-use-after-free` or `stack-buffer-overflow`
- `LSan: Detected memory leaks`
- `UBSan: undefined behavior` (null dereferences, integer overflow, misaligned pointer)
- `TSan: Data race detected`

### Phase 3: Root-Cause Remediation & Resolution Loop
If a sanitizer report is captured, treat it as a blocking defect:
1. **Parse Stack Trace:** Extract file path, line number, and function where violation originated.
2. **Apply Architectural Fixes:**
   - **For ASan/LSan:** Ensure strict RAII ownership (`std::unique_ptr`), correct array indexing/span bounds.
   - **For UBSan:** Add boundary validation checks, replace C-style casts with checked conversions.
   - **For TSan:** Eliminate uncoordinated shared mutation; use `std::atomic`, `std::scoped_lock`, or OpenMP reduction clauses.
3. **Re-Validate:** Re-execute the corresponding preset build and test run until all tests pass cleanly with zero sanitizer warnings.
