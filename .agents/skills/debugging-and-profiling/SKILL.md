---
name: debugging-and-profiling
description: Disciplined root-cause diagnosis, LLVM/GCC sanitizers (ASan, UBSan, TSan), and heap/memory profiling (Valgrind, Heaptrack) for C++23 bugs, races, and leaks. Use when diagnosing crashes, regressions, data races, or memory bloat.
---

# Debugging, Sanitizers & Memory Profiling

## 1. Disciplined Diagnosis Loop
1. **Reproduce deterministically** with the smallest test case (`ctest -R <TestName> --output-on-failure` or `correlation-cli` on a single frame `--min-frame 1 --max-frame 1`).
2. **Isolate & Form Hypothesis:** Inspect exact call stack, `.ninja_log` timings, or state transitions before editing code.
3. **Report Format:** `[File:Line] -> [Error Category] -> [Root Cause] -> [Fix Action]`.
4. **Lock with Regression Test:** Add a unit test reproducing the bug before/with the fix.

## 2. Sanitizers (ASan / UBSan / TSan)
Configure an isolated out-of-source sanitizer build directory (never pollute `build/` or `build-release/`):
```bash
# AddressSanitizer + UndefinedBehaviorSanitizer
cmake -B build-asan -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
  -DBUILD_GUI=OFF -DBUILD_TESTING=ON
cmake --build build-asan && ctest --test-dir build-asan --output-on-failure

# ThreadSanitizer (for TBB/jthread race checks)
cmake -B build-tsan -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_CXX_FLAGS="-fsanitize=thread -fno-omit-frame-pointer" \
  -DBUILD_GUI=OFF -DBUILD_TESTING=ON
```

## 3. Heap & Memory Profiling
- **Peak RSS / Allocation Hotspots:** `heaptrack ./build-release/src/correlation-cli <input>` → inspect with `heaptrack_print`.
- **Leak Check (when ASan unavailable):** `valgrind --leak-check=full --show-leak-kinds=definite ./build/tests/correlation_unit_tests --gtest_filter=<Suite>.*`.
- **Trajectory Memory Budget:** Watch for unintended deep copies of `Scene3D`, `Cell`, or `StructureAnalyzer::raw_histograms_` (pass `r_bin_width = 0.0` when only `NeighborGraph` is needed).
