---
name: cpp-performance-benchmark
description: Enforces data-oriented design (SoA vs AoS), SIMD loop friendliness, Google Benchmark integration, cache-line alignment, and microbenchmarking best practices.
---

# C++ Performance, Memory Locality, and Microbenchmarking Protocol

This skill governs high-performance C++ execution, memory-locality optimizations, SIMD vectorization contracts, and microbenchmarking protocols across `Correlation`.

## 1. Memory Locality & Data-Oriented Design (DOD)

1. **Contiguous Memory First:** Prefer `std::vector`, `std::array`, or contiguous `std::span` buffers over pointer-chasing structures (`std::list`, raw pointer arrays, node-based trees).
2. **SoA vs AoS Layout Selection:**
   - **Array of Structures (AoS):** Suitable for scalar access or non-vectorized operations (`struct Atom { Vector3D pos; Vector3D vel; double mass; };`).
   - **Structure of Arrays (SoA):** Required for high-throughput SIMD vector loops (`struct AtomBuffer { std::vector<float> x, y, z; };`) to ensure stride-1 sequential memory access and alignment.
3. **Cache Line Alignment (`alignas(64)`):**
   - Align thread-local accumulator structures to 64 bytes (`alignas(64)`) to eliminate false sharing in multi-threaded OpenMP/TBB reduction loops.

```cpp
// GOOD: Cache-line aligned thread-local accumulator
struct alignas(64) ThreadLocalHistogram {
    std::vector<uint64_t> bins;
};
```

---

## 2. SIMD Vectorization Contracts

1. **Vector-Friendly Loops:** Keep inner distance/histogram computation loops free of conditional branching (`if`/`else`), virtual function calls, or complex control flow.
2. **Contiguous Access:** Ensure pointer reads inside inner loops process contiguous memory blocks without strided pointer indirection.
3. **Use SIMD Intrinsics / Wrappers:** Prefer vectorized SIMD wrappers (`Vector3SIMD`, `dot_block`) or standard C++20 `std::span` over custom assembly instructions unless profiling explicitly demonstrates auto-vectorization failure.

---

## 3. Microbenchmarking Protocol (Google Benchmark)

1. **Dead-Code Elimination Guard:** Wrap benchmarked output variables in `benchmark::DoNotOptimize(...)` and call `benchmark::ClobberMemory()` to prevent compiler optimization passes from discarding calculated values.

```cpp
#include <benchmark/benchmark.h>

static void BM_RadialDistribution(benchmark::State& state) {
    auto trajectory = create_test_trajectory();
    RDFCalculator calculator;
    
    for (auto _ : state) {
        auto result = calculator.calculate(trajectory);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}
BENCHMARK(BM_RadialDistribution);
```

2. **Isolation & Statistical Rigor:**
   - Run benchmarks under Release builds with target architecture flags (`-O3 -march=native`).
   - Ensure benchmark targets are reproducible and isolated from external I/O operations.
