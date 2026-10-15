# Concurrency Architecture & Thread Governance

**Project:** `Correlation`  
**Standard:** Modern C++23, Intel oneAPI TBB, OpenMP SIMD  
**Document Revision:** 1.0.0 (Correlation v4.0.0)

---

## 1. Parallel Hierarchy Overview

`Correlation` employs a two-tier hybrid concurrency model designed to maximize hardware resource utilization across multi-core CPUs while preventing thread oversubscription and context-switching overhead:

1. **Outer Coarse-Grained Parallelism (Intel TBB):**
   - **Scope:** Trajectory frame decomposition, calculator pipeline task scheduling, and large batch transformations.
   - **Implementation:** `tbb::parallel_for`, `tbb::task_group`, and `tbb::concurrent_vector`.
   - **Work-Stealing:** Leverages Intel TBB's dynamic work-stealing scheduler to balance frame processing times across heterogeneous CPU cores and fluctuating frame load.

2. **Inner Fine-Grained Parallelism (OpenMP SIMD & Loop Kernels):**
   - **Scope:** Micro-level numerical kernels within individual calculator iterations (e.g., pair distance evaluation, histogram bin accumulations, SIMD vector reductions).
   - **Implementation:** `#pragma omp simd`, `#pragma omp parallel for reduction(...)`.
   - **Constraints:** Loop indices must be marked thread-private; shared reduction accumulators must maintain cache-line alignment (`alignas(64)`) to eliminate false sharing.

```
+-----------------------------------------------------------------------+
| Trajectory / Pipeline Orchestration (Intel TBB Work-Stealing Pool)     |
|   |--> Frame 0..N (tbb::parallel_for)                                 |
|   |--> Independent Calculator Tasks (tbb::task_group)                 |
+-----------------------------------------------------------------------+
                                  |
                                  v
+-----------------------------------------------------------------------+
| Numerical Calculation Kernel (OpenMP SIMD / Vectorized Inner Loops)   |
|   |--> Vectorized Pair Distance Evaluation (#pragma omp simd)         |
|   |--> Cache-Aligned Reductions (alignas(64) Histogram Accumulators)  |
+-----------------------------------------------------------------------+
```

---

## 2. Thread Oversubscription Prevention

Nested parallel regions (e.g., launching OpenMP worker teams inside TBB worker threads) cause severe thread oversubscription, excessive cache thrashing, and high latency.

### 2.1 Runtime Active Levels Guard

To strictly prevent nested thread pools from spawning multiplicative thread hierarchies, [`CorrelationEngine.cpp`](file:///home/isurwars/Projects/Correlation/src/analysis/CorrelationEngine.cpp) enforces active level restrictions at execution entry points:

```cpp
#if defined(_OPENMP)
  omp_set_max_active_levels(1);
#endif
```

This guarantees that:
- OpenMP worker teams are never instantiated recursively within existing parallel regions.
- When an OpenMP loop executes inside a TBB worker thread, OpenMP defaults to serialized execution with SIMD vectorization unless explicitly configured.

---

## 3. Cache Locality & False Sharing Prevention

1. **Alignment:** All shared parallel structures and thread-local histogram buffers must be aligned to 64 bytes (`alignas(64)`) to match modern x86_64 cache line boundaries.
2. **Padding:** Avoid packing multiple per-thread accumulators adjacently without padding.
3. **Thread Locality:** Favor thread-local accumulation buffers combined via parallel reductions over mutex-locked shared structures:

```cpp
struct alignas(64) ThreadAccumulator {
  std::vector<real_t> bins;
};
```

---

## 4. Environment Variables & Runtime Tuning

System administrators and high-performance computing (HPC) users should configure the following environment variables to tune thread allocation according to workload and core topology:

| Environment Variable | Recommended Setting | Description |
| :--- | :--- | :--- |
| `TBB_NUM_THREADS` | Physical core count (e.g., `nproc --all`) | Controls the maximum worker threads in Intel TBB's global task arena for trajectory frame processing. |
| `OMP_NUM_THREADS` | `1` (or `2` if running isolated OpenMP kernels) | Prevents thread oversubscription when running inside TBB outer parallelism. |
| `OMP_PROC_BIND` | `SPREAD` or `CLOSE` | Controls thread affinity. `SPREAD` is optimal for multi-socket NUMA systems to distribute memory bandwidth. |
| `OMP_PLACES` | `cores` | Binds OpenMP threads to physical CPU cores, avoiding hyperthread contention. |

---

## 5. Quality & Concurrency Gate Verification

All parallel algorithms and multi-threaded calculators must satisfy the following invariants before deployment:
- **Zero Data Races:** Validated using Clang ThreadSanitizer (`-fsanitize=thread`).
- **Deterministic Output:** Results across multiple thread counts must remain bit-exact or numerically identical within floating-point tolerance (`functional/analysis/ThreadDeterminismTests.cpp`).
- **Cognitive Complexity:** Functions coordinating parallel tasks must remain strictly $\le 25$ (`readability-function-cognitive-complexity`).
