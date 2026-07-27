---
name: cpp-concurrency-patterns
description: Enforces thread-safety rules for OpenMP and TBB reduction loops, std::scoped_lock synchronization, atomic state management, and cache alignment in parallel structural calculators.
---

# C++ Concurrency & Parallel Computing Patterns

This skill governs multi-threaded execution, OpenMP/TBB safety rules, and lock-free memory locality across Correlation calculators (`RDF`, `S(Q)`, `PAD`, `VACF`, `MSD`).

## 1. OpenMP Loop Safety Rules

### Rule 1: Thread-Private Iterators
Loop variables inside `#pragma omp parallel for` constructs must be thread-private. Never share loop indices across threads.

```cpp
// GOOD: OpenMP thread-private loop bounds
#pragma omp parallel for schedule(dynamic)
for (int i = 0; i < num_atoms; ++i) {
    for (int j = i + 1; j < num_atoms; ++j) {
        // Distance calculation
    }
}
```

### Rule 2: Cache-Line Alignment (`alignas(64)`)
Eliminate false sharing by padding thread-local accumulators to 64-byte cache lines:

```cpp
struct alignas(64) ThreadLocalHistogram {
    std::vector<uint64_t> bins;
};
```

---

## 2. Synchronization Primitives

### Rule 1: Safe Mutex Locking (`std::scoped_lock`)
Always use `std::scoped_lock` for locking one or more mutexes (C++17 RAII deadlock-free locking). Avoid raw `.lock()` / `.unlock()`.

```cpp
#include <mutex>

void update_shared_state(SharedData& data1, SharedData& data2) {
    std::scoped_lock lock(data1.mutex, data2.mutex);
    // Deadlock-safe multi-mutex mutation
}
```

### Rule 2: Thread-Safe Async Slint Dispatch
Background parallel worker threads MUST dispatch UI progress updates using `slint::invoke_from_event_loop()`:

```cpp
#include <slint/slint.h>

void calculate_async(slint::ComponentHandle<AppWindow> window) {
    std::thread([window]() {
        // Parallel computation loop
        #pragma omp parallel for
        for (int i = 0; i < total; ++i) { /* ... */ }

        // Async UI update
        slint::invoke_from_event_loop([window]() {
            window->set_progress(1.0f);
        });
    }).detach();
}
```

---

## 3. Concurrency Anti-Patterns

- ❌ Modifying shared `std::vector` inside OpenMP loops without thread-local buffers or mutexes (data race).
- ❌ Using `std::endl` or raw `std::cout` inside parallel loops (interleaved IO streams).
- ❌ Omitting `alignas(64)` on thread-local reduction arrays (false sharing cache degradation).
