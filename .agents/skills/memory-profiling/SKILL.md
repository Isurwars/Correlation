---
name: memory-profiling
description: Heap and memory footprint analysis using Valgrind Memcheck, Massif, and Heaptrack for large trajectory data processing. Use when diagnosing memory leaks or memory bloat.
---

# Memory Profiling & Heap Analysis

This skill provides protocols for inspecting, profiling, and optimizing heap memory usage during large trajectory analysis runs.

## 1. Tool Availability & Preflight

Before executing memory profiling tools, check if they are installed on the system:
```bash
command -v valgrind || echo "valgrind not found"
command -v heaptrack || echo "heaptrack not found"
```

If neither tool is installed, fall back to **AddressSanitizer / LeakSanitizer** via the repository preset:
```bash
cmake --preset asan && cmake --build --preset asan && ./build-asan/src/correlation <input_file> --r-max 10.0
```

---

## 2. Valgrind Memcheck Protocol

Detect uninitialized memory reads and heap leaks:

```bash
valgrind --tool=memcheck \
         --leak-check=full \
         --show-leak-kinds=all \
         --track-origins=yes \
         --log-file=memcheck.log \
         ./build/src/correlation trajectory.xyz -o output --r-max 10.0
```

---

## 3. Heaptrack Peak Memory Profiling

Track peak allocation hotspots during trajectory parsing:

```bash
heaptrack ./build/src/correlation trajectory.xyz -o output --r-max 10.0
heaptrack_print heaptrack.correlation.*.gz
```

---

## 4. High-Throughput Memory Guidelines

1. **Avoid Allocations in Loop Hotspots:** Pre-allocate `std::vector::reserve()` outside iteration loops.
2. **Buffer Reuse:** Pass reusable `std::vector` buffers by reference to avoid repetitive heap allocation/deallocation overhead.
3. **Move Semantics:** Use `std::move` when transferring large trajectory frames or distribution profile arrays.
