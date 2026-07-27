---
name: cpp-memory-profiler
description: Heap and memory footprint analysis using Valgrind Memcheck, Massif, and Heaptrack for large trajectory data processing.
---

# C++ Memory Profiling & Heap Analysis Skill

This skill provides protocols for inspecting, profiling, and optimizing heap memory usage during large trajectory analysis runs.

## 1. Valgrind Memcheck Protocol

Detect uninitialized memory reads and heap leaks:

```bash
valgrind --tool=memcheck \
         --leak-check=full \
         --show-leak-kinds=all \
         --track-origins=yes \
         --log-file=memcheck.log \
         ./build/src/correlation -i trajectory.xyz --rdf
```

---

## 2. Heaptrack Peak Memory Profiling

Track peak allocation hotspots during trajectory parsing:

```bash
heaptrack ./build/src/correlation -i trajectory.xyz --rdf
heaptrack_print heaptrack.correlation.*.gz
```

---

## 3. High-Throughput Memory Guidelines

1. **Avoid Allocations in Loop Hotspots:** Pre-allocate `std::vector::reserve()` outside iteration loops.
2. **Buffer Reuse:** Pass reusable `std::vector` buffers by reference to avoid repetitive heap allocation/deallocation overhead.
3. **Move Semantics:** Use `std::move` when transferring large trajectory frames or distribution profile arrays.
