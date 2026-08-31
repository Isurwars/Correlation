---
name: complexity-gate
description: Enforces and guides architectural refactoring to keep function cognitive complexity strictly below threshold (<= 25). Decomposes monolithic functions into modular subroutines and parameter structs without inline suppressions.
---

# Cognitive Complexity Gate & Refactoring Protocol

This skill directs the detection, diagnosis, and systematic refactoring of functions that trigger `readability-function-cognitive-complexity` violations in `clang-tidy`.

---

## 1. Quality Threshold & Guardrail
- **Cognitive Complexity Target:** $\le 25$.
- **Rule Enforcement:** Zero `NOLINT` or `NOLINTNEXTLINE` suppressions allowed.
- **Target Files:** All C++ header, implementation, and CUDA/HIP files (`*.hpp`, `*.cpp`, `*.cu`, `*.cuh`).

---

## 2. Refactoring Patterns

### Pattern A: Parameter Object Extraction
Replace functions taking 6+ related loose variables or returning multiple outputs with a typed struct:
```cpp
template <typename T>
struct SpatialPartitionData {
  GPUSearchGrid grid{};
  std::vector<T> wrapped_x;
  std::vector<T> wrapped_y;
  std::vector<T> wrapped_z;
  std::vector<int> element_ids;
  std::vector<int> atom_bin;
  std::vector<unsigned long long> bin_offsets;
  std::vector<unsigned long long> bin_indices;
  size_t num_bins{0};
  bool ignore_periodic_self_interactions{false};
};
```

### Pattern B: Pipeline Stage Decomposition
Extract discrete pipeline stages into dedicated anonymous namespace helpers:
1. `prepare_*()` or `build_*()` — Preprocessing and data structure construction.
2. `check_*()` or `has_*()` — State and condition inspection.
3. `execute_*()` — Core computation or kernel dispatch.
4. `unpack_*()` or `format_*()` — Result collation and tensor transformation.

### Pattern C: Flattening Deep Control Flow
- Replace nested ternary operators with clear structured branches or lookup tables.
- Use early returns/guard clauses to reduce indentation nesting levels.
- Extract nested loops inside conditional branches into dedicated helper functions.

---

## 3. Step-by-Step Workflow

1. **Diagnostic Scan:** Run `clang-tidy` to identify line numbers and exact nesting penalties:
   ```bash
   clang-tidy -p build path/to/file.cpp
   ```
2. **Decompose:** Identify the top 2-3 nesting clusters and extract them into single-responsibility functions.
3. **Internal Linkage:** Place helper functions in the anonymous `namespace { ... }` or mark them `static` to satisfy `misc-use-internal-linkage`.
4. **Const Correctness & Member Initialization:** Ensure all newly extracted structs have default member initializers and const-correct parameters.
5. **Verify:** Re-run `clang-tidy` and `correlation_unit_tests` to confirm clean compilation and zero warnings.
