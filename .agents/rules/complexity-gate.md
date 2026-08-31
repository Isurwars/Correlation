# Rule: Function Cognitive Complexity Gate

*Activation Mode: Glob (`**/*.{hpp,cpp,cxx,h,cc,c,cu,cuh}`)*

## 1. Core Principle & Threshold
- **Threshold Limit:** Every function must maintain a **Cognitive Complexity <= 25** as enforced by `clang-tidy` (`readability-function-cognitive-complexity`).
- **Zero Inline Suppressions:** Strictly prohibit `NOLINT`, `NOLINTNEXTLINE`, or any inline suppression comments to bypass the complexity check. Root causes must be resolved through architectural decomposition.

## 2. Refactoring & Decomposition Strategies
When a function exceeds the complexity threshold, apply the following structural patterns:

1. **Extract Data Structs (Parameter Objects):**
   - Bundle multi-variable setup, search grids, or execution states into typed structs (e.g., `SpatialPartitionData`, `TwoPassExecutionContext`).
   - Eliminate scattered state variables and deep conditional initialization blocks.

2. **Single Responsibility Sub-routines:**
   - Extract discrete operational stages into private helper functions (in anonymous `namespace { ... }` or `private:` methods):
     - Grid & coordinate partitioning (`build_spatial_partition()`)
     - Flag & condition checks (`has_active_bonds()`)
     - Multi-pass execution loops (`execute_two_pass_bonds()`)
     - Post-processing & tensor unpacking (`copy_and_unpack_histograms()`)

3. **Inversion of Control & Early Returns:**
   - Use guard clauses (`if (atom_count == 0) return;`) at the top of functions.
   - Avoid deep nesting (`if` inside `for` inside `for` inside `if`) which exponentially inflates cognitive complexity scores.

4. **Standard Algorithm Delegation:**
   - Replace manual nested loops for condition detection with standard library algorithms (`std::ranges::any_of`, `std::all_of`, `std::accumulate`).

## 3. Verification Protocol
- **Gate Check:** Before declaring any C++ implementation or refactoring complete, execute `clang-tidy` against `compile_commands.json` on the modified files:
  ```bash
  clang-tidy -p build <path/to/modified_file>
  ```
- Verify zero `readability-function-cognitive-complexity` warnings are emitted.
