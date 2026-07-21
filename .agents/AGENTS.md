# Global Agent Directives

## 1. Workspace Context
- **Project Domain:** `Correlation` — C++ atomic structural analysis suite (Pair Distribution Functions, Radial Distribution Functions, Planar Angle Distributions).
- **Architecture:** Target-based CMake build framework, modern C++ (C++20), OpenMP parallelization, and machine learning interatomic potential integrations (e.g., MACE, ORB-v3).

## 2. Rule Hierarchy & Discovery Strategy
1. **Rule Precedence:** Workspace-specific rules in `.agents/rules/` override general default behaviors.
2. **Context Economy (Caveman Discovery):** 
   - Never perform blind exploratory file reads or wide `grep` queries.
   - Always consult `graphify-out/GRAPH_REPORT.md` or `graphify-out/graph.json` first to identify exact file paths and module clusters.
3. **Execution Guardrails:** 
   - Never commit or edit generated build directories (`build/`, `graphify-out/`, `CMakeCache.txt`).
   - Validate modifications against `clang-format` and `clang-tidy` rules before task completion.

## 3. Delegation & Skill Invocation
- **C++ Refactoring / Code Generation:** Activate `skills/cpp-coding-standards.md` for RAII, modern C++20/23 type safety, and memory-locality guidelines.
- **Verification & Build:** Run project validation using CMake/CTest targets before declaring features complete.