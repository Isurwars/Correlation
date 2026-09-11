---
name: to-spec
description: Transform brainstorming notes, feature requests, or design discussions into an actionable, formal technical specification. Use when user asks to "spec this out", "write a spec", or before kicking off major architectural work.
---

# To-Spec

Turn unstructured discussions, user requests, or grilling sessions into a crisp, unambiguous technical specification for the `Correlation` suite.

## Phase 1: Information Gathering & Scope Extraction

Extract and establish the core technical boundaries:
1. **Problem Statement:** What physics or software deficiency is being solved?
2. **Domain Scope:** Which structural analysis calculators or modules are affected (`RDF`, `PDF`, `SQ`, `PAD`, `Rings`, `MSD`, `VACF`, `VDOS`, `LEF`, ML potentials)?
3. **Consumers & Interfaces:** Is this consumed by Slint UI, Python bindings (pybind11), C++ Core API, or CLI?
4. **Explicit Non-Goals:** What is intentionally omitted from this milestone?

## Phase 2: Specification Document Template

Generate a structured technical specification using the following schema:

```markdown
# Spec: [Feature / Module Name]

## 1. Overview & Architectural Motivation
[Concise executive summary of the addition or change.]

## 2. API Contract & Modern C++ Interfaces
- **Header:** `include/path/to/Header.hpp`
- **Class / Struct Definitions:** (Follow C++20 concepts, RAII, immutable by default)
- **Error Handling:** Use `std::expected<T, ErrorCode>` or clear exception boundaries.

## 3. Data Flow & Algorithmic Design
- **Input Data Structures:** Trajectory coordinates, simulation cell metrics, element types.
- **Parallelization Strategy:** OpenMP reduction / TBB parallel loop, memory alignment (`alignas(64)`).
- **Time & Space Complexity:** Bounds on execution time ($O(N^2)$, $O(N \log N)$ via neighbor lists) and RAM footprint.

## 4. UI / Slint Integration (If Applicable)
- **Slint Properties & Models:** VectorModel bindings, reactive callbacks, error states.
- **Thread Dispatch:** Safe worker-to-GUI dispatch mechanisms.

## 5. Test & Quality Plan
- **Unit Tests:** Atomic mathematical tests in `tests/unit/`.
- **Functional Tests:** Golden-sample regression tests with benchmark datasets in `tests/functional/`.
- **Complexity Gate:** Verify cognitive complexity $\le 25$.
```

## Phase 3: Review & Understanding Lock

1. Present the completed specification to the User (Lead Architect).
2. Prompt: *"Does this specification accurately capture all technical invariants and constraints? Please confirm to proceed to ticket breakdown or implementation."*
