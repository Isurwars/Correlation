# Rule: SOLID Architectural Principles for High-Performance C++23

*Activation Mode: Glob (`**/*.{hpp,cpp,cxx,h,cc}`)*

## 1. Architecture Baseline

All software components in `Correlation` must conform to the **SOLID** design principles. Because `Correlation` is a scientific computing and HPC application, SOLID principles must be applied without sacrificing computational throughput, cache locality, or zero-cost abstractions:

- **Inner Computational Kernels:** Use **compile-time static polymorphism** (C++20 concepts, templates, CRTP) to eliminate `vtable` indirection and inline critical loops.
- **Subsystem Boundaries & Services:** Use **pure abstract interfaces and dependency injection** (I/O parsers, UI view models, simulation engines, interatomic potentials).

---

## 2. Core Principles & C++23 Enforcement

### S — Single Responsibility Principle (SRP)
- **Constraint:** A class or module must have one, and only one, reason to change.
- **Separation of Concerns:**
  - **Data vs. Compute:** Data representation (e.g., `Trajectory`, `Cell`, `Atom` SoA/AoS) must not contain calculation algorithms or file parsing logic.
  - **Compute vs. I/O:** Calculators (`RDFCalculator`, `ChiralityCalculator`, etc.) must not read files, render graphics, or perform CLI output.
  - **Compute vs. UI:** ViewModels and controllers must orchestrate data flow and event dispatch without performing in-thread scientific computation.
- **Function Limit:** Keep Cognitive Complexity $\le 25$. Decompose multi-step algorithms into pure, cohesive helper functions.

### O — Open/Closed Principle (OCP)
- **Constraint:** Software entities must be open for extension, but closed for modification.
- **Compute Loops:** Define generic algorithms constrained by C++20 concepts (`template <AtomicStructure T>` or `template <CoordinateBuffer B>`). New structure representations must integrate without altering the core calculation kernels.
- **Boundary Services:** Define pure virtual interfaces (`ITrajectoryReader`, `IPotentialModel`, `IPlotRenderer`) or concept-constrained strategy types. Add support for new formats or ML models by implementing new classes, never by mutating existing parsing switches.

### L — Liskov Substitution Principle (LSP)
- **Constraint:** Subtypes and concept models must be substitutable for their base contracts without altering program correctness or violating performance assumptions.
- **Behavioral Invariants:**
  - Concept models or derived classes must fulfill all preconditions and postconditions specified by the interface/concept.
  - Derived classes must not throw unexpected exceptions not permitted by base contracts or `noexcept` specifications.
  - Avoid hidden side effects, state mutations, or unexpected allocations in methods marked `const`.

### I — Interface Segregation Principle (ISP)
- **Constraint:** Clients should not be forced to depend on interfaces or methods they do not use.
- **Lean Interfaces:** Avoid monolithic "God Interfaces".
  - Split broad interfaces into focused, cohesive roles (e.g., `IFrameReader`, `IFrameSeeker`, `ITrajectoryWriter`).
  - Compute engines must accept `std::span<const Vector3D>` or precise concept-constrained accessors rather than requiring an entire `Trajectory` or `AppBackend` reference.
  - For Slint UI, keep controller signals and properties segregated by panel rather than bundling all application state into a single bloated controller.

### D — Dependency Inversion Principle (DIP)
- **Constraint:** High-level modules must not depend on low-level modules; both must depend on abstractions. Abstractions must not depend on details; details must depend on abstractions.
- **Dependency Injection (DI):** Pass dependencies explicitly via constructor parameters (`std::shared_ptr<IInterface>`, `std::unique_ptr<IInterface>`, or reference) or template parameters.
- **Prohibit Global State & Hidden Singletons:** Do not instantiate concrete readers, disk loggers, or compute backends directly inside high-level controllers or analyzers. Inject them as abstractions.

---

## 3. High-Performance Guardrails (Zero-Cost Abstraction Rule)

1. **No Virtual Dispatch in Tight Loops:** Never call `virtual` methods inside nested loops or parallelized regions (`#pragma omp`, TBB `parallel_for`, SIMD loops). Use template parameters or concepts for hot loops.
2. **Value Semantics & Spans at Boundaries:** Prefer passing non-owning views (`std::span`, `std::string_view`) rather than exposing concrete containers (`std::vector`).
3. **No Dynamic Allocations in Inverted Dependencies:** Injected interfaces in performance-critical paths must operate on pre-allocated buffers or zero-copy spans.

---

## 4. Quality Gate Integration

- **Code Review Gate:** Every PR or code addition must be checked against the SOLID checklist in [code-review](file:///home/isurwars/Projects/Correlation/.agents/skills/code-review/SKILL.md).
- **Prohibit Workarounds:** Never bypass SOLID decoupling with friend classes, downcasts (`dynamic_cast`), or inline linter suppressions (`NOLINT`).

## Reference

- **See skill:** [solid-design](file:///home/isurwars/Projects/Correlation/.agents/skills/solid-design/SKILL.md) for C++23 patterns, anti-pattern catalogs, and refactoring guides.
- **See skill:** [cpp-standards](file:///home/isurwars/Projects/Correlation/.agents/skills/cpp-standards/SKILL.md) for RAII, immutability, and memory safety rules.
