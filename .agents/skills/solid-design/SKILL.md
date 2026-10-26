---
name: solid-design
description: Apply SOLID principles adapted for modern C++23 high-performance scientific computing. Use when designing classes, decoupling modules, or reviewing architectural boundaries.
---

# SOLID Architecture for Modern C++23

This skill guides the design, implementation, and refactoring of C++23 components to enforce **SOLID** architectural principles without incurring runtime performance penalties in scientific computing workflows.

---

## Architectural Duality: Static vs. Dynamic Polymorphism

In `Correlation`, SOLID design operates on two distinct tiers:

| Tier | Target Subsystems | Mechanism | Rationale |
| :--- | :--- | :--- | :--- |
| **Compute Tier** | Calculators, inner loops, math kernels, SIMD routines | **Static Polymorphism** (C++20/23 Concepts, Templates, CRTP) | Zero runtime overhead, full compiler inlining, no `vtable` pointer chasing, optimal cache locality. |
| **Service Tier** | File I/O, UI ViewModels, ML potentials, CLI controllers | **Dynamic Polymorphism** (Pure Virtual Interfaces, DI) | Modular decoupling, runtime pluggability, isolated unit test mockability. |

---

## 1. Single Responsibility Principle (SRP)

### Directives
- **Data vs Computation:** Classes holding state (e.g., `Trajectory`, `Cell`, `Atom`) should exclusively manage lifecycle, storage invariants, and efficient memory layout (SoA/AoS). They must not implement correlation algorithms or format-specific exporters.
- **Computation vs Presentation:** Calculators (`RDFCalculator`, `ChiralityCalculator`) accept input views (`std::span`) and return structured calculation results (`DistributionProfile`). They must never touch files, Slint UI components, or stdout.
- **Decompose Large Functions:** Cognitive complexity must stay $\le 25$. If a method coordinates parsing, validation, computing, and logging, decompose it into focused helper classes or private pure subroutines.

---

## 2. Open/Closed Principle (OCP)

### Directives
- **Generic Compute Kernels:** Constrain calculation algorithms with concepts:
  ```cpp
  template <typename T>
  concept CoordinateSource = requires(const T& source, size_t index) {
      { source.position(index) } -> std::convertible_to<Vector3D>;
      { source.size() } -> std::same_as<size_t>;
  };

  template <CoordinateSource Source>
  [[nodiscard]] DistributionProfile computeRDF(const Source& coords, double cutoff);
  ```
- **Pluggable Readers & ML Models:** Subsystem services must use abstract interfaces:
  ```cpp
  class ITrajectoryReader {
  public:
      virtual ~ITrajectoryReader() = default;
      [[nodiscard]] virtual std::expected<Frame, IOError> readNextFrame() = 0;
  };
  ```
  Adding a new file format parser or ML model must require adding a new class implementing the interface, with zero edits to the core analyzer orchestration logic.

---

## 3. Liskov Substitution Principle (LSP)

### Directives
- **Contract Fulfillment:** Subclasses must never weaken preconditions or strengthen postconditions.
- **Exception Uniformity:** Implementations of an abstract interface must not throw undocumented, implementation-specific exceptions. Prefer `std::expected<T, E>` for recoverable errors.
- **Const Correctness Invariants:** Methods marked `const` must not modify observable state via `mutable` or cast away `const`, nor cause hidden thread-safety hazards.

---

## 4. Interface Segregation Principle (ISP)

### Directives
- **Role-Based Interfaces:** Replace fat interfaces with small, focused interfaces.
  ```cpp
  // Anti-pattern: Fat Interface forcing consumers to implement unused methods
  class ITrajectoryService {
  public:
      virtual Frame readFrame(size_t index) = 0;
      virtual void writeFrame(const Frame& frame) = 0;
      virtual void computeHistogram() = 0;
      virtual void renderSvg() = 0;
  };

  // Segregated: Narrow, dedicated interfaces
  class IFrameReader {
  public:
      virtual ~IFrameReader() = default;
      [[nodiscard]] virtual std::expected<Frame, IOError> readFrame(size_t index) = 0;
  };

  class IFrameWriter {
  public:
      virtual ~IFrameWriter() = default;
      [[nodiscard]] virtual std::expected<void, IOError> writeFrame(const Frame& frame) = 0;
  };
  ```
- **Lean Boundaries with `std::span`:** Pass `std::span<const Vector3D>` instead of requiring the entire `Trajectory` or `AppBackend` class.

---

## 5. Dependency Inversion Principle (DIP)

### Directives
- **Constructor Injection:** Pass abstractions as constructor arguments.
  ```cpp
  class TrajectoryAnalyzer {
  public:
      explicit TrajectoryAnalyzer(std::unique_ptr<ITrajectoryReader> reader,
                                  std::shared_ptr<ICalculatorEngine> engine)
          : reader_(std::move(reader)), engine_(std::move(engine)) {}

      [[nodiscard]] AnalysisResult run();

  private:
      std::unique_ptr<ITrajectoryReader> reader_;
      std::shared_ptr<ICalculatorEngine> engine_;
  };
  ```
- **Prohibit Hidden Singletons:** Avoid global `getInstance()` calls within domain logic. High-level orchestrators must receive their dependencies explicitly.

---

## Anti-Pattern Catalog

| Anti-Pattern | Symptom | Remediation |
| :--- | :--- | :--- |
| **God Controller** | `AppBackend` directly reads files, computes histograms, formats SVG, and handles Slint events. | Extract dedicated services (`FileReaderService`, `CalculationService`, `PlotRendererService`). Inject into `AppBackend`. |
| **Virtual Call in Hot Loop** | Interface method invoked inside `#pragma omp parallel for` or SIMD reduction. | Convert the algorithm to a template constrained by a C++20 concept; instantiate statically. |
| **Fat Data Model** | `Trajectory` class contains calculation methods (`calculateRDF()`, `calculatePDF()`). | Move calculation methods into independent calculator classes/functions in `calculators/`. |
| **Concrete Tight Coupling** | Calculator creates a hardcoded `XYZReader` directly inside its constructor or calculation method. | Pass `std::span<const Vector3D>` into calculator; instantiate reader externally. |

---

## SOLID Pre-Commit Checklist

- [ ] Does every new class have a single, well-defined responsibility?
- [ ] Are calculation routines decoupled from I/O and visualization?
- [ ] Are performance-critical loops free of virtual function calls?
- [ ] Are interfaces minimal and focused on a single role?
- [ ] Are dependencies injected rather than constructed internally via raw `new` or singletons?
- [ ] Does the design satisfy Cognitive Complexity $\le 25$?
