# ADR 0002: Slint MVVM UI Architecture & Thread-Safe Event Dispatching

## Status
Accepted

## Context
A desktop application for materials science requires a responsive, high-framerate GUI capable of displaying rich plot previews, configuration matrices, and progress reporting without freezing during long-running atomic calculations. Legacy toolkits (Qt, wxWidgets, GTK) introduce complex build dependencies, high binary bloat, and manual thread synchronization pitfalls.

## Decision
1. **Framework Choice**: Standardize on **Slint** (Rust-backed declarative UI toolkit with C++ bindings) for the desktop GUI.
2. **MVVM Separation**:
   - **View (Slint)**: Pure declarative UI code (`ui/*.slint`), defining components, layout geometry, Material Design 3 tokens, and UI state models. Zero scientific computation logic in Slint files.
   - **Model / Backend (C++)**: Core analysis engines (`CorrelationEngine`, `AppBackend`, `DistributionFunctions`) operate purely in C++.
   - **ViewModel / Controller (C++)**: `AppController` coordinates bi-directional synchronization between the Slint event loop and backend computation threads.
3. **Thread Safety**: All heavy computational routines execute on worker threads (via `std::async` or TBB task groups). UI updates are strictly marshaled back to the Slint GUI thread via `slint::invoke_from_event_loop()`.

## Consequences
### Positive
- Sub-50MB lightweight standalone native binaries without heavy external runtime dependencies.
- Zero GUI thread freezing during large trajectory analysis.
- Clean separation of UI styling and scientific algorithms.

### Negative / Trade-offs
- Requires Rust toolchain (`cargo`, `rustc`) during build when compiling Slint from source.
- Property bindings must be carefully coordinated through `VectorModel` and model notifications to avoid redundant redraws.
