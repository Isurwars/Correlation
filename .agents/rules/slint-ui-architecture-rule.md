# Slint UI Architectural Directive

## 1. Core Architectural Separation (MVVM)
- **Declarative View Layer (`.slint`):** `.slint` files contain layout, visual styling, animations, state declarations (`in`, `out`, `in-out` properties), and user action callbacks (`callback`). No business logic or state computation belongs in `.slint`.
- **C++ Controller Layer (`.cpp`/`.hpp`):** C++ classes (`AppController`, `PlotController`, `FileIOHandler`) act as ViewModels/Controllers. They handle domain computation, file I/O, trajectory analysis, and expose signals/callbacks to Slint.
- **Generated Code Safety:** Never modify auto-generated CMake/Slint headers in `build/`. Always interact via public `slint::ComponentHandle<T>` interfaces.

## 2. Thread Safety & Event Loop Dispatch
- **Non-Blocking Main Loop:** Never execute blocking I/O, heavy trajectory analysis, or OpenMP computations directly on the main GUI thread.
- **Event Loop Dispatch:** All UI property mutations originating from background worker threads MUST be dispatched using `slint::invoke_from_event_loop()`:
  ```cpp
  slint::invoke_from_event_loop([app = app_handle_, data = std::move(results)] {
      app->set_analysis_progress(1.0f);
      app->set_status_message(slint::SharedString(data.summary));
  });
  ```

## 3. Design System & Component Hierarchy
- **Material Design Tokens:** Prefer tokenized styling definitions from `ui/material/` (`material_palette.slint`, spacing tokens, elevation, typography). Avoid arbitrary hardcoded color hex values inside isolated cards.
- **Component Decomposition:** Deconstruct large cards into modular components inside dedicated subdirectories (`ui/options/`, `ui/preview/`, `ui/run/`).
- **File Length Limits:** Keep individual `.slint` files modular and under 300 lines of code.

## 4. Property Bindings & Model Interop
- **Reactive Data Flow:** Use `in` properties for data pushed from C++ to UI. Use `out` properties for user state read by C++. Use `in-out` strictly when bi-directional binding is necessary.
- **Dynamic Collections:** Use `slint::VectorModel<T>` wrapped in `std::shared_ptr` for dynamic lists, dropdowns, and tabular views.
- **String Handling:** Pass strings using `slint::SharedString` to prevent memory copies.

## 5. Verification & Quality Gates
- **Auto-Reloading Preview:** Validate visual layout changes using `slint-viewer --auto-reload ui/AppWindow.slint`.
- **C++ GUI Targets:** Code changes must compile cleanly under `correlation_gui` and pass all GUI integration tests (`correlation_gui_tests`).
