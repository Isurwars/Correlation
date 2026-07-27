# Rule: Slint UI Architecture, Styling, and C++ Integration Standards

*Activation Mode: Glob (`**/*.slint`, `**/ui/**`, `**/gui/**`)*

## 1. Declarative Architecture & Separation of Concerns (MVVM)

- **Pure Declarative View Layer (`.slint`):** `.slint` files contain layout, visual styling, animations, state declarations (`in`, `out`, `in-out` properties), and user action callbacks (`callback`). No business logic or state computation belongs in `.slint`.
- **C++ Controller Layer (`.cpp`/`.hpp`):** C++ classes (`AppController`, `PlotController`, `FileIOHandler`) act as ViewModels/Controllers. They handle domain computation, file I/O, trajectory analysis, and expose signals/callbacks to Slint.
- **Generated Code Safety:** Never modify auto-generated CMake/Slint headers in `build/`. Always interact via public `slint::ComponentHandle<T>` interfaces.
- **Callback & Property Contracts:** Use typed `in` / `out` / `in-out` properties for state binding, and `callback` declarations for event emission.
- **Directory Layout:** Place all `.slint` source markup under `ui/` with logical subdirectories (`ui/options/`, `ui/preview/`, `ui/run/`, `ui/nav/`).

## 2. Modern UI/UX Design System & Styling

- **Centralized Design Tokens:** Maintain global design tokens (colors, font sizes, spacing, corner radii, transition durations) in `ui/material/` (`material_palette.slint`, spacing tokens, elevation, typography). Never hardcode hex colors in component files.
- **Modern Palette & Aesthetics:** Prefer rich dark modes, high-contrast accessible palettes, subtle gradients, and rounded container borders over default flat controls.
- **Responsive Layout Containers:** Mandate layout containers (`HorizontalBox`, `VerticalBox`, `GridBox`) with relative constraints (`min-width`, `preferred-width`, `max-width`). Avoid static pixel positioning.
- **Accessibility:** Annotate all interactive widgets with explicit `accessible-role` (e.g., `button`, `checkbox`, `slider`, `text-input`) and `accessible-label` properties.
- **Micro-Animations:** Animate property changes (hover, active states, expand/collapse) using `animate` directives with appropriate easing and `150ms`–`250ms` durations.

## 3. Component Decomposition & File Limits

- **Single Responsibility:** Each card or dialog gets its own `.slint` file inside logical subdirectories.
- **File Length Limit:** Keep individual `.slint` files under **300 lines**. Split complex components into child components.
- **Root Entry Point:** `ui/AppWindow.slint` serves strictly as the main window host importing and positioning sub-cards.

## 4. Thread Safety & Event Loop Dispatch

- **Non-Blocking Main Loop:** Never execute blocking I/O, heavy trajectory analysis, or OpenMP computations on the main GUI thread.
- **Event Loop Dispatch:** All UI property mutations from background threads MUST use `slint::invoke_from_event_loop()`.
- **Dynamic Collections:** Use `slint::VectorModel<T>` wrapped in `std::shared_ptr` for dynamic lists, dropdowns, and tabular views.
- **String Handling:** Pass strings using `slint::SharedString` to prevent memory copies.

## 5. C++ Target Integration & Build Safety

- **Target-Centric CMake:** Integrate `.slint` files using `slint_target_sources(target ui/AppWindow.slint)`.
- **C++20 Type Compatibility:** Bind C++ backends using generated header classes and type-safe `slint::` APIs.
- **Verification:** Code changes must compile cleanly and pass `correlation_gui_tests`.

## References
- **See skill:** [slint-component-architect](file:///home/isurwars/Projects/Correlation/.agents/skills/slint-component-architect/SKILL.md) for modular component patterns and widget design guidelines.
- **See skill:** [slint-cpp-integration](file:///home/isurwars/Projects/Correlation/.agents/skills/slint-cpp-integration/SKILL.md) for C++20 backend binding, async thread safety, and data model patterns.
- **See skill:** [slint-live-preview-loop](file:///home/isurwars/Projects/Correlation/.agents/skills/slint-live-preview-loop/SKILL.md) for live preview and hot-reload workflows.
- **See skill:** [slint-ui-development](file:///home/isurwars/Projects/Correlation/.agents/skills/slint-ui-development/SKILL.md) for declarative `.slint` layout design and property bindings.
