# Rule: Slint UI/UX Architecture, Styling, and C++ Integration Standards

*Activation Mode: Glob (`**/*.slint`, `**/ui/**`, `**/gui/**`)*

## 1. Declarative Architecture & Separation of Concerns
- **Pure Declarative UI:** Define all UI visual structure, styling, micro-animations, and layout in `.slint` files. Prohibit host business logic or computation inside UI callbacks/expressions.
- **Callback & Property Contracts:** Use typed `in` / `out` / `in-out` properties for state binding, and `callback` declarations for event emission.
- **Directory Layout:** Place all `.slint` source markup under `ui/` (or `src/ui/`), keeping host C++ domain code in `src/` or `include/`.

## 2. Modern UI/UX Design System & Styling
- **Centralized Design Tokens:** Maintain global design tokens (colors, font sizes, spacing, corner radii, transition durations) in a shared `Theme` struct/singleton in `.slint`.
- **Modern Palette & Glassmorphism:** Prefer rich dark modes, high-contrast accessible palettes, subtle gradients, and rounded container borders over default flat controls.
- **Responsive Layout Containers:** Mandate layout containers (`HorizontalBox`, `VerticalBox`, `GridBox`) with relative constraints (`min-width`, `preferred-width`, `max-width`, `flex`). Avoid static pixel positioning.
- **Accessibility:** Annotate all interactive widgets with explicit `accessible-role` (e.g., `button`, `checkbox`, `slider`) and `accessible-label` properties.

## 3. C++ Target Integration & Build Safety
- **Target-Centric CMake:** Integrate `.slint` files using modern CMake `slint_target_sources(target ui/main.slint)` commands.
- **C++20 Type Compatibility:** Bind C++ backends using `slint::VectorModel<T>`, `slint::SharedString`, and type-safe generated header classes (`ui_main.h`).
- **Thread Safety & Event Loop:** Background threads computing analytical data must dispatch state changes to the Slint UI thread strictly via `slint::invoke_from_event_loop(...)`.

## 4. Live Verification & Tooling
- **Hot Reloading:** Utilize `slint-viewer --auto-reload <file.slint>` during UI design iterations.
- **LSP Validation:** Ensure zero syntax errors or compiler warnings reported by the Slint Language Server prior to committing markup changes.

## References
- **See skill:** [slint-component-architect](file:///home/isurwars/Projects/Correlation/.agents/skills/slint-component-architect/SKILL.md) for modular component patterns and widget design guidelines.
- **See skill:** [slint-cpp-binding-generator](file:///home/isurwars/Projects/Correlation/.agents/skills/slint-cpp-binding-generator/SKILL.md) for C++20 backend binding and async thread safety patterns.
- **See skill:** [slint-live-preview-loop](file:///home/isurwars/Projects/Correlation/.agents/skills/slint-live-preview-loop/SKILL.md) for live preview and hot-reload workflows.
