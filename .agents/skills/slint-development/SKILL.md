---
name: slint-development
description: Unified guide for Slint UI components, Material Design 3 tokens, C++23 MVVM backend integration, and slint-lsp 1.19.0 / slint-compiler verification. Use when editing .slint files or Slint C++ ViewModels.
---

# Slint UI & C++23 MVVM Development

Single source of truth for designing, wiring, and validating `Correlation`'s Slint UI.

## 1. Material Design 3 Component & Token Rules
- **Design Tokens:** Always use `MaterialPalette`, `MaterialStyleMetrics`, `MaterialTypography`, and `MaterialAnimations` from `ui/material/material.slint`. Avoid hardcoded hex colors or magic spacings.
- **Standard Inputs:**
  - Use `NumberField` (or `ValidatedField { stepper: true; ... }`) for bounded numeric parameters (supports `−`/`+`, `UpArrow`/`DownArrow`, precision rounding).
  - Never squeeze 3 `NumberField`s into a narrow `< 260px` sidebar row—use labeled horizontal rows (`alignment: space-between`, fixed `width: 112px`).
- **Critical Layout Rule (`TouchArea` + `if`):** Conditional `if` children inside a `TouchArea` do **not** forward `preferred-height`/`min-height` to parent layouts. Wrap them in `content := VerticalLayout { ... }` and explicitly forward `min-height: content.min-height; preferred-height: content.preferred-height;` (see `ValidatedField.slint`).
- **Struct Bundling:** Pass cohesive state via `export struct` in `ui/Types.slint` (e.g., `PhysicsDiagnostics`, `CameraHudState`, `SlabClipState`) instead of prop-drilling 5+ loose primitives through host components.

## 2. C++23 Backend & Thread Safety
- **Event Loop Dispatch:** Never mutate `AppWindow` properties from worker threads. Always wrap cross-thread updates in `slint::invoke_from_event_loop` guarded by a shared `std::atomic<bool>` `alive_` token.
- **Build Speed:** Keep `COMPILATION_UNITS ${CORRELATION_SLINT_COMPILATION_UNITS}` (default `16`) and `target_precompile_headers(app_gui_obj PRIVATE AppWindow.h)` enabled in `src/CMakeLists.txt`.
- **Zero-Copy Models:** Use `slint::VectorModel<T>` for lists/tables and `slint::SharedPixelBuffer<slint::Rgba8Pixel>` for 3D viewport & heatmap images.

## 3. Mandatory Verification Gate
Always run **both** checks after any `.slint` edit:
1. **`slint-lsp` (1.19.0)** across `ui/**/*.slint` → must report **0 diagnostics** (catches redundant default-value `INFO`s like `width: 100%`, `spacing: 0px`, `single-line: true` that `slint-compiler` ignores).
2. **`slint-compiler`** on `ui/AppWindow.slint`:
   ```bash
   SLINT_BIN=$(find build-release build -name slint-compiler -type f 2>/dev/null | head -n 1)
   "$SLINT_BIN" ui/AppWindow.slint -L material=ui/material/material.slint --style material-dark -f cpp -o /dev/null
   ```
