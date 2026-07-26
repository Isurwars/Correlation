---
name: slint-ui-development
description: Directs declarative .slint UI design, layout component decomposition, Material design system usage, property bindings, reactive layouts, animations, and live visual preview with slint-viewer.
---

# Slint UI Declarative Development Protocol

This skill provides comprehensive instructions for designing, componentizing, styling, and debugging declarative Slint UI screens and components (`.slint`).

## 1. Syntax & Declaration Rules

### Component Declaration
Declare reusable UI elements as `export component Name inherits Parent` or local `component Name`.

```slint
import { Palette } from "../material/ui/styling/material_palette.slint";
import { Card } from "../material/ui/components/card.slint";

export component OptionsCard inherits Card {
    in property <string> card-title: "Settings";
    in-out property <bool> is-expanded: true;
    callback option-changed(string, float);

    // ... layout implementation ...
}
```

### Property Scope Guidelines
- **`in property`**: Read-only properties driven by external callers or C++ controllers.
- **`out property`**: Read-only internal state exposed for external inspection.
- **`in-out property`**: Bi-directional reactive state.
- **`private property`**: Internal component state encapsulation.

---

## 2. Layout Structure & Container Rules

### Recommended Layout Hierarchies
1. **`VerticalLayout` / `HorizontalLayout`**: Primary structural containers for sequential alignment. Use `spacing` and `padding` tokens rather than manual offsets.
2. **`GridLayout`**: Matrix layouts for form fields, parameter grids, and options cards.
3. **`ScrollView`**: Wrap dynamic or height-constrained panels to handle overflowing UI elements cleanly.

```slint
VerticalLayout {
    spacing: 12px;
    padding: 16px;

    HorizontalLayout {
        spacing: 8px;
        Text {
            text: root.card-title;
            font-size: 14px;
            font-weight: 600;
        }
    }

    ScrollView {
        VerticalLayout {
            // Dynamic form contents
        }
    }
}
```

---

## 3. Styling & Material Design Integration

- **Color Palettes**: Import `material_palette.slint` tokens (`Palette.primary`, `Palette.surface`, `Palette.on-surface`). Never hardcode hex colors (e.g. `#1a202c`) directly in component primitives.
- **Typography & Font Scaling**: Use standard point sizes (`12px` caption, `14px` body, `16px` subheader, `20px` title).
- **Interactive State Layers**: Standardize touch/hover feedback using state layers or `TouchArea`:
  ```slint
  TouchArea {
      pointer-event(event) => { /* handle hover/click */ }
      clicked => { root.option-changed("cutoff", 2.5); }
  }
  ```

---

## 4. Component Decomposition Guidelines

- **Single Responsibility**: Each card or dialog gets its own `.slint` file inside logical subdirectories (`ui/options/`, `ui/preview/`, `ui/run/`).
- **Maximum File Size**: Keep files under **300 lines**. Split nested tabs, complex menus, or custom dialogs into separate child components.
- **Root Entry Point**: `ui/AppWindow.slint` serves strictly as the main window host importing and positioning sub-cards.

---

## 5. Live Iteration & Preview Workflow

- **Auto-Reload Preview**: Test visual layout and responsiveness without triggering full C++ rebuilds:
  ```bash
  slint-viewer --auto-reload ui/AppWindow.slint
  ```
- **Syntax Check**: Run `slint-lsp` or `slint-viewer` to validate markup correctness prior to C++ compilation.
