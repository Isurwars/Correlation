---
name: slint-component-architect
description: Guides creation of modular Slint UI components, Material Design 3 design tokens, responsive layouts, property scoping, and modern visual aesthetics.
---

# Slint Component Architect & UX Guidelines

This skill provides comprehensive standards for constructing modular, accessible, Material 3 styled visual interfaces using the Slint GUI toolkit for the Correlation application.

## 1. Syntax, Property Scoping & Decomposition

### Declarative Property Scope
- **`in property`**: Read-only properties driven by external callers or C++ controllers.
- **`out property`**: Read-only internal state exposed for external inspection.
- **`in-out property`**: Bi-directional reactive state.
- **`private property`**: Internal component state encapsulation.

### Component Decomposition Rules
- **Single Responsibility**: Each card or dialog gets its own `.slint` file inside logical subdirectories (`ui/options/`, `ui/preview/`, `ui/run/`).
- **File Length Limit**: Keep individual `.slint` files under **300 lines**. Split complex components into child components.
- **Root Entry Point**: `ui/AppWindow.slint` serves strictly as the main window host importing and positioning sub-cards.

---

## 2. Design Tokens & Material 3 Styling (`ui/material/`)

### A. Material Palette (`MaterialPalette`)
Always import and use tokenized color definitions from `MaterialPalette`. Never hardcode raw hex values (`#1a202c`, `#000`) in `.slint` components.

```slint
import { MaterialPalette } from "material/material.slint";

// Recommended token pairings:
// Background/Surface: MaterialPalette.background, MaterialPalette.surface, MaterialPalette.surface_container
// Content/Text:       MaterialPalette.on_surface, MaterialPalette.on_surface_variant
// Primary Actions:    MaterialPalette.primary, MaterialPalette.on_primary
// Errors/Validation:  MaterialPalette.error, MaterialPalette.on_error
// Run State:          MaterialPalette.run, MaterialPalette.on_run
```

### B. Material Style Metrics (`MaterialStyleMetrics`)
Use standard metric tokens for spacing, padding, component heights, and corner radii:

- **Spacing Tokens:** `spacing_2` (2px), `spacing_4` (4px), `spacing_8` (8px), `spacing_12` (12px), `spacing_16` (16px), `spacing_24` (24px)
- **Size Metrics:** `size_40` (40px compact target), `size_56` (56px standard target)
- **Border Radii:** `border_radius_1` (4px), `border_radius_2` (8px), `border_radius_3` (12px)

### C. Material Typography & Dark/Light Themes
- **Typography Styles:** `display_large`, `headline_medium`, `title_medium`, `body_large`, `body_small`, `label_large`.
- **ColorScheme Binding:** Bind theme switching globally to `Palette.color-scheme` (`ColorScheme.dark` / `ColorScheme.light`).

---

## 3. Layout Structure & Container Rules

- **`VerticalBox` / `HorizontalBox`**: Primary structural containers for sequential alignment. Use `spacing` and `padding` tokens rather than manual pixel offsets.
- **`GridBox`**: Matrix layouts for form fields, parameter grids, and options cards.
- **`ScrollView`**: Wrap dynamic or height-constrained panels to handle overflowing UI elements cleanly.

---

## 4. Reusable Material Component Library

### A. Buttons & Cards
Use standard `FilledButton`, `IconButton`, `FilledCard`, and `Elevation` primitives:

```slint
import { FilledButton, FilledCard, Icons, MaterialPalette } from "material/material.slint";

FilledButton {
    text: "Run Analysis";
    icon: Icons.play_arrow;
    background: MaterialPalette.run;
    text_color: MaterialPalette.on_run;
    clicked => { root.run_analysis(); }
}
```

### B. Interactive Feedback & Accessibility Annotations
- **Accessibility:** Guarantee every interactive widget declares standard `accessible-role` (e.g., `button`, `checkbox`, `slider`, `text-input`) and descriptive `accessible-label` properties.
- **State Layers:** Provide visual hover and touch feedback on custom clickable areas using `StateLayer` and `TouchArea`.
- **Micro-Animations:** Animate property transitions using `animate` directives with `150ms`–`250ms` durations and easing curves.

---

## 5. Data Visualization Guidelines

- **Plotting Canvas Elements:** Use custom `Path` or `Image` elements to stream rendered data charts (e.g. Pair Distribution Functions g(r) or Radial Distribution Profiles).
- **Control Panels:** Group analysis parameters (cutoff radius, bin width, particle filters) in collapsable card widgets with active state indicators.

