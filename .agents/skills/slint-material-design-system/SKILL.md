---
name: slint-material-design-system
description: Directs Material Design 3 system usage in Slint UI, tokenized palettes, dark/light theme metrics, elevation, typography, and reusable Material components.
---

# Slint Material Design System Protocol

This skill governs the architecture, usage patterns, design tokens, and component library of the Material Design 3 system inside `ui/material/` for the Correlation application.

## 1. Design Tokens & Styling Metrics

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

```slint
import { MaterialStyleMetrics } from "material/material.slint";

Vertical {
    spacing: MaterialStyleMetrics.spacing_8;
    padding: MaterialStyleMetrics.spacing_12;
}
```

- **Spacing Tokens:** `spacing_2` (2px), `spacing_4` (4px), `spacing_8` (8px), `spacing_12` (12px), `spacing_16` (16px), `spacing_24` (24px)
- **Size Metrics:** `size_40` (40px compact target), `size_56` (56px standard target)
- **Border Radii:** `border_radius_1` (4px), `border_radius_2` (8px), `border_radius_3` (12px)

### C. Material Typography (`MaterialTypography`)
Use standardized typography tokens for consistent hierarchy:

```slint
import { MaterialText, MaterialTypography } from "material/material.slint";

MaterialText {
    text: "Section Header";
    style: MaterialTypography.label_large;
    color: MaterialPalette.on_surface_variant;
}
```

- **Styles:** `display_large`, `headline_medium`, `title_medium`, `body_large`, `body_small`, `label_large`

---

## 2. Material Component Library Usage

### A. Buttons (`FilledButton`, `IconButton`)
Use standard Material button variants with tokenized backgrounds and active state indicators:

```slint
import { FilledButton, IconButton, Icons } from "material/material.slint";

FilledButton {
    text: "Run Analysis";
    icon: Icons.play_arrow;
    background: MaterialPalette.run;
    text_color: MaterialPalette.on_run;
    enabled: root.file_loaded;
    clicked => { root.run_analysis(); }
}
```

### B. Cards & Containers (`FilledCard`, `Elevation`)
Wrap logical parameter groups inside `FilledCard` or `Elevation` primitives:

```slint
import { FilledCard, Vertical } from "material/material.slint";

export component ParameterCard inherits FilledCard {
    Vertical {
        padding: MaterialStyleMetrics.spacing_12;
        spacing: MaterialStyleMetrics.spacing_8;
        // Parameter controls
    }
}
```

### C. Text Inputs & Form Validation (`TextField`)
Form input fields support inline validation state, supporting error text, and label floats:

```slint
import { TextField } from "material/material.slint";

TextField {
    label: "Cutoff Radius (Å)";
    text <=> root.r_max;
    enabled: root.file_loaded;
    has_error: root.r_max_error != "";
    supporting_text: root.r_max_error;
}
```

### D. Dropdown Menus (`DropDownMenu`)
Custom compact dropdown menus for dataset selection:

```slint
import { DropDownMenu } from "material/material.slint";

DropDownMenu {
    label: "Select Plot";
    items: root.plot_items;
    current_index <=> root.selected_plot_index;
    selected(idx) => { root.plot_selected(idx); }
}
```

### E. Interactive State Layers (`StateLayer`, `TouchArea`)
Provide visual hover and touch feedback on custom clickable areas:

```slint
import { StateLayer } from "material/material.slint";

touch_area := TouchArea {
    accessible-role: button;
    accessible-label: "Clickable Region";
    
    StateLayer {
        width: 100%;
        height: 100%;
        background: MaterialPalette.on_surface;
        has_hover: touch_area.has_hover;
    }
}
```

---

## 3. Dark/Light Theme Switching

Theme state is bound globally to Slint's `Palette.color-scheme`:
- Light Theme: `ColorScheme.light`
- Dark Theme: `ColorScheme.dark`

```slint
import { Palette } from "std-widgets.slint";

IconButton {
    accessible-role: button;
    accessible-label: "Toggle Theme";
    icon: Palette.color-scheme == ColorScheme.dark ? Icons.light_mode : Icons.dark_mode;
    clicked => {
        Palette.color-scheme = Palette.color-scheme == ColorScheme.dark ? ColorScheme.light : ColorScheme.dark;
    }
}
```

---

## 4. Design Guidelines for Scientific Desktop UIs

1. **Compact Density:** Prefer `size_40` height targets over `size_56` for dense scientific parameter cards.
2. **Contrast Accessibility:** High-contrast text pairings (`on_surface` on `surface`, `on_primary` on `primary`).
3. **Responsive Flexibility:** Use relative layout wrappers (`Vertical`, `Horizontal`) with stretch metrics rather than static pixel absolute positions.
