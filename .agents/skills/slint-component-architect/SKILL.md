---
name: slint-component-architect
description: Guides creation of modular Slint UI components, design tokens, responsive layouts, and modern aesthetics.
---

# Slint Component Architect & UX Guidelines

This skill provides standards for constructing high-density, accessible, and responsive visual interfaces using the Slint GUI toolkit.

## Core Directives

1. **Design Token Centralization:** Define global UI design tokens (palette, spacing, typography, motion curves) in a unified global component or struct.
2. **Modular Widget Encapsulation:** Decompose complex UIs into clean, single-responsibility components with explicit `in`, `out`, and `in-out` properties and explicit `callback` declarations.
3. **Adaptive Auto-Layout:** Use `VerticalBox`, `HorizontalBox`, and `GridBox` wrappers with `spacing` and `padding` parameters to guarantee responsive window resizing.
4. **Interactive States & Micro-animations:** Animate property changes (e.g., hover scaling, background tint shifts, active button presses) using `animate` directives with custom easing functions.
5. **Accessibility Annotations:** Guarantee every interactive element declares standard `accessible-role` and descriptive `accessible-label` properties.

## Component Pattern Example

```slint
export global Theme {
    in-out property <color> background-dark: #0f172a;
    in-out property <color> surface-card: #1e293b;
    in-out property <color> accent-primary: #38bdf8;
    in-out property <color> text-primary: #f8fafc;
    in-out property <length> border-radius: 8px;
}

export component ActionButton inherits Rectangle {
    in property <string> text: "Submit";
    in property <bool> enabled: true;
    callback clicked();

    background: touch-area.has-hover ? Theme.accent-primary.darker(0.15) : Theme.accent-primary;
    border-radius: Theme.border-radius;
    animate background { duration: 150ms; easing: ease-in-out; }

    HorizontalBox {
        alignment: center;
        Text {
            text: root.text;
            color: Theme.text-primary;
            font-weight: 600;
        }
    }

    touch-area := TouchArea {
        enabled: root.enabled;
        clicked => { root.clicked(); }
        accessible-role: button;
        accessible-label: root.text;
    }
}
```

## Data Visualization Guidelines
- **Plotting Canvas Elements:** Use custom `Path` or `Image` elements to stream rendered data charts (e.g. Pair Distribution Functions g(r) or Radial Distribution Profiles).
- **Control Panels:** Group analysis parameters (cutoff radius, bin width, particle filters) in collapsable card widgets with active state indicators.
