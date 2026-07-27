---
name: slint-live-preview-loop
description: Executes hot-reloading preview sessions and Slint syntax verification for iterative UI development.
---

# Slint Live Preview & Verification Loop

This skill handles real-time live previewing, hot-reloading, and syntax diagnostics for Slint markup files (`.slint`).

## 1. Syntax & Markup Validation

Before submitting or compiling `.slint` changes:
1. Ensure all imported types exist in `Types.slint` or imported module files.
2. Check that component properties use proper `in`, `out`, or `in-out` annotations.
3. Validate layout nesting constraints (`Vertical`, `Horizontal`, `GridLayout`, `ScrollView`).

## 2. Preview Workflows

### Primary: `slint-viewer` (Hot-Reload)
When `slint-viewer` is available on the system:
```bash
slint-viewer --auto-reload -I ui/ ui/AppWindow.slint
```

### Fallback: CMake Target Verification
When `slint-viewer` is not installed, validate `.slint` markup by triggering a fast Slint-to-C++ header generation build:
```bash
cmake --build build --target correlation_gui_tests -j$(nproc)
```
If Slint code generator fails, inspect compiler output for syntax/line errors in the source `.slint` file.

## 3. Responsive Layout Testing

Validate layout behavior at key viewport dimensions:
- Minimum bounds: 1050 × 725 px
- Preferred bounds: 1180 × 795 px
- Max / High-DPI bounds: 1920 × 1080 px

Ensure scrolling wrappers (`ScrollView`) accommodate overflow elements cleanly at minimum bounds.
