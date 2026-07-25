---
name: slint-live-preview-loop
description: Executes hot-reloading preview sessions and Slint LSP syntax verification for iterative UI development.
---

# Slint Live Preview & Verification Loop

This skill handles real-time live previewing, hot-reloading, and visual/syntax diagnostics for Slint markup files.

## Core Directives

1. **Hot Reloading:** Launch `slint-viewer` in auto-reload mode to visualize `.slint` changes without recompiling host C++ code:
   ```bash
   slint-viewer --auto-reload -I ui/ ui/app_window.slint
   ```
2. **Syntax Validation:** Verify `.slint` file compilation and LSP warning output before integrating into C++ CMake build scripts.
3. **Responsive Previewing:** Test UI layout across multiple window dimensions (e.g. 800x600, 1920x1080) using `slint-viewer` rendering options.
4. **Subagent Feedback Loop:** Use live preview output to iteratively tune paddings, colors, and layout flexibility until pixel-perfect.
