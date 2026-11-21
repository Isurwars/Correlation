---
name: vulkan-3d-rendering
description: Guide Vulkan offscreen rendering, GLSL compute shaders, SPIR-V header regeneration, zero-copy Scene3D snapshots, and CPU SoftwareRasterizer parity in src/render3d/ and shaders/. Use when editing 3D rendering, halos, or compute shaders.
---

# Vulkan 3D Rendering & Compute Shaders

## 1. Shader Compilation Flow (GLSL → SPIR-V C++ Headers)
- GLSL compute & graphics shaders live in `shaders/` (e.g., `shaders/compute/*.comp`).
- Compiled SPIR-V bytecode is embedded as `inline constexpr uint32_t ...[]` headers in `include/render3d/vulkan/shaders/*Spv.hpp` so builds succeed even without `glslc`.
- Whenever modifying a `.comp`/`.vert`/`.frag` shader, recompile with `glslc -O` and regenerate the matching `*Spv.hpp` header, then verify `VulkanComputeServiceTests`.

## 2. Performance Invariants for Large Systems (16k–250k Atoms)
- **Zero-Copy Camera Interaction:** `Structure3DController` caches `published_scene_` and only clones `*scene_` when `scene_dirty_` is set. Never deep-copy `Scene3D` on camera orbit/zoom/pan.
- **Dirty-Tracked VRAM Uploads:** `VulkanRenderEngine` tracks `last_uploaded_scene_`. When geometry is unchanged, call `VulkanBufferManager::uploadUbo(camera_, dimensions_)` only—never re-upload atom/bond SSBOs on camera moves.
- **No $O(N^2)$ GPU Loops:** New or upgraded GPU neighbor kernels must use a 3D uniform spatial cell-list ($O(N)$), not `for (j = 0; j < num_atoms; ++j)`.
- **CPU Fallback Parity:** Keep `SoftwareRasterizer` fast (`sqrt`-free squared-distance halo falloff) and visually consistent with Vulkan impostor shaders.
- **Live Playback Sync:** Frame seek/step callbacks must call `syncPhaseSearchToWindow()` when phase search is active, and async descriptor completions must preserve/refresh the latest 2D heatmap without blanking on cache misses.
