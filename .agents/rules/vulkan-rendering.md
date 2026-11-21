---
trigger: glob
globs: "shaders/**, **/render3d/**"
---

# Rule: Vulkan 3D Rendering & Compute Shaders

*Activation Mode: Glob (`shaders/**`, `**/render3d/**`)*

## 1. Zero-Copy & Dirty-Tracked Uploads
- Never deep-copy `Scene3D` on camera orbit, pan, or zoom. Use dirty-tracked `published_scene_` snapshots in `Structure3DController`.
- In `VulkanRenderEngine`, skip `uploadScene()` when `scene_.get() == last_uploaded_scene_`; call `VulkanBufferManager::uploadUbo(camera_, dimensions_)` only.

## 2. Shader & SPIR-V Hygiene
- Keep GLSL sources in `shaders/` synchronized with precompiled C++ SPIR-V headers in `include/render3d/vulkan/shaders/*Spv.hpp` (`glslc -O`).
- Avoid $O(N^2)$ per-atom loops in GPU compute shaders; use 3D spatial cell-lists for large systems ($N \ge 10^4$).
- Maintain visual and behavioral parity with CPU `SoftwareRasterizer` (keep inner pixel loops free of `std::sqrt` and divisions where squared-distance math suffices).
