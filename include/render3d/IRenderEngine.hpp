/**
 * @file IRenderEngine.hpp
 * @brief Common interface for 3D molecular visualization render engines.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "app/viewmodel/PhaseSearchTypes.hpp"
#include "render3d/Camera3D.hpp"
#include "render3d/Scene3D.hpp"
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <string_view>

namespace correlation::render3d {

/**
 * @struct ViewportDimensions
 * @brief Width and height in pixels for the rendering viewport.
 */
struct ViewportDimensions {
  uint32_t width{800};  ///< Horizontal pixel resolution.
  uint32_t height{600}; ///< Vertical pixel resolution.
};

/**
 * @struct GpuHaloConfig
 * @brief Configuration parameters for GPU-driven halo visualization.
 */
struct GpuHaloConfig {
  bool enabled{false};
  app::PhaseDescriptor descriptor{app::PhaseDescriptor::Q6Bar};
  float min_threshold{0.0F};
  float max_threshold{1.0F};
  float halo_radius_ratio{2.70F};
  float halo_intensity{0.85F};
  std::array<float, 4> halo_color{1.0F, 0.84F, 0.0F, 1.0F};
  bool dim_non_matching{false};
  bool hide_non_matching{false};
};

/**
 * @struct RenderStyleConfig
 * @brief Visual shading and background styling configuration for 3D render engines.
 */
struct RenderStyleConfig {
  bool depth_cueing{true};
  RgbaColor background_color{.r = 25, .g = 25, .b = 31, .a = 255};
};

/**
 * @class IRenderEngine
 * @brief Abstract interface representing a headless offscreen 3D rendering backend.
 */
class IRenderEngine {
public:
  virtual ~IRenderEngine() = default;

  /**
   * @brief Initializes the rendering engine and framebuffer with dimensions.
   * @param dimensions Initial viewport dimensions.
   * @return True if initialization succeeded.
   */
  [[nodiscard]] virtual bool initialize(const ViewportDimensions &dimensions) = 0;

  /**
   * @brief Resizes the underlying framebuffer.
   * @param dimensions New viewport dimensions.
   */
  virtual void resize(const ViewportDimensions &dimensions) = 0;

  /**
   * @brief Assigns the active molecular scene to render.
   * @param scene Shared pointer to scene geometry.
   */
  virtual void setScene(std::shared_ptr<const Scene3D> scene) = 0;

  /**
   * @brief Updates camera parameters for viewpoint transformation.
   * @param camera 3D orbit camera.
   */
  virtual void setCamera(const Camera3D &camera) = 0;

  /**
   * @brief Sets visibility flags for optional features (bonds, wireframe).
   * @param render_bonds True to render covalent bonds.
   * @param render_wireframe True to render unit cell wireframe lines.
   */
  virtual void setFeatureVisibility(bool /*render_bonds*/, bool /*render_wireframe*/) {}

  /**
   * @brief Configures GPU-driven halo and phase search computation.
   * @param config Halo configuration descriptor and threshold state.
   */
  virtual void setGpuHaloConfig(const GpuHaloConfig & /*config*/) {}

  /**
   * @brief Configures visual shading parameters (depth cueing, background color).
   * @param style Visual style configuration.
   */
  virtual void setRenderStyle(const RenderStyleConfig & /*style*/) {}

  /**
   * @brief Executes a single frame render and returns RGBA8 pixel span.
   * @return Contiguous span of RGBA8 pixels (size = width * height * 4).
   */
  [[nodiscard]] virtual std::span<const uint8_t> renderFrame() = 0;

  /**
   * @brief Current viewport pixel dimensions.
   * @return Viewport dimensions.
   */
  [[nodiscard]] virtual ViewportDimensions dimensions() const noexcept = 0;

  /**
   * @brief Human-readable name of the backend.
   * @return Backend name string view.
   */
  [[nodiscard]] virtual std::string_view backendName() const noexcept = 0;

  /**
   * @brief Human-readable name of the active rendering device (e.g. GPU model or CPU).
   * @return Device name string.
   */
  [[nodiscard]] virtual std::string deviceName() const { return std::string(backendName()); }

  /**
   * @brief Indicates if the engine uses hardware GPU acceleration.
   * @return True for Vulkan/GPU, false for CPU fallback.
   */
  [[nodiscard]] virtual bool isHardwareAccelerated() const noexcept = 0;
};

} // namespace correlation::render3d
