/**
 * @file IRenderEngine.hpp
 * @brief Common interface for 3D molecular visualization render engines.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "render3d/Camera3D.hpp"
#include "render3d/Scene3D.hpp"
#include <cstdint>
#include <memory>
#include <span>
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
   * @brief Indicates if the engine uses hardware GPU acceleration.
   * @return True for Vulkan/GPU, false for CPU fallback.
   */
  [[nodiscard]] virtual bool isHardwareAccelerated() const noexcept = 0;
};

} // namespace correlation::render3d
