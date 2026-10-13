/**
 * @file RenderEngineFactory.hpp
 * @brief Factory for selecting and instantiating 3D rendering engines with automatic fallback.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "render3d/IRenderEngine.hpp"
#include <memory>

namespace correlation::render3d {

/**
 * @enum class BackendPreference
 * @brief Preferred rendering backend requested by application or user configuration.
 */
enum class BackendPreference {
  Auto,       ///< Try Vulkan GPU first; automatically fall back to CPU rasterizer.
  VulkanOnly, ///< Require Vulkan GPU; fail if unavailable.
  CpuOnly     ///< Explicitly use portable CPU software rasterizer.
};

/**
 * @class RenderEngineFactory
 * @brief Constructs initialized IRenderEngine instances with fallback handling.
 */
class RenderEngineFactory {
public:
  /**
   * @brief Creates and initializes an IRenderEngine instance.
   * @param preference Preferred backend (Auto, VulkanOnly, or CpuOnly).
   * @param dimensions Initial viewport pixel dimensions.
   * @return Initialized engine pointer, or nullptr if requested backend unavailable.
   */
  [[nodiscard]] static std::unique_ptr<IRenderEngine>
  createEngine(BackendPreference preference = BackendPreference::Auto,
               const ViewportDimensions &dimensions = ViewportDimensions{800, 600});
};

} // namespace correlation::render3d
