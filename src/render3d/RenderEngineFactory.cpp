/**
 * @file RenderEngineFactory.cpp
 * @brief Implementation of RenderEngineFactory with automatic fallback logic.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/RenderEngineFactory.hpp"
#include "render3d/cpu/CpuRenderEngine.hpp"

#ifdef CORRELATION_USE_VULKAN
#include "render3d/vulkan/VulkanRenderEngine.hpp"
#endif

namespace correlation::render3d {

std::unique_ptr<IRenderEngine>
RenderEngineFactory::createEngine(BackendPreference preference,
                                  const ViewportDimensions &dimensions) {
#ifdef CORRELATION_USE_VULKAN
  if (preference == BackendPreference::Auto || preference == BackendPreference::VulkanOnly) {
    auto vulkan_engine = std::make_unique<VulkanRenderEngine>();
    if (vulkan_engine->initialize(dimensions)) {
      return vulkan_engine;
    }
    if (preference == BackendPreference::VulkanOnly) {
      return nullptr;
    }
  }
#else
  if (preference == BackendPreference::VulkanOnly) {
    return nullptr;
  }
#endif

  auto cpu_engine = std::make_unique<CpuRenderEngine>();
  [[maybe_unused]] const bool initialized = cpu_engine->initialize(dimensions);
  return cpu_engine;
}

} // namespace correlation::render3d
