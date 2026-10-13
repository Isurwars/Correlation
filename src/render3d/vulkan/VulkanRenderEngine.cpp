/**
 * @file VulkanRenderEngine.cpp
 * @brief Implementation of VulkanRenderEngine headless pipeline.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/vulkan/VulkanRenderEngine.hpp"
#include "render3d/vulkan/VulkanContext.hpp"

namespace correlation::render3d {

VulkanRenderEngine::~VulkanRenderEngine() = default;

bool VulkanRenderEngine::initialize(const ViewportDimensions &dimensions) {
  if (!VulkanContext::isSupported()) {
    initialized_ = false;
    return false;
  }

  dimensions_ = dimensions;
  const std::size_t pixel_bytes =
      static_cast<std::size_t>(dimensions_.width) * dimensions_.height * 4;
  frame_pixels_.resize(pixel_bytes, 0);
  initialized_ = true;
  return true;
}

void VulkanRenderEngine::resize(const ViewportDimensions &dimensions) {
  dimensions_ = dimensions;
  const std::size_t pixel_bytes =
      static_cast<std::size_t>(dimensions_.width) * dimensions_.height * 4;
  frame_pixels_.resize(pixel_bytes, 0);
}

void VulkanRenderEngine::setScene(std::shared_ptr<const Scene3D> scene) {
  scene_ = std::move(scene);
}

void VulkanRenderEngine::setCamera(const Camera3D &camera) { camera_ = camera; }

void VulkanRenderEngine::setFeatureVisibility(bool /*render_bonds*/, bool /*render_wireframe*/) {}

std::span<const uint8_t> VulkanRenderEngine::renderFrame() { return frame_pixels_; }

ViewportDimensions VulkanRenderEngine::dimensions() const noexcept { return dimensions_; }

std::string_view VulkanRenderEngine::backendName() const noexcept {
  return "Vulkan (Headless Offscreen)";
}

bool VulkanRenderEngine::isHardwareAccelerated() const noexcept { return true; }

} // namespace correlation::render3d
