/**
 * @file CpuRenderEngine.cpp
 * @brief Implementation of CpuRenderEngine offscreen pipeline.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/cpu/CpuRenderEngine.hpp"

namespace correlation::render3d {

CpuRenderEngine::CpuRenderEngine(RasterizerOptions options) noexcept : rasterizer_(options) {}

bool CpuRenderEngine::initialize(const ViewportDimensions &dimensions) {
  resize(dimensions);
  return true;
}

void CpuRenderEngine::resize(const ViewportDimensions &dimensions) {
  dimensions_ = dimensions;
  framebuffer_.resize(dimensions_.width, dimensions_.height);
  camera_.setViewport(static_cast<float>(dimensions_.width),
                      static_cast<float>(dimensions_.height));
}

void CpuRenderEngine::setScene(std::shared_ptr<const Scene3D> scene) { scene_ = std::move(scene); }

void CpuRenderEngine::setCamera(const Camera3D &camera) {
  camera_ = camera;
  camera_.setViewport(static_cast<float>(dimensions_.width),
                      static_cast<float>(dimensions_.height));
}

void CpuRenderEngine::setFeatureVisibility(bool render_bonds, bool render_wireframe) {
  auto opts = rasterizer_.options();
  opts.render_bonds = render_bonds;
  opts.render_wireframe = render_wireframe;
  rasterizer_.setOptions(opts);
}

std::span<const uint8_t> CpuRenderEngine::renderFrame() {
  if (scene_ != nullptr) {
    rasterizer_.render(*scene_, camera_, framebuffer_);
  } else {
    framebuffer_.clear(rasterizer_.options().background_color);
  }
  return framebuffer_.colorBytes();
}

ViewportDimensions CpuRenderEngine::dimensions() const noexcept { return dimensions_; }

std::string_view CpuRenderEngine::backendName() const noexcept {
  return "CPU (Software Rasterizer)";
}

bool CpuRenderEngine::isHardwareAccelerated() const noexcept { return false; }

} // namespace correlation::render3d
