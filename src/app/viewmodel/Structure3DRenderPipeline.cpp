/**
 * @file Structure3DRenderPipeline.cpp
 * @brief Implementation of Structure3DRenderPipeline for background 3D offscreen rendering.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/viewmodel/Structure3DRenderPipeline.hpp"
#include "render3d/RenderEngineFactory.hpp"

#include <chrono>
#include <utility>

namespace correlation::app {

Structure3DRenderPipeline::Structure3DRenderPipeline(
    FrameDispatchedCallback callback, std::shared_ptr<std::atomic<bool>> alive,
    std::unique_ptr<render3d::IRenderEngine> engine)
    : callback_(std::move(callback)), alive_(std::move(alive)),
      engine_(
          engine ? std::move(engine)
                 : render3d::RenderEngineFactory::createEngine(render3d::BackendPreference::Auto)) {
  if (engine_ != nullptr) {
    static_cast<void>(engine_->initialize(viewport_dims_));
  }
  render_worker_ = std::jthread([this] { renderLoop(); });
}

Structure3DRenderPipeline::~Structure3DRenderPipeline() {
  {
    const std::scoped_lock lock(mutex_);
    stop_worker_.store(true, std::memory_order_release);
  }
  cv_.notify_all();
  if (render_worker_.joinable()) {
    render_worker_.join();
  }
}

void Structure3DRenderPipeline::initialize(const render3d::ViewportDimensions &dimensions) {
  const std::scoped_lock lock(mutex_);
  viewport_dims_ = dimensions;
  if (engine_ != nullptr) {
    static_cast<void>(engine_->initialize(viewport_dims_));
  }
}

void Structure3DRenderPipeline::resize(const render3d::ViewportDimensions &dimensions) {
  const std::scoped_lock lock(mutex_);
  viewport_dims_ = dimensions;
  if (engine_ != nullptr) {
    engine_->resize(viewport_dims_);
  }
  render_pending_ = true;
  cv_.notify_one();
}

void Structure3DRenderPipeline::updateSnapshot(const RenderPipelineSnapshot &snapshot,
                                               bool /*immediate*/) {
  {
    const std::scoped_lock lock(mutex_);
    current_snapshot_ = snapshot;
    render_pending_ = true;
  }
  cv_.notify_one();
}

void Structure3DRenderPipeline::requestRender(bool /*immediate*/) {
  {
    const std::scoped_lock lock(mutex_);
    render_pending_ = true;
  }
  cv_.notify_one();
}

render3d::ViewportDimensions Structure3DRenderPipeline::dimensions() const noexcept {
  const std::scoped_lock lock(mutex_);
  return viewport_dims_;
}

bool Structure3DRenderPipeline::isActive() const noexcept {
  const std::scoped_lock lock(mutex_);
  return current_snapshot_.is_active;
}

void Structure3DRenderPipeline::setActive(bool active) noexcept {
  {
    const std::scoped_lock lock(mutex_);
    current_snapshot_.is_active = active;
    if (active) {
      render_pending_ = true;
    }
  }
  if (active) {
    cv_.notify_one();
  }
}

void Structure3DRenderPipeline::dispatchFrame(const std::uint8_t *pixel_data, std::uint32_t width,
                                              std::uint32_t height) {
  if (pixel_data == nullptr || width == 0 || height == 0 || !callback_) {
    return;
  }

  const auto *rgba_pixels = reinterpret_cast<const slint::Rgba8Pixel *>(pixel_data);
  const slint::SharedPixelBuffer<slint::Rgba8Pixel> buffer(width, height, rgba_pixels);
  slint::Image img(buffer);

  slint::invoke_from_event_loop(
      [alive = alive_, callback = callback_, img = std::move(img)] mutable {
        if (alive != nullptr && !alive->load()) {
          return;
        }
        callback(img);
      });
}

void Structure3DRenderPipeline::renderLoop() {
  while (!stop_worker_.load(std::memory_order_relaxed)) {
    RenderPipelineSnapshot snapshot;

    {
      std::unique_lock lock(mutex_);
      cv_.wait_for(lock, std::chrono::milliseconds(50), [this] {
        return stop_worker_.load(std::memory_order_relaxed) ||
               (render_pending_ && current_snapshot_.is_active);
      });

      if (stop_worker_.load(std::memory_order_relaxed)) {
        break;
      }

      render_pending_ = false;
      snapshot = current_snapshot_;
    }

    if (!snapshot.is_active || snapshot.scene == nullptr || engine_ == nullptr) {
      continue;
    }

    const auto lod_policy =
        lod_controller_.evaluatePolicy(snapshot.scene->atoms().size(), snapshot.is_interacting);
    const bool show_bonds = snapshot.render_bonds && lod_policy.render_bonds;

    engine_->setFeatureVisibility(show_bonds, snapshot.render_wireframe);
    engine_->setCamera(snapshot.camera);
    engine_->setScene(snapshot.scene);

    const auto pixel_span = engine_->renderFrame();
    if (!pixel_span.empty()) {
      const auto dims = engine_->dimensions();
      dispatchFrame(pixel_span.data(), dims.width, dims.height);
    }
  }
}

} // namespace correlation::app
