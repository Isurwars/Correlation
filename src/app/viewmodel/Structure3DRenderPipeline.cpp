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
          engine
              ? [&engine, this] {
                  static_cast<void>(engine->initialize(viewport_dims_));
                  return std::move(engine);
                }()
              : render3d::RenderEngineFactory::createEngine(render3d::BackendPreference::Auto, viewport_dims_)) {
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
  {
    const std::scoped_lock lock(mutex_);
    viewport_dims_ = dimensions;
    reinit_pending_ = true;
    render_pending_ = true;
  }
  cv_.notify_one();
}

void Structure3DRenderPipeline::resize(const render3d::ViewportDimensions &dimensions) {
  {
    const std::scoped_lock lock(mutex_);
    viewport_dims_ = dimensions;
    render_pending_ = true;
  }
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

std::string Structure3DRenderPipeline::backendName() const {
  return engine_ != nullptr ? std::string(engine_->backendName()) : "None";
}

std::string Structure3DRenderPipeline::deviceName() const {
  return engine_ != nullptr ? engine_->deviceName() : "None";
}

bool Structure3DRenderPipeline::isHardwareAccelerated() const noexcept {
  return engine_ != nullptr && engine_->isHardwareAccelerated();
}

std::vector<std::uint8_t>
Structure3DRenderPipeline::renderSnapshotOffscreen(const RenderPipelineSnapshot &snapshot,
                                                   const render3d::ViewportDimensions &dims,
                                                   const render3d::RgbaColor &bg_color) {
  if (snapshot.scene == nullptr || dims.width == 0 || dims.height == 0) {
    return {};
  }
  auto temp_engine =
      render3d::RenderEngineFactory::createEngine(render3d::BackendPreference::Auto, dims);
  if (temp_engine == nullptr) {
    return {};
  }

  render3d::Camera3D snap_cam = snapshot.camera;
  snap_cam.setViewport(static_cast<float>(dims.width), static_cast<float>(dims.height));

  temp_engine->setFeatureVisibility(snapshot.render_bonds, snapshot.render_wireframe);
  temp_engine->setRenderStyle(render3d::RenderStyleConfig{
      .depth_cueing = snapshot.depth_cueing,
      .ssao_enabled = snapshot.ssao_enabled,
      .background_color = bg_color,
  });
  temp_engine->setCamera(snap_cam);
  temp_engine->setScene(snapshot.scene);

  const auto span = temp_engine->renderFrame();
  return {span.begin(), span.end()};
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
    render3d::ViewportDimensions target_dims{};
    bool need_reinit = false;

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
      need_reinit = std::exchange(reinit_pending_, false);
      target_dims = viewport_dims_;
      snapshot = current_snapshot_;
    }

    if (!snapshot.is_active || snapshot.scene == nullptr || engine_ == nullptr) {
      continue;
    }

    if (need_reinit) {
      static_cast<void>(engine_->initialize(target_dims));
    } else {
      const auto current_dims = engine_->dimensions();
      if (current_dims.width != target_dims.width || current_dims.height != target_dims.height) {
        engine_->resize(target_dims);
      }
    }

    const auto lod_policy =
        lod_controller_.evaluatePolicy(snapshot.scene->atoms().size(), snapshot.is_interacting);
    const bool show_bonds = snapshot.render_bonds && lod_policy.render_bonds;

    engine_->setFeatureVisibility(show_bonds, snapshot.render_wireframe);
    engine_->setRenderStyle(render3d::RenderStyleConfig{
        .depth_cueing = snapshot.depth_cueing,
        .ssao_enabled = snapshot.ssao_enabled,
    });
    engine_->setGpuHaloConfig(render3d::GpuHaloConfig{
        .enabled = snapshot.gpu_halos_enabled,
        .descriptor = snapshot.phase_descriptor,
        .min_threshold = snapshot.phase_min_threshold,
        .max_threshold = snapshot.phase_max_threshold,
        .halo_radius_ratio = snapshot.phase_halo_ratio,
        .halo_intensity = snapshot.phase_halo_intensity,
        .halo_color = snapshot.phase_halo_color,
        .dim_non_matching = snapshot.phase_dim_non_matching,
        .hide_non_matching = snapshot.phase_hide_non_matching,
    });
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
