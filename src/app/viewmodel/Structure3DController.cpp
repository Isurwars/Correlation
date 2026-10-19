/**
 * @file Structure3DController.cpp
 * @brief Implementation of Structure3DController for interactive Slint 3D visualization.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/viewmodel/Structure3DController.hpp"
#include "AppWindow.h"
#include "render3d/RenderEngineFactory.hpp"

#include <chrono>
#include <cmath>
#include <utility>

namespace correlation::app {

Structure3DController::Structure3DController(::AppWindow &window, const ProgramOptions &options)
    : window_(window), options_(options),
      engine_(render3d::RenderEngineFactory::createEngine(render3d::BackendPreference::Auto)),
      alive_(std::make_shared<std::atomic<bool>>(true)) {
  if (engine_ != nullptr) {
    static_cast<void>(engine_->initialize(viewport_dims_));
  }
  camera_.setViewport(static_cast<float>(viewport_dims_.width),
                      static_cast<float>(viewport_dims_.height));
  camera_.fitToBox(math::Vector3<float>{-10.0F, -10.0F, -10.0F},
                   math::Vector3<float>{10.0F, 10.0F, 10.0F});

  render_worker_ = std::jthread([this] { renderLoop(); });
}

Structure3DController::~Structure3DController() {
  if (alive_ != nullptr) {
    alive_->store(false);
  }
  {
    const std::scoped_lock lock(mutex_);
    stop_worker_.store(true);
  }
  cv_.notify_all();
  if (render_worker_.joinable()) {
    render_worker_.join();
  }
}

void Structure3DController::requestRender(bool /*immediate*/) {
  {
    const std::scoped_lock lock(mutex_);
    render_pending_ = true;
  }
  cv_.notify_one();
}

void Structure3DController::handleViewportResized(float width, float height) {
  if (width < 10.0F || height < 10.0F) {
    return;
  }
  const auto new_width = static_cast<uint32_t>(std::round(width));
  const auto new_height = static_cast<uint32_t>(std::round(height));

  {
    const std::scoped_lock lock(mutex_);
    if (viewport_dims_.width == new_width && viewport_dims_.height == new_height) {
      return;
    }
    viewport_dims_ = {.width = new_width, .height = new_height};
    camera_.setViewport(static_cast<float>(new_width), static_cast<float>(new_height));
    if (engine_ != nullptr) {
      engine_->resize(viewport_dims_);
    }
    render_pending_ = true;
  }
  cv_.notify_one();
}

void Structure3DController::handleRotate(float delta_x, float delta_y) {
  {
    const std::scoped_lock lock(mutex_);
    camera_.orbit(delta_x, delta_y);
    render_pending_ = true;
  }
  cv_.notify_one();
}

void Structure3DController::handlePan(float delta_x, float delta_y) {
  {
    const std::scoped_lock lock(mutex_);
    camera_.pan(delta_x, delta_y);
    render_pending_ = true;
  }
  cv_.notify_one();
}

void Structure3DController::handleZoom(float delta) {
  {
    const std::scoped_lock lock(mutex_);
    camera_.zoom(delta);
    render_pending_ = true;
  }
  cv_.notify_one();
}

void Structure3DController::handleResetView() {
  {
    const std::scoped_lock lock(mutex_);
    if (scene_ != nullptr) {
      camera_.fitToBox(scene_->boundsMin(), scene_->boundsMax());
    } else {
      camera_.fitToBox(math::Vector3<float>{-10.0F, -10.0F, -10.0F},
                       math::Vector3<float>{10.0F, 10.0F, 10.0F});
    }
    render_pending_ = true;
  }
  cv_.notify_one();
}

void Structure3DController::handleToggleBonds(bool enabled) {
  {
    const std::scoped_lock lock(mutex_);
    render_bonds_ = enabled;
    if (render_bonds_ && scene_ != nullptr) {
      scene_->recomputeBonds(bond_cutoff_scale_);
    }
    render_pending_ = true;
  }
  window_.set_structure_3d_bonds_enabled(enabled);
  cv_.notify_one();
}

void Structure3DController::handleToggleWireframe(bool enabled) {
  {
    const std::scoped_lock lock(mutex_);
    render_wireframe_ = enabled;
    render_pending_ = true;
  }
  window_.set_structure_3d_wireframe_enabled(enabled);
  cv_.notify_one();
}

void Structure3DController::handleSeekFrame(int frame_index) {
  if (frame_index < 0) {
    return;
  }
  int actual_frame = frame_index;
  {
    const std::scoped_lock lock(mutex_);
    if (streamer_.totalFrames() > 0 && scene_ != nullptr) {
      const auto target = std::min(static_cast<size_t>(frame_index), streamer_.totalFrames() - 1);
      static_cast<void>(streamer_.seekFrame(target, *scene_, render_bonds_, bond_cutoff_scale_));
      actual_frame = static_cast<int>(streamer_.currentFrameIndex());
      render_pending_ = true;
    }
  }
  window_.set_structure_3d_current_frame(actual_frame);
  cv_.notify_one();
}

void Structure3DController::handleStepFrame(bool forward) {
  int current_idx = 0;
  {
    const std::scoped_lock lock(mutex_);
    if (streamer_.totalFrames() > 0 && scene_ != nullptr) {
      if (forward) {
        static_cast<void>(streamer_.stepNext(*scene_, true, render_bonds_, bond_cutoff_scale_));
      } else {
        static_cast<void>(streamer_.stepPrev(*scene_, true, render_bonds_, bond_cutoff_scale_));
      }
      current_idx = static_cast<int>(streamer_.currentFrameIndex());
      render_pending_ = true;
    }
  }
  window_.set_structure_3d_current_frame(current_idx);
  cv_.notify_one();
}

void Structure3DController::updateBondCutoffScale(double cutoff_scale) {
  if (cutoff_scale <= 0.0) {
    return;
  }
  {
    const std::scoped_lock lock(mutex_);
    bond_cutoff_scale_ = cutoff_scale;
    if (scene_ != nullptr && render_bonds_) {
      scene_->recomputeBonds(bond_cutoff_scale_);
      render_pending_ = true;
    }
  }
  cv_.notify_one();
}

void Structure3DController::setCell(const core::Cell &cell) {
  auto new_scene = render3d::Scene3D::fromCell(cell, bond_cutoff_scale_, render_bonds_);
  {
    const std::scoped_lock lock(mutex_);
    scene_ = std::move(new_scene);
    if (scene_ != nullptr) {
      camera_.fitToBox(scene_->boundsMin(), scene_->boundsMax());
      if (engine_ != nullptr) {
        engine_->setScene(scene_);
      }
    }
    render_pending_ = true;
  }
  window_.set_structure_3d_current_frame(0);
  window_.set_structure_3d_total_frames(1);
  cv_.notify_one();
}

void Structure3DController::setTrajectory(std::shared_ptr<const core::Trajectory> trajectory) {
  if (!trajectory) {
    return;
  }
  const auto total_frames = static_cast<int>(trajectory->getFrameCount());
  {
    const std::scoped_lock lock(mutex_);
    if (total_frames > 0) {
      const auto &first_cell = trajectory->firstFrame();
      scene_ = render3d::Scene3D::fromCell(first_cell, bond_cutoff_scale_, render_bonds_);
      if (scene_ != nullptr) {
        camera_.fitToBox(scene_->boundsMin(), scene_->boundsMax());
        if (engine_ != nullptr) {
          engine_->setScene(scene_);
        }
      }
    }
    streamer_.setTrajectory(std::move(trajectory));
    render_pending_ = true;
  }
  window_.set_structure_3d_current_frame(0);
  window_.set_structure_3d_total_frames(total_frames);
  cv_.notify_one();
}

void Structure3DController::set3DModeActive(bool active) {
  {
    const std::scoped_lock lock(mutex_);
    if (is_active_ == active) {
      return;
    }
    is_active_ = active;
    if (is_active_) {
      render_pending_ = true;
    }
  }
  if (active) {
    cv_.notify_one();
  }
}

void Structure3DController::dispatchFrame(const uint8_t *pixel_data, uint32_t width,
                                          uint32_t height) {
  if (pixel_data == nullptr || width == 0 || height == 0) {
    return;
  }

  const auto *rgba_pixels = reinterpret_cast<const slint::Rgba8Pixel *>(pixel_data);
  const slint::SharedPixelBuffer<slint::Rgba8Pixel> buffer(width, height, rgba_pixels);
  slint::Image img(buffer);

  slint::invoke_from_event_loop([this, alive = alive_, img = std::move(img)] mutable {
    if (!alive->load()) {
      return;
    }
    window_.set_structure_3d_image(img);
  });
}

void Structure3DController::renderLoop() {
  while (!stop_worker_.load(std::memory_order_relaxed)) {
    render3d::Camera3D cam_copy;
    std::shared_ptr<render3d::Scene3D> scene_copy;
    bool active = false;
    bool bonds_vis = true;
    bool wireframe_vis = true;
    bool interacting = false;

    {
      std::unique_lock lock(mutex_);
      cv_.wait_for(lock, std::chrono::milliseconds(50), [this] {
        return stop_worker_.load(std::memory_order_relaxed) || (render_pending_ && is_active_);
      });

      if (stop_worker_.load(std::memory_order_relaxed)) {
        break;
      }

      render_pending_ = false;
      active = is_active_;
      cam_copy = camera_;
      scene_copy = scene_;
      bonds_vis = render_bonds_;
      wireframe_vis = render_wireframe_;
      interacting = is_interacting_;
    }

    if (!active || scene_copy == nullptr || engine_ == nullptr) {
      continue;
    }

    const auto lod_policy = lod_controller_.evaluatePolicy(scene_copy->atoms().size(), interacting);
    const bool show_bonds = bonds_vis && lod_policy.render_bonds;

    engine_->setFeatureVisibility(show_bonds, wireframe_vis);
    engine_->setCamera(cam_copy);
    engine_->setScene(scene_copy);

    const auto pixel_span = engine_->renderFrame();
    if (!pixel_span.empty()) {
      const auto dims = engine_->dimensions();
      dispatchFrame(pixel_span.data(), dims.width, dims.height);
    }
  }
}

} // namespace correlation::app
