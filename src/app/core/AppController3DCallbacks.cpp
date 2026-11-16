/**
 * @file AppController3DCallbacks.cpp
 * @brief 3D viewport, Camera, Phase Search, and Dual View callback registrations for AppController.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "AppWindow.h"
#include "app/core/AppController.hpp"
#include "app/services/FileIOHandler.hpp"
#include "app/viewmodel/PlotController.hpp"
#include "app/viewmodel/Structure3DController.hpp"

#include <algorithm>
#include <string_view>

namespace correlation::app {

void AppController::registerStructure3DCallbacks() {
  if (!structure_3d_controller_) {
    return;
  }
  registerCameraCallbacks();
  registerPhaseSearchCallbacks();
}

void AppController::registerCameraCallbacks() {
  window_.on_rotate_3d([this](float delta_x, float delta_y) {
    structure_3d_controller_->handleRotate(delta_x, delta_y);
  });

  window_.on_pan_3d([this](float delta_x, float delta_y) {
    structure_3d_controller_->handlePan(delta_x, delta_y);
  });

  window_.on_zoom_3d([this](float delta) { structure_3d_controller_->handleZoom(delta); });

  window_.on_reset_view_3d([this] { structure_3d_controller_->handleResetView(); });

  window_.on_set_camera_preset_3d(
      [this](int preset_idx) { structure_3d_controller_->handleSetCameraPreset(preset_idx); });

  window_.on_toggle_bonds_3d(
      [this](bool enabled) { structure_3d_controller_->handleToggleBonds(enabled); });

  window_.on_toggle_wireframe_3d(
      [this](bool enabled) { structure_3d_controller_->handleToggleWireframe(enabled); });

  window_.on_toggle_in_cell_3d(
      [this](bool enabled) { structure_3d_controller_->handleToggleInCell(enabled); });

  window_.on_seek_frame_3d([this](int frame_idx) {
    structure_3d_controller_->handleSeekFrame(frame_idx);
    if (file_io_handler_) {
      file_io_handler_->updateBoxDiagnosticsForFrame(
          static_cast<size_t>(std::max(0, window_.get_structure_3d_current_frame())));
    }
    if (window_.get_preview_active_tab() == 2 && plot_controller_ != nullptr) {
      plot_controller_->requestPlotUpdate(window_.get_selected_plot_index(), false);
    }
  });

  window_.on_step_frame_3d([this](bool forward) {
    structure_3d_controller_->handleStepFrame(forward);
    if (file_io_handler_) {
      file_io_handler_->updateBoxDiagnosticsForFrame(
          static_cast<size_t>(std::max(0, window_.get_structure_3d_current_frame())));
    }
    if (window_.get_preview_active_tab() == 2 && plot_controller_ != nullptr) {
      plot_controller_->requestPlotUpdate(window_.get_selected_plot_index(), false);
    }
  });

  window_.on_set_structure_3d_dist_metric([this](int metric_idx) {
    if (structure_3d_controller_) {
      structure_3d_controller_->handleSelectDistributionMetric(metric_idx);
    }
  });

  window_.on_toggle_structure_3d_dist_hud([this](bool visible) {
    if (structure_3d_controller_) {
      structure_3d_controller_->handleToggleDistributionHUD(visible);
    }
  });

  window_.on_viewport_3d_resized([this](float width, float height) {
    structure_3d_controller_->handleViewportResized(width, height);
  });
}

void AppController::handleSelectDualLeftPlot(int idx) {
  if (idx == 4) {
    window_.set_preview_dual_left_view(4);
    if (structure_3d_controller_) {
      structure_3d_controller_->setPhase2DMode(true);
      structure_3d_controller_->setPhaseSearchEnabled(true);
    }
    if (plot_controller_ != nullptr) {
      const auto keys = plot_controller_->availablePlotKeys();
      for (std::size_t i = 0; i < keys.size(); ++i) {
        if (keys[i] == "2D Phase Heatmap") {
          window_.set_selected_plot_index(static_cast<int>(i));
          break;
        }
      }
    }
    return;
  }

  window_.set_preview_dual_left_view(idx);
  if (plot_controller_ == nullptr) {
    return;
  }

  std::string_view target_prefix;
  if (idx == 0) {
    target_prefix = "g_r";
  } else if (idx == 1) {
    target_prefix = "S_q";
  } else if (idx == 2) {
    target_prefix = "PAD";
  } else if (idx == 3) {
    target_prefix = "XRD";
  }

  const auto keys = plot_controller_->availablePlotKeys();
  for (std::size_t i = 0; i < keys.size(); ++i) {
    if (keys[i] == target_prefix || keys[i].starts_with(target_prefix)) {
      window_.set_selected_plot_index(static_cast<int>(i));
      plot_controller_->requestPlotUpdate(static_cast<int>(i), true);
      break;
    }
  }
}

void AppController::registerPhaseSearchCallbacks() {
  window_.on_preview_tab_changed([this](int tab) {
    const bool is_3d = (tab == 1 || tab == 2);
    structure_3d_controller_->set3DModeActive(is_3d);
    if (tab == 2) {
      structure_3d_controller_->handleToggleDistributionHUD(false);
      if (plot_controller_ != nullptr) {
        plot_controller_->syncPlotListForActiveTab(2);
      }
    } else if (tab == 1) {
      structure_3d_controller_->handleToggleDistributionHUD(true);
    } else if (tab == 0 && plot_controller_ != nullptr) {
      plot_controller_->syncPlotListForActiveTab(0);
    }
    if ((tab == 0 || tab == 2) && plot_controller_ != nullptr) {
      plot_controller_->requestPlotUpdate(window_.get_selected_plot_index(), true);
    }
  });

  window_.on_select_dual_left_plot([this](int idx) { handleSelectDualLeftPlot(idx); });

  window_.on_toggle_phase_search_3d(
      [this](bool enabled) { structure_3d_controller_->setPhaseSearchEnabled(enabled); });

  window_.on_set_phase_preset_3d([this](int preset_idx) {
    structure_3d_controller_->setCrystalPreset(static_cast<CrystalPreset>(preset_idx));
  });

  window_.on_set_phase_descriptor_3d([this](int desc_idx) {
    structure_3d_controller_->setPhaseDescriptor(static_cast<PhaseDescriptor>(desc_idx));
  });

  window_.on_set_phase_colormap_3d([this](int cmap_idx) {
    structure_3d_controller_->setPhaseColormap(static_cast<plotters::ColormapType>(cmap_idx));
  });

  window_.on_set_phase_thresholds_3d([this](float min_val, float max_val) {
    structure_3d_controller_->setPhaseThresholds(min_val, max_val);
  });

  window_.on_set_phase_intensity_3d(
      [this](float intensity) { structure_3d_controller_->setPhaseHaloIntensity(intensity); });

  window_.on_set_phase_radius_ratio_3d(
      [this](float ratio) { structure_3d_controller_->setPhaseHaloRadiusRatio(ratio); });

  window_.on_set_atom_scale_3d(
      [this](float scale) { structure_3d_controller_->setAtomScale(scale); });

  window_.on_set_phase_dim_non_matching_3d(
      [this](bool dim) { structure_3d_controller_->setPhaseDimNonMatching(dim); });

  window_.on_set_phase_hide_non_matching_3d(
      [this](bool hide) { structure_3d_controller_->setPhaseHideNonMatching(hide); });

  window_.on_set_phase_2d_mode_3d(
      [this](bool is_2d) { structure_3d_controller_->setPhase2DMode(is_2d); });

  window_.on_set_phase_descriptor_y_3d([this](int desc_idx) {
    structure_3d_controller_->setPhaseDescriptorY(static_cast<PhaseDescriptor>(desc_idx));
  });

  window_.on_set_phase_thresholds_2d_3d([this](float min_x, float max_x, float min_y, float max_y) {
    structure_3d_controller_->setPhaseThresholds2D(min_x, max_x, min_y, max_y);
  });

  window_.on_export_phase_heatmap_png_3d(
      [this] { structure_3d_controller_->exportHeatmapPng("phase_heatmap_2d.png"); });
}

} // namespace correlation::app
