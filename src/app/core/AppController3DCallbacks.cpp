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
#include <array>
#include <format>
#include <nfd.h>
#include <string>
#include <string_view>

namespace correlation::app {

namespace {

void handleExportTrajectoryAnimation(const AppWindow &window, Structure3DController &controller,
                                     bool as_gif, int bg_mode) {
  const std::array<nfdfilteritem_t, 1> filters =
      as_gif ? std::array<nfdfilteritem_t, 1>{{{.name = "Animated GIF", .spec = "gif"}}}
             : std::array<nfdfilteritem_t, 1>{{{.name = "PNG Sequence ZIP", .spec = "zip"}}};
  nfdu8char_t *out_path = nullptr;
  const std::string base = std::string(window.get_file_basename());
  const std::string ext = as_gif ? ".gif" : "_frames.zip";
  const std::string default_name = (base.empty() ? "trajectory" : base) + ext;
  const nfdresult_t res =
      NFD_SaveDialogU8(&out_path, filters.data(), filters.size(), nullptr, default_name.c_str());
  if (res != NFD_OKAY || out_path == nullptr) {
    return;
  }
  std::string target_path(out_path);
  NFD_FreePathU8(out_path);
  const std::string suffix = as_gif ? ".gif" : ".zip";
  if (!target_path.ends_with(suffix)) {
    target_path += suffix;
  }
  std::size_t exported_frames = 0;
  const bool exported =
      controller.exportTrajectoryAnimation(target_path, as_gif, 15, bg_mode, &exported_frames);
  if (exported) {
    window.invoke_show_toast(
        slint::SharedString(
            std::format("Exported {} frames ({})", exported_frames, as_gif ? "GIF" : "PNG ZIP")),
        1);
  } else {
    window.invoke_show_toast(slint::SharedString("Failed to export trajectory animation"), 2);
  }
}

void handleExportSnapshot(const AppWindow &window, Structure3DController &controller, int scale_idx,
                          int bg_mode) {
  constexpr std::array<nfdfilteritem_t, 1> K_FILTERS{{{.name = "PNG Image", .spec = "png"}}};
  nfdu8char_t *out_path = nullptr;
  const std::string base = std::string(window.get_file_basename());
  const std::string default_name = base.empty() ? "structure_3d.png" : (base + "_3d.png");
  const nfdresult_t res = NFD_SaveDialogU8(&out_path, K_FILTERS.data(), K_FILTERS.size(), nullptr,
                                           default_name.c_str());
  if (res != NFD_OKAY || out_path == nullptr) {
    return;
  }
  std::string target_path(out_path);
  NFD_FreePathU8(out_path);
  if (!target_path.ends_with(".png") && !target_path.ends_with(".PNG")) {
    target_path += ".png";
  }
  std::uint32_t width = 0;
  std::uint32_t height = 0;
  const bool exported =
      controller.exportSnapshotPng(target_path, scale_idx, bg_mode, &width, &height);
  if (exported) {
    window.invoke_show_toast(
        slint::SharedString(std::format("Saved 3D snapshot ({}×{})", width, height)), 1);
  } else {
    window.invoke_show_toast(slint::SharedString("Failed to export 3D snapshot"), 2);
  }
}

} // namespace

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

  window_.on_toggle_projection_3d(
      [this](bool enabled) { structure_3d_controller_->handleToggleProjection(enabled); });

  window_.on_toggle_depth_cueing_3d(
      [this](bool enabled) { structure_3d_controller_->handleToggleDepthCueing(enabled); });

  window_.on_toggle_ssao_3d(
      [this](bool enabled) { structure_3d_controller_->handleToggleSsao(enabled); });

  window_.on_set_visual_supercell_3d([this](int rep_a, int rep_b, int rep_c) {
    structure_3d_controller_->handleSetVisualSupercell(rep_a, rep_b, rep_c);
  });

  window_.on_reset_visual_supercell_3d(
      [this] { structure_3d_controller_->handleResetVisualSupercell(); });

  window_.on_export_trajectory_animation_3d([this](bool as_gif, int bg_mode) {
    handleExportTrajectoryAnimation(window_, *structure_3d_controller_, as_gif, bg_mode);
  });

  window_.on_focus_picked_atom_3d([this] { structure_3d_controller_->handleFocusPickedAtom(); });

  window_.on_tick_camera_animation_3d([this] { structure_3d_controller_->tickCameraAnimation(); });

  window_.on_toggle_element_visibility_3d([this](const slint::SharedString &symbol) {
    structure_3d_controller_->handleToggleElementVisibility(std::string(symbol));
  });

  window_.on_show_all_elements_3d([this] { structure_3d_controller_->handleShowAllElements(); });

  window_.on_export_snapshot_3d([this](int scale_idx, int bg_mode) {
    handleExportSnapshot(window_, *structure_3d_controller_, scale_idx, bg_mode);
  });

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

  window_.on_pick_atom_3d([this](float norm_x, float norm_y) {
    structure_3d_controller_->handlePickAtom(norm_x, norm_y);
  });

  window_.on_clear_picked_atom_3d([this] { structure_3d_controller_->clearPickedAtom(); });

  window_.on_measure_atom_3d([this](float norm_x, float norm_y) {
    structure_3d_controller_->handleMeasureAtom(norm_x, norm_y);
  });

  window_.on_clear_measurement_3d([this] { structure_3d_controller_->clearMeasurement(); });

  window_.on_set_slab_clip_3d([this](const SlabClipState &state) {
    structure_3d_controller_->setSlabClipConfig(SlabClipConfig{
        .enabled = state.enabled,
        .x_min = state.x_min,
        .x_max = state.x_max,
        .y_min = state.y_min,
        .y_max = state.y_max,
        .z_min = state.z_min,
        .z_max = state.z_max,
    });
  });

  window_.on_reset_slab_clip_3d([this] { structure_3d_controller_->resetSlabClip(); });
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
