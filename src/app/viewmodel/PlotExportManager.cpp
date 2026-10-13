/**
 * @file PlotExportManager.cpp
 * @brief Implementation of PlotExportManager.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/viewmodel/PlotExportManager.hpp"
#include "app/core/AppOptions.hpp"
#include "app/formatters/PlotExportService.hpp"
#include "plotters/SvgPlotter.hpp"

#include <nfd.h>
#include <slint.h>

#include <array>
#include <filesystem>
#include <utility>
#include <vector>

namespace correlation::app {

PlotExportManager::PlotExportManager() = default;

PlotExportManager::~PlotExportManager() {
  if (dialog_thread_.joinable()) {
    dialog_thread_.join();
  }
}

bool PlotExportManager::isDialogActive() const noexcept { return dialog_active_.load(); }

std::expected<void, std::string> PlotExportManager::savePlot(
    const std::string &filepath, const std::string &plot_name,
    const correlation::analysis::Histogram *hist, correlation::plotters::PlotConfig config,
    std::span<const PinnedRun> pinned_runs, const std::map<std::string, real_t> &ashcroft_weights) {
  if (hist == nullptr) {
    return std::unexpected("Histogram is null");
  }

  config.use_native_text = true;
  if (config.preset_size == correlation::plotters::PlotConfig::PresetSize::Default) {
    config.width = 1200.0;
    config.height = 900.0;
  }

  if (pinned_runs.empty()) {
    return PlotExportService::exportHistogram(filepath, *hist, config, ashcroft_weights);
  }

  std::vector<correlation::plotters::LabeledHistogram> datasets;
  datasets.reserve(pinned_runs.size() + 1);
  datasets.push_back({.label = "Current", .hist = hist});
  for (const auto &pinned_run : pinned_runs) {
    const auto hist_it = pinned_run.histograms.find(plot_name);
    if (hist_it != pinned_run.histograms.end()) {
      datasets.push_back({.label = pinned_run.label, .hist = &hist_it->second});
    }
  }
  const std::string comp_key = PlotExportService::getComparisonKey(hist);
  return PlotExportService::exportComparison(filepath, datasets, comp_key, config);
}

void PlotExportManager::exportPlotAsync(const std::string &default_base_path,
                                        const std::string &plot_name,
                                        const correlation::analysis::Histogram *hist,
                                        correlation::plotters::PlotConfig config,
                                        std::vector<PinnedRun> pinned_runs,
                                        std::map<std::string, real_t> ashcroft_weights,
                                        StatusCallback on_status) {
  if (dialog_active_.exchange(true)) {
    return;
  }

  if (hist == nullptr) {
    dialog_active_.store(false);
    return;
  }

  std::string default_path = default_base_path;
  if (!default_path.empty()) {
    default_path += "_" + plot_name;
  } else {
    default_path = plot_name;
  }

  const std::filesystem::path def_path(default_path);
  std::string default_dir = def_path.parent_path().string();
  std::string default_name = def_path.filename().string();

  if (dialog_thread_.joinable()) {
    dialog_thread_.join();
  }

  dialog_thread_ = std::jthread([this, default_dir = std::move(default_dir),
                                 default_name = std::move(default_name), plot_name, hist,
                                 config = std::move(config), pinned_runs = std::move(pinned_runs),
                                 ashcroft_weights = std::move(ashcroft_weights),
                                 on_status = std::move(on_status)] mutable {
    std::array<nfdfilteritem_t, 2> filter_list = {
        {
            {
                .name = "SVG Image",
                .spec = "svg",
            },
            {
                .name = "PDF Document",
                .spec = "pdf",
            },
        },
    };
    const nfdfiltersize_t filter_count = filter_list.size();

    nfdchar_t *out_path = nullptr;
    const nfdresult_t result =
        NFD_SaveDialogU8(&out_path, filter_list.data(), filter_count,
                         default_dir.empty() ? nullptr : default_dir.c_str(),
                         default_name.empty() ? nullptr : default_name.c_str());

    if (result == NFD_OKAY) {
      std::string filepath(out_path);
      NFD_FreePathU8(out_path);
      slint::invoke_from_event_loop([this, filepath = std::move(filepath), plot_name, hist, config,
                                     pinned_runs = std::move(pinned_runs),
                                     ashcroft_weights = std::move(ashcroft_weights), on_status] {
        dialog_active_.store(false);
        const auto save_res =
            savePlot(filepath, plot_name, hist, config, pinned_runs, ashcroft_weights);
        if (save_res.has_value()) {
          on_status(plot_name + " plot saved successfully.");
        } else {
          on_status("Failed to save plot: " + save_res.error());
        }
      });
    } else if (result == NFD_CANCEL) {
      slint::invoke_from_event_loop([this, on_status] {
        dialog_active_.store(false);
        on_status(AppDefaults::MSG_SAVE_CANCELLED);
      });
    } else {
      std::string error_msg = "Error: ";
      error_msg += NFD_GetError();
      slint::invoke_from_event_loop([this, error_msg = std::move(error_msg), on_status] {
        dialog_active_.store(false);
        on_status(error_msg);
      });
    }
  });
}

} // namespace correlation::app
