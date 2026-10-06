/**
 * @file FileIOHandler.cpp
 * @brief Implementation of FileIOHandler.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/services/FileIOHandler.hpp"
#include "AppWindow.h"
#include "app/core/AppController.hpp"
#include "app/services/InputValidator.hpp"
#include "physics/PhysicalData.hpp"
#include <filesystem>
#include <format>
#include <nfd.h>

namespace correlation::app {

FileIOHandler::FileIOHandler(::AppWindow &window, TrajectoryLoader &loader,
                             AnalysisDispatcher &dispatcher, ProgramOptions &options,
                             AppController &controller)
    : window_(window), loader_(loader), dispatcher_(dispatcher), options_(options),
      controller_(controller) {}

FileIOHandler::~FileIOHandler() {
  if (dialog_thread_.joinable()) {
    dialog_thread_.join();
  }
  if (load_thread_.joinable()) {
    load_thread_.join();
  }
}

void FileIOHandler::executeWriteFiles(const std::string &filepath) {
  ProgramOptions opts = controller_.handleOptionsfromUI();
  opts.output_file_base = filepath;

  const auto format_opts = window_.get_export_format_options();
  opts.use_zip = format_opts.use_zip;
  opts.use_csv = format_opts.export_csv;
  opts.use_hdf5 = format_opts.export_hdf5;
  opts.use_parquet = format_opts.export_arrow;
  opts.export_images = format_opts.export_images;

  opts.export_algorithms.clear();
  const auto export_algos_model = window_.get_export_analysis_items();
  if (export_algos_model) {
    for (size_t i = 0; i < export_algos_model->row_count(); ++i) {
      const auto maybe_item = export_algos_model->row_data(i);
      if (maybe_item && maybe_item->enabled) {
        opts.export_algorithms.emplace_back(maybe_item->id.data());
      }
    }
  }

  options_ = opts;

  const auto write_res = dispatcher_.writeFiles(options_);
  if (write_res) {
    window_.set_analysis_status_text(slint::SharedString(AppDefaults::MSG_FILES_WRITTEN));
  } else {
    window_.set_analysis_status_text(slint::SharedString(write_res.error()));
  }
}

void FileIOHandler::handleWriteFiles() {
  if (dialog_active_.exchange(true)) {
    return;
  }

  if (dialog_thread_.joinable()) {
    dialog_thread_.join();
  }
  const auto format_opts = window_.get_export_format_options();

#ifndef CORRELATION_USE_HDF5
  if (format_opts.export_hdf5) {
    dialog_active_.store(false);
    window_.set_analysis_status_text("Error: HDF5 format is not available.");
    return;
  }
#endif

#ifndef CORRELATION_USE_ARROW
  if (format_opts.export_arrow) {
    dialog_active_.store(false);
    window_.set_analysis_status_text("Error: Parquet format is not available.");
    return;
  }
#endif

  const bool is_zip = format_opts.use_zip;

  dialog_thread_ = std::jthread([this, is_zip]() {
    nfdchar_t *out_path = nullptr;
    nfdresult_t result = NFD_CANCEL;

    if (is_zip) {
      std::array<nfdfilteritem_t, 1> filter_list = {{{
          .name = "Consolidated ZIP Bundle",
          .spec = "zip",
      }}};
      result = NFD_SaveDialogU8(&out_path, filter_list.data(), filter_list.size(), nullptr,
                                "correlation_export.zip");
    } else {
      result = NFD_PickFolderU8(&out_path, nullptr);
    }

    if (result == NFD_OKAY) {
      std::string filepath(out_path);
      NFD_FreePathU8(out_path);

      slint::invoke_from_event_loop([this, filepath = std::move(filepath)]() {
        dialog_active_.store(false);
        executeWriteFiles(filepath);
      });
    } else if (result == NFD_CANCEL) {
      slint::invoke_from_event_loop([this]() {
        dialog_active_.store(false);
        window_.set_analysis_status_text(slint::SharedString(AppDefaults::MSG_SAVE_CANCELLED));
      });
    } else {
      std::string error_msg = "Error: ";
      error_msg += NFD_GetError();
      slint::invoke_from_event_loop([this, error_msg = std::move(error_msg)]() {
        dialog_active_.store(false);
        window_.set_analysis_status_text(slint::SharedString(error_msg));
      });
    }
  });
}

void FileIOHandler::startLoadingTrajectory(const std::string &filepath) {
  window_.set_in_file_text(slint::SharedString(filepath));
  window_.set_file_status_text(slint::SharedString("Loading file..."));
  window_.set_file_loading(true);
  window_.set_timer_running(true);
  window_.set_text_opacity(true);
  window_.set_progress(0.0F);

  if (load_thread_.joinable()) {
    load_thread_.join();
  }

  load_thread_ = std::jthread([this, filepath]() {
    auto progress_cb = [this](float progress, const std::string &msg) {
      slint::invoke_from_event_loop([progress, msg, this]() {
        window_.set_progress(progress);
        if (!msg.empty()) {
          window_.set_file_status_text(slint::SharedString(msg));
        }
      });
    };

    std::string message;
    bool success = false;
    const auto load_res = loader_.loadFile(filepath, progress_cb);
    if (load_res) {
      message = *load_res;
      success = true;
      options_.input_file = filepath;
    } else {
      message = std::string(AppDefaults::MSG_ERROR_LOADING) + load_res.error();
    }

    slint::invoke_from_event_loop([this, filepath, message, success]() {
      window_.set_file_status_text(slint::SharedString(message));
      window_.set_file_loading(false);
      window_.set_timer_running(false);
      window_.set_text_opacity(false);

      if (success && loader_.cell() != nullptr) {
        window_.set_file_loaded(true);

        // Atom Counts
        auto atom_counts_map = loader_.getAtomCounts();
        auto slint_atom_counts = std::make_shared<slint::VectorModel<AtomCount>>();
        for (const auto &[symbol, count] : atom_counts_map) {
          slint_atom_counts->push_back({
              .symbol = slint::SharedString(symbol),
              .count = count,
          });
        }
        window_.set_atom_counts(slint_atom_counts);

        // Bond Cutoffs
        controller_.setBondCutoffs();

        // File Info & Basename
        updateFileMetadata(filepath);

        // Box Diagnostics
        updateBoxDiagnostics(loader_.cell());

        {
          auto opts = window_.get_analysis_options();
          opts.time_step =
              slint::SharedString(std::format("{:.2f}", loader_.getRecommendedTimeStep()));
          window_.set_analysis_options(opts);
        }

        // Update Run Analysis Card Frame Info
        {
          auto opts = window_.get_analysis_options();
          opts.min_frame = "1";
          opts.frame_stride = "1";
          window_.set_analysis_options(opts);
        }
        {
          auto opts = window_.get_analysis_options();
          opts.max_frame = slint::SharedString(std::to_string(loader_.getFrameCount()));
          window_.set_analysis_options(opts);
        }
        static_cast<void>(controller_.getInputValidator()->validateInputs());
      }
    });
  });
}

void FileIOHandler::handleReloadFile() {
  const std::string &path = options_.input_file;
  if (!path.empty()) {
    startLoadingTrajectory(path);
  }
}

void FileIOHandler::updateBoxDiagnostics(const core::Cell *cell) {
  if (cell == nullptr) {
    return;
  }
  const real_t vol = cell->volume();
  const std::string vol_str = (vol > 0.0) ? std::format("{:.2f} Å³", vol) : "N/A";
  window_.set_box_volume(slint::SharedString(vol_str));

  const auto &params = cell->latticeParameters();
  const std::string dims_str =
      std::format("a: {:.2f}  b: {:.2f}  c: {:.2f} Å", params[0], params[1], params[2]);
  window_.set_box_dimensions(slint::SharedString(dims_str));

  real_t total_mass = 0.0;
  for (const auto &atom : cell->atoms()) {
    const auto *elem_data = physics::detail::find(atom.element().symbol);
    if (elem_data != nullptr) {
      total_mass += elem_data->mass;
    }
  }
  if (vol > 0.0 && total_mass > 0.0) {
    constexpr auto DA_PER_A3_TO_G_CM3 = static_cast<real_t>(1.66053906660);
    const real_t density_g_cm3 = (total_mass / vol) * DA_PER_A3_TO_G_CM3;
    window_.set_box_density(slint::SharedString(std::format("{:.2f} g/cm³", density_g_cm3)));
  } else {
    window_.set_box_density(slint::SharedString("N/A"));
  }
}

void FileIOHandler::updateFileMetadata(const std::string &filepath) {
  const std::string basename = std::filesystem::path(filepath).filename().string();
  window_.set_file_basename(slint::SharedString(basename));
  window_.set_num_frames(static_cast<int>(loader_.getFrameCount()));
  window_.set_total_atoms(static_cast<int>(loader_.getTotalAtomCount()));
  window_.set_removed_frames_count(static_cast<int>(loader_.getRemovedFrameCount()));
}

void FileIOHandler::handleBrowseFile() {
  if (dialog_active_.exchange(true)) {
    return;
  }

  if (dialog_thread_.joinable()) {
    dialog_thread_.join();
  }

  dialog_thread_ = std::jthread([this]() {
    std::array<nfdfilteritem_t, 10> filter_list = {
        {{
             .name = "Supported Structure Files",
             .spec = "arc,car,cell,cif,dat,md,outmol,poscar,contcar,vasp,xdatcar",
         },
         {
             .name = "Materials Studio CAR",
             .spec = "car",
         },
         {
             .name = "Materials Studio ARC",
             .spec = "arc",
         },
         {
             .name = "CASTEP CELL",
             .spec = "cell",
         },
         {
             .name = "CASTEP MD",
             .spec = "md",
         },
         {
             .name = "CIF files",
             .spec = "cif",
         },
         {
             .name = "ONETEP DAT",
             .spec = "dat",
         },
         {
             .name = "DMol3 Outmol",
             .spec = "outmol",
         },
         {
             .name = "VASP POSCAR/CONTCAR",
             .spec = "poscar,contcar,vasp",
         },
         {
             .name = "VASP XDATCAR",
             .spec = "xdatcar",
         }}};
    const nfdfiltersize_t filter_count = filter_list.size();

    nfdchar_t *out_path = nullptr;
    const nfdresult_t result =
        NFD_OpenDialogU8(&out_path, filter_list.data(), filter_count, nullptr);

    if (result == NFD_OKAY) {
      std::string filepath(out_path);
      NFD_FreePathU8(out_path);

      slint::invoke_from_event_loop([this, filepath = std::move(filepath)]() {
        dialog_active_.store(false);
        startLoadingTrajectory(filepath);
      });
    } else if (result == NFD_CANCEL) {
      slint::invoke_from_event_loop([this]() {
        dialog_active_.store(false);
        const std::string message = AppDefaults::MSG_FILE_SELECTION_CANCELLED;
        window_.set_file_status_text(slint::SharedString(message));
        window_.set_timer_running(false);
        window_.set_text_opacity(false);
      });
    } else {
      std::string error_msg = "Error opening file dialog: ";
      error_msg += NFD_GetError();
      slint::invoke_from_event_loop([this, error_msg = std::move(error_msg)]() {
        dialog_active_.store(false);
        window_.set_file_status_text(slint::SharedString(error_msg));
        window_.set_timer_running(false);
        window_.set_text_opacity(false);
      });
    }
  });
}

} // namespace correlation::app
