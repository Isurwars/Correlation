/**
 * @file FileIOHandler.cpp
 * @brief Implementation of FileIOHandler.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/FileIOHandler.hpp"
#include "AppWindow.h"
#include "app/AppController.hpp"
#include "app/InputValidator.hpp"
#include <filesystem>
#include <format>
#include <nfd.h>

namespace correlation::app {

FileIOHandler::FileIOHandler(::AppWindow &window, AppBackend &backend, AppController &controller)
    : window_(window), backend_(backend), controller_(controller) {}

FileIOHandler::~FileIOHandler() {
  if (dialog_thread_.joinable()) {
    dialog_thread_.join();
  }
  if (load_thread_.joinable()) {
    load_thread_.join();
  }
}

void FileIOHandler::executeWriteFiles(const std::string &filepath) {
  std::filesystem::path file_path_obj(filepath);
  const std::string ext = file_path_obj.extension().string();

#ifndef CORRELATION_USE_HDF5
  if (ext == ".h5" || ext == ".hdf5") {
    window_.set_analysis_status_text("Error: HDF5 format is not available.");
    return;
  }
#endif

#ifndef CORRELATION_USE_ARROW
  if (ext == ".parquet") {
    window_.set_analysis_status_text("Error: Parquet format is not available.");
    return;
  }
#endif

  const bool use_hdf5 = (ext == ".h5" || ext == ".hdf5");
  const bool use_parquet = (ext == ".parquet");
  const bool use_csv = (!use_hdf5 && !use_parquet);

  if (use_csv && ext != ".csv") {
    file_path_obj.replace_extension(".csv");
  }

  if (file_path_obj.has_extension()) {
    file_path_obj.replace_extension("");
  }

  ProgramOptions opts = controller_.handleOptionsfromUI();
  opts.output_file_base = file_path_obj.string();
  opts.use_csv = use_csv;
  opts.use_hdf5 = use_hdf5;
  opts.use_parquet = use_parquet;
  backend_.setOptions(opts);

  const auto write_res = backend_.write_files();
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

  dialog_thread_ = std::thread([this]() {
    std::array<nfdfilteritem_t, 3> filter_list = {{{
                                                       .name = "Comma Separated Values",
                                                       .spec = "csv",
                                                   },
                                                   {
                                                       .name = "Hierarchical Data Format",
                                                       .spec = "h5,hdf5",
                                                   },
                                                   {
                                                       .name = "Apache Parquet",
                                                       .spec = "parquet",
                                                   }}};
    const nfdfiltersize_t filter_count = filter_list.size();

    nfdchar_t *out_path = nullptr;
    const nfdresult_t result =
        NFD_SaveDialogU8(&out_path, filter_list.data(), filter_count, nullptr, nullptr);

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

  load_thread_ = std::thread([this, filepath]() {
    backend_.setProgressCallback([this](float progress, const std::string &msg) {
      slint::invoke_from_event_loop([progress, msg, this]() {
        window_.set_progress(progress);
        if (!msg.empty()) {
          window_.set_file_status_text(slint::SharedString(msg));
        }
      });
    });

    std::string message;
    bool success = false;
    try {
      message = backend_.load_file(filepath);
      success = true;
    } catch (const std::exception &e) {
      message = std::string(AppDefaults::MSG_ERROR_LOADING) + std::string(e.what());
    }

    slint::invoke_from_event_loop([this, message, success]() {
      window_.set_file_status_text(slint::SharedString(message));
      window_.set_file_loading(false);
      window_.set_timer_running(false);
      window_.set_text_opacity(false);

      if (success && backend_.cell() != nullptr) {
        window_.set_file_loaded(true);

        // Atom Counts
        auto atom_counts_map = backend_.getAtomCounts();
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

        // File Info
        window_.set_num_frames(static_cast<int>(backend_.getFrameCount()));
        window_.set_total_atoms(static_cast<int>(backend_.getTotalAtomCount()));
        window_.set_removed_frames_count(static_cast<int>(backend_.getRemovedFrameCount()));
        {
          auto opts = window_.get_analysis_options();
          opts.time_step =
              slint::SharedString(std::format("{:.2f}", backend_.getRecommendedTimeStep()));
          window_.set_analysis_options(opts);
        }

        // Update Run Analysis Card Frame Info
        {
          auto opts = window_.get_analysis_options();
          opts.min_frame = "1";
          window_.set_analysis_options(opts);
        }
        {
          auto opts = window_.get_analysis_options();
          opts.max_frame = slint::SharedString(std::to_string(backend_.getFrameCount()));
          window_.set_analysis_options(opts);
        }
        static_cast<void>(controller_.getInputValidator()->validateInputs());
      }
    });
  });
}

void FileIOHandler::handleBrowseFile() {
  if (dialog_active_.exchange(true)) {
    return;
  }

  if (dialog_thread_.joinable()) {
    dialog_thread_.join();
  }

  dialog_thread_ = std::thread([this]() {
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
