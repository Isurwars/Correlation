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
#include "app/viewmodel/PlotController.hpp"
#include "app/viewmodel/Structure3DController.hpp"
#include "physics/PhysicalData.hpp"
#include "readers/NaturalSort.hpp"
#include "render3d/ElementColors.hpp"
#include "writers/structure/BaseStructureWriter.hpp"
#include "writers/structure/StructureWriterFactory.hpp"
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
  alive_->store(false, std::memory_order_release);
  if (dialog_thread_.joinable()) {
    dialog_thread_.join();
  }
  if (load_thread_.joinable()) {
    load_thread_.join();
  }
}

namespace {

[[nodiscard]] ProgramOptions buildExportOptionsFromUI(const AppWindow &window,
                                                      AppController &controller) {
  ProgramOptions opts = controller.handleOptionsfromUI();
  const auto format_opts = window.get_export_format_options();
  opts.use_zip = format_opts.use_zip;
  opts.use_csv = format_opts.export_csv;
  opts.use_hdf5 = format_opts.export_hdf5;
  opts.use_parquet = format_opts.export_arrow;
  opts.export_images = format_opts.export_images;
  opts.export_algorithms.clear();
  if (!window.get_export_all_analyses()) {
    if (const auto model = window.get_export_analysis_items(); model != nullptr) {
      for (size_t i = 0; i < model->row_count(); ++i) {
        if (const auto item = model->row_data(i); item.has_value() && item->enabled) {
          opts.export_algorithms.emplace_back(item->id.data());
        }
      }
    }
  }
  return opts;
}

void dispatchStatusUpdate(AppWindow &window, std::shared_ptr<std::atomic<bool>> alive,
                          std::string status_msg) {
  slint::invoke_from_event_loop(
      [&window, alive = std::move(alive), status_msg = std::move(status_msg)] {
        if (!alive->load(std::memory_order_acquire)) {
          return;
        }
        window.set_analysis_status_text(slint::SharedString(status_msg));
      });
}

} // namespace

void FileIOHandler::executeWriteFiles(const std::string &filepath) {
  ProgramOptions opts = buildExportOptionsFromUI(window_, controller_);
  opts.output_file_base = filepath;
  options_ = std::move(opts);

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
  ProgramOptions export_opts = buildExportOptionsFromUI(window_, controller_);

  dialog_thread_ =
      std::jthread([this, alive = alive_, is_zip, export_opts = std::move(export_opts)] mutable {
        nfdchar_t *out_path = nullptr;
        nfdresult_t result = NFD_CANCEL;

        if (is_zip) {
          std::array<nfdfilteritem_t, 1> filter_list = {
              {
                  {
                      .name = "Consolidated ZIP Bundle",
                      .spec = "zip",
                  },
              },
          };
          result = NFD_SaveDialogU8(&out_path, filter_list.data(), filter_list.size(), nullptr,
                                    "correlation_export.zip");
        } else {
          result = NFD_PickFolderU8(&out_path, nullptr);
        }

        if (result == NFD_OKAY) {
          export_opts.output_file_base = std::string(out_path);
          NFD_FreePathU8(out_path);

          // Perform heavy disk/archive I/O on the background worker thread (off UI thread)
          const auto write_res = dispatcher_.writeFiles(export_opts);
          const std::string status_msg =
              write_res ? std::string(AppDefaults::MSG_FILES_WRITTEN) : write_res.error();
          dialog_active_.store(false);

          slint::invoke_from_event_loop(
              [this, alive, export_opts = std::move(export_opts), status_msg] mutable {
                if (!alive->load(std::memory_order_acquire)) {
                  return;
                }
                options_ = std::move(export_opts);
                window_.set_analysis_status_text(slint::SharedString(status_msg));
              });
        } else if (result == NFD_CANCEL) {
          dialog_active_.store(false);
          dispatchStatusUpdate(window_, alive, std::string(AppDefaults::MSG_SAVE_CANCELLED));
        } else {
          std::string error_msg = "Error: ";
          error_msg += NFD_GetError();
          dialog_active_.store(false);
          dispatchStatusUpdate(window_, alive, std::move(error_msg));
        }
      });
}

void FileIOHandler::handleExportStructure() {
  const auto *cell = loader_.cell();
  if (cell == nullptr) {
    window_.set_analysis_status_text("Error: No structure loaded to export or convert.");
    return;
  }

  if (dialog_active_.exchange(true)) {
    return;
  }

  if (dialog_thread_.joinable()) {
    dialog_thread_.join();
  }

  const correlation::core::Cell cell_copy = *cell;

  dialog_thread_ = std::jthread([this, alive = alive_, cell_copy] {
    nfdchar_t *out_path = nullptr;
    std::array<nfdfilteritem_t, 4> const filter_list = {
        {
            {.name = "VASP POSCAR", .spec = "vasp,poscar"},
            {.name = "CASTEP Cell", .spec = "cell"},
            {.name = "Crystallographic Information File", .spec = "cif"},
            {.name = "Accelrys CAR", .spec = "car"},
        },
    };

    nfdresult_t const result = NFD_SaveDialogU8(&out_path, filter_list.data(), filter_list.size(),
                                                nullptr, "structure.vasp");

    if (result == NFD_OKAY) {
      std::string const filepath(out_path);
      NFD_FreePathU8(out_path);

      std::string status_msg;
      try {
        const std::filesystem::path path(filepath);
        std::string ext = path.extension().string();
        if (ext.empty()) {
          ext = path.filename().string();
        }
        const auto *const writer =
            writers::structure::StructureWriterFactory::instance().getWriterForExtension(ext);
        if (writer == nullptr) {
          status_msg = "Error: Unsupported export format for " + filepath;
        } else {
          writer->write(filepath, cell_copy);
          status_msg = "Successfully exported structure to " + filepath;
        }
      } catch (const std::exception &e) {
        status_msg = std::string("Error exporting structure: ") + e.what();
      }
      dialog_active_.store(false);

      slint::invoke_from_event_loop([this, alive, status_msg = std::move(status_msg)] {
        if (!alive->load(std::memory_order_acquire)) {
          return;
        }
        window_.set_analysis_status_text(slint::SharedString(status_msg));
      });
    } else {
      dialog_active_.store(false);
    }
  });
}

void FileIOHandler::startLoadingTrajectory(const std::string &filepath,
                                           std::function<void(bool)> on_complete) {
  loaded_filepaths_ = {filepath};
  window_.set_in_file_text(slint::SharedString(filepath));
  window_.set_file_status_text(slint::SharedString("Loading file..."));
  window_.set_file_loading(true);
  window_.set_timer_running(true);
  window_.set_text_opacity(true);
  window_.set_progress(0.0F);

  if (load_thread_.joinable()) {
    load_thread_.join();
  }

  load_thread_ =
      std::jthread([this, alive = alive_, filepath, on_complete = std::move(on_complete)] {
        const auto progress_cb = [this, alive](float progress, const std::string &msg) {
          slint::invoke_from_event_loop([progress, msg, this, alive] {
            if (!alive->load(std::memory_order_acquire)) {
              return;
            }
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
        } else {
          message = std::string(AppDefaults::MSG_ERROR_LOADING) + load_res.error();
        }

        slint::invoke_from_event_loop([this, alive, filepath, message, success, on_complete] {
          if (!alive->load(std::memory_order_acquire)) {
            return;
          }
          window_.set_file_status_text(slint::SharedString(message));
          window_.set_file_loading(false);
          window_.set_timer_running(false);
          window_.set_text_opacity(false);

          if (success) {
            options_.input_file = filepath;
            onTrajectoryLoadSuccess(filepath);
          }

          if (on_complete) {
            on_complete(success);
          }
        });
      });
}

void FileIOHandler::startLoadingTrajectories(const std::vector<std::string> &filepaths,
                                             std::function<void(bool)> on_complete) {
  if (filepaths.empty()) {
    return;
  }
  if (filepaths.size() == 1) {
    startLoadingTrajectory(filepaths.front(), std::move(on_complete));
    return;
  }

  std::vector<std::string> sorted_paths = filepaths;
  readers::naturalSort(sorted_paths);
  loaded_filepaths_ = sorted_paths;

  const std::string display_name =
      sorted_paths.front() + " (+" + std::to_string(sorted_paths.size() - 1) + " files)";
  window_.set_in_file_text(slint::SharedString(display_name));
  window_.set_file_status_text(slint::SharedString("Loading multi-file trajectory..."));
  window_.set_file_loading(true);
  window_.set_timer_running(true);
  window_.set_text_opacity(true);
  window_.set_progress(0.0F);

  if (load_thread_.joinable()) {
    load_thread_.join();
  }

  load_thread_ = std::jthread([this, alive = alive_, sorted_paths = std::move(sorted_paths),
                               on_complete = std::move(on_complete)] {
    const auto progress_cb = [this, alive](float progress, const std::string &msg) {
      slint::invoke_from_event_loop([progress, msg, this, alive] {
        if (!alive->load(std::memory_order_acquire)) {
          return;
        }
        window_.set_progress(progress);
        if (!msg.empty()) {
          window_.set_file_status_text(slint::SharedString(msg));
        }
      });
    };

    std::string message;
    bool success = false;
    const auto load_res = loader_.loadFiles(sorted_paths, progress_cb);
    if (load_res) {
      message = *load_res;
      success = true;
    } else {
      message = std::string(AppDefaults::MSG_ERROR_LOADING) + load_res.error();
    }

    slint::invoke_from_event_loop([this, alive, sorted_paths, message, success, on_complete] {
      if (!alive->load(std::memory_order_acquire)) {
        return;
      }
      window_.set_file_status_text(slint::SharedString(message));
      window_.set_file_loading(false);
      window_.set_timer_running(false);
      window_.set_text_opacity(false);

      if (success) {
        options_.input_file = sorted_paths.front();
        onTrajectoryLoadSuccess(sorted_paths.front());
      }

      if (on_complete) {
        on_complete(success);
      }
    });
  });
}

void FileIOHandler::refreshAtomCountsUI() {
  const auto ordered_counts = loader_.getOrderedAtomCounts();
  const auto slint_atom_counts = std::make_shared<slint::VectorModel<AtomCount>>();
  for (const auto &[symbol, count] : ordered_counts) {
    const auto rgba = render3d::getElementColor(symbol);
    const bool is_unknown = (physics::detail::find(symbol) == nullptr);
    slint_atom_counts->push_back({
        .symbol = slint::SharedString(symbol),
        .count = count,
        .color = slint::Color::from_rgb_uint8(rgba.r, rgba.g, rgba.b),
        .is_unknown = is_unknown,
    });
  }
  window_.set_atom_counts(slint_atom_counts);
}

void FileIOHandler::handleRenameElement(const std::string &old_symbol,
                                        const std::string &new_symbol) {
  if (!loader_.renameElement(old_symbol, new_symbol)) {
    refreshAtomCountsUI();
    window_.set_analysis_status_text(slint::SharedString(
        std::format("Error: Cannot rename species '{}' to '{}'.", old_symbol, new_symbol)));
    return;
  }

  refreshAtomCountsUI();
  controller_.setBondCutoffs();

  const int current_frame = std::max(0, window_.get_structure_3d_current_frame());
  updateBoxDiagnosticsForFrame(static_cast<size_t>(current_frame));

  {
    auto opts = window_.get_analysis_options();
    opts.time_step = slint::SharedString(std::format("{:.2f}", loader_.getRecommendedTimeStep()));
    window_.set_analysis_options(opts);
  }

  updateStructure3DView();
  if (current_frame > 0 && controller_.getStructure3DController() != nullptr) {
    controller_.getStructure3DController()->handleSeekFrame(current_frame);
  }

  if (controller_.getPlotController() != nullptr) {
    controller_.getPlotController()->populatePlotList();
  }

  static_cast<void>(controller_.getInputValidator()->validateInputs());
  window_.set_analysis_status_text(
      slint::SharedString(std::format("Mapped species '{}' → '{}'.", old_symbol, new_symbol)));
}

void FileIOHandler::onTrajectoryLoadSuccess(const std::string &filepath) {
  if (loader_.cell() == nullptr) {
    return;
  }

  window_.set_file_loaded(true);
  window_.set_analysis_done(false);
  window_.set_results_stale(false);

  // Atom Counts
  refreshAtomCountsUI();

  // Bond Cutoffs
  controller_.setBondCutoffs();

  // MRU Recent Files
  controller_.addRecentFile(filepath);

  // File Info & Basename
  updateFileMetadata(filepath);

  // Box Diagnostics
  updateBoxDiagnostics(loader_.cell());

  {
    auto opts = window_.get_analysis_options();
    opts.time_step = slint::SharedString(std::format("{:.2f}", loader_.getRecommendedTimeStep()));
    opts.min_frame = "1";
    opts.frame_stride = "1";
    opts.max_frame = slint::SharedString(std::to_string(loader_.getFrameCount()));
    window_.set_analysis_options(opts);
  }

  window_.set_active_panel(1);
  if (window_.get_left_col_width() < 180.0F) {
    const float last_w = window_.get_last_col_width();
    window_.set_left_col_width(last_w >= 180.0F ? last_w : 220.0F);
  }

  updateStructure3DView();
  window_.set_preview_active_tab(1);

  if (controller_.getPlotController() != nullptr) {
    controller_.getPlotController()->populatePlotList();
  }

  static_cast<void>(controller_.getInputValidator()->validateInputs());
}

void FileIOHandler::handleReloadFile() {
  if (loaded_filepaths_.size() > 1) {
    startLoadingTrajectories(loaded_filepaths_);
    return;
  }
  const std::string &path = options_.input_file;
  if (!path.empty()) {
    startLoadingTrajectory(path);
  }
}

void FileIOHandler::updateStructure3DView() {
  if (controller_.getStructure3DController() == nullptr) {
    return;
  }
  if (loader_.trajectoryShared() != nullptr) {
    controller_.getStructure3DController()->setTrajectory(loader_.trajectoryShared());
  } else if (loader_.cell() != nullptr) {
    controller_.getStructure3DController()->setCell(*loader_.cell());
  }
  controller_.getStructure3DController()->set3DModeActive(true);
}

void FileIOHandler::updateBoxDiagnosticsForFrame(size_t frame_index) {
  const auto *traj = loader_.trajectory();
  if (traj == nullptr || traj->getFrameCount() == 0) {
    return;
  }
  const size_t clamped_index = std::min(frame_index, traj->getFrameCount() - 1);
  const core::Cell frame_cell = traj->getFrame(clamped_index);
  updateBoxDiagnostics(&frame_cell);
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
  std::string basename = std::filesystem::path(filepath).filename().string();
  if (loaded_filepaths_.size() > 1) {
    basename += " (+" + std::to_string(loaded_filepaths_.size() - 1) + " files)";
  }
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

  dialog_thread_ = std::jthread([this] {
    std::array<nfdfilteritem_t, 11> filter_list = {
        {
            {
                .name = "Supported Structure Files",
                .spec = "arc,car,cell,cif,dat,dump,lammpstrj,md,outmol,poscar,contcar,vasp,xdatcar,"
                        "xyz,exyz,extxyz,mace,chgnet,gap,quip,nequip,gro,pdb,ent,*",
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
                .name = "LAMMPS Dump",
                .spec = "dump,lammpstrj,*",
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
            },
        },
    };
    const auto filter_count = static_cast<nfdfiltersize_t>(filter_list.size());

    const nfdpathset_t *path_set = nullptr;
    const nfdresult_t result =
        NFD_OpenDialogMultipleU8(&path_set, filter_list.data(), filter_count, nullptr);

    if (result == NFD_OKAY) {
      nfdpathsetsize_t num_paths = 0;
      NFD_PathSet_GetCount(path_set, &num_paths);
      std::vector<std::string> filepaths;
      filepaths.reserve(num_paths);
      for (nfdpathsetsize_t i = 0; i < num_paths; ++i) {
        nfdchar_t *out_path = nullptr;
        if (NFD_PathSet_GetPathU8(path_set, i, &out_path) == NFD_OKAY && out_path != nullptr) {
          filepaths.emplace_back(out_path);
          NFD_PathSet_FreePathU8(out_path);
        }
      }
      NFD_PathSet_Free(path_set);

      slint::invoke_from_event_loop([this, filepaths = std::move(filepaths)] {
        dialog_active_.store(false);
        if (filepaths.size() == 1) {
          startLoadingTrajectory(filepaths.front());
        } else if (!filepaths.empty()) {
          startLoadingTrajectories(filepaths);
        }
      });
    } else if (result == NFD_CANCEL) {
      slint::invoke_from_event_loop([this] {
        dialog_active_.store(false);
        const std::string message = AppDefaults::MSG_FILE_SELECTION_CANCELLED;
        window_.set_file_status_text(slint::SharedString(message));
        window_.set_timer_running(false);
        window_.set_text_opacity(false);
      });
    } else {
      std::string error_msg = "Error opening file dialog: ";
      error_msg += NFD_GetError();
      slint::invoke_from_event_loop([this, error_msg = std::move(error_msg)] {
        dialog_active_.store(false);
        window_.set_file_status_text(slint::SharedString(error_msg));
        window_.set_timer_running(false);
        window_.set_text_opacity(false);
      });
    }
  });
}

} // namespace correlation::app
