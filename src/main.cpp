// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "../include/AppBackend.hpp"
#include "../include/PortableFileDialogs.hpp"
#include "app_window.h"

std::unique_ptr<pfd::open_file> current_file_dialog;

int main() {
  // Create the main window from the Slint UI definition
  auto ui = AppWindow::create();
  // Create an instance of our application backend
  AppBackend backend;

  // Set up the callbacks from the UI to our C++ backend
  ui->on_load_file([&](const slint::SharedString &path) {
    backend.load_file(static_cast<std::string>(path));
  });

  ui->on_run_analysis([&]() {
    // Create a ProgramOptions object from the UI properties
    ProgramOptions options;
    options.input_file = ui->get_in_file_text().data();
    options.output_file_base = options.input_file;
    options.r_cut = ui->get_r_cut();
    options.r_bin_width = ui->get_r_bin_width();
    options.angle_bin_width = ui->get_angle_bin_width();
    options.bond_factor = ui->get_bond_factor();
    options.smoothing = ui->get_smoothing();
    options.smoothing_sigma = ui->get_smoothing_sigma();
    options.smoothing_kernel =
        static_cast<KernelType>(ui->get_smoothing_kernel());

    backend.run_analysis(options);
  });

  ui->on_browse_file([&]() {
    // Define the file filters
    std::vector<std::string> filters = {"Supported Structure Files",
                                        "*.car *.cell *.cif *.dat",
                                        "Materials Studio CAR",
                                        "*.car",
                                        "CASTEP CELL",
                                        "*.cell",
                                        "CIF files",
                                        "*.cif",
                                        "ONETEP DAT",
                                        "*.dat",
                                        "All Files",
                                        "*"};

    // Create and show the file dialog, but don't call result() immediately
    current_file_dialog = std::make_unique<pfd::open_file>(
        "Select a structure file", "", filters, pfd::opt::multiselect);

    ui->set_timer_running(true);
  });

  ui->on_check_file_dialog_status([&]() {
    // Check if the dialog is ready without blocking
    if (current_file_dialog && current_file_dialog->ready(0)) {
      // Dialog is ready, get the result
      auto selection = current_file_dialog->result();
      if (!selection.empty()) {
        ui->set_in_file_text(slint::SharedString(selection[0]));
      }
      // The dialog is finished, stop the timer and reset the unique_ptr
      ui->set_timer_running(false);
      current_file_dialog.reset();
    }
  });

  // Run the Slint event loop
  ui->run();

  return 0;
}
