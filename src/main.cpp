// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "../include/AppBackend.hpp"
#include "app_window.h"

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

  // Run the Slint event loop
  ui->run();

  return 0;
}
