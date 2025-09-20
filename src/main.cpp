// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#if defined(_WIN32)
#define NOMINMAX // Prevents Windows.h from defining min and max macros
#include <Windows.h>
#endif

#include <string>

#include "../include/AppBackend.hpp"
#include "../include/PortableFileDialogs.hpp"

#include "app_window.h"

std::unique_ptr<pfd::open_file> current_file_dialog;

int main() {
  // Create the main window from the Slint UI definition
  auto ui = AppWindow::create();
  // Create an instance of our application backend
  AppBackend backend;

  ui->on_run_analysis([&backend, &ui]() {
    try {
      // Create a ProgramOptions object from the UI properties
      ProgramOptions options;
      const std::string input_path = ui->get_in_file_text().data();
      options.input_file = input_path;
      options.output_file_base = input_path;
      options.r_cut = ui->get_r_cut();
      options.r_bin_width = ui->get_r_bin_width();
      options.angle_bin_width = ui->get_angle_bin_width();
      options.bond_factor = ui->get_bond_factor();
      options.smoothing = ui->get_smoothing();
      options.smoothing_sigma = ui->get_smoothing_sigma();
      options.smoothing_kernel =
          static_cast<KernelType>(ui->get_smoothing_kernel());

      backend.run_analysis(options);
      ui->set_status_text("Analysis ended successfully.");
    } catch (const std::exception &e) {
      ui->set_status_text(
          slint::SharedString("Error: " + std::string(e.what())));
    }
  });

  ui->on_browse_file([&ui]() {
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

  ui->on_check_file_dialog_status([&backend, &ui]() {
    // Check if the dialog is ready without blocking
    if (current_file_dialog && current_file_dialog->ready(0)) {
      // Dialog is ready, get the result
      auto selection = current_file_dialog->result();
      if (!selection.empty()) {
        try {
          // Use std::string constructor for safety
          std::string selected_path = selection[0];
          ui->set_in_file_text(slint::SharedString(selected_path));
          backend.load_file(selected_path);
          ui->set_status_text("File: " + slint::SharedString(selected_path) +
                              " loaded successfully.");
        } catch (const std::exception &e) {
          ui->set_status_text(slint::SharedString("Error loading file: " +
                                                  std::string(e.what())));
        }
      } else {
        // Handle the case where the user canceled the dialog
        ui->set_status_text("File selection cancelled.");
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

#if _WIN32
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
                   LPSTR lpCmdLine, int nCmdShow) {
  return main();
}
#endif
