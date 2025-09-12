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

  ui->on_run_analysis([&backend]() { backend.run_analysis(); });

  // Run the Slint event loop
  ui->run();

  return 0;
}
