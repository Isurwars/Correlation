/**
 * @file main.cpp
 * @brief Application entry point.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#if defined(_WIN32)
#define NOMINMAX
#define IDI_ICON1 101
#include <Windows.h>
#endif

#include "AppWindow.h"
#include "app/AppBackend.hpp"
#include "app/AppController.hpp"

#include <stdlib.h>

/**
 * @brief Main entry point of the application.
 *
 * Initializes the Slint UI, backend components, and the AppController.
 * Starts the Slint event loop.
 *
 * @return Exit code (0 for success).
 */
int main() {

#if defined(_WIN32)
  // Load the icon from the application's resource file
  HICON hIcon = LoadIcon(GetModuleHandle(nullptr), MAKEINTRESOURCE(IDI_ICON1));
#endif

#if !defined(_WIN32)
  setenv("SLINT_BACKEND", "winit-skia", 1);
#else
  _putenv_s("SLINT_BACKEND", "winit-skia");
#endif

  auto ui = AppWindow::create();
  correlation::app::AppBackend backend;

  // Instantiate the controller, which sets up all callbacks
  correlation::app::AppController controller(*ui, backend);

  ui->run();

  return 0;
}

#if _WIN32
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
                   LPSTR lpCmdLine, int nCmdShow) {
  return main();
}
#endif
