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

#include "app/AppBackend.hpp"
#include "app/AppController.hpp"
#include "AppWindow.h"

int main() {

#if defined(_WIN32)
  // Load the icon from the application's resource file
  HICON hIcon = LoadIcon(GetModuleHandle(nullptr), MAKEINTRESOURCE(IDI_ICON1));
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
