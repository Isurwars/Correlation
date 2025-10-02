// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#if defined(_WIN32)
#define NOMINMAX
#define IDI_ICON1 101
#include <Windows.h>
#endif

#include "../include/AppBackend.hpp"
#include "../include/AppController.hpp"

#include "app_window.h"

int main() {

#if defined(_WIN32)
  // Load the icon from the application's resource file
  HICON hIcon = LoadIcon(GetModuleHandle(nullptr), MAKEINTRESOURCE(IDI_ICON1));
#endif

  auto ui = AppWindow::create();
  AppBackend backend;

  // Instantiate the controller, which sets up all callbacks
  AppController controller(*ui, backend);

  ui->run();

  return 0;
}

#if _WIN32
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
                   LPSTR lpCmdLine, int nCmdShow) {
  return main();
}
#endif
