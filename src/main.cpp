/**
 * @file main.cpp
 * @brief Application entry point.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#ifdef _WIN32
#define NOMINMAX
#define IDI_ICON1 101
#include <Windows.h>
#include <cstdio>
#include <cstdlib>
#endif

#include "AppWindow.h"
#include "app/AppBackend.hpp"
#include "app/AppController.hpp"

#include <exception>
#include <iostream>

namespace {

/**
 * @brief Attempts to create the Slint AppWindow with fallback to software rendering.
 * @return Constructed Slint component handle for AppWindow.
 */
slint::ComponentHandle<AppWindow> createAppWindow() {
  try {
    return AppWindow::create();
  } catch (const std::exception &err) {
    std::cerr << "Primary Slint backend initialization failed: " << err.what() << "\n";
    std::cerr << "Falling back to software rendering backend...\n";
#ifdef _WIN32
    _putenv_s("SLINT_BACKEND", "software");
#else
    setenv("SLINT_BACKEND", "software", 1);
#endif
    return AppWindow::create();
  } catch (...) {
    std::cerr << "Primary Slint backend failed with unknown error. Falling back to software renderer...\n";
#ifdef _WIN32
    _putenv_s("SLINT_BACKEND", "software");
#else
    setenv("SLINT_BACKEND", "software", 1);
#endif
    return AppWindow::create();
  }
}

} // namespace

/**
 * @brief Main entry point of the application.
 *
 * Initializes the Slint UI, backend components, and the AppController.
 * Starts the Slint event loop with fallback error handling for graphics initialization.
 *
 * @return Exit code (0 for success, 1 on fatal error).
 */
int main() {
#ifdef _WIN32
  // Attach to parent console if launched from terminal (CMD/PowerShell)
  if (AttachConsole(ATTACH_PARENT_PROCESS) != 0) {
    FILE *dummy_file = nullptr;
    freopen_s(&dummy_file, "CONOUT$", "w", stdout);
    freopen_s(&dummy_file, "CONOUT$", "w", stderr);
    freopen_s(&dummy_file, "CONIN$", "r", stdin);
  }

  // Load the icon from the application's resource file
  LoadIcon(GetModuleHandle(nullptr), MAKEINTRESOURCE(IDI_ICON1));
#endif

  try {
    auto window = createAppWindow();
    correlation::app::AppBackend backend;
    correlation::app::AppController const controller(*window, backend);

    window->run();
    return 0;
  } catch (const std::exception &err) {
    std::cerr << "Fatal startup error: " << err.what() << "\n";
#ifdef _WIN32
    MessageBoxA(nullptr, err.what(), "Correlation - Startup Error", MB_ICONERROR | MB_OK);
#endif
    return 1;
  } catch (...) {
    std::cerr << "An unknown fatal error occurred during application startup.\n";
#ifdef _WIN32
    MessageBoxA(nullptr, "An unknown fatal error occurred during application startup.",
                "Correlation - Startup Error", MB_ICONERROR | MB_OK);
#endif
    return 1;
  }
}

#if _WIN32
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow) { return main(); }
#endif


