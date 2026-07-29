/**
 * @file main.cpp
 * @brief Application entry point.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include <cstdio>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>

#ifdef _WIN32
#define NOMINMAX
#define IDI_ICON1 101
#include <Windows.h>

namespace {
void setupWindowsDebugEnvironment() {
  // Ensure Rust/Slint outputs detailed backtraces on panics
  const char *existing_backtrace = std::getenv("RUST_BACKTRACE");
  if (existing_backtrace == nullptr || std::string(existing_backtrace)[0] == '\0') {
    _putenv_s("RUST_BACKTRACE", "full");
  }

  // Attach to parent console if launched from terminal (CMD/PowerShell)
  if (AttachConsole(ATTACH_PARENT_PROCESS) != 0) {
    FILE *dummy_file = nullptr;
    freopen_s(&dummy_file, "CONOUT$", "w", stdout);
    freopen_s(&dummy_file, "CONOUT$", "w", stderr);
    freopen_s(&dummy_file, "CONIN$", "r", stdin);
  } else {
    // If launched from GUI Explorer without a console, mirror stderr to a launch debug log file
    std::error_code ec;
    auto temp_dir = std::filesystem::temp_directory_path(ec);
    if (!ec) {
      auto log_path = temp_dir / "correlation_debug_launch.log";
      FILE *log_file = nullptr;
      if (freopen_s(&log_file, log_path.string().c_str(), "w", stderr) == 0) {
        std::cerr << "=== Correlation Debug Launch Log ===\n";
      }
    }
  }

  // Set top-level unhandled exception filter
  SetUnhandledExceptionFilter([](EXCEPTION_POINTERS *info) -> LONG {
    char msg[512];
    std::snprintf(msg, sizeof(msg), "Unhandled Windows Exception 0x%08X at Address %p.\nLog written to temp directory.",
                  static_cast<unsigned int>(info->ExceptionRecord->ExceptionCode),
                  info->ExceptionRecord->ExceptionAddress);
    std::cerr << "FATAL: " << msg << "\n";
    std::cerr.flush();
    MessageBoxA(nullptr, msg, "Correlation - Unhandled Crash", MB_ICONERROR | MB_OK);
    return EXCEPTION_EXECUTE_HANDLER;
  });
}
} // namespace
#endif

#include "AppWindow.h"
#include "app/AppBackend.hpp"
#include "app/AppController.hpp"

namespace {

/**
 * @brief Custom terminate handler to log uncaught exceptions or panics.
 */
void handleTerminate() {
  std::cerr << "FATAL: std::terminate called due to unhandled exception or panic.\n";
  if (auto exc_ptr = std::current_exception()) {
    try {
      std::rethrow_exception(exc_ptr);
    } catch (const std::exception &e) {
      std::cerr << "Exception detail: " << e.what() << "\n";
#ifdef _WIN32
      MessageBoxA(nullptr, e.what(), "Correlation - Uncaught Exception", MB_ICONERROR | MB_OK);
#endif
    } catch (...) {
      std::cerr << "Unknown exception type caught during terminate.\n";
    }
  }
  std::cerr.flush();
  std::abort();
}

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
  std::set_terminate(handleTerminate);

#ifdef _WIN32
  setupWindowsDebugEnvironment();
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
    MessageBoxA(nullptr, "An unknown fatal error occurred during application startup.", "Correlation - Startup Error",
                MB_ICONERROR | MB_OK);
#endif
    return 1;
  }
}

#if _WIN32
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow) { return main(); }
#endif
