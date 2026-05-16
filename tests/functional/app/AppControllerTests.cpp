// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "AppWindow.h"
#include "app/AppBackend.hpp"
#include "app/AppController.hpp"

#include <gtest/gtest.h>
#include <stdlib.h>

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#endif

// Test fixture for the AppController class
class AppControllerTests : public ::testing::Test {
protected:
  void SetUp() override {
    // Ensure Slint uses the headless backend during tests to avoid UI/X11
    // requirements
#if !defined(_WIN32)
    setenv("SLINT_BACKEND", "winit-headless", 1);
    setenv("WINIT_UNIX_BACKEND", "x11",
           1); // sometimes winit headless needs this or just use
               // SLINT_BACKEND=software or testing
#else
    _putenv_s("SLINT_BACKEND", "winit-headless");
#endif
  }
};

TEST_F(AppControllerTests, ConstructorInitializesCorrectly) {
  // Slint testing backend setup if headless env var is not enough.
  // slint::testing::backend::init(); // available in slint-testing

  // Arrange
  auto ui = AppWindow::create();

  correlation::app::AppBackend backend;

  // Act & Assert
  EXPECT_NO_THROW(
      { correlation::app::AppController controller(*ui, backend); });
}
