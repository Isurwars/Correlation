// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
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
    // Ensure Slint uses the software backend during tests to avoid OpenGL
    // requirements in headless CI/Xvfb environments
#if !defined(_WIN32)
    setenv("SLINT_BACKEND", "software", 1);
#else
    _putenv_s("SLINT_BACKEND", "software");
#endif
  }

  correlation::app::ProgramOptions callHandleOptionsfromUI(correlation::app::AppController &controller, AppWindow &ui) {
    return controller.handleOptionsfromUI(ui);
  }
};

TEST_F(AppControllerTests, ConstructorInitializesCorrectly) {
  // Slint testing backend setup if headless env var is not enough.
  // slint::testing::backend::init(); // available in slint-testing

  // Arrange
  auto ui = AppWindow::create();

  correlation::app::AppBackend backend;

  // Act & Assert
  EXPECT_NO_THROW({ correlation::app::AppController controller(*ui, backend); });
}

TEST_F(AppControllerTests, HandlesBondCutoffsCorrectly) {
  // Arrange
  auto ui = AppWindow::create();
  correlation::app::AppBackend backend;
  
  std::string file_path = "../../examples/a-PdSi/a-PdSi.car";
  if (!std::filesystem::exists(file_path)) {
    file_path = "../examples/a-PdSi/a-PdSi.car";
  }
  if (!std::filesystem::exists(file_path)) {
    file_path = "examples/a-PdSi/a-PdSi.car";
  }
  backend.load_file(file_path);

  correlation::app::AppController controller(*ui, backend);

  // Set up elements in UI that match the elements in loaded cell (Pd, Si)
  auto atom_counts = std::make_shared<slint::VectorModel<AtomCount>>();
  atom_counts->push_back({"Pd", 80});
  atom_counts->push_back({"Si", 20});
  ui->set_atom_counts(atom_counts);

  // Set up cutoffs: Pd-Pd, Pd-Si, Si-Si
  // In upper-triangular order: Pd-Pd, Pd-Si, Si-Si
  auto cutoffs = std::make_shared<slint::VectorModel<BondCutoff>>();
  cutoffs->push_back({"Pd", "Pd", "2.80"});
  cutoffs->push_back({"Pd", "Si", "3.10"});
  cutoffs->push_back({"Si", "Si", "2.50"});
  ui->set_bond_cutoffs(cutoffs);

  // Act
  auto opts = callHandleOptionsfromUI(controller, *ui);

  // Assert: verify that the flat cutoffs are mapped to symmetric 2D matrix
  ASSERT_EQ(opts.bond_cutoffs_sq.size(), 2);
  ASSERT_EQ(opts.bond_cutoffs_sq[0].size(), 2);
  ASSERT_EQ(opts.bond_cutoffs_sq[1].size(), 2);

  EXPECT_NEAR(opts.bond_cutoffs_sq[0][0], 2.50 * 2.50, 1e-5); // Si-Si
  EXPECT_NEAR(opts.bond_cutoffs_sq[0][1], 3.10 * 3.10, 1e-5); // Si-Pd
  EXPECT_NEAR(opts.bond_cutoffs_sq[1][0], 3.10 * 3.10, 1e-5); // Pd-Si
  EXPECT_NEAR(opts.bond_cutoffs_sq[1][1], 2.80 * 2.80, 1e-5); // Pd-Pd
}

