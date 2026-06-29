// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "AppWindow.h"
#include "app/AppBackend.hpp"
#include "app/AppController.hpp"
#include <filesystem>

#include <gtest/gtest.h>

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

  static correlation::app::ProgramOptions callHandleOptionsfromUI(correlation::app::AppController &controller) {
    return controller.handleOptionsfromUI();
  }
};

TEST_F(AppControllerTests, ConstructorInitializesCorrectly) {
  // Slint testing backend setup if headless env var is not enough.
  // slint::testing::backend::init(); // available in slint-testing

  // Arrange
  auto window = AppWindow::create();

  correlation::app::AppBackend backend;

  // Act & Assert
  EXPECT_NO_THROW({ correlation::app::AppController controller(*window, backend); });
}

TEST_F(AppControllerTests, HandlesBondCutoffsCorrectly) {
  // Arrange
  auto window = AppWindow::create();
  correlation::app::AppBackend backend;

  std::string file_path = "../../examples/a-PdSi/a-PdSi.car";
  if (!std::filesystem::exists(file_path)) {
    file_path = "../examples/a-PdSi/a-PdSi.car";
  }
  if (!std::filesystem::exists(file_path)) {
    file_path = "examples/a-PdSi/a-PdSi.car";
  }
  backend.load_file(file_path);

  correlation::app::AppController controller(*window, backend);

  // Set up elements in UI that match the elements in loaded cell (Pd, Si)
  auto atom_counts = std::make_shared<slint::VectorModel<AtomCount>>();
  atom_counts->push_back({.symbol = "Pd", .count = 80});
  atom_counts->push_back({.symbol = "Si", .count = 20});
  window->set_atom_counts(atom_counts);

  // Set up cutoffs: Pd-Pd, Pd-Si, Si-Si
  // In upper-triangular order: Pd-Pd, Pd-Si, Si-Si
  auto cutoffs = std::make_shared<slint::VectorModel<BondCutoff>>();
  cutoffs->push_back({.element1 = "Pd", .element2 = "Pd", .distance = "2.80"});
  cutoffs->push_back({.element1 = "Pd", .element2 = "Si", .distance = "3.10"});
  cutoffs->push_back({.element1 = "Si", .element2 = "Si", .distance = "2.50"});
  window->set_bond_cutoffs(cutoffs);

  // Act
  auto opts = callHandleOptionsfromUI(controller);

  // Assert: verify that the flat cutoffs are mapped to symmetric 2D matrix
  ASSERT_EQ(opts.bond_cutoffs_sq.size(), 2);
  ASSERT_EQ(opts.bond_cutoffs_sq[0].size(), 2);
  ASSERT_EQ(opts.bond_cutoffs_sq[1].size(), 2);

  EXPECT_NEAR(opts.bond_cutoffs_sq[0][0], 2.50 * 2.50, 1e-5); // Si-Si
  EXPECT_NEAR(opts.bond_cutoffs_sq[0][1], 3.10 * 3.10, 1e-5); // Si-Pd
  EXPECT_NEAR(opts.bond_cutoffs_sq[1][0], 3.10 * 3.10, 1e-5); // Pd-Si
  EXPECT_NEAR(opts.bond_cutoffs_sq[1][1], 2.80 * 2.80, 1e-5); // Pd-Pd
}

TEST_F(AppControllerTests, SynchronizesOptionsToAndFromUI) {
  auto window = AppWindow::create();
  correlation::app::AppBackend backend;

  correlation::app::ProgramOptions backend_opts;
  backend_opts.input_file = "dummy_input.arc";
  backend_opts.smoothing = true;
  backend_opts.r_max = 12.5;
  backend_opts.r_bin_width = 0.05;
  backend_opts.q_max = 15.0;
  backend_opts.q_bin_width = 0.1;
  backend_opts.r_int_max = 8.0;
  backend_opts.angle_bin_width = 2.0;
  backend_opts.dihedral_bin_width = 3.0;
  backend_opts.max_ring_size = 8;
  backend_opts.smoothing_sigma = 0.03;
  backend_opts.smoothing_kernel = correlation::math::KernelType::Gaussian;
  backend_opts.material_type = 1;
  backend_opts.min_frame = 2; // 0-based index 2 -> 3rd frame
  backend_opts.max_frame = 10;
  backend_opts.time_step = 2.0;

  backend.setOptions(backend_opts);

  correlation::app::AppController controller(*window, backend);

  // 1. Verify synchronization TO the UI
  controller.handleOptionstoUI();

  EXPECT_EQ(window->get_in_file_text(), "dummy_input.arc");
  EXPECT_TRUE(window->get_analysis_options().smoothing_enabled);
  EXPECT_EQ(window->get_analysis_options().r_max, "12.50");
  EXPECT_EQ(window->get_analysis_options().r_bin_width, "0.05");
  EXPECT_EQ(window->get_analysis_options().q_max, "15.00");
  EXPECT_EQ(window->get_analysis_options().q_bin_width, "0.10");
  EXPECT_EQ(window->get_analysis_options().r_int_max, "8.00");
  EXPECT_EQ(window->get_analysis_options().angle_bin_width, "2.00");
  EXPECT_EQ(window->get_analysis_options().dihedral_bin_width, "3.00");
  EXPECT_EQ(window->get_analysis_options().max_ring_size, "8");
  EXPECT_EQ(window->get_analysis_options().smoothing_sigma, "0.03");
  EXPECT_EQ(window->get_analysis_options().smoothing_kernel, static_cast<int>(correlation::math::KernelType::Gaussian));
  EXPECT_EQ(window->get_analysis_options().material_type, 1);
  EXPECT_EQ(window->get_analysis_options().min_frame, "3"); // UI is 1-based (2 + 1 = 3)
  EXPECT_EQ(window->get_analysis_options().max_frame, "10");
  EXPECT_EQ(window->get_analysis_options().time_step, "2.00");

  // 2. Modify properties in UI and verify synchronization FROM the UI
  auto ui_opts = window->get_analysis_options();
  ui_opts.r_max = "14.20";
  ui_opts.min_frame = "Start";
  ui_opts.max_frame = "End";
  window->set_analysis_options(ui_opts);

  auto parsed_opts = controller.handleOptionsfromUI();
  EXPECT_NEAR(parsed_opts.r_max, 14.20, 1e-5);
  EXPECT_EQ(parsed_opts.min_frame, 0);
  EXPECT_EQ(parsed_opts.max_frame, -1);
}

TEST_F(AppControllerTests, PopulatesRecommendedBondCutoffs) {
  auto window = AppWindow::create();
  correlation::app::AppBackend backend;

  std::string file_path = "../../examples/a-PdSi/a-PdSi.car";
  if (!std::filesystem::exists(file_path)) {
    file_path = "../examples/a-PdSi/a-PdSi.car";
  }
  if (!std::filesystem::exists(file_path)) {
    file_path = "examples/a-PdSi/a-PdSi.car";
  }
  backend.load_file(file_path);

  correlation::app::AppController controller(*window, backend);

  controller.setBondCutoffs();

  auto cutoffs = window->get_bond_cutoffs();
  ASSERT_EQ(cutoffs->row_count(), 3); // Pd-Pd, Pd-Si, Si-Si

  bool found_pd_pd = false;
  bool found_pd_si = false;
  bool found_si_si = false;

  for (size_t i = 0; i < cutoffs->row_count(); ++i) {
    auto item = cutoffs->row_data(i).value();
    std::string el1 = item.element1.data();
    std::string el2 = item.element2.data();
    if (el1 == "Pd" && el2 == "Pd") {
      found_pd_pd = true;
    }
    if ((el1 == "Pd" && el2 == "Si") || (el1 == "Si" && el2 == "Pd")) {
      found_pd_si = true;
    }
    if (el1 == "Si" && el2 == "Si") {
      found_si_si = true;
    }
  }

  EXPECT_TRUE(found_pd_pd);
  EXPECT_TRUE(found_pd_si);
  EXPECT_TRUE(found_si_si);
}

TEST_F(AppControllerTests, UpdatesActiveGroupFlagsWhenCalculatorsToggled) {
  auto window = AppWindow::create();
  correlation::app::AppBackend backend;

  correlation::app::ProgramOptions opts = backend.options();
  opts.active_calculators["g(r), J(r), G(r)"] = true;
  backend.setOptions(opts);

  correlation::app::AppController controller(*window, backend);

  controller.updateActiveGroupFlags();
  EXPECT_TRUE(window->get_has_radial_active());

  backend.setCalculatorActive("g(r), J(r), G(r)", false);
  controller.updateActiveGroupFlags();
  EXPECT_FALSE(window->get_has_radial_active());
}

TEST_F(AppControllerTests, HandlesCalculatorToggleSignal) {
  auto window = AppWindow::create();
  correlation::app::AppBackend backend;

  correlation::app::AppController controller(*window, backend);

  // Trigger toggle calculator signal from the UI
  window->invoke_toggle_calculator("g(r), J(r), G(r)", false);

  // Verify backend active state has been changed
  EXPECT_FALSE(backend.options().active_calculators.at("g(r), J(r), G(r)"));

  window->invoke_toggle_calculator("g(r), J(r), G(r)", true);
  EXPECT_TRUE(backend.options().active_calculators.at("g(r), J(r), G(r)"));
}
