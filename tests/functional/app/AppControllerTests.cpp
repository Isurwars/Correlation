// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "AppWindow.h"
#include "analysis/AnalysisTypes.hpp"
#include "app/core/AppController.hpp"
#include "app/core/AppOptions.hpp"
#include "app/services/AnalysisDispatcher.hpp"
#include "app/services/TrajectoryLoader.hpp"
#include "app/viewmodel/PlotController.hpp"
#include "app/viewmodel/Structure3DController.hpp"
#include <csignal>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string_view>

#include <gtest/gtest.h>

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#endif

namespace {

[[nodiscard]] auto resolveExamplePath(std::string_view relative_subpath) -> std::string {
  for (const auto &prefix : {"../../", "../", "./"}) {
    const std::string candidate = prefix + std::string(relative_subpath);
    if (std::filesystem::exists(candidate)) {
      return candidate;
    }
  }
  return "";
}

template <typename T>
[[nodiscard]] auto getRow(const std::shared_ptr<slint::Model<T>> &model, size_t index) -> T {
  const auto opt = model->row_data(index);
  if (!opt.has_value()) {
    throw std::runtime_error("Row index out of bounds or missing");
  }
  return *opt;
}

void verifyCutoffCell(const correlation::analysis::BondCutoffRange &range, double expected_min,
                      double expected_max) {
  EXPECT_NEAR(range.max_sq, expected_max, 1e-5);
  EXPECT_NEAR(range.min_sq, expected_min, 1e-5);
}

void verifyBondCutoffMatrix(const correlation::analysis::BondCutoffMatrix &cutoffs) {
  ASSERT_EQ(cutoffs.rows(), 2);
  ASSERT_EQ(cutoffs.cols(), 2);

  verifyCutoffCell(cutoffs[0][0], 0.65 * 0.65, 2.50 * 2.50);
  verifyCutoffCell(cutoffs[0][1], 0.70 * 0.70, 3.10 * 3.10);
  verifyCutoffCell(cutoffs[1][0], 0.70 * 0.70, 3.10 * 3.10);
  verifyCutoffCell(cutoffs[1][1], 0.60 * 0.60, 2.80 * 2.80);
}

auto findCutoffPair(const std::shared_ptr<slint::Model<BondCutoff>> &cutoffs, std::string_view el_a,
                    std::string_view el_b) -> std::optional<BondCutoff> {
  for (size_t i = 0; i < cutoffs->row_count(); ++i) {
    auto item = getRow(cutoffs, i);
    const std::string_view elem1(item.element1.data());
    const std::string_view elem2(item.element2.data());
    const bool direct_match = (elem1 == el_a && elem2 == el_b);
    const bool reverse_match = (elem1 == el_b && elem2 == el_a);
    if (direct_match || reverse_match) {
      return item;
    }
  }
  return std::nullopt;
}

void verifyPairDistances(const std::optional<BondCutoff> &item, std::string_view expected_min,
                         std::string_view expected_max) {
  ASSERT_TRUE(item.has_value());
  EXPECT_EQ(item->min_distance.data(), expected_min);
  EXPECT_EQ(item->max_distance.data(), expected_max);
}

void verifyRecommendedCutoffs(const std::shared_ptr<slint::Model<BondCutoff>> &cutoffs) {
  ASSERT_EQ(cutoffs->row_count(), 3);
  verifyPairDistances(findCutoffPair(cutoffs, "Pd", "Pd"), "1.44", "2.88");
  verifyPairDistances(findCutoffPair(cutoffs, "Pd", "Si"), "1.42", "2.83");
  verifyPairDistances(findCutoffPair(cutoffs, "Si", "Si"), "1.39", "2.78");
}

void verifyRestoredBondCutoffs(const std::shared_ptr<slint::Model<BondCutoff>> &restored) {
  ASSERT_EQ(restored->row_count(), 3);
  for (size_t i = 0; i < restored->row_count(); ++i) {
    const auto item = getRow(restored, i);
    if (item.element1 == "Pd" && item.element2 == "Pd") {
      EXPECT_EQ(item.min_distance.data(), std::string("1.44"));
      EXPECT_EQ(item.max_distance.data(), std::string("2.88"));
    }
  }
}

void verifyUniformCutoffs(const std::shared_ptr<slint::Model<BondCutoff>> &cutoffs,
                          std::string_view expected_max) {
  ASSERT_GT(cutoffs->row_count(), 0);
  for (size_t i = 0; i < cutoffs->row_count(); ++i) {
    const auto opt = getRow(cutoffs, i);
    EXPECT_EQ(opt.min_distance.data(), std::string("0.00"));
    EXPECT_EQ(opt.max_distance.data(), expected_max);
  }
}

void verifyCurveFilterMatches(const std::shared_ptr<slint::Model<CurveToggleItem>> &curves,
                              std::string_view pattern) {
  for (std::size_t i = 0; i < curves->row_count(); ++i) {
    const auto item = getRow(curves, i);
    const std::string label = std::string(item.label);
    if (label.contains(pattern)) {
      EXPECT_TRUE(item.matches_filter);
    } else {
      EXPECT_FALSE(item.matches_filter);
    }
  }
}

void verifyAllCurvesMatchFilter(const std::shared_ptr<slint::Model<CurveToggleItem>> &curves) {
  for (std::size_t i = 0; i < curves->row_count(); ++i) {
    const auto item = getRow(curves, i);
    EXPECT_TRUE(item.matches_filter);
  }
}

void verifyPlotStyling(const correlation::plotters::PlotConfig &config) {
  EXPECT_EQ(config.preset_size, correlation::plotters::PlotConfig::PresetSize::SingleColumn);
  EXPECT_EQ(config.palette, correlation::plotters::PlotConfig::Palette::Viridis);
  EXPECT_TRUE(config.show_markers);
  EXPECT_TRUE(config.fill_area);
}

void verifyPlotMetrics(const correlation::plotters::PlotConfig &config) {
  EXPECT_FLOAT_EQ(config.font_scale, 1.5F);
  EXPECT_FLOAT_EQ(config.line_width, 4.0F);
  EXPECT_FLOAT_EQ(config.marker_size, 5.0F);
}

void verifyActivePlotConfig(const correlation::plotters::PlotConfig &config) {
  verifyPlotStyling(config);
  verifyPlotMetrics(config);
}

} // namespace

// Test fixture for the AppController class
class AppControllerTests : public ::testing::Test {
public:
  [[nodiscard]] correlation::app::TrajectoryLoader &loader() noexcept { return loader_; }
  [[nodiscard]] const correlation::app::TrajectoryLoader &loader() const noexcept {
    return loader_;
  }

  [[nodiscard]] correlation::app::AnalysisDispatcher &dispatcher() noexcept { return dispatcher_; }
  [[nodiscard]] const correlation::app::AnalysisDispatcher &dispatcher() const noexcept {
    return dispatcher_;
  }

  [[nodiscard]] correlation::app::ProgramOptions &options() noexcept { return options_; }
  [[nodiscard]] const correlation::app::ProgramOptions &options() const noexcept {
    return options_;
  }

protected:
  void SetUp() override {
    // Ensure Slint uses the software backend during tests to avoid OpenGL
    // requirements in headless CI/Xvfb environments
#ifndef _WIN32
    [[maybe_unused]] const int setenv_res = setenv("SLINT_BACKEND", "software", 1);
    [[maybe_unused]] const bool sig_res = (std::signal(SIGPIPE, SIG_IGN) != SIG_ERR);
#else
    [[maybe_unused]] const errno_t putenv_res = _putenv_s("SLINT_BACKEND", "software");
#endif
  }

  static correlation::app::ProgramOptions
  callHandleOptionsfromUI(correlation::app::AppController &controller) {
    return controller.handleOptionsfromUI();
  }

private:
  correlation::app::TrajectoryLoader loader_;
  correlation::app::AnalysisDispatcher dispatcher_;
  correlation::app::ProgramOptions options_;
};

TEST_F(AppControllerTests, ConstructorInitializesCorrectly) {
  auto window = AppWindow::create();
  EXPECT_NO_THROW({
    const correlation::app::AppController controller(*window, loader(), dispatcher(), options());
  });
}

TEST_F(AppControllerTests, HandlesBondCutoffsCorrectly) {
  auto window = AppWindow::create();

  const std::string file_path = resolveExamplePath("examples/a-PdSi/a-PdSi.car");
  ASSERT_FALSE(file_path.empty());
  ASSERT_TRUE(loader().loadFile(file_path).has_value());

  correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  const auto atom_counts = std::make_shared<slint::VectorModel<AtomCount>>();
  atom_counts->push_back({
      .symbol = "Pd",
      .count = 80,
      .color = {},
      .is_unknown = false,
  });
  atom_counts->push_back({
      .symbol = "Si",
      .count = 20,
      .color = {},
      .is_unknown = false,
  });
  window->set_atom_counts(atom_counts);

  const auto cutoffs = std::make_shared<slint::VectorModel<BondCutoff>>();
  cutoffs->push_back({
      .element1 = "Pd",
      .element2 = "Pd",
      .min_distance = "0.60",
      .max_distance = "2.80",
  });
  cutoffs->push_back({
      .element1 = "Pd",
      .element2 = "Si",
      .min_distance = "0.70",
      .max_distance = "3.10",
  });
  cutoffs->push_back({
      .element1 = "Si",
      .element2 = "Si",
      .min_distance = "0.65",
      .max_distance = "2.50",
  });
  window->set_bond_cutoffs(cutoffs);

  const auto opts = callHandleOptionsfromUI(controller);
  verifyBondCutoffMatrix(opts.bond_cutoffs);
}

TEST_F(AppControllerTests, SynchronizesOptionsToAndFromUI) {
  auto window = AppWindow::create();

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

  options() = backend_opts;

  correlation::app::AppController controller(*window, loader(), dispatcher(), options());

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
  EXPECT_EQ(window->get_analysis_options().smoothing_kernel,
            static_cast<int>(correlation::math::KernelType::Gaussian));
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

  const auto parsed_opts = controller.handleOptionsfromUI();
  EXPECT_NEAR(parsed_opts.r_max, 14.20, 1e-5);
  EXPECT_EQ(parsed_opts.min_frame, 0);
  EXPECT_EQ(parsed_opts.max_frame, -1);
}

TEST_F(AppControllerTests, PopulatesRecommendedBondCutoffs) {
  auto window = AppWindow::create();

  const std::string file_path = resolveExamplePath("examples/a-PdSi/a-PdSi.car");
  ASSERT_FALSE(file_path.empty());
  ASSERT_TRUE(loader().loadFile(file_path).has_value());

  correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  controller.setBondCutoffs();

  const auto cutoffs = window->get_bond_cutoffs();
  verifyRecommendedCutoffs(cutoffs);
}

TEST_F(AppControllerTests, ResetsBondCutoffsToDefaultsWhenInvoked) {
  auto window = AppWindow::create();

  const std::string file_path = resolveExamplePath("examples/a-PdSi/a-PdSi.car");
  ASSERT_FALSE(file_path.empty());
  ASSERT_TRUE(loader().loadFile(file_path).has_value());

  const correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  // Set modified/custom cutoffs
  const auto cutoffs = std::make_shared<slint::VectorModel<BondCutoff>>();
  cutoffs->push_back(
      {.element1 = "Pd", .element2 = "Pd", .min_distance = "9.99", .max_distance = "9.99"});
  cutoffs->push_back(
      {.element1 = "Pd", .element2 = "Si", .min_distance = "9.99", .max_distance = "9.99"});
  cutoffs->push_back(
      {.element1 = "Si", .element2 = "Si", .min_distance = "9.99", .max_distance = "9.99"});
  window->set_bond_cutoffs(cutoffs);

  // Invoke reset callback
  window->invoke_reset_bond_cutoffs();

  // Verify recommended defaults were restored
  const auto restored = window->get_bond_cutoffs();
  verifyRestoredBondCutoffs(restored);
}

TEST_F(AppControllerTests, ResetsOptionBlocksToDefaultsWhenInvoked) {
  auto window = AppWindow::create();
  const correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  // Set modified/custom values for all option blocks
  auto opts = window->get_analysis_options();
  opts.r_max = "99.00";
  opts.r_bin_width = "1.23";
  opts.angle_bin_width = "45.00";
  opts.dihedral_bin_width = "60.00";
  opts.q_max = "50.00";
  opts.q_bin_width = "0.99";
  opts.r_int_max = "25.00";
  opts.max_ring_size = "20";
  opts.smoothing_enabled = false;
  opts.smoothing_sigma = "0.85";
  opts.smoothing_kernel = 2;
  opts.lef_cutoff = "12.00";
  opts.lef_sigma = "0.95";
  opts.hyper_samples = "500";
  opts.time_step = "5.00";
  opts.min_frame = "10";
  opts.max_frame = "20";
  window->set_analysis_options(opts);

  // 1. Reset RDF options
  window->invoke_reset_rdf_options();
  opts = window->get_analysis_options();
  EXPECT_EQ(opts.r_max, "20.00");
  EXPECT_EQ(opts.r_bin_width, "0.02");

  // 2. Reset Angle options
  window->invoke_reset_angle_options();
  opts = window->get_analysis_options();
  EXPECT_EQ(opts.angle_bin_width, "0.25");
  EXPECT_EQ(opts.dihedral_bin_width, "0.25");

  // 3. Reset SQ options
  window->invoke_reset_sq_options();
  opts = window->get_analysis_options();
  EXPECT_EQ(opts.q_max, "20.00");
  EXPECT_EQ(opts.q_bin_width, "0.02");
  EXPECT_EQ(opts.r_int_max, "10.00");

  // 4. Reset Rings options
  window->invoke_reset_rings_options();
  opts = window->get_analysis_options();
  EXPECT_EQ(opts.max_ring_size, "12");

  // 5. Reset Smoothing options
  window->invoke_reset_smoothing_options();
  opts = window->get_analysis_options();
  EXPECT_TRUE(opts.smoothing_enabled);
  EXPECT_EQ(opts.smoothing_sigma, "0.10");
  EXPECT_EQ(opts.smoothing_kernel, 0);

  // 6. Reset Advanced options
  window->invoke_reset_advanced_options();
  opts = window->get_analysis_options();
  EXPECT_EQ(opts.lef_cutoff, "5.00");
  EXPECT_EQ(opts.lef_sigma, "0.20");
  EXPECT_EQ(opts.hyper_samples, "10000");

  // 7. Reset Trajectory options
  window->invoke_reset_trajectory_options();
  opts = window->get_analysis_options();
  EXPECT_EQ(opts.time_step, "1.00");
  EXPECT_EQ(opts.min_frame, "1");
  EXPECT_EQ(opts.max_frame, "End");

  // 8. Reset Export options
  ExportConfig cfg;
  cfg.size_preset = 3;
  cfg.palette = 2;
  cfg.font_scale = "2.5";
  cfg.line_width = "5.0";
  cfg.marker_size = "8.0";
  cfg.show_legend = false;
  cfg.show_grid = false;
  cfg.show_markers = true;
  cfg.fill_area = true;
  cfg.manual_ticks = true;
  cfg.x_major_step = "2.0";
  cfg.y_major_step = "1.0";
  cfg.minor_subdivisions = 5;
  window->set_export_config(cfg);
  window->invoke_reset_export_settings();
  const auto restored_cfg = window->get_export_config();
  EXPECT_EQ(restored_cfg.size_preset, 0);
  EXPECT_EQ(restored_cfg.palette, 0);
  EXPECT_EQ(restored_cfg.font_scale, "1.0");
  EXPECT_EQ(restored_cfg.line_width, "3.0");
  EXPECT_EQ(restored_cfg.marker_size, "3.5");
  EXPECT_TRUE(restored_cfg.show_legend);
  EXPECT_TRUE(restored_cfg.show_grid);
  EXPECT_FALSE(restored_cfg.show_markers);
  EXPECT_FALSE(restored_cfg.fill_area);
  EXPECT_FALSE(restored_cfg.manual_ticks);
  EXPECT_EQ(restored_cfg.x_major_step, "1.0");
  EXPECT_EQ(restored_cfg.y_major_step, "0.5");
  EXPECT_EQ(restored_cfg.minor_subdivisions, 0);

  // 9. Reset Material Type
  opts.material_type = 2;
  window->set_analysis_options(opts);
  window->invoke_reset_material_type();
  EXPECT_EQ(window->get_analysis_options().material_type, 0);

  // 10. Reset Analyses Selection
  window->invoke_reset_analyses_selection();
  EXPECT_GT(window->get_calculator_groups()->row_count(), 0);

  // 11. Clear Comparison Curves
  const auto comp_curves = std::make_shared<slint::VectorModel<ComparisonCurveData>>();
  comp_curves->push_back({
      .id = 1,
      .label = "Run 1",
      .color_hex = {},
      .dash_style = 0,
      .visible = true,
      .file_path = "",
  });
  window->set_comparison_curves(comp_curves);
  EXPECT_EQ(window->get_comparison_curves()->row_count(), 1);
  window->invoke_clear_comparison_curves();
  EXPECT_EQ(window->get_comparison_curves()->row_count(), 0);

  // 12. Reset XRD options
  opts = window->get_analysis_options();
  opts.xrd_lambda = "0.7107";
  opts.xrd_theta_min = "5.0";
  opts.xrd_theta_max = "90.0";
  opts.xrd_bin_width = "0.01";
  opts.xrd_radiation_preset = 1;
  window->set_analysis_options(opts);
  window->invoke_reset_xrd_options();
  opts = window->get_analysis_options();
  EXPECT_EQ(opts.xrd_lambda, "1.5406");
  EXPECT_EQ(opts.xrd_theta_min, "10.0");
  EXPECT_EQ(opts.xrd_theta_max, "140.0");
  EXPECT_EQ(opts.xrd_bin_width, "0.05");
  EXPECT_EQ(opts.xrd_radiation_preset, 0);
}

TEST_F(AppControllerTests, HandlesXRDOptionsAndCutoffHeuristics) {
  auto window = AppWindow::create();

  const std::string file_path = resolveExamplePath("examples/a-PdSi/a-PdSi.car");
  ASSERT_FALSE(file_path.empty());
  ASSERT_TRUE(loader().loadFile(file_path).has_value());

  const correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  // Test XRD preset changed callback
  window->invoke_xrd_radiation_preset_changed(1); // Mo-Kalpha (0.7107)
  auto opts = window->get_analysis_options();
  EXPECT_EQ(opts.xrd_radiation_preset, 1);
  EXPECT_EQ(opts.xrd_lambda, "0.7107");

  window->invoke_xrd_radiation_preset_changed(2); // Co-Kalpha (1.7890)
  opts = window->get_analysis_options();
  EXPECT_EQ(opts.xrd_radiation_preset, 2);
  EXPECT_EQ(opts.xrd_lambda, "1.7890");

  window->invoke_xrd_radiation_preset_changed(3); // Cr-Kalpha (2.2897)
  opts = window->get_analysis_options();
  EXPECT_EQ(opts.xrd_radiation_preset, 3);
  EXPECT_EQ(opts.xrd_lambda, "2.2897");

  // Test Cutoff Heuristics
  window->invoke_apply_scaled_cutoffs(1.15F);
  auto cutoffs = window->get_bond_cutoffs();
  ASSERT_GT(cutoffs->row_count(), 0);

  window->invoke_set_uniform_cutoff(3.4F);
  cutoffs = window->get_bond_cutoffs();
  verifyUniformCutoffs(cutoffs, "3.40");

  // Test Min, Max, and Global Cutoff signals
  window->invoke_apply_min_factor(0.6F);
  cutoffs = window->get_bond_cutoffs();
  ASSERT_GT(cutoffs->row_count(), 0);

  window->invoke_apply_max_factor(1.2F);
  cutoffs = window->get_bond_cutoffs();
  ASSERT_GT(cutoffs->row_count(), 0);

  window->invoke_apply_global_cutoff(3.5F);
  cutoffs = window->get_bond_cutoffs();
  verifyUniformCutoffs(cutoffs, "3.50");
}

TEST_F(AppControllerTests, UpdatesActiveGroupFlagsWhenCalculatorsToggled) {
  auto window = AppWindow::create();

  options().active_calculators["g(r), J(r), G(r)"] = true;

  correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  controller.updateActiveGroupFlags();
  EXPECT_TRUE(window->get_has_radial_active());

  options().active_calculators["g(r), J(r), G(r)"] = false;
  controller.updateActiveGroupFlags();
  EXPECT_FALSE(window->get_has_radial_active());
}

TEST_F(AppControllerTests, HandlesCalculatorToggleSignal) {
  auto window = AppWindow::create();

  const correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  // Trigger toggle calculator signal from the UI
  window->invoke_toggle_calculator("g(r), J(r), G(r)", false);

  // Verify options active state has been changed
  ASSERT_TRUE(options().active_calculators.contains("g(r), J(r), G(r)"));
  EXPECT_FALSE(options().active_calculators.at("g(r), J(r), G(r)"));

  window->invoke_toggle_calculator("g(r), J(r), G(r)", true);
  ASSERT_TRUE(options().active_calculators.contains("g(r), J(r), G(r)"));
  EXPECT_TRUE(options().active_calculators.at("g(r), J(r), G(r)"));
}

TEST_F(AppControllerTests, PopulatesTableAndDynamicProperties) {
  auto window = AppWindow::create();

  std::string file_path = "../../examples/a-PdSi/a-PdSi.car";
  if (!std::filesystem::exists(file_path)) {
    file_path = "../examples/a-PdSi/a-PdSi.car";
  }
  if (!std::filesystem::exists(file_path)) {
    file_path = "examples/a-PdSi/a-PdSi.car";
  }
  ASSERT_TRUE(loader().loadFile(file_path).has_value());

  // Set the RDF calculator active
  options().active_calculators["g(r), J(r), G(r)"] = true;

  correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  // Run analysis
  EXPECT_TRUE(dispatcher().runAnalysis(loader().trajectoryMut(), options()).has_value());

  // Populate plot list (which also sets dynamic properties)
  controller.getPlotController()->populatePlotList();

  // Verify plot items is not empty
  const auto plot_items = window->get_plot_items();
  ASSERT_GT(plot_items->row_count(), 0);

  // Trigger plot update for index 0
  controller.getPlotController()->requestPlotUpdate(0, true);

  // Verify table headers are set
  const auto headers = window->get_table_headers();
  ASSERT_GT(headers->row_count(), 0);

  // Verify table rows are set
  const auto rows = window->get_table_rows();
  ASSERT_GT(rows->row_count(), 0);
}

TEST_F(AppControllerTests, GuiLaunchAndEventLoopSmokeTest) {
  auto window = AppWindow::create();
  const correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  slint::Timer quit_timer;
  quit_timer.start(slint::TimerMode::SingleShot, std::chrono::milliseconds(50),
                   [] { slint::quit_event_loop(); });

  EXPECT_NO_THROW(window->run());
}

TEST_F(AppControllerTests, MiniSideBarPanelNavigationAndCollapse) {
  auto window = AppWindow::create();
  const correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  // Initial state defaults to panel 0 (Operations)
  EXPECT_EQ(window->get_active_panel(), 0);

  // Switch across all 6 panels
  for (int panel = 0; panel <= 5; ++panel) {
    window->set_active_panel(panel);
    EXPECT_EQ(window->get_active_panel(), panel);
  }

  // Collapse drawer (-1)
  window->set_active_panel(-1);
  EXPECT_EQ(window->get_active_panel(), -1);
}

namespace {
[[nodiscard]] auto findPdSiExample() -> std::string {
  for (const auto *path : {
           "../../examples/a-PdSi/a-PdSi.car",
           "../examples/a-PdSi/a-PdSi.car",
           "examples/a-PdSi/a-PdSi.car",
       }) {
    if (std::filesystem::exists(path)) {
      return path;
    }
  }
  return "";
}
} // namespace

TEST_F(AppControllerTests, ResetRdfAndAngleCallbacksIntegration) {
  auto window = AppWindow::create();
  const correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  auto opts = window->get_analysis_options();
  opts.r_max = "99.00";
  opts.r_bin_width = "1.00";
  opts.angle_bin_width = "99.00";
  window->set_analysis_options(opts);

  window->invoke_reset_rdf_options();
  EXPECT_EQ(window->get_analysis_options().r_max, "20.00");
  EXPECT_EQ(window->get_analysis_options().r_bin_width, "0.02");

  window->invoke_reset_angle_options();
  EXPECT_EQ(window->get_analysis_options().angle_bin_width, "0.25");
}

TEST_F(AppControllerTests, ResetSqRingsAndSmoothingCallbacksIntegration) {
  auto window = AppWindow::create();
  const correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  auto opts = window->get_analysis_options();
  opts.q_max = "99.00";
  opts.max_ring_size = "99";
  opts.smoothing_sigma = "99.00";
  opts.material_type = 2;
  window->set_analysis_options(opts);

  window->invoke_reset_sq_options();
  EXPECT_EQ(window->get_analysis_options().q_max, "20.00");

  window->invoke_reset_rings_options();
  EXPECT_EQ(window->get_analysis_options().max_ring_size, "12");

  window->invoke_reset_material_type();
  EXPECT_EQ(window->get_analysis_options().material_type, 0);

  window->invoke_reset_smoothing_options();
  EXPECT_EQ(window->get_analysis_options().smoothing_sigma, "0.10");
}

TEST_F(AppControllerTests, XRDRadiationPresetCallbackIntegration) {
  auto window = AppWindow::create();
  const correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  window->invoke_xrd_radiation_preset_changed(1);
  EXPECT_EQ(window->get_analysis_options().xrd_lambda, "0.7107");
}

TEST_F(AppControllerTests, BondCutoffsUniformCallbackIntegration) {
  auto window = AppWindow::create();
  const auto file_path = findPdSiExample();
  ASSERT_FALSE(file_path.empty());
  ASSERT_TRUE(loader().loadFile(file_path).has_value());

  correlation::app::AppController controller(*window, loader(), dispatcher(), options());
  controller.setBondCutoffs();

  window->invoke_set_uniform_cutoff(3.5F);
  const auto cutoffs_uniform = window->get_bond_cutoffs();
  ASSERT_GT(cutoffs_uniform->row_count(), 0);
  for (size_t idx = 0; idx < cutoffs_uniform->row_count(); ++idx) {
    const auto item = cutoffs_uniform->row_data(idx);
    ASSERT_TRUE(item.has_value());
    EXPECT_EQ(item->max_distance, "3.50");
  }
}

TEST_F(AppControllerTests, BondCutoffsGlobalCallbackIntegration) {
  auto window = AppWindow::create();
  const auto file_path = findPdSiExample();
  ASSERT_FALSE(file_path.empty());
  ASSERT_TRUE(loader().loadFile(file_path).has_value());

  correlation::app::AppController controller(*window, loader(), dispatcher(), options());
  controller.setBondCutoffs();

  window->invoke_apply_global_cutoff(4.0F);
  const auto cutoffs_global = window->get_bond_cutoffs();
  ASSERT_GT(cutoffs_global->row_count(), 0);
  for (size_t idx = 0; idx < cutoffs_global->row_count(); ++idx) {
    const auto item = cutoffs_global->row_data(idx);
    ASSERT_TRUE(item.has_value());
    EXPECT_EQ(item->max_distance, "4.00");
  }
}

TEST_F(AppControllerTests, ExportAnalysesToggleCallbackIntegration) {
  auto window = AppWindow::create();
  const correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  const auto model = std::make_shared<slint::VectorModel<ExportAnalysisItem>>();
  ExportAnalysisItem item;
  item.id = "rdf";
  item.name = "RDF";
  item.enabled = true;
  model->push_back(item);
  window->set_export_analysis_items(model);

  window->invoke_toggle_all_export_analyses(false);
  const auto export_model = window->get_export_analysis_items();
  ASSERT_TRUE(export_model != nullptr);
  ASSERT_GT(export_model->row_count(), 0);
  auto row = export_model->row_data(0);
  ASSERT_TRUE(row.has_value());
  EXPECT_FALSE(row->enabled);

  window->invoke_toggle_export_analysis("rdf", true);
  row = export_model->row_data(0);
  ASSERT_TRUE(row.has_value());
  EXPECT_TRUE(row->enabled);
}

TEST_F(AppControllerTests, ExportFormatAndLayoutCallbacksIntegration) {
  auto window = AppWindow::create();
  const correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  window->invoke_reset_export_options();
  const auto export_config = window->get_export_format_options();
  EXPECT_TRUE(export_config.export_csv);
  EXPECT_TRUE(export_config.export_images);
  EXPECT_FALSE(export_config.export_hdf5);

  EXPECT_NO_THROW(window->invoke_layout_geometry_changed(280.0F, 380.0F));
}

namespace {
/// Finds the row index of a calculator group by name in the window's group model.
[[nodiscard]] std::optional<size_t> findGroupRow(const AppWindow &window, std::string_view name) {
  const auto groups = window.get_calculator_groups();
  for (size_t idx = 0; groups != nullptr && idx < groups->row_count(); ++idx) {
    const auto row = groups->row_data(idx);
    if (row.has_value() && std::string_view(row->name.data()) == name) {
      return idx;
    }
  }
  return std::nullopt;
}

[[nodiscard]] CalculatorGroup groupAt(const AppWindow &window, size_t idx) {
  const auto row = window.get_calculator_groups()->row_data(idx);
  if (!row.has_value()) {
    throw std::runtime_error("Calculator group row missing");
  }
  return *row;
}
} // namespace

TEST_F(AppControllerTests, CalculatorToggleDoesNotCollapseExpandedSubTree) {
  auto window = AppWindow::create();
  const correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  const auto scattering = findGroupRow(*window, "Scattering");
  ASSERT_TRUE(scattering.has_value());
  const auto model_before = window->get_calculator_groups();

  window->invoke_set_calculator_group_expanded("Scattering", true);
  ASSERT_TRUE(groupAt(*window, *scattering).expanded);

  const auto group = groupAt(*window, *scattering);
  ASSERT_GT(group.calculators->row_count(), 0U);
  const auto calc = getRow(group.calculators, 0);
  window->invoke_toggle_calculator(calc.id, !calc.enabled);

  // Sub-tree stays open and the root model is updated in place (not replaced).
  EXPECT_TRUE(groupAt(*window, *scattering).expanded);
  EXPECT_EQ(window->get_calculator_groups(), model_before);
}

TEST_F(AppControllerTests, GroupTriStateToggleKeepsExpansionAndUpdatesCounts) {
  auto window = AppWindow::create();
  const correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  const auto radial = findGroupRow(*window, "Radial");
  ASSERT_TRUE(radial.has_value());
  window->invoke_set_calculator_group_expanded("Radial", true);

  window->invoke_toggle_calculator_group("Radial", false);
  auto group = groupAt(*window, *radial);
  EXPECT_TRUE(group.expanded);
  EXPECT_EQ(group.enabled_count, 0);

  window->invoke_toggle_calculator_group("Radial", true);
  group = groupAt(*window, *radial);
  EXPECT_TRUE(group.expanded);
  EXPECT_EQ(group.enabled_count, group.total_count);
}

TEST_F(AppControllerTests, SelectAllAndSelectNoneCalculators) {
  auto window = AppWindow::create();
  const correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  window->invoke_toggle_all_calculators(false);
  EXPECT_EQ(window->get_enabled_calculator_count(), 0);

  window->invoke_toggle_all_calculators(true);
  EXPECT_EQ(window->get_enabled_calculator_count(), window->get_total_calculator_count());
  EXPECT_GT(window->get_total_calculator_count(), 0);
}

TEST_F(AppControllerTests, ExpandAndCollapseAllCalculatorGroups) {
  auto window = AppWindow::create();
  const correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  const auto groups = window->get_calculator_groups();
  ASSERT_GT(groups->row_count(), 0U);

  window->invoke_set_all_calculator_groups_expanded(true);
  for (size_t idx = 0; idx < groups->row_count(); ++idx) {
    EXPECT_TRUE(groupAt(*window, idx).expanded);
  }

  window->invoke_set_all_calculator_groups_expanded(false);
  for (size_t idx = 0; idx < groups->row_count(); ++idx) {
    EXPECT_FALSE(groupAt(*window, idx).expanded);
  }
}

TEST_F(AppControllerTests, PlotZoomLimitsSynchronizeToUIAndReset) {
  auto window = AppWindow::create();

  const std::string sample_file = resolveExamplePath("examples/a-PdSi/a-PdSi.car");
  ASSERT_FALSE(sample_file.empty());
  ASSERT_TRUE(loader().loadFile(sample_file).has_value());

  options().active_calculators["g(r), J(r), G(r)"] = true;
  correlation::app::AppController controller(*window, loader(), dispatcher(), options());
  EXPECT_TRUE(dispatcher().runAnalysis(loader().trajectoryMut(), options()).has_value());

  auto *plot_ctrl = controller.getPlotController();
  ASSERT_NE(plot_ctrl, nullptr);
  plot_ctrl->populatePlotList();
  window->set_selected_plot_index(0);
  plot_ctrl->requestPlotUpdate(0, true);

  // Initial axis settings are auto
  const auto initial_settings = window->get_axis_settings();
  EXPECT_TRUE(initial_settings.auto_x);
  EXPECT_TRUE(initial_settings.auto_y);

  // Apply a zoom rect (x: 100->300, y: 100->300)
  plot_ctrl->handleZoomRect(100.0F, 100.0F, 300.0F, 300.0F);

  // Viewport should have synced into axis_settings
  const auto zoomed_settings = window->get_axis_settings();
  EXPECT_FALSE(zoomed_settings.auto_x);
  EXPECT_FALSE(zoomed_settings.auto_y);
  EXPECT_FALSE(zoomed_settings.x_min.empty());
  EXPECT_FALSE(zoomed_settings.x_max.empty());

  // Reset zoom
  plot_ctrl->handleResetZoom();
  const auto reset_settings = window->get_axis_settings();
  EXPECT_TRUE(reset_settings.auto_x);
  EXPECT_TRUE(reset_settings.auto_y);
  EXPECT_TRUE(reset_settings.x_min.empty());
  EXPECT_TRUE(reset_settings.x_max.empty());
}

TEST_F(AppControllerTests, PlotCopyToClipboardSetsPayloadAndStatus) {
  auto window = AppWindow::create();

  const std::string sample_file = resolveExamplePath("examples/a-PdSi/a-PdSi.car");
  ASSERT_FALSE(sample_file.empty());
  ASSERT_TRUE(loader().loadFile(sample_file).has_value());

  options().active_calculators["g(r), J(r), G(r)"] = true;
  correlation::app::AppController controller(*window, loader(), dispatcher(), options());
  EXPECT_TRUE(dispatcher().runAnalysis(loader().trajectoryMut(), options()).has_value());

  auto *plot_ctrl = controller.getPlotController();
  ASSERT_NE(plot_ctrl, nullptr);
  plot_ctrl->populatePlotList();
  window->set_selected_plot_index(0);
  plot_ctrl->requestPlotUpdate(0, true);

  plot_ctrl->handleCopyPlotToClipboard();

  const std::string payload = std::string(window->get_clipboard_payload());
  EXPECT_FALSE(payload.empty());
  EXPECT_NE(payload.find("<svg"), std::string::npos);
  EXPECT_NE(std::string(window->get_analysis_status_text()).find("copied to clipboard"),
            std::string::npos);
}

TEST_F(AppControllerTests, PlotPeakTaggingAndCoordinationNumberIntegration) {
  auto window = AppWindow::create();

  const std::string sample_file = resolveExamplePath("examples/a-PdSi/a-PdSi.car");
  ASSERT_FALSE(sample_file.empty());
  ASSERT_TRUE(loader().loadFile(sample_file).has_value());

  options().active_calculators["g(r), J(r), G(r)"] = true;
  correlation::app::AppController controller(*window, loader(), dispatcher(), options());
  EXPECT_TRUE(dispatcher().runAnalysis(loader().trajectoryMut(), options()).has_value());

  auto *plot_ctrl = controller.getPlotController();
  ASSERT_NE(plot_ctrl, nullptr);
  plot_ctrl->populatePlotList();
  window->set_selected_plot_index(0);
  plot_ctrl->requestPlotUpdate(0, true);

  // Simulate range tagging across the first coordination shell
  window->invoke_tag_peak_range(120.0F, 240.0F);

  const auto peaks_model = window->get_tagged_peaks();
  ASSERT_EQ(peaks_model->row_count(), 1U);
  const auto peak_item = getRow(peaks_model, 0);
  EXPECT_FALSE(std::string(peak_item.r_peak).empty());
  EXPECT_FALSE(std::string(peak_item.coordination_num).empty());
  EXPECT_GE(peak_item.screen_r_peak, 0.0F);

  // Clear peaks
  window->invoke_clear_tagged_peaks();
  EXPECT_EQ(window->get_tagged_peaks()->row_count(), 0U);
}

TEST_F(AppControllerTests, CurveFilterAndInvertSelectionIntegration) {
  auto window = AppWindow::create();

  const std::string sample_file = resolveExamplePath("examples/a-PdSi/a-PdSi.car");
  ASSERT_FALSE(sample_file.empty());
  ASSERT_TRUE(loader().loadFile(sample_file).has_value());

  options().active_calculators["g(r), J(r), G(r)"] = true;
  correlation::app::AppController controller(*window, loader(), dispatcher(), options());
  EXPECT_TRUE(dispatcher().runAnalysis(loader().trajectoryMut(), options()).has_value());

  auto *plot_ctrl = controller.getPlotController();
  ASSERT_NE(plot_ctrl, nullptr);
  plot_ctrl->populatePlotList();
  window->set_selected_plot_index(0);
  plot_ctrl->requestPlotUpdate(0, true);

  auto curves = window->get_curve_toggle_items();
  ASSERT_GT(curves->row_count(), 1U);

  // Apply search filter "Si"
  window->invoke_filter_curves(slint::SharedString("Si"));
  curves = window->get_curve_toggle_items();
  verifyCurveFilterMatches(curves, "Si");

  // Invert matching curves
  window->invoke_invert_curves();
  curves = window->get_curve_toggle_items();
  ASSERT_TRUE(curves->row_data(0).has_value());

  // Clear filter
  window->invoke_filter_curves(slint::SharedString(""));
  curves = window->get_curve_toggle_items();
  verifyAllCurvesMatchFilter(curves);
}

TEST_F(AppControllerTests, ConsolidatedPlotSettingsSynchronizationIntegration) {
  auto window = AppWindow::create();

  const std::string sample_file = resolveExamplePath("examples/a-PdSi/a-PdSi.car");
  ASSERT_FALSE(sample_file.empty());
  ASSERT_TRUE(loader().loadFile(sample_file).has_value());

  options().active_calculators["g(r), J(r), G(r)"] = true;
  correlation::app::AppController controller(*window, loader(), dispatcher(), options());
  EXPECT_TRUE(dispatcher().runAnalysis(loader().trajectoryMut(), options()).has_value());

  auto *plot_ctrl = controller.getPlotController();
  ASSERT_NE(plot_ctrl, nullptr);
  plot_ctrl->populatePlotList();
  window->set_selected_plot_index(0);
  plot_ctrl->requestPlotUpdate(0, true);

  // Modify publication export styling (as set via CurveSidebar Tab 1)
  auto cfg = window->get_export_config();
  cfg.size_preset = 1; // Single Column
  cfg.palette = 2;     // Viridis
  cfg.font_scale = slint::SharedString("1.5");
  cfg.line_width = slint::SharedString("4.0");
  cfg.marker_size = slint::SharedString("5.0");
  cfg.show_markers = true;
  cfg.fill_area = true;
  window->set_export_config(cfg);

  // Apply settings through UI callback
  window->invoke_apply_axis_settings();

  // Verify plot config reflects consolidated publication settings
  const auto active_config = plot_ctrl->buildPlotConfigFromUI();
  verifyActivePlotConfig(active_config);

  // Reset export settings through UI callback
  window->invoke_reset_export_settings();
  const auto reset_cfg = window->get_export_config();
  EXPECT_EQ(reset_cfg.size_preset, 0);
  EXPECT_EQ(reset_cfg.palette, 0);
  EXPECT_FALSE(reset_cfg.show_markers);
}

TEST_F(AppControllerTests, DualViewAndTabOrderLifecycleIntegration) {
  auto window = AppWindow::create();
  correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  auto *structure_ctrl = controller.getStructure3DController();
  ASSERT_NE(structure_ctrl, nullptr);

  // Tab 0: 2D Plot (is3DModeActive == false)
  window->set_preview_active_tab(0);
  window->invoke_preview_tab_changed(0);
  EXPECT_EQ(window->get_preview_active_tab(), 0);
  EXPECT_FALSE(structure_ctrl->is3DModeActive());

  // Tab 1: 3D View (is3DModeActive == true)
  window->set_preview_active_tab(1);
  window->invoke_preview_tab_changed(1);
  EXPECT_EQ(window->get_preview_active_tab(), 1);
  EXPECT_TRUE(structure_ctrl->is3DModeActive());

  // Tab 2: Dual View (is3DModeActive == true)
  window->set_preview_active_tab(2);
  window->invoke_preview_tab_changed(2);
  EXPECT_EQ(window->get_preview_active_tab(), 2);
  EXPECT_TRUE(structure_ctrl->is3DModeActive());
  EXPECT_FLOAT_EQ(window->get_preview_dual_split_ratio(), 0.5F);

  // Modify dual split ratio (e.g. user drags SplitterHandle)
  window->set_preview_dual_split_ratio(0.65F);
  EXPECT_FLOAT_EQ(window->get_preview_dual_split_ratio(), 0.65F);

  // Tab 3: Data Explorer (is3DModeActive == false)
  window->set_preview_active_tab(3);
  window->invoke_preview_tab_changed(3);
  EXPECT_EQ(window->get_preview_active_tab(), 3);
  EXPECT_FALSE(structure_ctrl->is3DModeActive());
}

TEST_F(AppControllerTests, CameraPresetsOrientationIntegration) {
  auto window = AppWindow::create();
  correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  auto *structure_ctrl = controller.getStructure3DController();
  ASSERT_NE(structure_ctrl, nullptr);

  // Preset 0: Top (XY)
  window->invoke_set_camera_preset_3d(0);
  EXPECT_FLOAT_EQ(structure_ctrl->camera().azimuth(), 0.0F);
  EXPECT_FLOAT_EQ(structure_ctrl->camera().elevation(), 89.0F);

  // Preset 1: Front (XZ)
  window->invoke_set_camera_preset_3d(1);
  EXPECT_FLOAT_EQ(structure_ctrl->camera().azimuth(), 0.0F);
  EXPECT_FLOAT_EQ(structure_ctrl->camera().elevation(), 0.0F);

  // Preset 2: Side (YZ)
  window->invoke_set_camera_preset_3d(2);
  EXPECT_FLOAT_EQ(structure_ctrl->camera().azimuth(), 90.0F);
  EXPECT_FLOAT_EQ(structure_ctrl->camera().elevation(), 0.0F);

  // Preset 3: Isometric
  window->invoke_set_camera_preset_3d(3);
  EXPECT_FLOAT_EQ(structure_ctrl->camera().azimuth(), 45.0F);
  EXPECT_NEAR(structure_ctrl->camera().elevation(), 35.26F, 0.01F);
}

TEST_F(AppControllerTests, SampleDatasetAndRecentFilesIntegration) {
  auto window = AppWindow::create();
  correlation::app::AppController controller(*window, loader(), dispatcher(), options());

  // Verify recent files recording and UI sync
  controller.addRecentFile("test_path_1.car");
  controller.addRecentFile("test_path_2.car");
  auto recents = window->get_recent_files();
  ASSERT_EQ(recents->row_count(), 2);
  EXPECT_EQ(getRow(recents, 0), slint::SharedString("test_path_2.car"));
  EXPECT_EQ(getRow(recents, 1), slint::SharedString("test_path_1.car"));

  // Verify adding duplicate file promotes it to front without duplicating
  controller.addRecentFile("test_path_1.car");
  recents = window->get_recent_files();
  ASSERT_EQ(recents->row_count(), 2);
  EXPECT_EQ(getRow(recents, 0), slint::SharedString("test_path_1.car"));
  EXPECT_EQ(getRow(recents, 1), slint::SharedString("test_path_2.car"));
}

TEST_F(AppControllerTests, DualViewPerFrameOverlayIntegration) {
  auto window = AppWindow::create();

  const std::string file_path = resolveExamplePath("examples/a-PdSi/a-PdSi.car");
  ASSERT_FALSE(file_path.empty());
  ASSERT_TRUE(loader().loadFile(file_path).has_value());

  options().active_calculators["g(r), J(r), G(r)"] = true;

  correlation::app::AppController controller(*window, loader(), dispatcher(), options());
  EXPECT_TRUE(dispatcher().runAnalysis(loader().trajectoryMut(), options()).has_value());

  auto *plot_ctrl = controller.getPlotController();
  auto *structure_ctrl = controller.getStructure3DController();
  ASSERT_NE(plot_ctrl, nullptr);
  ASSERT_NE(structure_ctrl, nullptr);

  const auto *dist_funcs = dispatcher().getDistributionFunctions();
  ASSERT_NE(dist_funcs, nullptr);

  // Set up 2-frame distribution store
  const auto store = std::make_shared<correlation::analysis::PerFrameDistributionStore>(2);
  const auto &bins = dist_funcs->getAllHistograms().at("g_r").bins;
  const auto &mean_curve = dist_funcs->getAllHistograms().at("g_r").partials.at("Total");
  store->setBins(correlation::analysis::DistributionMetric::RDF, bins);
  store->setMeanCurve(correlation::analysis::DistributionMetric::RDF, mean_curve);
  store->setFrameCurve(correlation::analysis::DistributionMetric::RDF, 0, mean_curve);
  store->setFrameCurve(correlation::analysis::DistributionMetric::RDF, 1, mean_curve);
  structure_ctrl->setDistributionStore(store);

  plot_ctrl->populatePlotList();

  // Tab 0: Standard 2D Plot view with Total and partials
  window->set_preview_active_tab(0);
  plot_ctrl->requestPlotUpdate(0, true);
  const auto toggle_items_tab0 = window->get_curve_toggle_items();
  ASSERT_GT(toggle_items_tab0->row_count(), 1);
  EXPECT_EQ(getRow(toggle_items_tab0, 0).label, slint::SharedString("Total"));

  // Tab 2: Dual View with Ensemble Mean (accented) and Frame 1 (faded)
  window->set_preview_active_tab(2);
  window->set_structure_3d_current_frame(0);
  plot_ctrl->requestPlotUpdate(0, true);
  const auto toggle_items_tab2 = window->get_curve_toggle_items();
  ASSERT_EQ(toggle_items_tab2->row_count(), 2);
  EXPECT_EQ(getRow(toggle_items_tab2, 0).label, slint::SharedString("Ensemble Mean"));
  EXPECT_EQ(getRow(toggle_items_tab2, 1).label, slint::SharedString("Frame 1"));

  // Seek to frame 1 in 3D viewport -> updates Dual View to Frame 2
  window->invoke_seek_frame_3d(1);
  plot_ctrl->requestPlotUpdate(0, true);
  const auto toggle_items_frame1 = window->get_curve_toggle_items();
  ASSERT_EQ(toggle_items_frame1->row_count(), 2);
  EXPECT_EQ(getRow(toggle_items_frame1, 0).label, slint::SharedString("Ensemble Mean"));
  EXPECT_EQ(getRow(toggle_items_frame1, 1).label, slint::SharedString("Frame 2"));

  // Switching back to Tab 0 restores standard partials
  window->set_preview_active_tab(0);
  window->invoke_preview_tab_changed(0);
  const auto restored_items = window->get_curve_toggle_items();
  ASSERT_GT(restored_items->row_count(), 1);
  EXPECT_EQ(getRow(restored_items, 0).label, slint::SharedString("Total"));
}
