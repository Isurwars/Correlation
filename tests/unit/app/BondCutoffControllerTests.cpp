/**
 * @file BondCutoffControllerTests.cpp
 * @brief Unit tests for BondCutoffController.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "AppWindow.h"
#include "analysis/DistributionFunctions.hpp"
#include "app/core/AppOptions.hpp"
#include "app/services/AnalysisDispatcher.hpp"
#include "app/services/TrajectoryLoader.hpp"
#include "app/viewmodel/BondCutoffController.hpp"

#include <filesystem>
#include <gtest/gtest.h>
#include <optional>
#include <stdexcept>
#include <string>

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#endif

namespace correlation::app {
namespace {

std::string getTestDataDir() {
  const std::vector<std::string> candidates = {
      "../tests/data/",
      "tests/data/",
      "data/",
      "../../tests/data/",
  };
  for (const auto &dir : candidates) {
    if (std::filesystem::exists(dir + "xyz/clean.xyz")) {
      return dir;
    }
  }
  return "../tests/data/";
}

template <typename T> [[nodiscard]] T requireValue(const std::optional<T> &opt) {
  if (!opt.has_value()) {
    throw std::runtime_error("Optional value is unexpectedly empty");
  }
  return *opt;
}

class BondCutoffControllerTests : public ::testing::Test {
public:
  [[nodiscard]] AppWindow &window() {
    if (!window_.has_value()) {
      throw std::runtime_error("Window is not initialized");
    }
    return **window_;
  }

protected:
  void SetUp() override {
#ifndef _WIN32
    setenv("SLINT_BACKEND", "software", 1);
#else
    _putenv_s("SLINT_BACKEND", "software");
#endif
    window_.emplace(AppWindow::create());
  }

private:
  std::optional<slint::ComponentHandle<AppWindow>> window_;
};

TEST_F(BondCutoffControllerTests, HandlesNullCellGracefully) {
  auto &win = window();
  TrajectoryLoader loader;
  ProgramOptions options;
  BondCutoffController controller(win, loader, options);

  EXPECT_NO_THROW(controller.setBondCutoffs());
  const auto matrix = controller.getBondCutoffs();
  EXPECT_TRUE(matrix.empty());

  EXPECT_NO_THROW(controller.applyScaledCutoffs(1.2F));
  EXPECT_NO_THROW(controller.setUniformCutoff(2.5F));
  EXPECT_NO_THROW(controller.applyCovalentFactor(1.1F, FactorBound::Min));
  EXPECT_NO_THROW(controller.applyCovalentFactor(1.1F, FactorBound::Max));
  EXPECT_NO_THROW(controller.applyGlobalCutoff(3.0F));
}

TEST_F(BondCutoffControllerTests, IgnoresInvalidScaleAndCutoffFactors) {
  auto &win = window();
  TrajectoryLoader loader;
  ProgramOptions options;
  const std::string filepath = getTestDataDir() + "xyz/clean.xyz";
  const auto res = loader.loadFile(filepath);
  ASSERT_TRUE(res.has_value());
  ASSERT_NE(loader.cell(), nullptr);

  BondCutoffController controller(win, loader, options);
  controller.setBondCutoffs();
  const int initial_trigger = win.get_bond_cutoffs_reset_trigger();

  controller.applyScaledCutoffs(-1.0F);
  controller.applyScaledCutoffs(0.0F);
  controller.setUniformCutoff(-2.0F);
  controller.applyCovalentFactor(0.0F, FactorBound::Min);
  controller.applyCovalentFactor(-0.5F, FactorBound::Max);

  EXPECT_EQ(win.get_bond_cutoffs_reset_trigger(), initial_trigger);
}

TEST_F(BondCutoffControllerTests, PopulatesAndParsesCutoffsWithLoadedCell) {
  auto &win = window();
  TrajectoryLoader loader;
  ProgramOptions options;
  const std::string filepath = getTestDataDir() + "xyz/clean.xyz";
  const auto res = loader.loadFile(filepath);
  ASSERT_TRUE(res.has_value());
  ASSERT_NE(loader.cell(), nullptr);

  BondCutoffController controller(win, loader, options);
  const int initial_trigger = win.get_bond_cutoffs_reset_trigger();

  controller.setBondCutoffs();
  EXPECT_GT(win.get_bond_cutoffs_reset_trigger(), initial_trigger);

  const auto cutoffs_model = win.get_bond_cutoffs();
  ASSERT_NE(cutoffs_model, nullptr);
  EXPECT_GT(cutoffs_model->row_count(), 0);

  const auto matrix = controller.getBondCutoffs();
  EXPECT_FALSE(matrix.empty());
  EXPECT_EQ(matrix.size(), loader.cell()->elements().size());
}

TEST_F(BondCutoffControllerTests, SetUniformCutoffAndGlobalCutoff) {
  auto &win = window();
  TrajectoryLoader loader;
  ProgramOptions options;
  const std::string filepath = getTestDataDir() + "xyz/clean.xyz";
  const auto res = loader.loadFile(filepath);
  ASSERT_TRUE(res.has_value());
  ASSERT_NE(loader.cell(), nullptr);

  BondCutoffController controller(win, loader, options);
  controller.setUniformCutoff(2.80F);

  auto cutoffs_model = win.get_bond_cutoffs();
  ASSERT_NE(cutoffs_model, nullptr);
  ASSERT_GT(cutoffs_model->row_count(), 0);

  const auto first_row = requireValue(cutoffs_model->row_data(0));
  EXPECT_EQ(first_row.min_distance, "0.00");
  EXPECT_EQ(first_row.max_distance, "2.80");

  controller.applyGlobalCutoff(3.50F);
  cutoffs_model = win.get_bond_cutoffs();
  const auto updated_row = requireValue(cutoffs_model->row_data(0));
  EXPECT_EQ(updated_row.min_distance, "0.00");
  EXPECT_EQ(updated_row.max_distance, "3.50");
}

TEST_F(BondCutoffControllerTests, ApplyCovalentFactorMinAndMax) {
  auto &win = window();
  TrajectoryLoader loader;
  ProgramOptions options;
  const std::string filepath = getTestDataDir() + "xyz/clean.xyz";
  const auto res = loader.loadFile(filepath);
  ASSERT_TRUE(res.has_value());
  ASSERT_NE(loader.cell(), nullptr);

  BondCutoffController controller(win, loader, options);
  controller.setBondCutoffs();

  // Apply Min Factor
  controller.applyCovalentFactor(0.85F, FactorBound::Min);
  auto cutoffs_model = win.get_bond_cutoffs();
  ASSERT_NE(cutoffs_model, nullptr);
  auto row = requireValue(cutoffs_model->row_data(0));
  EXPECT_FALSE(row.min_distance.empty());

  // Apply Max Factor
  controller.applyCovalentFactor(1.25F, FactorBound::Max);
  cutoffs_model = win.get_bond_cutoffs();
  row = requireValue(cutoffs_model->row_data(0));
  EXPECT_FALSE(row.max_distance.empty());
}

TEST_F(BondCutoffControllerTests, InitializesCalibrationPairsOnSetBondCutoffs) {
  auto &win = window();
  TrajectoryLoader loader;
  ProgramOptions options;
  const std::string filepath = getTestDataDir() + "xyz/clean.xyz";
  const auto res = loader.loadFile(filepath);
  ASSERT_TRUE(res.has_value());

  BondCutoffController controller(win, loader, options);
  controller.setBondCutoffs();

  const auto calib_pairs = win.get_calibration_pairs();
  ASSERT_NE(calib_pairs, nullptr);
  EXPECT_GT(calib_pairs->row_count(), 0);

  const auto pair0 = requireValue(calib_pairs->row_data(0));
  EXPECT_EQ(pair0.pair_index, 0);
  EXPECT_FALSE(pair0.label.empty());
  EXPECT_FALSE(pair0.r_min.empty());
  EXPECT_FALSE(pair0.r_max.empty());
  EXPECT_TRUE(pair0.is_valid);
}

TEST_F(BondCutoffControllerTests, AutoDetectCutoffsUpdatesBoundsAndPeak) {
  auto &win = window();
  TrajectoryLoader loader;
  ProgramOptions options;
  const std::string filepath = getTestDataDir() + "xyz/clean.xyz";
  const auto res = loader.loadFile(filepath);
  ASSERT_TRUE(res.has_value());

  AnalysisDispatcher dispatcher;
  auto distribution_funcs =
      std::make_unique<correlation::analysis::DistributionFunctions>(*loader.cell());
  correlation::analysis::Histogram hist;
  hist.title = "Radial Distribution Function";
  hist.bins = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
  hist.partials["Total"] = {0.0, 0.2, 1.0, 3.5, 1.8, 0.4, 0.6, 1.2, 1.1};
  hist.partials["Si-Si"] = {0.0, 0.2, 1.0, 3.5, 1.8, 0.4, 0.6, 1.2, 1.1};
  const std::string title = hist.title;
  distribution_funcs->addHistogram(title, std::move(hist));
  dispatcher.setDistributionFunctions(std::move(distribution_funcs));

  BondCutoffController controller(win, loader, options, &dispatcher);
  controller.setBondCutoffs();

  controller.autoDetectCutoffs(0);

  const auto calib_pairs = win.get_calibration_pairs();
  ASSERT_NE(calib_pairs, nullptr);
  const auto updated_pair = requireValue(calib_pairs->row_data(0));
  EXPECT_NE(updated_pair.r_peak, "—");
  EXPECT_NE(updated_pair.r_min, updated_pair.r_max);

  // Check sync with bond_cutoffs table
  const auto cutoffs_model = win.get_bond_cutoffs();
  ASSERT_NE(cutoffs_model, nullptr);
  const auto cutoff_row = requireValue(cutoffs_model->row_data(0));
  EXPECT_EQ(cutoff_row.min_distance, updated_pair.r_min);
  EXPECT_EQ(cutoff_row.max_distance, updated_pair.r_max);
}

TEST_F(BondCutoffControllerTests, ResetCalibrationPairRestoresDefaults) {
  auto &win = window();
  TrajectoryLoader loader;
  ProgramOptions options;
  const std::string filepath = getTestDataDir() + "xyz/clean.xyz";
  const auto res = loader.loadFile(filepath);
  ASSERT_TRUE(res.has_value());

  BondCutoffController controller(win, loader, options);
  controller.setBondCutoffs();

  controller.setCutoffRange(0, 1.80F, 4.20F);
  auto calib_pairs = win.get_calibration_pairs();
  const auto modified = requireValue(calib_pairs->row_data(0));
  EXPECT_EQ(modified.r_min, "1.80");
  EXPECT_EQ(modified.r_max, "4.20");

  controller.resetCalibrationPair(0);
  calib_pairs = win.get_calibration_pairs();
  const auto restored = requireValue(calib_pairs->row_data(0));
  EXPECT_NE(restored.r_min, "1.80");
  EXPECT_NE(restored.r_max, "4.20");
  EXPECT_EQ(restored.r_peak, "—");
}

} // namespace
} // namespace correlation::app
