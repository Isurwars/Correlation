// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "app/AppBackend.hpp"

#include <algorithm>
#include <filesystem>
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string getTestDataDir() {
  std::vector<std::string> const candidates = {
      "../../tests/data/",
      "../tests/data/",
      "tests/data/",
      "data/",
  };
  for (const auto &dir : candidates) {
    if (std::filesystem::exists(dir + "xyz/clean.xyz")) {
      return dir;
    }
  }
  return "../../tests/data/";
}

// Test fixture for the AppBackend class
class AppBackendTests : public ::testing::Test {};

TEST_F(AppBackendTests, DefaultConstructorInitializesCorrectly) {
  // Arrange & Act
  correlation::app::AppBackend const backend;
  correlation::app::ProgramOptions const opts = backend.options();

  // Assert
  EXPECT_DOUBLE_EQ(opts.r_max, correlation::app::AppDefaults::R_MAX);
  EXPECT_DOUBLE_EQ(opts.r_bin_width, correlation::app::AppDefaults::R_BIN_WIDTH);
  EXPECT_DOUBLE_EQ(opts.q_max, correlation::app::AppDefaults::Q_MAX);
  EXPECT_DOUBLE_EQ(opts.q_bin_width, correlation::app::AppDefaults::Q_BIN_WIDTH);
  EXPECT_DOUBLE_EQ(opts.r_int_max, correlation::app::AppDefaults::R_INT_MAX);
  EXPECT_DOUBLE_EQ(opts.angle_bin_width, correlation::app::AppDefaults::ANGLE_BIN_WIDTH);
  EXPECT_DOUBLE_EQ(opts.smoothing_sigma, correlation::app::AppDefaults::SMOOTHING_SIGMA);
  EXPECT_EQ(opts.smoothing_kernel, correlation::app::AppDefaults::SMOOTHING_KERNEL);
  EXPECT_DOUBLE_EQ(opts.time_step, correlation::app::AppDefaults::TIME_STEP);

  EXPECT_EQ(backend.cell(), nullptr);
  EXPECT_EQ(backend.getFrameCount(), 0);
  EXPECT_EQ(backend.getTotalAtomCount(), 0);
}

TEST_F(AppBackendTests, SetOptionsModifiesState) {
  // Arrange
  correlation::app::AppBackend backend;
  correlation::app::ProgramOptions opts;
  opts.r_max = 50.0;
  opts.smoothing_sigma = 0.5;
  opts.min_frame = 10;
  opts.max_frame = 20;

  // Act
  backend.setOptions(opts);

  // Assert
  correlation::app::ProgramOptions const new_opts = backend.options();
  EXPECT_DOUBLE_EQ(new_opts.r_max, 50.0);
  EXPECT_DOUBLE_EQ(new_opts.smoothing_sigma, 0.5);
  EXPECT_EQ(new_opts.min_frame, 10);
  EXPECT_EQ(new_opts.max_frame, 20);
}

TEST_F(AppBackendTests, LoadInvalidFileThrowsException) {
  // Arrange
  correlation::app::AppBackend backend;

  // Act & Assert
  EXPECT_THROW(backend.load_file("nonexistent_file_that_should_not_exist.xyz"), std::runtime_error);
}

TEST_F(AppBackendTests, RecommendedTimeStepWithNoCellReturnsDefault) {
  // Arrange
  correlation::app::AppBackend const backend;

  // Act
  double const time_step = backend.getRecommendedTimeStep();

  // Assert
  EXPECT_DOUBLE_EQ(time_step, correlation::app::AppDefaults::TIME_STEP);
}

TEST_F(AppBackendTests, LoadValidXYZFile) {
  // Arrange
  correlation::app::AppBackend backend;
  std::string const data_dir = getTestDataDir();
  std::string const filepath = data_dir + "xyz/clean.xyz";

  // Act
  std::string const load_status = backend.load_file(filepath);

  // Assert
  EXPECT_NE(load_status.find("File loaded:"), std::string::npos);
  EXPECT_EQ(backend.getFrameCount(), 1);
  EXPECT_EQ(backend.getTotalAtomCount(), 2);

  auto const counts = backend.getAtomCounts();
  EXPECT_EQ(counts.at("C"), 1);
  EXPECT_EQ(counts.at("O"), 1);

  const correlation::core::Cell *cell = backend.cell();
  ASSERT_NE(cell, nullptr);
  EXPECT_EQ(cell->atomCount(), 2);
}

TEST_F(AppBackendTests, LoadValidCarFileAndRunAnalysisAndWriteFiles) {
  // Arrange
  correlation::app::AppBackend backend;
  std::string const data_dir = getTestDataDir();
  std::string const filepath = data_dir + "car/clean.car";

  // Act & Assert 1: Load file
  std::string const load_status = backend.load_file(filepath);
  EXPECT_NE(load_status.find("File loaded:"), std::string::npos);
  EXPECT_EQ(backend.getFrameCount(), 1);
  EXPECT_EQ(backend.getTotalAtomCount(), 2);

  const correlation::core::Cell *cell = backend.cell();
  ASSERT_NE(cell, nullptr);
  EXPECT_EQ(cell->atomCount(), 2);

  // Act & Assert 2: Recommended values
  double const rec_time_step = backend.getRecommendedTimeStep();
  EXPECT_NEAR(rec_time_step, 1.3470, 1e-3);

  auto const cutoffs = backend.getRecommendedBondCutoffs();
  ASSERT_FALSE(cutoffs.empty());
  EXPECT_EQ(cutoffs.size(), cell->elements().size());

  // Act & Assert 3: Setup options and run analysis
  correlation::app::ProgramOptions opts = backend.options();
  opts.r_max = 5.0;
  opts.r_bin_width = 0.1;
  opts.smoothing = true;
  opts.smoothing_sigma = 0.1;
  // Enable RDF calculator
  opts.active_calculators["RDF"] = true;
  backend.setOptions(opts);

  std::string const run_status = backend.run_analysis();
  EXPECT_TRUE(run_status.empty()) << "Analysis failed: " << run_status;

  // Verify histograms
  auto const hist_names = backend.getAvailableHistogramNames();
  EXPECT_FALSE(hist_names.empty());
  EXPECT_NE(std::ranges::find(hist_names, "g_r"), hist_names.end());

  const auto *rdf_hist = backend.getHistogram("g_r");
  ASSERT_NE(rdf_hist, nullptr);
  EXPECT_FALSE(rdf_hist->bins.empty());

  const auto weights = backend.getAshcroftWeights();
  EXPECT_FALSE(weights.empty());

  // Act & Assert 4: Write output files to a temp directory
  auto tmp_dir = std::filesystem::temp_directory_path() / "correlation_backend_test";
  std::filesystem::create_directories(tmp_dir);
  struct TempDirGuard {
    std::filesystem::path path;
    explicit TempDirGuard(std::filesystem::path tmp_path) : path(std::move(tmp_path)) {}
    TempDirGuard(const TempDirGuard &) = delete;
    TempDirGuard &operator=(const TempDirGuard &) = delete;
    TempDirGuard(TempDirGuard &&) = delete;
    TempDirGuard &operator=(TempDirGuard &&) = delete;
    ~TempDirGuard() { std::filesystem::remove_all(path); }
  };
  const TempDirGuard guard{tmp_dir};

  std::string const out_base = (tmp_dir / "result").string();

  opts.output_file_base = out_base;
  opts.use_csv = true;
  opts.use_hdf5 = false;
  opts.use_parquet = false;
  backend.setOptions(opts);

  std::string const write_status = backend.write_files();
  EXPECT_TRUE(write_status.empty()) << "Write files failed: " << write_status;

  // Verify file was written
  EXPECT_TRUE(std::filesystem::exists(out_base + "_g.csv"));
}

} // namespace