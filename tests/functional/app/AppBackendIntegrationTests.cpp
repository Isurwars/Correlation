// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "app/AppBackend.hpp"

#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace correlation::testing {
namespace {

[[nodiscard]] std::string findTestDataPath(const std::string &relative_path) {
  const std::vector<std::string> search_roots = {
      "../../tests/data/",
      "../tests/data/",
      "tests/data/",
      "data/",
  };

  for (const auto &root : search_roots) {
    auto full_path = std::filesystem::path(root) / relative_path;
    if (std::filesystem::exists(full_path)) {
      return full_path.string();
    }
  }
  return "../../tests/data/" + relative_path;
}

class TempDirectoryGuard {
public:
  explicit TempDirectoryGuard(std::filesystem::path path) : path_(std::move(path)) {
    std::filesystem::create_directories(path_);
  }

  ~TempDirectoryGuard() {
    std::error_code error_code;
    std::filesystem::remove_all(path_, error_code);
  }

  TempDirectoryGuard(const TempDirectoryGuard &) = delete;
  TempDirectoryGuard &operator=(const TempDirectoryGuard &) = delete;
  TempDirectoryGuard(TempDirectoryGuard &&) = delete;
  TempDirectoryGuard &operator=(TempDirectoryGuard &&) = delete;

  [[nodiscard]] const std::filesystem::path &path() const noexcept { return path_; }

private:
  std::filesystem::path path_;
};

TEST(AppBackendIntegrationTests, EndToEndMultiFrameTrajectoryPipeline) {
  const std::string trajectory_path = findTestDataPath("xdatcar/Si.xdatcar");
  ASSERT_TRUE(std::filesystem::exists(trajectory_path)) << "Trajectory test fixture not found: " << trajectory_path;

  correlation::app::AppBackend backend;
  const std::string load_status = backend.load_file(trajectory_path);
  EXPECT_FALSE(load_status.empty());
  EXPECT_EQ(backend.getFrameCount(), 3);
  ASSERT_NE(backend.cell(), nullptr);
  EXPECT_EQ(backend.cell()->atomCount(), 4);

  const auto temp_dir = std::filesystem::temp_directory_path() / "correlation_e2e_backend_test";
  const TempDirectoryGuard temp_guard{temp_dir};
  const std::string out_base = (temp_guard.path() / "analysis_run").string();

  correlation::app::ProgramOptions options = backend.options();
  options.output_file_base = out_base;
  options.r_max = 5.0;
  options.r_bin_width = 0.25;
  options.q_max = 5.0;
  options.q_bin_width = 0.5;
  options.time_step = 1.0;
  options.use_csv = true;
  options.use_hdf5 = false;
  options.use_parquet = false;

  options.active_calculators["RDF"] = true;
  options.active_calculators["SQ"] = true;
  options.active_calculators["VACF"] = true;
  options.active_calculators["MSD"] = true;
  backend.setOptions(options);

  bool progress_invoked = false;
  backend.setProgressCallback([&progress_invoked](float progress, const std::string &) {
    if (progress > 0.0F) {
      progress_invoked = true;
    }
  });

  const auto run_result = backend.run_analysis();
  ASSERT_TRUE(run_result.has_value()) << "run_analysis failed: " << (run_result ? "" : run_result.error());
  EXPECT_TRUE(progress_invoked);
  EXPECT_FALSE(backend.is_cancelled());

  const auto available_hists = backend.getAvailableHistogramNames();
  EXPECT_FALSE(available_hists.empty());
  EXPECT_NE(std::ranges::find(available_hists, "g_r"), available_hists.end());
  EXPECT_NE(std::ranges::find(available_hists, "VACF"), available_hists.end());
  EXPECT_NE(std::ranges::find(available_hists, "MSD"), available_hists.end());

  const auto write_result = backend.write_files();
  ASSERT_TRUE(write_result.has_value()) << "write_files failed: " << (write_result ? "" : write_result.error());

  EXPECT_TRUE(std::filesystem::exists(out_base + "_g.csv"));
  EXPECT_TRUE(std::filesystem::exists(out_base + "_J.csv"));
  EXPECT_TRUE(std::filesystem::exists(out_base + "_VACF.csv"));
  EXPECT_TRUE(std::filesystem::exists(out_base + "_MSD.csv"));
}

TEST(AppBackendIntegrationTests, AbortsGracefullyOnMissingTrajectory) {
  correlation::app::AppBackend backend;
  const auto run_result = backend.run_analysis();
  ASSERT_FALSE(run_result.has_value());
  EXPECT_EQ(run_result.error(), correlation::app::AppDefaults::MSG_ANALYSIS_ABORTED);
}

} // namespace
} // namespace correlation::testing
