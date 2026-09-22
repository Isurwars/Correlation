// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "app/PlotExportService.hpp"
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>

namespace {

using correlation::app::PlotExportService;

TEST(PlotExportServiceTests, GetComparisonKeyDefaultsToTotal) {
  EXPECT_EQ(PlotExportService::getComparisonKey(nullptr), "Total");

  correlation::analysis::Histogram hist;
  hist.partials["Total"] = {};
  EXPECT_EQ(PlotExportService::getComparisonKey(&hist), "Total");
}

TEST(PlotExportServiceTests, GetComparisonKeySelectsFirstPartialWhenTotalMissing) {
  correlation::analysis::Histogram hist;
  hist.partials["Si-O"] = {};
  hist.partials["Si-Si"] = {};
  EXPECT_EQ(PlotExportService::getComparisonKey(&hist), "Si-O");
}

TEST(PlotExportServiceTests, ExportHistogramToSvgFile) {
  correlation::analysis::Histogram hist;
  hist.title = "Test Export";
  hist.bins = {1.0, 2.0, 3.0};
  hist.partials["Total"] = {0.5, 1.5, 0.2};

  correlation::plotters::PlotConfig config;

  const auto temp_path = std::filesystem::temp_directory_path() / "test_export_hist.svg";
  std::error_code ec;
  std::filesystem::remove(temp_path, ec);

  const auto result = PlotExportService::exportHistogram(temp_path.string(), hist, config);
  ASSERT_TRUE(result.has_value());
  EXPECT_TRUE(std::filesystem::exists(temp_path));

  std::ifstream in(temp_path);
  std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  EXPECT_NE(content.find("<svg"), std::string::npos);

  std::filesystem::remove(temp_path, ec);
}

TEST(PlotExportServiceTests, ExportComparisonToSvgFile) {
  correlation::analysis::Histogram hist1;
  hist1.title = "Run 1";
  hist1.bins = {1.0, 2.0};
  hist1.partials["Total"] = {0.5, 1.0};

  correlation::analysis::Histogram hist2;
  hist2.title = "Run 2";
  hist2.bins = {1.0, 2.0};
  hist2.partials["Total"] = {0.6, 1.1};

  std::vector<correlation::plotters::LabeledHistogram> datasets = {
      {"Run 1", &hist1},
      {"Run 2", &hist2},
  };

  correlation::plotters::PlotConfig config;

  const auto temp_path = std::filesystem::temp_directory_path() / "test_export_comp.svg";
  std::error_code ec;
  std::filesystem::remove(temp_path, ec);

  const auto result = PlotExportService::exportComparison(temp_path.string(), datasets, "Total", config);
  ASSERT_TRUE(result.has_value());
  EXPECT_TRUE(std::filesystem::exists(temp_path));

  std::ifstream in(temp_path);
  std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  EXPECT_NE(content.find("<svg"), std::string::npos);

  std::filesystem::remove(temp_path, ec);
}

TEST(PlotExportServiceTests, ExportFailsGracefullyOnInvalidDirectory) {
  correlation::analysis::Histogram hist;
  correlation::plotters::PlotConfig config;

  const std::string invalid_path = "/non_existent_dir_12345/sub/test.svg";
  const auto result = PlotExportService::exportHistogram(invalid_path, hist, config);
  EXPECT_FALSE(result.has_value());
}

} // namespace
