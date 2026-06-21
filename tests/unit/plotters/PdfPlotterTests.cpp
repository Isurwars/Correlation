// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DistributionFunctions.hpp"
#include "plotters/PdfPlotter.hpp"

#include <gtest/gtest.h>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>

namespace correlation::testing {

using namespace correlation::plotters;
using namespace correlation::analysis;

class PdfPlotterTests : public ::testing::Test {
protected:
  std::string temp_pdf_path;

  void SetUp() override {
    auto temp_dir = std::filesystem::temp_directory_path();
    temp_pdf_path = (temp_dir / "test_output.pdf").string();
  }

  void TearDown() override {
    if (std::filesystem::exists(temp_pdf_path)) {
      std::filesystem::remove(temp_pdf_path);
    }
  }

  bool fileExistsAndIsNotEmpty(const std::string &path) {
    if (!std::filesystem::exists(path)) {
      return false;
    }
    return std::filesystem::file_size(path) > 0;
  }
};

TEST_F(PdfPlotterTests, RendersEmptyHistogramAsPdfGracefully) {
  Histogram empty_hist;
  empty_hist.title = "Empty Plot";
  empty_hist.x_label = "Distance";
  empty_hist.y_label = "g(r)";

  PlotConfig config;
  renderHistogramAsPdf(empty_hist, temp_pdf_path, config);

  EXPECT_TRUE(fileExistsAndIsNotEmpty(temp_pdf_path));
}

TEST_F(PdfPlotterTests, RendersValidHistogramAsPdfCorrectly) {
  Histogram hist;
  hist.title = "Radial Distribution Function g(r)";
  hist.x_label = "r";
  hist.y_label = "g";
  hist.x_unit = "A";
  hist.y_unit = "A^-1";

  hist.bins = {1.0, 2.0, 3.0, 4.0, 5.0};
  hist.partials["Total"] = {0.1, 0.5, 1.2, 0.8, 0.2};
  hist.partials["Si-O"] = {0.0, 0.3, 0.9, 0.4, 0.1};

  PlotConfig config;
  config.theme = PlotConfig::Theme::Light;
  config.show_grid = true;
  config.fill_area = true;
  config.show_markers = true;

  renderHistogramAsPdf(hist, temp_pdf_path, config);

  EXPECT_TRUE(fileExistsAndIsNotEmpty(temp_pdf_path));
}

TEST_F(PdfPlotterTests, RendersComparisonPdfCorrectly) {
  Histogram hist1;
  hist1.title = "Comparison";
  hist1.bins = {1.0, 2.0, 3.0};
  hist1.partials["Total"] = {0.5, 1.0, 1.5};

  Histogram hist2;
  hist2.title = "Comparison";
  hist2.bins = {1.0, 2.0, 3.0};
  hist2.partials["Total"] = {0.6, 1.1, 1.6};

  std::vector<LabeledHistogram> datasets = {{"Run A", &hist1}, {"Run B", &hist2}};

  PlotConfig config;
  config.theme = PlotConfig::Theme::Dark;
  config.fill_area = true;
  renderComparisonPdf(datasets, {"Total", temp_pdf_path}, config);

  EXPECT_TRUE(fileExistsAndIsNotEmpty(temp_pdf_path));
}

} // namespace correlation::testing
