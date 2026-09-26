/**
 * @file PlotTableFormatterTests.cpp
 * @brief Unit tests for PlotTableFormatter service.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/formatters/PlotTableFormatter.hpp"
#include <gtest/gtest.h>
#include <optional>
#include <stdexcept>

namespace correlation::app {
namespace {

template <typename T> [[nodiscard]] T requireValue(const std::optional<T> &opt) {
  if (!opt.has_value()) {
    throw std::runtime_error("Optional value is unexpectedly empty");
  }
  return *opt;
}

TEST(PlotTableFormatterTests, NullHistogramReturnsEmptyModels) {
  const auto model = PlotTableFormatter::formatTable(nullptr);
  ASSERT_NE(model.headers, nullptr);
  ASSERT_NE(model.rows, nullptr);
  EXPECT_EQ(model.headers->row_count(), 0);
  EXPECT_EQ(model.rows->row_count(), 0);

  const auto keys = PlotTableFormatter::extractSortedPartialKeys(nullptr);
  EXPECT_TRUE(keys.empty());
}

TEST(PlotTableFormatterTests, EmptyBinsReturnsHeaderOnly) {
  correlation::analysis::Histogram hist;
  hist.x_label = "r";
  hist.x_unit = "Å";

  const auto model = PlotTableFormatter::formatTable(&hist);
  ASSERT_NE(model.headers, nullptr);
  ASSERT_NE(model.rows, nullptr);
  EXPECT_EQ(model.headers->row_count(), 1);
  EXPECT_EQ(requireValue(model.headers->row_data(0)), slint::SharedString("r (Å)"));
  EXPECT_EQ(model.rows->row_count(), 0);
}

TEST(PlotTableFormatterTests, PartialsOrderingPlacesTotalFirst) {
  correlation::analysis::Histogram hist;
  hist.x_label = "q";
  hist.x_unit = "1/Å";
  hist.partials["Si-O"] = {1.0, 2.0};
  hist.partials["O-O"] = {3.0, 4.0};
  hist.partials["Total"] = {5.0, 6.0};
  hist.partials["Si-Si"] = {7.0, 8.0};

  const auto keys = PlotTableFormatter::extractSortedPartialKeys(&hist);
  ASSERT_EQ(keys.size(), 4);
  EXPECT_EQ(keys[0], "Total");
  EXPECT_EQ(keys[1], "O-O");
  EXPECT_EQ(keys[2], "Si-O");
  EXPECT_EQ(keys[3], "Si-Si");

  const auto model = PlotTableFormatter::formatTable(&hist);
  ASSERT_NE(model.headers, nullptr);
  EXPECT_EQ(model.headers->row_count(), 5);
  EXPECT_EQ(requireValue(model.headers->row_data(0)), slint::SharedString("q (1/Å)"));
  EXPECT_EQ(requireValue(model.headers->row_data(1)), slint::SharedString("Total"));
  EXPECT_EQ(requireValue(model.headers->row_data(2)), slint::SharedString("O-O"));
  EXPECT_EQ(requireValue(model.headers->row_data(3)), slint::SharedString("Si-O"));
  EXPECT_EQ(requireValue(model.headers->row_data(4)), slint::SharedString("Si-Si"));
}

TEST(PlotTableFormatterTests, FormatsRowsAndNumericPrecision) {
  correlation::analysis::Histogram hist;
  hist.x_label = "Distance";
  hist.bins = {1.25, 2.50};
  hist.partials["Total"] = {0.1234567, 10.5};

  const auto model = PlotTableFormatter::formatTable(&hist);
  ASSERT_NE(model.rows, nullptr);
  ASSERT_EQ(model.rows->row_count(), 2);

  const auto row0 = requireValue(model.rows->row_data(0));
  ASSERT_NE(row0.values, nullptr);
  ASSERT_EQ(row0.values->row_count(), 2);
  EXPECT_EQ(requireValue(row0.values->row_data(0)), slint::SharedString("1.2500"));
  EXPECT_EQ(requireValue(row0.values->row_data(1)), slint::SharedString("0.123457"));

  const auto row1 = requireValue(model.rows->row_data(1));
  ASSERT_NE(row1.values, nullptr);
  ASSERT_EQ(row1.values->row_count(), 2);
  EXPECT_EQ(requireValue(row1.values->row_data(0)), slint::SharedString("2.5000"));
  EXPECT_EQ(requireValue(row1.values->row_data(1)), slint::SharedString("10.5"));
}

TEST(PlotTableFormatterTests, SmoothedPartialsTakesPrecedence) {
  correlation::analysis::Histogram hist;
  hist.x_label = "r";
  hist.bins = {1.0};
  hist.partials["Total"] = {10.0};
  hist.smoothed_partials["Total"] = {20.0};

  const auto model = PlotTableFormatter::formatTable(&hist);
  ASSERT_NE(model.rows, nullptr);
  ASSERT_EQ(model.rows->row_count(), 1);

  const auto row = requireValue(model.rows->row_data(0));
  ASSERT_NE(row.values, nullptr);
  EXPECT_EQ(requireValue(row.values->row_data(1)), slint::SharedString("20"));
}

} // namespace
} // namespace correlation::app
