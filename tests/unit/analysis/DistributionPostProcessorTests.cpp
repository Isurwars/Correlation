// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DistributionFunctions.hpp"
#include "analysis/DistributionPostProcessor.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>

namespace correlation::testing {
namespace {

using namespace correlation::analysis;
using correlation::real_t;

TEST(DistributionPostProcessorTests, SmoothSingleHistogramGeneratesSmoothedPartials) {
  Histogram hist;
  hist.bins = {
      static_cast<real_t>(0.0), static_cast<real_t>(1.0), static_cast<real_t>(2.0),
      static_cast<real_t>(3.0), static_cast<real_t>(4.0),
  };
  hist.partials["Total"] = {
      static_cast<real_t>(0.0),  static_cast<real_t>(10.0), static_cast<real_t>(20.0),
      static_cast<real_t>(10.0), static_cast<real_t>(0.0),
  };

  DistributionPostProcessor::smooth(hist, "RDF", static_cast<real_t>(0.5));

  ASSERT_TRUE(hist.smoothed_partials.contains("Total"));
  EXPECT_EQ(hist.smoothed_partials["Total"].size(), hist.partials["Total"].size());
  // Peak should still be around the center
  EXPECT_GT(hist.smoothed_partials["Total"][2], hist.smoothed_partials["Total"][0]);
}

TEST(DistributionPostProcessorTests, AngularDistributionsScaleEffectiveSigma) {
  Histogram hist_spatial;
  hist_spatial.bins = {
      static_cast<real_t>(0.0), static_cast<real_t>(1.0), static_cast<real_t>(2.0),
      static_cast<real_t>(3.0), static_cast<real_t>(4.0),
  };
  hist_spatial.partials["Total"] = {
      static_cast<real_t>(0.0),  static_cast<real_t>(10.0), static_cast<real_t>(20.0),
      static_cast<real_t>(10.0), static_cast<real_t>(0.0),
  };

  Histogram hist_angular = hist_spatial;

  DistributionPostProcessor::smooth(hist_spatial, "RDF", static_cast<real_t>(0.1));
  DistributionPostProcessor::smooth(hist_angular, "PAD", static_cast<real_t>(0.1));

  // Angular smoothing with 20x domain scale flattens the peak more than spatial smoothing
  EXPECT_LT(hist_angular.smoothed_partials["Total"][2], hist_spatial.smoothed_partials["Total"][2]);
}

TEST(DistributionPostProcessorTests, StructureFactorScalesEffectiveSigma) {
  Histogram hist_spatial;
  hist_spatial.bins = {
      static_cast<real_t>(0.00), static_cast<real_t>(0.05), static_cast<real_t>(0.10),
      static_cast<real_t>(0.15), static_cast<real_t>(0.20),
  };
  hist_spatial.partials["Total"] = {
      static_cast<real_t>(0.0),  static_cast<real_t>(10.0), static_cast<real_t>(20.0),
      static_cast<real_t>(10.0), static_cast<real_t>(0.0),
  };

  Histogram hist_sq = hist_spatial;

  DistributionPostProcessor::smooth(hist_spatial, "RDF", static_cast<real_t>(0.1));
  DistributionPostProcessor::smooth(hist_sq, "S_q", static_cast<real_t>(0.1));

  // S(Q) smoothing with 2.0x domain scale flattens the peak more than spatial smoothing
  EXPECT_LT(hist_sq.smoothed_partials["Total"][2], hist_spatial.smoothed_partials["Total"][2]);
}

TEST(DistributionPostProcessorTests, MatchesAliasedAndFormattedNamesForScaling) {
  Histogram base_hist;
  base_hist.bins = {
      static_cast<real_t>(0.00), static_cast<real_t>(0.05), static_cast<real_t>(0.10),
      static_cast<real_t>(0.15), static_cast<real_t>(0.20),
  };
  base_hist.partials["Total"] = {
      static_cast<real_t>(0.0),  static_cast<real_t>(10.0), static_cast<real_t>(20.0),
      static_cast<real_t>(10.0), static_cast<real_t>(0.0),
  };

  Histogram h_pad_raw = base_hist;
  Histogram h_pad_norm = base_hist;
  Histogram h_sq_alias = base_hist;

  DistributionPostProcessor::smooth(h_pad_raw, "PAD(θ) RAW", static_cast<real_t>(0.1));
  DistributionPostProcessor::smooth(h_pad_norm, "PAD(θ) Normalized", static_cast<real_t>(0.1));
  DistributionPostProcessor::smooth(h_sq_alias, "S(Q)", static_cast<real_t>(0.1));

  EXPECT_LT(h_pad_raw.smoothed_partials["Total"][2], base_hist.partials["Total"][2]);
  EXPECT_LT(h_pad_norm.smoothed_partials["Total"][2], base_hist.partials["Total"][2]);
  EXPECT_LT(h_sq_alias.smoothed_partials["Total"][2], base_hist.partials["Total"][2]);
}

TEST(DistributionPostProcessorTests, SmoothAllOperatesOnHistogramMap) {
  std::map<std::string, Histogram> map;

  Histogram hist1;
  hist1.bins = {static_cast<real_t>(0.0), static_cast<real_t>(1.0), static_cast<real_t>(2.0)};
  hist1.partials["P1"] = {
      static_cast<real_t>(1.0),
      static_cast<real_t>(2.0),
      static_cast<real_t>(1.0),
  };
  map["H1"] = hist1;

  Histogram hist2;
  hist2.bins = {static_cast<real_t>(0.0), static_cast<real_t>(1.0), static_cast<real_t>(2.0)};
  hist2.partials["P2"] = {
      static_cast<real_t>(3.0),
      static_cast<real_t>(4.0),
      static_cast<real_t>(3.0),
  };
  map["H2"] = hist2;

  DistributionPostProcessor::smoothAll(map, static_cast<real_t>(0.5));

  EXPECT_TRUE(map["H1"].smoothed_partials.contains("P1"));
  EXPECT_TRUE(map["H2"].smoothed_partials.contains("P2"));
}

TEST(DistributionPostProcessorTests, SmoothOnDistributionFunctionsThrowsOnMissingHistogram) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  DistributionFunctions dists(cell);

  EXPECT_THROW(DistributionPostProcessor::smooth(dists, "NonExistent", static_cast<real_t>(0.5)),
               std::runtime_error);
}

} // namespace
} // namespace correlation::testing
