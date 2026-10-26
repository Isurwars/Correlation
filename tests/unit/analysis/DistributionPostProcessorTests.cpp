// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DistributionFunctions.hpp"
#include "analysis/DistributionPostProcessor.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>

namespace correlation::testing {

using namespace correlation::analysis;

TEST(DistributionPostProcessorTests, SmoothSingleHistogramGeneratesSmoothedPartials) {
  Histogram hist;
  hist.bins = {0.0, 1.0, 2.0, 3.0, 4.0};
  hist.partials["Total"] = {0.0, 10.0, 20.0, 10.0, 0.0};

  DistributionPostProcessor::smooth(hist, "RDF", 0.5);

  ASSERT_TRUE(hist.smoothed_partials.contains("Total"));
  EXPECT_EQ(hist.smoothed_partials["Total"].size(), hist.partials["Total"].size());
  // Peak should still be around the center
  EXPECT_GT(hist.smoothed_partials["Total"][2], hist.smoothed_partials["Total"][0]);
}

TEST(DistributionPostProcessorTests, AngularDistributionsScaleEffectiveSigma) {
  Histogram hist_spatial;
  hist_spatial.bins = {0.0, 1.0, 2.0, 3.0, 4.0};
  hist_spatial.partials["Total"] = {0.0, 10.0, 20.0, 10.0, 0.0};

  Histogram hist_angular = hist_spatial;

  DistributionPostProcessor::smooth(hist_spatial, "RDF", 0.1);
  DistributionPostProcessor::smooth(hist_angular, "PAD", 0.1);

  // Angular smoothing with 15x domain scale flattens the peak more than spatial smoothing
  EXPECT_LT(hist_angular.smoothed_partials["Total"][2], hist_spatial.smoothed_partials["Total"][2]);
}

TEST(DistributionPostProcessorTests, SmoothAllOperatesOnHistogramMap) {
  std::map<std::string, Histogram> map;

  Histogram hist1;
  hist1.bins = {0.0, 1.0, 2.0};
  hist1.partials["P1"] = {1.0, 2.0, 1.0};
  map["H1"] = hist1;

  Histogram hist2;
  hist2.bins = {0.0, 1.0, 2.0};
  hist2.partials["P2"] = {3.0, 4.0, 3.0};
  map["H2"] = hist2;

  DistributionPostProcessor::smoothAll(map, 0.5);

  EXPECT_TRUE(map["H1"].smoothed_partials.contains("P1"));
  EXPECT_TRUE(map["H2"].smoothed_partials.contains("P2"));
}

TEST(DistributionPostProcessorTests, SmoothOnDistributionFunctionsThrowsOnMissingHistogram) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  DistributionFunctions dists(cell);

  EXPECT_THROW(DistributionPostProcessor::smooth(dists, "NonExistent", 0.5), std::runtime_error);
}

} // namespace correlation::testing
