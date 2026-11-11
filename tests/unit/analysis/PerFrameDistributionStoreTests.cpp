/**
 * @file PerFrameDistributionStoreTests.cpp
 * @brief Unit tests for PerFrameDistributionStore curve cache.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "analysis/PerFrameDistributionStore.hpp"

#include <gtest/gtest.h>

namespace {

using namespace correlation::analysis;
using correlation::real_t;

class PerFrameDistributionStoreTests : public ::testing::Test {
protected:
  void SetUp() override { store_ = std::make_unique<PerFrameDistributionStore>(10); }

  [[nodiscard]] auto &store() noexcept { return *store_; }
  [[nodiscard]] const auto &store() const noexcept { return *store_; }

private:
  std::unique_ptr<PerFrameDistributionStore> store_;
};

TEST_F(PerFrameDistributionStoreTests, InitializesWithCorrectFrameCount) {
  EXPECT_EQ(store().getTotalFrames(), 10U);
  EXPECT_FALSE(store().hasData(DistributionMetric::RDF));
  EXPECT_FALSE(store().hasData(DistributionMetric::PAD));
  EXPECT_FALSE(store().hasData(DistributionMetric::SQ));
}

TEST_F(PerFrameDistributionStoreTests, StoresAndRetrievesBinsAndMean) {
  const std::vector<real_t> bins = {
      static_cast<real_t>(0.1), static_cast<real_t>(0.2), static_cast<real_t>(0.3),
      static_cast<real_t>(0.4), static_cast<real_t>(0.5),
  };
  const std::vector<real_t> mean = {
      static_cast<real_t>(1.0), static_cast<real_t>(2.5), static_cast<real_t>(3.2),
      static_cast<real_t>(1.8), static_cast<real_t>(0.4),
  };

  store().setBins(DistributionMetric::RDF, bins);
  store().setMeanCurve(DistributionMetric::RDF, mean);

  EXPECT_TRUE(store().hasData(DistributionMetric::RDF));
  EXPECT_EQ(store().getBins(DistributionMetric::RDF), bins);
  EXPECT_EQ(store().getMeanCurve(DistributionMetric::RDF), mean);
}

TEST_F(PerFrameDistributionStoreTests, StoresAndRetrievesPerFrameCurves) {
  const std::vector<real_t> bins = {
      static_cast<real_t>(1.0),
      static_cast<real_t>(2.0),
      static_cast<real_t>(3.0),
  };
  const std::vector<real_t> frame0 = {
      static_cast<real_t>(0.5),
      static_cast<real_t>(1.2),
      static_cast<real_t>(0.8),
  };
  const std::vector<real_t> frame5 = {
      static_cast<real_t>(0.9),
      static_cast<real_t>(2.1),
      static_cast<real_t>(1.1),
  };

  store().setBins(DistributionMetric::PAD, bins);
  store().setFrameCurve(DistributionMetric::PAD, 0, frame0);
  store().setFrameCurve(DistributionMetric::PAD, 5, frame5);

  EXPECT_TRUE(store().hasFrame(DistributionMetric::PAD, 0));
  EXPECT_FALSE(store().hasFrame(DistributionMetric::PAD, 1));
  EXPECT_TRUE(store().hasFrame(DistributionMetric::PAD, 5));

  EXPECT_EQ(store().getFrameCurve(DistributionMetric::PAD, 0), frame0);
  EXPECT_EQ(store().getFrameCurve(DistributionMetric::PAD, 5), frame5);
  EXPECT_TRUE(store().getFrameCurve(DistributionMetric::PAD, 1).empty());
}

TEST_F(PerFrameDistributionStoreTests, HandlesOutOfBoundsGracefully) {
  const std::vector<real_t> curve = {static_cast<real_t>(1.0), static_cast<real_t>(2.0)};
  // Should not crash or throw
  store().setFrameCurve(DistributionMetric::SQ, 15, curve);
  EXPECT_FALSE(store().hasFrame(DistributionMetric::SQ, 15));
  EXPECT_TRUE(store().getFrameCurve(DistributionMetric::SQ, 15).empty());
}

TEST_F(PerFrameDistributionStoreTests, ValidatesMetadata) {
  const auto &rdf_meta = store().getMetricData(DistributionMetric::RDF);
  EXPECT_EQ(rdf_meta.name, "g(r)");
  EXPECT_EQ(rdf_meta.x_unit, "Å");

  const auto &pad_meta = store().getMetricData(DistributionMetric::PAD);
  EXPECT_EQ(pad_meta.name, "PAD");
  EXPECT_EQ(pad_meta.x_unit, "deg");

  const auto &sq_meta = store().getMetricData(DistributionMetric::SQ);
  EXPECT_EQ(sq_meta.name, "S(Q)");
  EXPECT_EQ(sq_meta.x_unit, "Å⁻¹");
}

} // namespace
