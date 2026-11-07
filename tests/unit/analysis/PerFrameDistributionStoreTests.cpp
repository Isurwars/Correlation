/**
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
*/

#include "analysis/PerFrameDistributionStore.hpp"

#include <gtest/gtest.h>

using namespace correlation::analysis;

class PerFrameDistributionStoreTests : public ::testing::Test {
protected:
  void SetUp() override {
    store_ = std::make_unique<PerFrameDistributionStore>(10);
  }

  std::unique_ptr<PerFrameDistributionStore> store_;
};

TEST_F(PerFrameDistributionStoreTests, InitializesWithCorrectFrameCount) {
  EXPECT_EQ(store_->getTotalFrames(), 10U);
  EXPECT_FALSE(store_->hasData(DistributionMetric::RDF));
  EXPECT_FALSE(store_->hasData(DistributionMetric::PAD));
  EXPECT_FALSE(store_->hasData(DistributionMetric::SQ));
}

TEST_F(PerFrameDistributionStoreTests, StoresAndRetrievesBinsAndMean) {
  const std::vector<real_t> bins = {0.1, 0.2, 0.3, 0.4, 0.5};
  const std::vector<real_t> mean = {1.0, 2.5, 3.2, 1.8, 0.4};

  store_->setBins(DistributionMetric::RDF, bins);
  store_->setMeanCurve(DistributionMetric::RDF, mean);

  EXPECT_TRUE(store_->hasData(DistributionMetric::RDF));
  EXPECT_EQ(store_->getBins(DistributionMetric::RDF), bins);
  EXPECT_EQ(store_->getMeanCurve(DistributionMetric::RDF), mean);
}

TEST_F(PerFrameDistributionStoreTests, StoresAndRetrievesPerFrameCurves) {
  const std::vector<real_t> bins = {1.0, 2.0, 3.0};
  const std::vector<real_t> frame0 = {0.5, 1.2, 0.8};
  const std::vector<real_t> frame5 = {0.9, 2.1, 1.1};

  store_->setBins(DistributionMetric::PAD, bins);
  store_->setFrameCurve(DistributionMetric::PAD, 0, frame0);
  store_->setFrameCurve(DistributionMetric::PAD, 5, frame5);

  EXPECT_TRUE(store_->hasFrame(DistributionMetric::PAD, 0));
  EXPECT_FALSE(store_->hasFrame(DistributionMetric::PAD, 1));
  EXPECT_TRUE(store_->hasFrame(DistributionMetric::PAD, 5));

  EXPECT_EQ(store_->getFrameCurve(DistributionMetric::PAD, 0), frame0);
  EXPECT_EQ(store_->getFrameCurve(DistributionMetric::PAD, 5), frame5);
  EXPECT_TRUE(store_->getFrameCurve(DistributionMetric::PAD, 1).empty());
}

TEST_F(PerFrameDistributionStoreTests, HandlesOutOfBoundsGracefully) {
  const std::vector<real_t> curve = {1.0, 2.0};
  // Should not crash or throw
  store_->setFrameCurve(DistributionMetric::SQ, 15, curve);
  EXPECT_FALSE(store_->hasFrame(DistributionMetric::SQ, 15));
  EXPECT_TRUE(store_->getFrameCurve(DistributionMetric::SQ, 15).empty());
}

TEST_F(PerFrameDistributionStoreTests, ValidatesMetadata) {
  const auto &rdf_meta = store_->getMetricData(DistributionMetric::RDF);
  EXPECT_EQ(rdf_meta.name, "g(r)");
  EXPECT_EQ(rdf_meta.x_unit, "Å");

  const auto &pad_meta = store_->getMetricData(DistributionMetric::PAD);
  EXPECT_EQ(pad_meta.name, "PAD");
  EXPECT_EQ(pad_meta.x_unit, "deg");

  const auto &sq_meta = store_->getMetricData(DistributionMetric::SQ);
  EXPECT_EQ(sq_meta.name, "S(Q)");
  EXPECT_EQ(sq_meta.x_unit, "Å⁻¹");
}

