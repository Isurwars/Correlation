/**
 * @file PerFrameDistributionStoreTests.cpp
 * @brief Unit tests for PerFrameDistributionStore curve cache.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "analysis/DistributionFunctions.hpp"
#include "analysis/PerFrameDistributionStore.hpp"
#include "analysis/TrajectoryAnalyzer.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

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
  EXPECT_FALSE(store().hasData(DistributionMetric::XRD));
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

  const auto &xrd_meta = store().getMetricData(DistributionMetric::XRD);
  EXPECT_EQ(xrd_meta.name, "XRD");
  EXPECT_EQ(xrd_meta.x_unit, "deg");
  EXPECT_EQ(xrd_meta.y_unit, "a.u.");
}

TEST_F(PerFrameDistributionStoreTests, StoresAndRetrievesXRDPerFrameAndMean) {
  const std::vector<real_t> bins = {static_cast<real_t>(20.0), static_cast<real_t>(40.0),
                                    static_cast<real_t>(60.0)};
  const std::vector<real_t> mean = {static_cast<real_t>(100.0), static_cast<real_t>(500.0),
                                    static_cast<real_t>(250.0)};
  const std::vector<real_t> frame2 = {static_cast<real_t>(90.0), static_cast<real_t>(520.0),
                                      static_cast<real_t>(240.0)};

  store().setBins(DistributionMetric::XRD, bins);
  store().setMeanCurve(DistributionMetric::XRD, mean);
  store().setFrameCurve(DistributionMetric::XRD, 2, frame2);

  EXPECT_TRUE(store().hasData(DistributionMetric::XRD));
  EXPECT_TRUE(store().hasFrame(DistributionMetric::XRD, 2));
  EXPECT_FALSE(store().hasFrame(DistributionMetric::XRD, 1));
  EXPECT_EQ(store().getBins(DistributionMetric::XRD), bins);
  EXPECT_EQ(store().getMeanCurve(DistributionMetric::XRD), mean);
  EXPECT_EQ(store().getFrameCurve(DistributionMetric::XRD, 2), frame2);
}

correlation::core::Trajectory makeSyntheticTrajectory(size_t frame_count) {
  correlation::core::Trajectory traj;
  for (size_t frame_idx = 0; frame_idx < frame_count; ++frame_idx) {
    correlation::core::Cell cell({static_cast<real_t>(10.0), static_cast<real_t>(10.0),
                                  static_cast<real_t>(10.0), static_cast<real_t>(90.0),
                                  static_cast<real_t>(90.0), static_cast<real_t>(90.0)});
    const auto shift = static_cast<real_t>(frame_idx) * static_cast<real_t>(0.1);
    cell.addAtom("Ar", {static_cast<real_t>(1.0) + shift, static_cast<real_t>(1.0),
                        static_cast<real_t>(1.0)});
    cell.addAtom("Ar", {static_cast<real_t>(4.0), static_cast<real_t>(1.0),
                        static_cast<real_t>(1.0)});
    cell.addAtom("Ar", {static_cast<real_t>(1.0), static_cast<real_t>(4.5),
                        static_cast<real_t>(1.0)});
    cell.addAtom("Ar", {static_cast<real_t>(4.5), static_cast<real_t>(4.5),
                        static_cast<real_t>(4.5)});
    traj.addFrame(cell);
  }
  return traj;
}

size_t countPopulatedFrames(const PerFrameDistributionStore &store, DistributionMetric metric) {
  size_t populated_frames = 0;
  for (size_t frame_idx = 0; frame_idx < store.getTotalFrames(); ++frame_idx) {
    const bool has_curve =
        store.hasFrame(metric, frame_idx) && !store.getFrameCurve(metric, frame_idx).empty();
    populated_frames += static_cast<size_t>(has_curve);
  }
  return populated_frames;
}

TEST_F(PerFrameDistributionStoreTests, ComputeMeanPopulatesAllFramesWhenStrided) {
  auto traj = makeSyntheticTrajectory(5);

  AnalysisSettings settings;
  settings.r_max = static_cast<real_t>(4.5);
  settings.r_bin_width = static_cast<real_t>(0.5);
  settings.frame_stride = 2;
  settings.smoothing = false;
  settings.active_calculators = {
      {"Radial Distribution Function", true},
      {"Plane Angle Distribution", false},
      {"Static Structure Factor", false},
      {"XRD", false},
  };

  const size_t start_frame = 1; // Analyzed frames: 1, 3; skipped frames: 0, 2, 4
  const TrajectoryAnalyzer analyzer(traj, settings.r_max, BondCutoffMatrix{},
                                    StartFrame{start_frame}, EndFrame{5});
  const auto dists = DistributionFunctions::computeMean(traj, analyzer, start_frame, settings);

  ASSERT_NE(dists, nullptr);
  const auto per_frame = dists->getPerFrameStore();
  ASSERT_NE(per_frame, nullptr);
  EXPECT_EQ(per_frame->getTotalFrames(), 5U);
  EXPECT_TRUE(per_frame->hasData(DistributionMetric::RDF));
  EXPECT_FALSE(per_frame->hasData(DistributionMetric::SQ));
  EXPECT_EQ(countPopulatedFrames(*per_frame, DistributionMetric::RDF), 5U);
}

} // namespace
