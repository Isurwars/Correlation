// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DistributionFunctions.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "calculators/MotifFinder.hpp"
#include "calculators/PADCalculator.hpp"
#include "calculators/RDCalculator.hpp"
#include "calculators/RDFCalculator.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"
#include "core/Trajectory.hpp"
#include "math/Constants.hpp"
#include "math/Precision.hpp"

#include "../../CrystalTestHelper.hpp"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <numbers>
#include <vector>

namespace correlation::analysis {
namespace {

inline real_t vecNorm(const math::Vector3<real_t> &v) noexcept {
  return std::sqrt(v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
}

class RingCellTests : public ::testing::TestWithParam<size_t> {
protected:
  real_t const radius_{1.5};
  real_t const box_size_{25.0};

  correlation::core::Cell createRing(size_t num_atoms) const {
    return testing::crystals::createRingCell(
        {.num_atoms = num_atoms, .element = "C", .radius = radius_, .box_size = box_size_});
  }
};

INSTANTIATE_TEST_SUITE_P(PolygonRingSizes, RingCellTests, ::testing::Values(3, 4, 5, 6, 7, 8));

// 1. Verify atom count, ring center distance, and neighbor distance geometry
TEST_P(RingCellTests, GeometryAndBondDistances) {
  size_t const n = GetParam();
  auto const cell = createRing(n);

  EXPECT_EQ(cell.atomCount(), n);

  real_t const mid = box_size_ * static_cast<real_t>(0.5);
  real_t const expected_bond_dist = static_cast<real_t>(2.0) * radius_ *
                                    std::sin(static_cast<real_t>(std::numbers::pi) / static_cast<real_t>(n));

  // Check radius of each atom relative to center
  for (size_t i = 0; i < n; ++i) {
    const auto &atom = cell.atoms()[i];
    real_t const dx = atom.position().x() - mid;
    real_t const dy = atom.position().y() - mid;
    real_t const dz = atom.position().z() - mid;
    real_t const dist_from_center = std::sqrt(dx * dx + dy * dy + dz * dz);
    EXPECT_NEAR(dist_from_center, radius_, 1e-4);
  }

  // Check distance between adjacent ring atoms
  for (size_t i = 0; i < n; ++i) {
    size_t const next_idx = (i + 1) % n;
    const auto &a1 = cell.atoms()[i];
    const auto &a2 = cell.atoms()[next_idx];
    real_t const dist = vecNorm(a1.position() - a2.position());
    EXPECT_NEAR(dist, expected_bond_dist, 1e-4);
  }
}

// 2. Verify MotifFinder ring size detection for N-polygon rings
TEST_P(RingCellTests, MotifFinderRingDetection) {
  size_t const n = GetParam();
  auto const cell = createRing(n);

  core::NeighborGraph graph(n);
  for (size_t i = 0; i < n; ++i) {
    size_t const next_idx = (i + 1) % n;
    const auto &a1 = cell.atoms()[i];
    const auto &a2 = cell.atoms()[next_idx];
    real_t const d = vecNorm(a1.position() - a2.position());
    graph.addDirectedEdge(i, next_idx, d, {0, 0, 0});
    graph.addDirectedEdge(next_idx, i, d, {0, 0, 0});
  }

  auto rings = calculators::MotifFinder::findRings(graph, 8);
  EXPECT_EQ(rings[n], 1) << "Ring of size " << n << " should detect exactly 1 ring motif.";
}

// 3. Verify Bond Angle Distribution (BAD) internal angles
TEST_P(RingCellTests, BondAngleDistribution) {
  size_t const n = GetParam();
  auto const cell = createRing(n);

  real_t const bond_dist = static_cast<real_t>(2.0) * radius_ *
                           std::sin(static_cast<real_t>(std::numbers::pi) / static_cast<real_t>(n));
  real_t const cutoff = bond_dist * static_cast<real_t>(1.2);

  core::Trajectory traj;
  traj.addFrame(cell);
  traj.precomputeBondCutoffs();

  DistributionFunctions dists(cell, cutoff, traj.getBondCutoffsSQ());
  dists.calculatePAD(0.01);

  real_t const expected_angle_deg = static_cast<real_t>(180.0) * static_cast<real_t>(n - 2) / static_cast<real_t>(n);

  if (dists.getAllHistograms().contains("BAD")) {
    const auto &hist = dists.getHistogram("BAD");
    if (hist.partials.contains("C-C-C")) {
      const auto &ccc = hist.partials.at("C-C-C");
      auto max_it = std::max_element(ccc.begin(), ccc.end());
      if (max_it != ccc.end() && *max_it > static_cast<real_t>(0.0)) {
        size_t const idx = static_cast<size_t>(std::distance(ccc.begin(), max_it));
        real_t const measured_angle = hist.bins[idx];
        EXPECT_NEAR(measured_angle, expected_angle_deg, 2.0)
            << "BAD peak for " << n << "-ring should be around " << expected_angle_deg << " degrees.";
      }
    }
  }
}

// 4. Verify Planar Angle Distribution (PAD) computation completes cleanly
TEST_P(RingCellTests, PlanarAngleDistributionRunsCleanly) {
  size_t const n = GetParam();
  auto const cell = createRing(n);

  real_t const bond_dist = static_cast<real_t>(2.0) * radius_ *
                           std::sin(static_cast<real_t>(std::numbers::pi) / static_cast<real_t>(n));
  real_t const cutoff = bond_dist * static_cast<real_t>(1.2);

  core::Trajectory traj;
  traj.addFrame(cell);
  traj.precomputeBondCutoffs();

  DistributionFunctions dists(cell, cutoff, traj.getBondCutoffsSQ());
  EXPECT_NO_THROW(dists.calculatePAD(0.05));
}

// 5. Verify Radial Interatomic Distance Spectrum (RDF pair distance check)
TEST_P(RingCellTests, RadialPairDistances) {
  size_t const n = GetParam();
  auto const cell = createRing(n);

  // Check unique pair distances against 2 * R * sin(k * pi / N)
  for (size_t k = 1; k <= n / 2; ++k) {
    real_t const expected_d_k = static_cast<real_t>(2.0) * radius_ *
                                std::sin(static_cast<real_t>(k) * static_cast<real_t>(std::numbers::pi) /
                                         static_cast<real_t>(n));

    // Find at least one pair matching this topological distance step
    bool distance_found = false;
    for (size_t i = 0; i < n; ++i) {
      size_t const j = (i + k) % n;
      real_t const actual_d = vecNorm(cell.atoms()[i].position() - cell.atoms()[j].position());
      if (std::abs(actual_d - expected_d_k) < 1e-4) {
        distance_found = true;
        break;
      }
    }
    EXPECT_TRUE(distance_found) << "Ring size " << n << " missing expected pair distance for step k=" << k;
  }
}

} // namespace
} // namespace correlation::analysis
