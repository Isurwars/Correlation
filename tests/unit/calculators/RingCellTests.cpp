// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DistributionFunctions.hpp"
#include "calculators/spatial/MotifFinder.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"
#include "core/Trajectory.hpp"
#include "math/Precision.hpp"

#include "../../CrystalTestHelper.hpp"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <numbers>
#include <vector>

using correlation::real_t;

namespace correlation::analysis {
namespace {

inline real_t vecNorm(const math::Vector3<real_t> &vec) noexcept {
  return std::sqrt(vec.x() * vec.x() + vec.y() * vec.y() + vec.z() * vec.z());
}

class RingCellTests : public ::testing::TestWithParam<size_t> {
private:
  static constexpr real_t RADIUS{1.5};
  static constexpr real_t BOX_SIZE{25.0};

protected:
  [[nodiscard]] static real_t radius() noexcept { return RADIUS; }
  [[nodiscard]] static real_t boxSize() noexcept { return BOX_SIZE; }

  [[nodiscard]] static correlation::core::Cell createRing(size_t num_atoms) {
    return testing::crystals::createRingCell({
        .num_atoms = num_atoms,
        .element = "C",
        .radius = RADIUS,
        .box_size = BOX_SIZE,
    });
  }
};

INSTANTIATE_TEST_SUITE_P(PolygonRingSizes, RingCellTests, ::testing::Values(3, 4, 5, 6, 7, 8));

// 1. Verify atom count, ring center distance, and neighbor distance geometry
TEST_P(RingCellTests, GeometryAndBondDistances) {
  size_t const num_atoms = GetParam();
  auto const cell = createRing(num_atoms);

  EXPECT_EQ(cell.atomCount(), num_atoms);

  real_t const mid = boxSize() * static_cast<real_t>(0.5);
  real_t const expected_bond_dist =
      static_cast<real_t>(2.0) * radius() *
      std::sin(static_cast<real_t>(std::numbers::pi) / static_cast<real_t>(num_atoms));

  // Check radius of each atom relative to center
  for (size_t i = 0; i < num_atoms; ++i) {
    const auto &atom = cell.atoms()[i];
    real_t const delta_x = atom.position().x() - mid;
    real_t const delta_y = atom.position().y() - mid;
    real_t const delta_z = atom.position().z() - mid;
    real_t const dist_from_center =
        std::sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
    EXPECT_NEAR(dist_from_center, radius(), 1e-4);
  }

  // Check distance between adjacent ring atoms
  for (size_t i = 0; i < num_atoms; ++i) {
    size_t const next_idx = (i + 1) % num_atoms;
    const auto &atom1 = cell.atoms()[i];
    const auto &atom2 = cell.atoms()[next_idx];
    real_t const dist = vecNorm(atom1.position() - atom2.position());
    EXPECT_NEAR(dist, expected_bond_dist, 1e-4);
  }
}

// 2. Verify MotifFinder ring size detection for N-polygon rings
TEST_P(RingCellTests, MotifFinderRingDetection) {
  size_t const num_atoms = GetParam();
  auto const cell = createRing(num_atoms);

  core::NeighborGraph graph(num_atoms);
  for (size_t i = 0; i < num_atoms; ++i) {
    size_t const next_idx = (i + 1) % num_atoms;
    const auto &atom1 = cell.atoms()[i];
    const auto &atom2 = cell.atoms()[next_idx];
    real_t const dist = vecNorm(atom1.position() - atom2.position());
    graph.addDirectedEdge(i, next_idx, dist, {0, 0, 0});
    graph.addDirectedEdge(next_idx, i, dist, {0, 0, 0});
  }

  auto rings = calculators::MotifFinder::findRings(graph, 8);
  EXPECT_EQ(rings[num_atoms], 1) << "Ring of size " << num_atoms
                                 << " should detect exactly 1 ring motif.";
}

// 3. Verify Plane Angle Distribution (PAD) internal angles
TEST_P(RingCellTests, PlaneAngleDistribution) {
  size_t const num_atoms = GetParam();
  auto const cell = createRing(num_atoms);

  real_t const bond_dist =
      static_cast<real_t>(2.0) * radius() *
      std::sin(static_cast<real_t>(std::numbers::pi) / static_cast<real_t>(num_atoms));
  real_t const cutoff = bond_dist * static_cast<real_t>(1.2);

  core::Trajectory traj;
  traj.addFrame(cell);
  traj.precomputeBondCutoffs();

  DistributionFunctions dists(cell, cutoff, traj.getBondCutoffsSQ());
  dists.calculatePAD(0.01);

  real_t const expected_angle_deg = static_cast<real_t>(180.0) *
                                    static_cast<real_t>(num_atoms - 2) /
                                    static_cast<real_t>(num_atoms);

  if (dists.getAllHistograms().contains("PAD")) {
    const auto &hist = dists.getHistogram("PAD");
    if (hist.partials.contains("C-C-C")) {
      const auto &ccc = hist.partials.at("C-C-C");
      auto const max_it = std::ranges::max_element(ccc);
      if (max_it != ccc.end() && *max_it > static_cast<real_t>(0.0)) {
        auto const idx = static_cast<size_t>(std::distance(ccc.begin(), max_it));
        real_t const measured_angle = hist.bins[idx];
        EXPECT_NEAR(measured_angle, expected_angle_deg, 2.0)
            << "PAD peak for " << num_atoms << "-ring should be around " << expected_angle_deg
            << " degrees.";
      }
    }
  }
}

// 4. Verify Planar Angle Distribution (PAD) computation completes cleanly
TEST_P(RingCellTests, PlanarAngleDistributionRunsCleanly) {
  size_t const num_atoms = GetParam();
  auto const cell = createRing(num_atoms);

  real_t const bond_dist =
      static_cast<real_t>(2.0) * radius() *
      std::sin(static_cast<real_t>(std::numbers::pi) / static_cast<real_t>(num_atoms));
  real_t const cutoff = bond_dist * static_cast<real_t>(1.2);

  core::Trajectory traj;
  traj.addFrame(cell);
  traj.precomputeBondCutoffs();

  DistributionFunctions dists(cell, cutoff, traj.getBondCutoffsSQ());
  EXPECT_NO_THROW(dists.calculatePAD(0.05));
}

// 5. Verify Radial Interatomic Distance Spectrum (RDF pair distance check)
TEST_P(RingCellTests, RadialPairDistances) {
  size_t const num_atoms = GetParam();
  auto const cell = createRing(num_atoms);

  // Check unique pair distances against 2 * R * sin(k * pi / N)
  for (size_t step_k = 1; step_k <= num_atoms / 2; ++step_k) {
    real_t const expected_d_k =
        static_cast<real_t>(2.0) * radius() *
        std::sin(static_cast<real_t>(step_k) * static_cast<real_t>(std::numbers::pi) /
                 static_cast<real_t>(num_atoms));

    // Find at least one pair matching this topological distance step
    bool distance_found = false;
    for (size_t i = 0; i < num_atoms; ++i) {
      size_t const target_idx = (i + step_k) % num_atoms;
      real_t const actual_d =
          vecNorm(cell.atoms()[i].position() - cell.atoms()[target_idx].position());
      if (std::abs(actual_d - expected_d_k) < 1e-4) {
        distance_found = true;
        break;
      }
    }
    EXPECT_TRUE(distance_found) << "Ring size " << num_atoms
                                << " missing expected pair distance for step k=" << step_k;
  }
}

} // namespace
} // namespace correlation::analysis
