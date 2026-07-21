// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Cell.hpp"
#include "math/Constants.hpp"
#include "math/LinearAlgebra.hpp"
#include "math/Precision.hpp"

#include <cmath>
#include <gtest/gtest.h>

namespace correlation::testing {

using namespace correlation::core;
namespace {
class CellTests : public ::testing::Test {
protected:
  Cell orthogonal_cell_{{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}};
};
} // namespace

// --- Unitary Tests: Constructors & Lattice ---

TEST_F(CellTests, DefaultConstructorInitializesEmpty) {
  const Cell cell{};
  EXPECT_TRUE(cell.isEmpty());
  EXPECT_EQ(cell.atomCount(), 0);
  EXPECT_EQ(cell.volume(), 0.0);
}

TEST_F(CellTests, ParameterConstructorSetsCorrectVolume) {
  const std::array<real_t, 6> params = {4.0, 4.0, 4.0, 90.0, 90.0, 90.0};
  const Cell cell(params);
  EXPECT_NEAR(cell.volume(), 64.0, 1e-9);
}

TEST_F(CellTests, VectorConstructorCalculatesParameters) {
  const Cell cell({2.0, 0.0, 0.0}, {0.0, 3.0, 0.0}, {0.0, 0.0, 4.0});
  const auto &params = cell.lattice_parameters();

  EXPECT_NEAR(params[0], 2.0, 1e-9);
  EXPECT_NEAR(params[1], 3.0, 1e-9);
  EXPECT_NEAR(params[2], 4.0, 1e-9);
  EXPECT_NEAR(params[3], 90.0, 1e-9);
  EXPECT_NEAR(params[4], 90.0, 1e-9);
  EXPECT_NEAR(params[5], 90.0, 1e-9);
  EXPECT_NEAR(cell.volume(), 24.0, 1e-9);
}

TEST_F(CellTests, NonOrthogonalVolumeIsCorrect) {
  const std::array<real_t, 6> params = {5.0, 6.0, 7.0, 80.0, 90.0, 100.0};
  const Cell cell(params);

  // Expected volume calculation
  const real_t cos_a = static_cast<real_t>(std::cos(80.0 * correlation::math::deg_to_rad));
  const real_t cos_b = static_cast<real_t>(std::cos(90.0 * correlation::math::deg_to_rad));
  const real_t cos_g = static_cast<real_t>(std::cos(100.0 * correlation::math::deg_to_rad));
  const real_t vol_sqrt = static_cast<real_t>(1.0) - (cos_a * cos_a) - (cos_b * cos_b) - (cos_g * cos_g) +
                          static_cast<real_t>(2.0) * (cos_a * cos_b * cos_g);
  const real_t expected_volume =
      static_cast<real_t>(5.0) * static_cast<real_t>(6.0) * static_cast<real_t>(7.0) * std::sqrt(vol_sqrt);

  EXPECT_NEAR(cell.volume(), expected_volume, correlation::is_single_precision ? 1e-5 : 1e-9);
}

TEST_F(CellTests, RuleOfFiveWorksCorrectly) {
  const std::array<real_t, 6> params = {4.0, 4.0, 4.0, 90.0, 90.0, 90.0};
  Cell cell(params);
  cell.addAtom("H", {0.5, 0.5, 0.5});
  cell.setEnergy(1.23);

  // Copy Constructor
  Cell const cell_copy(cell);
  EXPECT_EQ(cell_copy.volume(), cell.volume());
  EXPECT_EQ(cell_copy.atomCount(), cell.atomCount());
  EXPECT_NEAR(cell_copy.getEnergy(), 1.23, correlation::is_single_precision ? 1e-6 : 1e-9);

  // Copy Assignment
  Cell cell_assigned{};
  cell_assigned = cell;
  EXPECT_EQ(cell_assigned.volume(), cell.volume());
  EXPECT_EQ(cell_assigned.atomCount(), cell.atomCount());
}

TEST_F(CellTests, MoveSemanticsLeavesMovedFromStateEmpty) {
  Cell cell({{4.0, 4.0, 4.0, 90.0, 90.0, 90.0}});
  cell.addAtom("H", {0.5, 0.5, 0.5});

  Cell const cell_moved(std::move(cell));
  EXPECT_TRUE(cell.isEmpty());
  EXPECT_EQ(cell_moved.atomCount(), 1);
}

TEST_F(CellTests, AccessorsWorkCorrectly) {
  Cell cell{};
  cell.setLatticeParameters({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  EXPECT_NEAR(cell.volume(), 1000.0, correlation::is_single_precision ? 1e-5 : 1e-9);

  cell.setEnergy(-13.6);
  EXPECT_NEAR(cell.getEnergy(), -13.6, correlation::is_single_precision ? 1e-6 : 1e-9);
}

// --- Unitary Tests: Limit Cases (Lattice) ---

TEST_F(CellTests, ThrowsOnInvalidLatticeParameters) {
  EXPECT_THROW(Cell({-5.0, 1.0, 1.0, 90.0, 90.0, 90.0}), std::invalid_argument);
  EXPECT_THROW(Cell({0.0, 1.0, 1.0, 90.0, 90.0, 90.0}), std::invalid_argument);
  EXPECT_THROW(Cell({1.0, 1.0, 1.0, 0.0, 90.0, 90.0}), std::invalid_argument);
  EXPECT_THROW(Cell({1.0, 1.0, 1.0, 180.0, 90.0, 90.0}), std::invalid_argument);
}

// --- Unitary Tests: PBC Methods ---

TEST_F(CellTests, WrapPositionsCorrectlyMapsAtoms) {
  Cell cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  cell.addAtom("H", {12.0, -3.0, 25.0});
  cell.wrapPositions();
  const auto &pos = cell.atoms().front().position();
  EXPECT_NEAR(pos.x(), 2.0, correlation::is_single_precision ? 1e-5 : 1e-9);
  EXPECT_NEAR(pos.y(), 7.0, correlation::is_single_precision ? 1e-5 : 1e-9);
  EXPECT_NEAR(pos.z(), 5.0, correlation::is_single_precision ? 1e-5 : 1e-9);
}

TEST_F(CellTests, WrapPositionsHandlesLargeDisplacements) {
  Cell cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  cell.addAtom("H", {1002.5, -997.5, 0.0});
  cell.wrapPositions();
  const auto &pos = cell.atoms().front().position();
  EXPECT_NEAR(pos.x(), 2.5, correlation::is_single_precision ? 1e-3 : 1e-9);
  EXPECT_NEAR(pos.y(), 2.5, correlation::is_single_precision ? 1e-3 : 1e-9);
}

TEST_F(CellTests, MinimumImageCalculatesShortestVector) {
  Cell const cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  auto image = cell.minimumImage({7.0, -8.0, 2.0});
  EXPECT_NEAR(image.x(), -3.0, correlation::is_single_precision ? 1e-5 : 1e-9);
  EXPECT_NEAR(image.y(), 2.0, correlation::is_single_precision ? 1e-5 : 1e-9);
  EXPECT_NEAR(image.z(), 2.0, correlation::is_single_precision ? 1e-5 : 1e-9);
}

TEST_F(CellTests, MinimumImageHandlesDistancesLargerThanBox) {
  Cell const cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  auto image = cell.minimumImage({15.0, 0.0, 0.0});
  EXPECT_NEAR(std::abs(image.x()), 5.0, correlation::is_single_precision ? 1e-5 : 1e-9);
}

TEST_F(CellTests, InverseLatticeVectorsAreCorrect) {
  const Cell cell({{2.0, 4.0, 5.0, 90.0, 90.0, 90.0}});
  const auto &inv = cell.inverseLatticeVectors();
  EXPECT_NEAR(inv[0][0], 0.5, correlation::is_single_precision ? 1e-5 : 1e-9);
  EXPECT_NEAR(inv[1][1], 0.25, correlation::is_single_precision ? 1e-5 : 1e-9);
  EXPECT_NEAR(inv[2][2], 0.2, correlation::is_single_precision ? 1e-5 : 1e-9);
}

// --- Unitary Tests: Element Management ---

TEST_F(CellTests, AddAtomRegistersNewElements) {
  Cell cell{};
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("O", {1.0, 1.0, 1.0});
  cell.addAtom("Si", {2.0, 2.0, 2.0});
  EXPECT_EQ(cell.elements().size(), 2);
  EXPECT_EQ(cell.atomCount(), 3);
}

TEST_F(CellTests, FindElementWorksCorrectly) {
  Cell cell{};
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  auto si_opt = cell.findElement("Si");
  ASSERT_TRUE(si_opt.has_value());
  EXPECT_EQ(si_opt->symbol, "Si");
  EXPECT_FALSE(cell.findElement("Au").has_value());
}

TEST_F(CellTests, ConstructorThrowsOnZeroOrSingularVolume) {
  // Linearly dependent lattice vectors: zero volume
  EXPECT_THROW(Cell({1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {0.0, 0.0, 1.0}), std::logic_error);
  EXPECT_THROW(Cell({0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}), std::logic_error);

  // Degenerate parameters yielding volume near 0
  EXPECT_THROW(Cell({10.0, 10.0, 10.0, 1.0, 1.0, 179.0}), std::logic_error);
}

TEST_F(CellTests, SetLatticeParametersBoundaryAngles) {
  // Extreme but valid angles
  EXPECT_NO_THROW(Cell({10.0, 10.0, 10.0, 0.1, 90.0, 90.0}));
  EXPECT_NO_THROW(Cell({10.0, 10.0, 10.0, 179.9, 90.0, 90.0}));

  // Invalid angles throwing invalid_argument
  EXPECT_THROW(Cell({10.0, 10.0, 10.0, -0.1, 90.0, 90.0}), std::invalid_argument);
  EXPECT_THROW(Cell({10.0, 10.0, 10.0, 180.0, 90.0, 90.0}), std::invalid_argument);
}

TEST_F(CellTests, MinimumImageHandlesInfiniteAndNaNDistance) {
  Cell const cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});

  // Checking that passing Inf or NaN distances doesn't crash but propagates or returns predictably
  auto mi_nan = cell.minimumImage({std::numeric_limits<real_t>::quiet_NaN(), 0.0, 0.0});
  EXPECT_TRUE(std::isnan(mi_nan.x()));

  auto mi_inf = cell.minimumImage({std::numeric_limits<real_t>::infinity(), 0.0, 0.0});
  EXPECT_TRUE(std::isnan(mi_inf.x())); // std::round(inf) yields nan/indefinite, which makes x() NaN
}

// --- Extreme / Edge-Case Tests ---

TEST_F(CellTests, TriclinicCellWrapPositions) {
  // A triclinic cell: a=10, b=10, c=10, alpha=80, beta=85, gamma=75
  Cell cell({{10.0, 10.0, 10.0, 80.0, 85.0, 75.0}});
  // Place an atom well outside the cell in Cartesian coordinates
  cell.addAtom("H", {25.0, -15.0, 30.0});
  cell.wrapPositions();

  const auto &pos = cell.atoms().front().position();
  // After wrapping, fractional coordinates must be in [0, 1).
  // Convert back to fractional and verify
  auto frac = cell.inverseLatticeVectors() * pos;
  EXPECT_GE(frac.x(), 0.0 - 1e-9);
  EXPECT_LT(frac.x(), 1.0 + 1e-9);
  EXPECT_GE(frac.y(), 0.0 - 1e-9);
  EXPECT_LT(frac.y(), 1.0 + 1e-9);
  EXPECT_GE(frac.z(), 0.0 - 1e-9);
  EXPECT_LT(frac.z(), 1.0 + 1e-9);
}

TEST_F(CellTests, TriclinicCellMinimumImage) {
  // Non-orthogonal cell
  Cell const cell({{5.0, 5.0, 5.0, 60.0, 60.0, 60.0}});

  // Distance vector that spans more than half the cell in some direction
  auto image = cell.minimumImage({4.0, 4.0, 4.0});
  real_t const image_length = correlation::math::norm(image);

  // The minimum image vector must be shorter than or equal to half the max box extent
  // For a cell with a=5, the maximum half-diagonal is bounded
  real_t const half_diagonal =
      static_cast<real_t>(0.5) * std::sqrt(static_cast<real_t>(5.0) * static_cast<real_t>(5.0) *
                                           static_cast<real_t>(3)); // conservative upper bound
  EXPECT_LE(image_length, half_diagonal + 1e-6);

  // The zero vector should map to zero
  auto zero_image = cell.minimumImage({0.0, 0.0, 0.0});
  EXPECT_NEAR(zero_image.x(), 0.0, 1e-9);
  EXPECT_NEAR(zero_image.y(), 0.0, 1e-9);
  EXPECT_NEAR(zero_image.z(), 0.0, 1e-9);
}

TEST_F(CellTests, ExtremelySmallCell) {
  // Very small cell — should not cause underflow or precision issues
  // Note: Cell::updateLattice rejects volume <= 1e-9, so 0.01^3 = 1e-6 is valid
  const std::array<real_t, 6> params = {0.01, 0.01, 0.01, 90.0, 90.0, 90.0};
  Cell cell(params);
  EXPECT_NEAR(cell.volume(), 1e-6, correlation::is_single_precision ? 1e-10 : 1e-12);

  cell.addAtom("H", {0.005, 0.005, 0.005});
  cell.wrapPositions();
  const auto &pos = cell.atoms().front().position();
  EXPECT_NEAR(pos.x(), 0.005, correlation::is_single_precision ? 1e-9 : 1e-12);
}

TEST_F(CellTests, ExtremelyLargeCell) {
  // Very large cell — should not overflow
  const std::array<real_t, 6> params = {1e6, 1e6, 1e6, 90.0, 90.0, 90.0};
  Cell cell(params);
  EXPECT_NEAR(cell.volume(), 1e18, correlation::is_single_precision ? 1.5e11 : 1e9);

  cell.addAtom("H", {5e5, 5e5, 5e5});
  auto image = cell.minimumImage({3e5, 0.0, 0.0});
  EXPECT_NEAR(image.x(), 3e5, 1e-3);
}

TEST_F(CellTests, HighAtomCount) {
  // Stress test: add many atoms without crashing
  Cell cell({{100.0, 100.0, 100.0, 90.0, 90.0, 90.0}});
  const size_t N_ATOMS = 10000;
  for (size_t i = 0; i < N_ATOMS; ++i) {
    real_t const pos = static_cast<real_t>(i) * static_cast<real_t>(0.01);
    cell.addAtom("H", {pos, pos, pos});
  }
  EXPECT_EQ(cell.atomCount(), N_ATOMS);
  EXPECT_EQ(cell.elements().size(), 1);

  // Wrap all positions — should not crash or take excessively long
  EXPECT_NO_THROW(cell.wrapPositions());
}

TEST_F(CellTests, AcosNumericalNoiseClamping) {
  Cell const cell({1.0, 0.0, 0.0}, {-1.0, 1.01e-8, 0.0}, {0.0, 0.0, 1.0});
  auto params = cell.lattice_parameters();
  EXPECT_FALSE(std::isnan(params[3]));
  EXPECT_FALSE(std::isnan(params[4]));
  EXPECT_FALSE(std::isnan(params[5]));
}

TEST_F(CellTests, FractionalCartesianRoundTripPrecision) {
  const std::array<real_t, 6> params = {12.0, 15.0, 18.0, 85.0, 95.0, 70.0};
  const Cell cell(params);
  const auto &lat = cell.latticeVectors();
  const auto &inv = cell.inverseLatticeVectors();

  // Test several fractional points (s_x, s_y, s_z) in [0, 1)
  const std::vector<correlation::math::Vector3<real_t>> test_fracs = {
      {0.1, 0.2, 0.3}, {0.0, 0.0, 0.0}, {0.99, 0.5, 0.25}, {0.333, 0.666, 0.123}};

  for (const auto &orig_frac : test_fracs) {
    // r = lattice * s
    const correlation::math::Vector3<real_t> cart = lat * orig_frac;
    // s_calc = inv_lattice * r
    const correlation::math::Vector3<real_t> calc_frac = inv * cart;

    EXPECT_NEAR(calc_frac.x(), orig_frac.x(), correlation::is_single_precision ? 1e-5 : 1e-9);
    EXPECT_NEAR(calc_frac.y(), orig_frac.y(), correlation::is_single_precision ? 1e-5 : 1e-9);
    EXPECT_NEAR(calc_frac.z(), orig_frac.z(), correlation::is_single_precision ? 1e-5 : 1e-9);
  }
}

} // namespace correlation::testing
