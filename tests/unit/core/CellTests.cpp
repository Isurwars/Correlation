// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Cell.hpp"
#include "math/Constants.hpp"
#include "math/LinearAlgebra.hpp"

#include <gtest/gtest.h>
#include <cmath>

namespace correlation::testing {

using namespace correlation::core;

class CellTests : public ::testing::Test {
protected:
  Cell orthogonal_cell_{{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}};
};

// --- Unitary Tests: Constructors & Lattice ---

TEST_F(CellTests, DefaultConstructorInitializesEmpty) {
  const Cell cell{};
  EXPECT_TRUE(cell.isEmpty());
  EXPECT_EQ(cell.atomCount(), 0);
  EXPECT_EQ(cell.volume(), 0.0);
}

TEST_F(CellTests, ParameterConstructorSetsCorrectVolume) {
  const std::array<double, 6> params = {4.0, 4.0, 4.0, 90.0, 90.0, 90.0};
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
  const std::array<double, 6> params = {5.0, 6.0, 7.0, 80.0, 90.0, 100.0};
  const Cell cell(params);

  // Expected volume calculation
  const double cos_a = std::cos(80.0 * correlation::math::deg_to_rad);
  const double cos_b = std::cos(90.0 * correlation::math::deg_to_rad);
  const double cos_g = std::cos(100.0 * correlation::math::deg_to_rad);
  const double vol_sqrt = 1.0 - cos_a * cos_a - cos_b * cos_b - cos_g * cos_g +
                          2 * cos_a * cos_b * cos_g;
  const double expected_volume = 5.0 * 6.0 * 7.0 * std::sqrt(vol_sqrt);

  EXPECT_NEAR(cell.volume(), expected_volume, 1e-9);
}

TEST_F(CellTests, RuleOfFiveWorksCorrectly) {
  const std::array<double, 6> params = {4.0, 4.0, 4.0, 90.0, 90.0, 90.0};
  Cell cell(params);
  cell.addAtom("H", {0.5, 0.5, 0.5});
  cell.setEnergy(1.23);

  // Copy Constructor
  Cell cell_copy(cell);
  EXPECT_EQ(cell_copy.volume(), cell.volume());
  EXPECT_EQ(cell_copy.atomCount(), cell.atomCount());
  EXPECT_DOUBLE_EQ(cell_copy.getEnergy(), 1.23);

  // Copy Assignment
  Cell cell_assigned{};
  cell_assigned = cell;
  EXPECT_EQ(cell_assigned.volume(), cell.volume());
  EXPECT_EQ(cell_assigned.atomCount(), cell.atomCount());
}

TEST_F(CellTests, MoveSemanticsLeavesMovedFromStateEmpty) {
  Cell cell({{4.0, 4.0, 4.0, 90.0, 90.0, 90.0}});
  cell.addAtom("H", {0.5, 0.5, 0.5});
  
  Cell cell_moved(std::move(cell));
  EXPECT_TRUE(cell.isEmpty());
  EXPECT_EQ(cell_moved.atomCount(), 1);
}

TEST_F(CellTests, AccessorsWorkCorrectly) {
  Cell cell{};
  cell.setLatticeParameters({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  EXPECT_NEAR(cell.volume(), 1000.0, 1e-9);
  
  cell.setEnergy(-13.6);
  EXPECT_DOUBLE_EQ(cell.getEnergy(), -13.6);
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
  EXPECT_NEAR(pos.x(), 2.0, 1e-9);
  EXPECT_NEAR(pos.y(), 7.0, 1e-9);
  EXPECT_NEAR(pos.z(), 5.0, 1e-9);
}

TEST_F(CellTests, WrapPositionsHandlesLargeDisplacements) {
  Cell cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  cell.addAtom("H", {1002.5, -997.5, 0.0});
  cell.wrapPositions();
  const auto &pos = cell.atoms().front().position();
  EXPECT_NEAR(pos.x(), 2.5, 1e-9);
  EXPECT_NEAR(pos.y(), 2.5, 1e-9);
}

TEST_F(CellTests, MinimumImageCalculatesShortestVector) {
  Cell cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  auto mi = cell.minimumImage({7.0, -8.0, 2.0});
  EXPECT_NEAR(mi.x(), -3.0, 1e-9);
  EXPECT_NEAR(mi.y(), 2.0, 1e-9);
  EXPECT_NEAR(mi.z(), 2.0, 1e-9);
}

TEST_F(CellTests, MinimumImageHandlesDistancesLargerThanBox) {
  Cell cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  auto mi = cell.minimumImage({15.0, 0.0, 0.0});
  EXPECT_NEAR(std::abs(mi.x()), 5.0, 1e-9);
}

TEST_F(CellTests, InverseLatticeVectorsAreCorrect) {
  const Cell cell({{2.0, 4.0, 5.0, 90.0, 90.0, 90.0}});
  const auto &inv = cell.inverseLatticeVectors();
  EXPECT_NEAR(inv[0][0], 0.5, 1e-9);
  EXPECT_NEAR(inv[1][1], 0.25, 1e-9);
  EXPECT_NEAR(inv[2][2], 0.2, 1e-9);
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

} // namespace correlation::testing
