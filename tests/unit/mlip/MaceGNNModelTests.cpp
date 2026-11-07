// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Cell.hpp"
#include "mlip/MaceGNNModel.hpp"

#include <gtest/gtest.h>

namespace correlation::mlip {
namespace {

TEST(MaceGNNModelTests, InterfaceAndAccessors) {
  MaceGNNModel model("non_existent_mace.pt", "cpu", 6.0);

  EXPECT_EQ(model.cutoff(), static_cast<real_t>(6.0));
  model.setCutoff(5.0);
  EXPECT_EQ(model.cutoff(), static_cast<real_t>(5.0));
  EXPECT_TRUE(model.getModelName().contains("non_existent_mace.pt"));
  EXPECT_TRUE(model.getModelName().starts_with("MACE ["));

#ifndef CORRELATION_HAS_LIBTORCH
  EXPECT_FALSE(model.isLoaded());
  correlation::core::Cell cell(std::array<real_t, 6>{10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  auto const output = model.evaluate(cell);
  EXPECT_EQ(output.forces.size(), 1U);
#endif
}

TEST(MaceGNNModelTests, MoveSemantics) {
  MaceGNNModel model1("sample_mace.model", "cpu", 5.5);
  EXPECT_EQ(model1.cutoff(), static_cast<real_t>(5.5));

  MaceGNNModel model2 = std::move(model1);
  EXPECT_EQ(model2.cutoff(), static_cast<real_t>(5.5));
  EXPECT_TRUE(model2.getModelName().contains("sample_mace.model"));

  MaceGNNModel model3("another_mace.pt", "cpu", 3.0);
  model3 = std::move(model2);
  EXPECT_EQ(model3.cutoff(), static_cast<real_t>(5.5));
}

TEST(MaceGNNModelTests, FallbackOnEmptyCell) {
  MaceGNNModel const model("dummy_mace.pt", "cpu", 5.0);
  correlation::core::Cell const empty_cell(
      std::array<real_t, 6>{10.0, 10.0, 10.0, 90.0, 90.0, 90.0});

  auto const output = model.evaluate(empty_cell);
  EXPECT_TRUE(output.forces.empty());
  EXPECT_DOUBLE_EQ(output.total_energy, 0.0);
}

TEST(MaceGNNModelTests, FallbackMultiAtomCellForces) {
  MaceGNNModel const model("dummy_mace.pt", "cpu", 5.0);
  correlation::core::Cell cell(std::array<real_t, 6>{10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("C", {0.0, 0.0, 0.0});
  cell.addAtom("O", {1.2, 0.0, 0.0});

  auto const output = model.evaluate(cell);
  EXPECT_EQ(output.forces.size(), 2U);
}

} // namespace
} // namespace correlation::mlip

