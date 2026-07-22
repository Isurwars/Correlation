/**
 * @file MLIPCalculatorTests.cpp
 * @brief Unit tests for MLIPCalculator and MLIPInterface.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/CalculatorFactory.hpp"
#include "calculators/MLIPCalculator.hpp"
#include "core/Cell.hpp"
#include "mlip/MLIPInterface.hpp"

#include <gtest/gtest.h>

namespace correlation::calculators {

namespace {
class MockORBv3Model : public correlation::mlip::MLIPInterface {
public:
  [[nodiscard]] std::string getModelName() const override { return "ORB-v3"; }

  [[nodiscard]] correlation::mlip::MLIPOutput evaluate(const correlation::core::Cell &cell) const override {
    correlation::mlip::MLIPOutput out;
    const size_t n_atoms = cell.atoms().size();
    out.total_energy = static_cast<real_t>(-42.5);
    out.per_atom_energy.resize(n_atoms, static_cast<real_t>(-4.25));
    out.forces.resize(n_atoms, correlation::math::Vector3<real_t>{0.1, -0.2, 0.3});
    out.stress = correlation::math::Matrix3<real_t>();
    return out;
  }
};
} // namespace

TEST(MLIPCalculatorTests, DiscoveryInCalculatorFactory) {
  const auto &factory = CalculatorFactory::instance();
  const auto *calc = factory.getCalculator("ML Interatomic Potential (ORB-v3)");

  ASSERT_NE(calc, nullptr);
  EXPECT_EQ(calc->getShortName(), "MLIP");
  EXPECT_EQ(calc->getGroup(), "Machine Learning");
  EXPECT_TRUE(calc->isFrameCalculator());
  EXPECT_FALSE(calc->isTrajectoryCalculator());
}

TEST(MLIPCalculatorTests, CalculateWithMockORBv3Model) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {1.5, 1.5, 1.5});

  MockORBv3Model const orb_model;
  auto const output = MLIPCalculator::calculate(cell, &orb_model);

  EXPECT_NEAR(output.total_energy, static_cast<real_t>(-42.5), 1e-5);
  ASSERT_EQ(output.per_atom_energy.size(), 2U);
  EXPECT_NEAR(output.per_atom_energy[0], static_cast<real_t>(-4.25), 1e-5);
  ASSERT_EQ(output.forces.size(), 2U);
  EXPECT_NEAR(output.forces[0].x(), static_cast<real_t>(0.1), 1e-5);
}

TEST(MLIPCalculatorTests, FallbackExecutionWithoutModel) {
  correlation::core::Cell cell({5.0, 5.0, 5.0, 90.0, 90.0, 90.0});
  cell.addAtom("C", {0.0, 0.0, 0.0});

  auto const output = MLIPCalculator::calculate(cell, nullptr);

  EXPECT_NEAR(output.total_energy, static_cast<real_t>(0.0), 1e-5);
  ASSERT_EQ(output.forces.size(), 1U);
}

} // namespace correlation::calculators
