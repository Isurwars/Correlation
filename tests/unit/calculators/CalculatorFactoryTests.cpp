// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/BaseCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"

#include <gtest/gtest.h>
#include <memory>

namespace correlation::testing {

using namespace correlation::calculators;
namespace {
class MockCalculator : public BaseCalculator {
public:
  [[nodiscard]] std::string_view getName() const override { return "MockCalculator"; }
  [[nodiscard]] std::string_view getShortName() const override { return "Mock"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Test"; }
  [[nodiscard]] std::string_view getDescription() const override {
    return "A mock calculator for testing.";
  }
  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }
};
} // namespace
TEST(CalculatorFactoryTests, SingletonInstanceIsUnique) {
  const auto &factory1 = CalculatorFactory::instance();
  const auto &factory2 = CalculatorFactory::instance();
  EXPECT_EQ(&factory1, &factory2);
}

TEST(CalculatorFactoryTests, GetRegisteredCalculators) {
  const auto &factory = CalculatorFactory::instance();
  const auto &calculators = factory.getCalculators();
  // There should be some default calculators registered automatically
  EXPECT_GT(calculators.size(), 0);
}

TEST(CalculatorFactoryTests, LookupStandardCalculators) {
  const auto &factory = CalculatorFactory::instance();
  const auto &calculators = factory.getCalculators();
  ASSERT_FALSE(calculators.empty());
  std::string_view const first_calc_name = calculators[0]->getName();

  const BaseCalculator *retrieved = factory.getCalculator(first_calc_name);
  ASSERT_NE(retrieved, nullptr);
  EXPECT_EQ(retrieved->getName(), first_calc_name);
}

TEST(CalculatorFactoryTests, RegisterAndLookupCustomCalculator) {
  auto &factory = CalculatorFactory::instance();

  // Register a custom mock calculator
  auto mock = std::make_unique<MockCalculator>();
  bool const result = factory.registerCalculator(std::move(mock));
  EXPECT_TRUE(result);

  // Verify it can be retrieved
  const BaseCalculator *retrieved = factory.getCalculator("MockCalculator");
  ASSERT_NE(retrieved, nullptr);
  EXPECT_EQ(retrieved->getName(), "MockCalculator");
  EXPECT_EQ(retrieved->getShortName(), "Mock");
  EXPECT_EQ(retrieved->getGroup(), "Test");
}

TEST(CalculatorFactoryTests, LookupNonExistentReturnsNullptr) {
  const auto &factory = CalculatorFactory::instance();
  const BaseCalculator *retrieved = factory.getCalculator("NonExistentCalculatorName");
  EXPECT_EQ(retrieved, nullptr);
}

TEST(CalculatorFactoryTests, EmptyNameLookupReturnsNullptr) {
  const auto &factory = CalculatorFactory::instance();
  const BaseCalculator *retrieved = factory.getCalculator("");
  EXPECT_EQ(retrieved, nullptr);
}

TEST(CalculatorFactoryTests, BaseCalculatorSatisfiesSegregatedContracts) {
  const MockCalculator mock;

  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  // Polymorphic binding to IFrameCalculator role interface
  const IFrameCalculator &frame_calc = mock;
  EXPECT_NO_THROW({
    correlation::analysis::DistributionFunctions dists(cell);
    correlation::analysis::AnalysisSettings const settings;
    frame_calc.calculateFrame(dists, settings);
  });

  // Polymorphic binding to ITrajectoryCalculator role interface
  const ITrajectoryCalculator &traj_calc = mock;
  EXPECT_NO_THROW({
    correlation::analysis::DistributionFunctions dists(cell);
    correlation::core::Trajectory const traj;
    correlation::analysis::AnalysisSettings const settings;
    traj_calc.calculateTrajectory(dists, traj, settings);
  });
}

} // namespace correlation::testing
