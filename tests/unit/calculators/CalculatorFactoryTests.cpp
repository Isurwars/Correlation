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
  [[nodiscard]] std::string getName() const override { return "MockCalculator"; }
  [[nodiscard]] std::string getShortName() const override { return "Mock"; }
  [[nodiscard]] std::string getGroup() const override { return "Test"; }
  [[nodiscard]] std::string getDescription() const override { return "A mock calculator for testing."; }
  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }
};
} // namespace
TEST(CalculatorFactoryTests, SingletonInstanceIsUnique) {
  auto &factory1 = CalculatorFactory::instance();
  auto &factory2 = CalculatorFactory::instance();
  EXPECT_EQ(&factory1, &factory2);
}

TEST(CalculatorFactoryTests, GetRegisteredCalculators) {
  auto &factory = CalculatorFactory::instance();
  const auto &calculators = factory.getCalculators();
  // There should be some default calculators registered automatically
  EXPECT_GT(calculators.size(), 0);
}

TEST(CalculatorFactoryTests, LookupStandardCalculators) {
  auto &factory = CalculatorFactory::instance();
  // "RDF" is a standard calculator registered under calculators_obj
  // Wait, let's make sure the name matches. RDFCalculator name is "RDF" or "Radial Distribution Function"?
  // Let's verify RDFCalculator::getName() or short name.
  // Actually, let's look up by the name.
  // Let's find one standard calculator. We know RDFCalculator is registered.
  // Let's check getCalculator with different potential names or just use any from getCalculators()
  const auto &calculators = factory.getCalculators();
  ASSERT_FALSE(calculators.empty());
  std::string const first_calc_name = calculators[0]->getName();

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
  auto &factory = CalculatorFactory::instance();
  const BaseCalculator *retrieved = factory.getCalculator("NonExistentCalculatorName");
  EXPECT_EQ(retrieved, nullptr);
}

TEST(CalculatorFactoryTests, EmptyNameLookupReturnsNullptr) {
  auto &factory = CalculatorFactory::instance();
  const BaseCalculator *retrieved = factory.getCalculator("");
  EXPECT_EQ(retrieved, nullptr);
}

} // namespace correlation::testing
