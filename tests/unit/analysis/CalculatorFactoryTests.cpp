// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/CalculatorFactory.hpp"
#include "calculators/BaseCalculator.hpp"

#include <gtest/gtest.h>
#include <memory>

namespace correlation::testing {

using namespace correlation::calculators;

class MockCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "MockCalculator"; }
  std::string getShortName() const override { return "Mock"; }
  std::string getGroup() const override { return "Test"; }
  std::string getDescription() const override { return "A mock calculator for testing."; }
  bool isFrameCalculator() const override { return true; }
  bool isTrajectoryCalculator() const override { return false; }
};

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
  std::string first_calc_name = calculators[0]->getName();
  
  const BaseCalculator *retrieved = factory.getCalculator(first_calc_name);
  ASSERT_NE(retrieved, nullptr);
  EXPECT_EQ(retrieved->getName(), first_calc_name);
}

TEST(CalculatorFactoryTests, RegisterAndLookupCustomCalculator) {
  auto &factory = CalculatorFactory::instance();
  
  // Register a custom mock calculator
  auto mock = std::make_unique<MockCalculator>();
  bool result = factory.registerCalculator(std::move(mock));
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

} // namespace correlation::testing
