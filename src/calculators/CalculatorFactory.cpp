/**
 * @file CalculatorFactory.cpp
 * @brief Implementation of the calculator factory.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "calculators/CalculatorFactory.hpp"

namespace correlation::calculators {

CalculatorFactory &CalculatorFactory::instance() {
  static CalculatorFactory instance_val;
  return instance_val;
}

bool CalculatorFactory::registerCalculator(
    std::unique_ptr<BaseCalculator> calculator) {
  if (!calculator)
    return false;
  calculators_.push_back(std::move(calculator));
  return true;
}

const std::vector<std::unique_ptr<BaseCalculator>> &
CalculatorFactory::getCalculators() const {
  return calculators_;
}

const BaseCalculator *
CalculatorFactory::getCalculator(const std::string &name) const {
  for (const auto &calc : calculators_) {
    if (calc->getName() == name) {
      return calc.get();
    }
  }
  return nullptr;
}

} // namespace correlation::calculators
