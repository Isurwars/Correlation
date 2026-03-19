// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/CalculatorFactory.hpp"

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

const BaseCalculator* CalculatorFactory::getCalculator(const std::string& name) const {
  for (const auto& calc : calculators_) {
    if (calc->getName() == name) {
      return calc.get();
    }
  }
  return nullptr;
}
