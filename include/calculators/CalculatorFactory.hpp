/**
 * @file CalculatorFactory.hpp
 * @brief Factory for instantiating analysis calculator objects.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "BaseCalculator.hpp"
#include <memory>
#include <string>
#include <vector>

/**
 * @brief Registry for all calculators, enabling automatic discovery.
 */
class CalculatorFactory {
public:
  static CalculatorFactory &instance();

  /**
   * @brief Registers a new calculator.
   * @param calculator Unique pointer to the calculator instance.
   * @return true if registration was successful.
   */
  bool registerCalculator(std::unique_ptr<BaseCalculator> calculator);

  /**
   * @brief Returns all registered calculators.
   */
  const std::vector<std::unique_ptr<BaseCalculator>> &getCalculators() const;

  /**
   * @brief Gets a calculator by its name.
   * @param name The name of the calculator.
   * @return Pointer to the calculator, or nullptr if not found.
   */
  const BaseCalculator* getCalculator(const std::string& name) const;

private:
  CalculatorFactory() = default;
  std::vector<std::unique_ptr<BaseCalculator>> calculators_;
};
