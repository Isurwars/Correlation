/**
 * @file ICalculatorRegistry.hpp
 * @brief Abstract interface for querying and discovering analysis calculators.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include <memory>
#include <string_view>
#include <vector>

namespace correlation::calculators {

class BaseCalculator;

/**
 * @brief Abstract interface defining registry lookup capabilities for analysis calculators.
 */
class ICalculatorRegistry {
public:
  ICalculatorRegistry() = default;
  virtual ~ICalculatorRegistry() = default;

  ICalculatorRegistry(const ICalculatorRegistry &) = delete;
  ICalculatorRegistry &operator=(const ICalculatorRegistry &) = delete;
  ICalculatorRegistry(ICalculatorRegistry &&) = delete;
  ICalculatorRegistry &operator=(ICalculatorRegistry &&) = delete;

  /**
   * @brief Returns all registered calculators.
   * @return Constant reference to vector of unique pointers to BaseCalculator.
   */
  [[nodiscard]] virtual const std::vector<std::unique_ptr<BaseCalculator>> &
  getCalculators() const = 0;

  /**
   * @brief Gets a calculator by its full name or recognized abbreviation/alias.
   * @param[in] name The identifier or name of the calculator.
   * @return Pointer to the calculator instance, or nullptr if not found.
   */
  [[nodiscard]] virtual const BaseCalculator *getCalculator(std::string_view name) const = 0;
};

} // namespace correlation::calculators
