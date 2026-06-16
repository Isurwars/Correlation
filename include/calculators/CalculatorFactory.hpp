/**
 * @file CalculatorFactory.hpp
 * @brief Factory for instantiating analysis calculator objects.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseCalculator.hpp"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace correlation::calculators {

/**
 * @brief Registry for all calculators, enabling automatic discovery.
 */
class CalculatorFactory {
public:
  /** @return Singleton instance of the CalculatorFactory. */
  static CalculatorFactory &instance();

  /**
   * @brief Registers a new calculator.
   * @param calculator Unique pointer to the calculator instance.
   * @return true if registration was successful.
   */
  bool registerCalculator(std::unique_ptr<BaseCalculator> calculator);

  /**
   * @brief Exception-safe template helper for static auto-registration.
   * @tparam T The specific Calculator class type to instantiate.
   * @param name The human-readable name of the calculator (for error logging).
   * @return true if registration succeeded, false if an exception occurred.
   */
  template <typename T> static bool registerTypeSafe(const char *name) noexcept {
    try {
      return instance().registerCalculator(std::make_unique<T>());
    } catch (const std::exception &e) {
      std::cerr << "[FATAL] Failed to statically register calculator '" << (name != nullptr ? name : "unknown")
                << "': " << e.what() << '\n';
      return false;
    } catch (...) {
      std::cerr << "[FATAL] Failed to statically register calculator '" << (name != nullptr ? name : "unknown")
                << "' due to an unknown exception." << '\n';
      return false;
    }
  }

  /**
   * @brief Returns all registered calculators.
   * @return Constant reference to the internal vector of calculators.
   */
  [[nodiscard]] const std::vector<std::unique_ptr<BaseCalculator>> &getCalculators() const;

  /**
   * @brief Gets a calculator by its name.
   * @param name The name of the calculator.
   * @return Pointer to the calculator, or nullptr if not found.
   */
  [[nodiscard]] const BaseCalculator *getCalculator(const std::string &name) const;

private:
  CalculatorFactory() = default;
  std::vector<std::unique_ptr<BaseCalculator>> calculators_; ///< Storage for calculator instances.
};

} // namespace correlation::calculators
