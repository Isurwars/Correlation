/**
 * @file CalculatorFactory.cpp
 * @brief Implementation of the calculator factory.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/CalculatorFactory.hpp"

namespace correlation::calculators {

CalculatorFactory &CalculatorFactory::instance() {
  static CalculatorFactory instance_val;
  return instance_val;
}

bool CalculatorFactory::registerCalculator(std::unique_ptr<BaseCalculator> calculator) {
  if (!calculator) {
    return false;
  }
  calculators_.push_back(std::move(calculator));
  return true;
}

const std::vector<std::unique_ptr<BaseCalculator>> &CalculatorFactory::getCalculators() const {
  return calculators_;
}

const BaseCalculator *CalculatorFactory::getCalculator(std::string_view name) const {
  for (const auto &calc : calculators_) {
    if (calc->getName() == name || calc->getShortName() == name) {
      return calc.get();
    }
  }

  static const std::unordered_map<std::string_view, std::string_view> LEGACY_ALIASES = {
      {"Steinhardt Parameter — GPU Accelerated", "Steinhardt Order Parameters (GPU)"},
      {"S(Q) — GPU Accelerated", "Static Structure Factor (GPU)"},
      {"XRD — GPU Accelerated", "X-Ray Diffraction (GPU)"},
      {"Steinhardt Parameter", "Steinhardt Order Parameters"},
      {"Radial Distribution Function (RDF)", "Radial Distribution Function"},
      {"g(r), J(r), G(r)", "Radial Distribution Function"},
      {"BAD", "Plane Angle Distribution"},
      {"PAD", "Plane Angle Distribution"},
      {"DAD", "Dihedral Angle Distribution"},
      {"CN", "Coordination Number"},
      {"CNA", "Common Neighbor Analysis"},
      {"RD", "Ring Size Distribution"},
      {"Cluster Analysis", "Cluster Size Distribution"},
      {"Hydrogen Bond", "Hydrogen Bond Distribution"},
      {"S(K)", "Static Structure Factor"},
      {"S_K", "Static Structure Factor"},
      {"XRD", "X-Ray Diffraction"},
      {"Local Entropy", "Local Entropy Fingerprint"},
      {"σ²_N(R), χ_H(R)", "Hyperuniformity Diagnostics"},
      {"Dynamic", "Dynamical"},
      {"vDoS", "Vibrational Density of States"}};

  if (auto it = LEGACY_ALIASES.find(name); it != LEGACY_ALIASES.end()) {
    return getCalculator(it->second);
  }

  return nullptr;
}

} // namespace correlation::calculators
