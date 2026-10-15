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
  const auto *calc_ptr = calculator.get();
  const auto &name = calc_ptr->getName();
  const auto &short_name = calc_ptr->getShortName();

  calculators_.push_back(std::move(calculator));

  if (!name.empty()) {
    index_.try_emplace(name, calc_ptr);
  }
  if (!short_name.empty()) {
    index_.try_emplace(short_name, calc_ptr);
  }
  return true;
}

const std::vector<std::unique_ptr<BaseCalculator>> &CalculatorFactory::getCalculators() const {
  return calculators_;
}

namespace {

const std::unordered_map<std::string_view, std::string_view> &getLegacyAliases() {
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
      {"vDoS", "Vibrational Density of States"},
  };
  return LEGACY_ALIASES;
}

} // namespace

const BaseCalculator *CalculatorFactory::getCalculator(std::string_view name) const {
  if (const auto iter = index_.find(name); iter != index_.end()) {
    return iter->second;
  }

  const auto &aliases = getLegacyAliases();
  if (const auto alias_iter = aliases.find(name); alias_iter != aliases.end()) {
    if (const auto resolved_iter = index_.find(alias_iter->second); resolved_iter != index_.end()) {
      return resolved_iter->second;
    }
  }

  return nullptr;
}

} // namespace correlation::calculators
