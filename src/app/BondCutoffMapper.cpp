/**
 * @file BondCutoffMapper.cpp
 * @brief Implementation of decoupled atomic bond cutoff transformations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/BondCutoffMapper.hpp"
#include "physics/PhysicalData.hpp"

#include <format>
#include <stdexcept>

namespace correlation::app {

real_t BondCutoffMapper::parseDistanceSafe(const std::string &str) {
  try {
    return static_cast<real_t>(std::stod(str));
  } catch (const std::exception &) {
    return static_cast<real_t>(0.0);
  }
}

int BondCutoffMapper::findElementIndex(std::span<const correlation::core::Element> elements,
                                       const std::string &symbol) {
  for (size_t idx = 0; idx < elements.size(); ++idx) {
    if (elements[idx].symbol == symbol) {
      return static_cast<int>(idx);
    }
  }
  return -1;
}

std::vector<CutoffEntry>
BondCutoffMapper::createDefaultCutoffEntries(std::span<const correlation::core::Element> elements) {
  std::vector<CutoffEntry> entries;
  const size_t num_elements = elements.size();
  entries.reserve((num_elements * (num_elements + 1)) / 2);

  auto safe_get_radius = [](const std::string &symbol) -> real_t {
    try {
      return physics::getCovalentRadius(symbol);
    } catch (const std::out_of_range &) {
      return static_cast<real_t>(1.5);
    }
  };

  for (size_t i = 0; i < num_elements; ++i) {
    const real_t radius_a = safe_get_radius(elements[i].symbol);
    for (size_t j = i; j < num_elements; ++j) {
      const real_t radius_b = safe_get_radius(elements[j].symbol);
      const real_t sum_radii = radius_a + radius_b;
      const real_t min_d = sum_radii * static_cast<real_t>(0.6);
      const real_t max_d = sum_radii * static_cast<real_t>(1.2);

      entries.push_back(CutoffEntry{
          .element1 = elements[i].symbol,
          .element2 = elements[j].symbol,
          .min_distance = std::format("{:.2f}", min_d),
          .max_distance = std::format("{:.2f}", max_d),
      });
    }
  }

  return entries;
}

correlation::analysis::BondCutoffMatrix
BondCutoffMapper::parseCutoffMatrix(std::span<const CutoffEntry> entries,
                                    std::span<const correlation::core::Element> elements) {
  const size_t num_elements = elements.size();
  correlation::analysis::BondCutoffMatrix cutoffs(
      num_elements, std::vector<correlation::analysis::BondCutoffRange>(
                        num_elements, correlation::analysis::BondCutoffRange{
                                          .min_sq = static_cast<real_t>(0.0),
                                          .max_sq = static_cast<real_t>(0.0)}));

  for (const auto &entry : entries) {
    const real_t min_dist = parseDistanceSafe(entry.min_distance);
    const real_t max_dist = parseDistanceSafe(entry.max_distance);

    const int idx1 = findElementIndex(elements, entry.element1);
    const int idx2 = findElementIndex(elements, entry.element2);

    if (idx1 != -1 && idx2 != -1) {
      const correlation::analysis::BondCutoffRange range{
          .min_sq = min_dist * min_dist,
          .max_sq = max_dist * max_dist,
      };
      cutoffs[idx1][idx2] = range;
      cutoffs[idx2][idx1] = range;
    }
  }

  return cutoffs;
}

} // namespace correlation::app
