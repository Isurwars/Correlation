/**
 * @file RDCalculator.cpp
 * @brief Implementation of the reduced distribution function calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/RDCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "calculators/MotifFinder.hpp"

#include <stdexcept>
#include <utility>

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
// NOLINTNEXTLINE(cert-err58-cpp)
const bool registered = CalculatorFactory::registerTypeSafe<RDCalculator>("RDCalculator");
} // namespace

void RDCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                  const correlation::analysis::AnalysisSettings &settings) const {
  if (dists.neighbors() == nullptr) {
    return;
  }
  dists.addHistogram("RD", calculate(dists.neighbors()->neighborGraph(), settings.max_ring_size));
}

correlation::analysis::Histogram RDCalculator::calculate(const correlation::core::NeighborGraph &graph,
                                                         size_t max_ring_size) {
  if (max_ring_size < 3) {
    throw std::invalid_argument("Max ring size must be at least 3");
  }

  auto rings = correlation::calculators::MotifFinder::findRings(graph, max_ring_size);

  correlation::analysis::Histogram f_motif;
  f_motif.x_label = "Ring Size";
  f_motif.title = "Ring Distribution";
  f_motif.y_label = "Frequency";
  f_motif.x_unit = "atoms";
  f_motif.y_unit = "counts";
  f_motif.description = "Ring Distribution";
  f_motif.file_suffix = "_RD";

  size_t const num_bins = max_ring_size >= 3 ? (max_ring_size - 2) : 0;
  f_motif.bins.resize(num_bins);
  for (size_t i = 0; i < num_bins; ++i) {
    f_motif.bins[i] = static_cast<real_t>(i + 3);
  }

  auto &partial_hist = f_motif.partials["Rings"];
  partial_hist.assign(num_bins, 0.0);

  size_t total_counts = 0;
  for (const auto &[size, count] : rings) {
    if (size >= 3 && std::cmp_less_equal(size, max_ring_size)) {
      auto const bin = static_cast<size_t>(size - 3);
      partial_hist[bin] = static_cast<real_t>(count);
      total_counts += count;
    }
  }

  if (total_counts > 0) {
    auto const total_counts_d = static_cast<real_t>(total_counts);
    for (size_t i = 0; i < num_bins; ++i) {
      partial_hist[i] /= total_counts_d;
    }
  }

  f_motif.partials["Total"] = partial_hist;

  return f_motif;
}

} // namespace correlation::calculators
