// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/MDCalculator.hpp"
#include "calculators/MotifFinder.hpp"
#include <stdexcept>

Histogram MDCalculator::calculate(const NeighborGraph &graph,
                                  size_t max_ring_size) {
  if (max_ring_size < 3) {
    throw std::invalid_argument("Max ring size must be at least 3");
  }

  auto rings = calculators::MotifFinder::findRings(graph, max_ring_size);

  Histogram f_motif;
  f_motif.bin_label = "Ring Size";

  // Create bins from 3 to max_ring_size
  size_t num_bins = max_ring_size - 2;
  f_motif.bins.resize(num_bins);
  for (size_t i = 0; i < num_bins; ++i) {
    f_motif.bins[i] = static_cast<double>(i + 3);
  }

  auto &partial_hist = f_motif.partials["Rings"];
  partial_hist.assign(num_bins, 0.0);

  double total_counts = 0;
  for (const auto &[size, count] : rings) {
    if (size >= 3 && size <= static_cast<int>(max_ring_size)) {
      size_t bin = static_cast<size_t>(size - 3);
      partial_hist[bin] = static_cast<double>(count);
      total_counts += count;
    }
  }

  // Not strictly normalized the same way as angular distributions,
  // since this is discrete occurrences. The "Total" can just be
  // a copy of the "Rings" partial if desired, or we just leave it
  // with "Rings". For consistency with FileWriter, we might just
  // provide a Total if we expect one, but FileWriter works on
  // any partial keys. We will provide "Total" as identical for now.

  f_motif.partials["Total"] = partial_hist;

  // Discrete rings often preferred as raw counts rather than a probability
  // density with bin_width=1. We'll leave them as counts to allow users
  // to see exactly how many 5-membered rings were found.

  return f_motif;
}
