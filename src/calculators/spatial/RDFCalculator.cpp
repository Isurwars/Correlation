/**
 * @file RDFCalculator.cpp
 * @brief Implementation of the radial distribution function calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/spatial/RDFCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "calculators/spatial/DistanceCalculator.hpp"
#include "math/Constants.hpp"
#include "math/Precision.hpp"
#include "math/SIMDUtils.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>

namespace correlation::calculators {

namespace {
std::string getPartialKey(const correlation::core::Cell &cell, size_t type1, size_t type2) {
  const auto &elements = cell.elements();
  if (type1 > type2) {
    std::swap(type1, type2);
  }
  return elements[type1].symbol + "-" + elements[type2].symbol;
}

std::string getInversePartialKey(const correlation::core::Cell &cell, size_t type1, size_t type2) {
  const auto &elements = cell.elements();
  if (type1 < type2) {
    std::swap(type1, type2);
  }
  return elements[type1].symbol + "-" + elements[type2].symbol;
}

// Static registration of the calculator in the factory
const bool REGISTERED = CalculatorFactory::registerTypeSafe<RDFCalculator>("RDFCalculator");

struct ActiveCellContext {
  correlation::core::Cell supercell_holder;
  bool is_expanded{false};
  real_t replication_factor{1.0};

  void initialize(const correlation::core::Cell &original, real_t r_max, real_t max_radius) {
    const auto widths = original.perpendicularWidths();
    const real_t eps = std::numeric_limits<real_t>::epsilon();
    if (widths[0] > eps && widths[1] > eps && widths[2] > eps) {
      const real_t target_width = static_cast<real_t>(2.0) * r_max;
      if (widths[0] < target_width || widths[1] < target_width || widths[2] < target_width) {
        supercell_holder = original.autoSupercell(r_max, 10, max_radius);
        if (supercell_holder.atomCount() > original.atomCount() && original.atomCount() > 0) {
          is_expanded = true;
          replication_factor = static_cast<real_t>(supercell_holder.atomCount()) /
                               static_cast<real_t>(original.atomCount());
        }
      }
    }
  }
};

struct RDFSettings {
  real_t r_max;
  real_t r_bin_width;
  size_t num_bins;
};

void accumulateRawCounts(const correlation::core::Cell &cell,
                         const correlation::analysis::StructureAnalyzer * /*neighbors*/,
                         RDFSettings settings, correlation::analysis::Histogram &h_r) {
  const auto &elements = cell.elements();
  const size_t num_elements = elements.size();

  RawHistogramTensor standalone_histograms;
  correlation::core::NeighborGraph dummy_graph;

  real_t const cutoff_sq = settings.r_max * settings.r_max;
  correlation::analysis::BondCutoffMatrix const empty_bonds;
  DistanceCalculationConfig const hist_config{
      .r_max = settings.r_max,
      .r_bin_width = settings.r_bin_width,
      .num_bins = settings.num_bins,
  };
  DistanceCalculator::compute(cell, cutoff_sq, empty_bonds, false, dummy_graph,
                              &standalone_histograms, hist_config);

  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      std::string const key = getPartialKey(cell, i, j);
      auto &partial_hist = h_r.partials[key];
      partial_hist.assign(settings.num_bins, 0.0);

      if (i < standalone_histograms.size() && j < standalone_histograms[i].size()) {
        const auto &src_bins = standalone_histograms[i][j];
        size_t const copy_count = std::min(settings.num_bins, src_bins.size());
        std::copy_n(src_bins.begin(), copy_count, partial_hist.begin());
      }

      // For self-pairs (A-A), each pair is counted once in the upper triangular
      // raw_histograms. We multiply by 2 to account for both A_1 -> A_2 and A_2
      // -> A_1 interactions.
      if (i == j) {
        correlation::math::scale_bins(partial_hist.data(), static_cast<real_t>(2.0),
                                      settings.num_bins);
      }
    }
  }
}

struct RDFNormalizationSettings {
  real_t volume;
  real_t bin_width;
  size_t num_bins;
};

void normalizeDistributions(const correlation::core::Cell &cell,
                            std::span<const real_t> element_counts,
                            RDFNormalizationSettings settings,
                            const correlation::analysis::Histogram &h_r,
                            correlation::analysis::Histogram &g_r,
                            correlation::analysis::Histogram &g_r_reduced,
                            correlation::analysis::Histogram &j_r) {
  const auto &elements = cell.elements();
  const size_t num_elements = elements.size();

  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      std::string const key = getPartialKey(cell, i, j);
      std::string const inversekey = getInversePartialKey(cell, i, j);

      const real_t n_i = element_counts[i];
      const real_t n_j = element_counts[j];

      const auto &h_ij = h_r.partials.at(key);

      g_r.partials[key].assign(settings.num_bins, 0.0);
      g_r_reduced.partials[key].assign(settings.num_bins, 0.0);
      j_r.partials[key].assign(settings.num_bins, 0.0);
      j_r.partials[inversekey].assign(settings.num_bins, 0.0);

      // Second Pass: Normalize the raw counts H(r) into target distribution
      // functions. g(r) normalization constant: V / (4 * pi * dr * n_i * n_j).
      // The r^2 term is applied per-bin inside the SIMD kernel.
      const real_t g_norm_constant =
          settings.volume / (correlation::math::FOUR_PI * settings.bin_width * n_i * n_j);
      const real_t rho_j = n_j / settings.volume;
      const real_t inv_ni_dr = static_cast<real_t>(1.0) / (n_i * settings.bin_width);
      const real_t inv_nj_dr = static_cast<real_t>(1.0) / (n_j * settings.bin_width);
      const real_t pi4_rho_j = correlation::math::FOUR_PI * rho_j;

      correlation::math::RDFNormalizationParams<real_t> const params{
          .hist_data = h_ij.data(),
          .radial_bins = g_r.bins.data(),
          .g_norm = g_norm_constant,
          .inv_Ni_dr = inv_ni_dr,
          .inv_Nj_dr = inv_nj_dr,
          .pi4_rho_j = pi4_rho_j,
          .g_out = g_r.partials[key].data(),
          .G_out = g_r_reduced.partials[key].data(),
          .J_out = j_r.partials[key].data(),
          .Jinv_out = j_r.partials[inversekey].data(),
          .count = settings.num_bins,
      };
      correlation::math::normalize_rdf_bins(params);
    }
  }
}

struct RDFWeightingSettings {
  real_t rho_0;
  size_t num_bins;
};

void weightPartials(const correlation::core::Cell &cell,
                    const std::map<std::string, real_t> &ashcroft_weights,
                    RDFWeightingSettings settings, correlation::analysis::Histogram &g_r,
                    correlation::analysis::Histogram &g_r_reduced) {
  const auto &elements = cell.elements();
  const size_t num_elements = elements.size();

  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      std::string const key = getPartialKey(cell, i, j);
      real_t const weight = ashcroft_weights.at(key);

      // 1. Weight g_r partial: g_ij_weighted(r) = w_ij * g_ij(r)
      auto &g_part = g_r.partials.at(key);
      for (size_t k = 0; k < settings.num_bins; ++k) {
        g_part[k] *= weight;
      }

      // 2. Weight G_r partial: G_ij_weighted(r) = w_ij * 4 * pi * rho_0 * r * (g_ij(r) - 1)
      // Since we already weighted g_part, we have g_part[k] = w_ij * g_ij(r).
      // Thus, G_ij_weighted(r) = 4 * pi * rho_0 * r * (g_part[k] - weight)
      auto &g_part_reduced = g_r_reduced.partials.at(key);
      for (size_t k = 0; k < settings.num_bins; ++k) {
        const real_t r_k = g_r.bins[k];
        if (r_k < 1e-9) {
          g_part_reduced[k] = 0.0;
        } else {
          g_part_reduced[k] =
              correlation::math::FOUR_PI * settings.rho_0 * r_k * (g_part[k] - weight);
        }
      }
    }
  }
}

std::vector<real_t> countCellElements(const correlation::core::Cell &cell) {
  const size_t num_elements = cell.elements().size();
  std::vector<real_t> element_counts(num_elements, static_cast<real_t>(0.0));
  for (const auto &atom : cell.atoms()) {
    const int elem_id = atom.elementId();
    if (std::cmp_greater_equal(elem_id, 0) && std::cmp_less(elem_id, num_elements)) {
      element_counts[static_cast<size_t>(elem_id)] += static_cast<real_t>(1.0);
    }
  }
  return element_counts;
}

void computeTotalDistributions(correlation::analysis::Histogram &g_r,
                               correlation::analysis::Histogram &g_r_reduced,
                               correlation::analysis::Histogram &j_r,
                               const std::map<std::string, real_t> &ashcroft_weights, real_t rho_0,
                               size_t num_bins) {
  auto &total_g = g_r.partials["Total"];
  total_g.assign(num_bins, 0.0);
  for (const auto &[key, g_partial] : g_r.partials) {
    if (key == "Total") {
      continue;
    }

    real_t const weight = ashcroft_weights.at(key);

    for (size_t k = 0; k < num_bins; ++k) {
      total_g[k] += g_partial[k] * weight;
    }
  }

  auto &total_j = j_r.partials["Total"];
  auto &total_g_reduced = g_r_reduced.partials["Total"];
  total_j.assign(num_bins, 0.0);
  total_g_reduced.assign(num_bins, 0.0);

  for (size_t k = 0; k < num_bins; ++k) {
    const real_t r_k = g_r.bins[k];
    if (r_k < 1e-9) {
      continue;
    }

    total_j[k] = correlation::math::FOUR_PI * r_k * r_k * rho_0 * total_g[k];
    total_g_reduced[k] =
        correlation::math::FOUR_PI * rho_0 * r_k * (total_g[k] - static_cast<real_t>(1.0));
  }
}

void rescaleSupercellCounts(correlation::analysis::Histogram &h_r,
                            const ActiveCellContext &cell_context) {
  if (!cell_context.is_expanded || cell_context.replication_factor <= static_cast<real_t>(1.0)) {
    return;
  }
  const real_t inv_scale = static_cast<real_t>(1.0) / cell_context.replication_factor;
  for (auto &[key, counts] : h_r.partials) {
    for (auto &val : counts) {
      val *= inv_scale;
    }
  }
}
} // namespace

void RDFCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                   const correlation::analysis::AnalysisSettings &settings) const {
  auto results = calculate(dists.cell(), dists.neighbors(), dists.getAshcroftWeights(),
                           settings.r_max, settings.r_bin_width);
  for (auto &[name, histogram] : results) {
    dists.addHistogram(name, std::move(histogram));
  }
}

std::map<std::string, correlation::analysis::Histogram>
RDFCalculator::calculate(const correlation::core::Cell &cell,
                         const correlation::analysis::StructureAnalyzer *neighbors,
                         const std::map<std::string, real_t> &ashcroft_weights, real_t r_max,
                         real_t r_bin_width, real_t max_radius) {
  if (r_bin_width <= 0) {
    throw std::invalid_argument("Bin width must be positive, got: " + std::to_string(r_bin_width));
  }
  if (r_max <= 0) {
    throw std::invalid_argument("Cutoff radius must be positive, got: " + std::to_string(r_max));
  }
  if (r_max > max_radius) {
    throw std::invalid_argument("Cutoff radius r_max (" + std::to_string(r_max) +
                                " Å) exceeds safe maximum radius of " + std::to_string(max_radius) +
                                " Å");
  }

  ActiveCellContext cell_context;
  cell_context.initialize(cell, r_max, max_radius);
  const auto &active_cell = cell_context.is_expanded ? cell_context.supercell_holder : cell;

  const real_t volume = active_cell.volume();
  if (volume <= std::numeric_limits<real_t>::epsilon()) {
    throw std::logic_error("Cell volume must be positive, got: " + std::to_string(volume));
  }

  const auto num_atoms = static_cast<real_t>(active_cell.atomCount());
  if (num_atoms == 0.0) {
    return {};
  }

  const std::vector<real_t> element_counts = countCellElements(active_cell);

  const auto num_bins = static_cast<size_t>(std::floor(r_max / r_bin_width));
  const real_t rho_0 = num_atoms / volume;

  correlation::analysis::Histogram h_r;
  correlation::analysis::Histogram g_r;
  correlation::analysis::Histogram g_r_reduced;
  correlation::analysis::Histogram j_r;
  h_r.bins.resize(num_bins);
  g_r.bins.resize(num_bins);
  g_r_reduced.bins.resize(num_bins);
  j_r.bins.resize(num_bins);
  h_r.x_label = "r";
  h_r.title = "H(r) — Distance Histogram";
  h_r.y_label = "H(r)";
  h_r.x_unit = "Å";
  h_r.y_unit = "counts";
  h_r.description = "Radial Distribution Function";
  h_r.file_suffix = "_H";

  g_r.x_label = "r";
  g_r.title = "g(r) — Radial Distribution Function";
  g_r.y_label = "g(r)";
  g_r.x_unit = "Å";
  g_r.y_unit = "";
  g_r.description = "Radial Distribution Function";
  g_r.file_suffix = "_g";

  g_r_reduced.x_label = "r";
  g_r_reduced.title = "G(r) — Reduced Pair Distribution Function";
  g_r_reduced.y_label = "G(r)";
  g_r_reduced.x_unit = "Å";
  g_r_reduced.y_unit = "Å⁻²";
  g_r_reduced.description = "Radial Distribution Function";
  g_r_reduced.file_suffix = "_G_reduced";

  j_r.x_label = "r";
  j_r.title = "J(r) — Reduced Pair Distribution";
  j_r.y_label = "J(r)";
  j_r.x_unit = "Å";
  j_r.y_unit = "Å⁻¹";
  j_r.description = "Radial Distribution Function";
  j_r.file_suffix = "_J";

  for (size_t i = 0; i < num_bins; ++i) {
    const real_t r_i = (static_cast<real_t>(i) + static_cast<real_t>(0.5)) * r_bin_width;
    h_r.bins[i] = r_i;
    g_r.bins[i] = r_i;
    g_r_reduced.bins[i] = r_i;
    j_r.bins[i] = r_i;
  }

  accumulateRawCounts(active_cell, neighbors,
                      {
                          .r_max = r_max,
                          .r_bin_width = r_bin_width,
                          .num_bins = num_bins,
                      },
                      h_r);

  normalizeDistributions(active_cell, element_counts,
                         {
                             .volume = volume,
                             .bin_width = r_bin_width,
                             .num_bins = num_bins,
                         },
                         h_r, g_r, g_r_reduced, j_r);

  computeTotalDistributions(g_r, g_r_reduced, j_r, ashcroft_weights, rho_0, num_bins);

  correlation::analysis::Histogram g_r_unweighted = g_r;
  g_r_unweighted.title = "g(r) — Unweighted Radial Distribution Function";
  g_r_unweighted.file_suffix = "_g_unweighted";

  weightPartials(active_cell, ashcroft_weights,
                 {
                     .rho_0 = rho_0,
                     .num_bins = num_bins,
                 },
                 g_r, g_r_reduced);

  rescaleSupercellCounts(h_r, cell_context);

  std::map<std::string, correlation::analysis::Histogram> results;
  results["H_r"] = std::move(h_r);
  results["g_r_unweighted"] = std::move(g_r_unweighted);
  results["J_r"] = std::move(j_r);
  results["g_r"] = std::move(g_r);
  results["G_r"] = std::move(g_r_reduced);

  return results;
}

} // namespace correlation::calculators
