/**
 * @file TetrahedralCalculator.cpp
 * @brief Chau-Hardwick / Errington-Debenedetti tetrahedral order parameter calculator
 * implementation.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/order/TetrahedralCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "core/NeighborGraph.hpp"
#include "math/LinearAlgebra.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <map>
#include <sstream>
#include <tbb/parallel_for.h>

namespace correlation::calculators {

namespace {

const bool REGISTERED =
    CalculatorFactory::registerTypeSafe<TetrahedralCalculator>("TetrahedralCalculator");

/**
 * @brief Computes Chau-Hardwick tetrahedrality from 4 neighbor displacement vectors.
 *
 * Evaluates q_tet = 1 - 3/8 * sum_{j=0}^2 sum_{k=j+1}^3 (cos(psi_jk) + 1/3)^2.
 */
[[nodiscard]] real_t
computeTetrahedrality(const std::array<math::Vector3<real_t>, 4> &displacements) noexcept {
  std::array<math::Vector3<real_t>, 4> unit_vecs{};
  for (size_t i = 0; i < 4; ++i) {
    const real_t length = math::norm(displacements.at(i));
    if (length < static_cast<real_t>(1e-8)) {
      return static_cast<real_t>(0.0);
    }
    unit_vecs.at(i) = displacements.at(i) / length;
  }

  auto term_sum = static_cast<real_t>(0.0);
  constexpr real_t K_ONE_THIRD = static_cast<real_t>(1.0) / static_cast<real_t>(3.0);

  for (size_t j = 0; j < 3; ++j) {
    for (size_t k = j + 1; k < 4; ++k) {
      const real_t cos_psi = std::clamp(unit_vecs.at(j) * unit_vecs.at(k), static_cast<real_t>(-1.0),
                                        static_cast<real_t>(1.0));
      const real_t dev = cos_psi + K_ONE_THIRD;
      term_sum += dev * dev;
    }
  }

  const real_t q_tet =
      static_cast<real_t>(1.0) - (static_cast<real_t>(3.0) / static_cast<real_t>(8.0)) * term_sum;
  return std::clamp(q_tet, static_cast<real_t>(0.0), static_cast<real_t>(1.0));
}

/**
 * @brief Extracts 4 nearest neighbor vectors from NeighborGraph.
 */
[[nodiscard]] std::optional<std::array<math::Vector3<real_t>, 4>>
extractNearestFromGraph(std::span<const core::Neighbor> neighbors) {
  if (neighbors.size() < 4) {
    return std::nullopt;
  }

  std::vector<core::Neighbor> sorted_neighbors(neighbors.begin(), neighbors.end());
  std::partial_sort(sorted_neighbors.begin(), sorted_neighbors.begin() + 4, sorted_neighbors.end(),
                    [](const core::Neighbor &lhs, const core::Neighbor &rhs) noexcept {
                      return lhs.distance < rhs.distance;
                    });

  std::array<math::Vector3<real_t>, 4> displacements{};
  for (size_t k = 0; k < 4; ++k) {
    displacements.at(k) = sorted_neighbors.at(k).r_ij;
  }
  return displacements;
}

/**
 * @brief Extracts 4 nearest neighbor vectors using Cell minimum-image distance search.
 */
[[nodiscard]] std::optional<std::array<math::Vector3<real_t>, 4>>
extractNearestFromCell(const correlation::core::Cell &cell, size_t target_idx) {
  const auto &atoms = cell.atoms();
  const size_t num_atoms = atoms.size();
  if (num_atoms < 5) {
    return std::nullopt;
  }

  const auto &target_pos = atoms[target_idx].position();

  struct DistVec {
    real_t dist_sq{0.0};
    math::Vector3<real_t> r_ij;
  };

  std::vector<DistVec> candidates;
  candidates.reserve(num_atoms - 1);

  for (size_t j = 0; j < num_atoms; ++j) {
    if (j == target_idx) {
      continue;
    }
    const math::Vector3<real_t> disp = cell.minimumImage(atoms[j].position() - target_pos);
    const real_t dist_sq_val = math::norm_sq(disp);
    candidates.push_back(DistVec{.dist_sq = dist_sq_val, .r_ij = disp});
  }

  std::partial_sort(
      candidates.begin(), candidates.begin() + 4, candidates.end(),
      [](const DistVec &lhs, const DistVec &rhs) noexcept { return lhs.dist_sq < rhs.dist_sq; });

  std::array<math::Vector3<real_t>, 4> displacements{};
  for (size_t k = 0; k < 4; ++k) {
    displacements.at(k) = candidates.at(k).r_ij;
  }
  return displacements;
}

/**
 * @brief Retrieves 4 nearest neighbor displacement vectors for atom.
 */
[[nodiscard]] std::optional<std::array<math::Vector3<real_t>, 4>>
findNearest4Neighbors(const correlation::core::Cell &cell,
                      const correlation::analysis::StructureAnalyzer *neighbors, size_t atom_idx) {
  if (neighbors != nullptr) {
    const auto &graph = neighbors->neighborGraph();
    if (atom_idx < graph.nodeCount()) {
      auto result = extractNearestFromGraph(graph.getNeighbors(atom_idx));
      if (result.has_value()) {
        return result;
      }
    }
  }
  return extractNearestFromCell(cell, atom_idx);
}

} // namespace

void TetrahedralCalculator::calculateFrame(
    correlation::analysis::DistributionFunctions &dists,
    const correlation::analysis::AnalysisSettings & /*settings*/) const {
  auto hist = calculate(dists.cell(), dists.neighbors());
  dists.addHistogram("Q_tet", std::move(hist));
}

std::vector<real_t>
TetrahedralCalculator::calculatePerAtom(const correlation::core::Cell &cell,
                                        const correlation::analysis::StructureAnalyzer *neighbors,
                                        TetrahedralParams /*params*/) {
  const auto &atoms = cell.atoms();
  const size_t num_atoms = atoms.size();
  if (num_atoms < 5) {
    return std::vector<real_t>(num_atoms, static_cast<real_t>(0.0));
  }

  std::vector<real_t> q_tet_values(num_atoms, static_cast<real_t>(0.0));

  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_atoms),
                    [&](const tbb::blocked_range<size_t> &range) {
                      for (size_t i = range.begin(); i < range.end(); ++i) {
                        const auto nearest = findNearest4Neighbors(cell, neighbors, i);
                        if (nearest.has_value()) {
                          q_tet_values[i] = computeTetrahedrality(*nearest);
                        } else {
                          q_tet_values[i] = static_cast<real_t>(0.0);
                        }
                      }
                    });

  return q_tet_values;
}

correlation::analysis::Histogram
TetrahedralCalculator::calculate(const correlation::core::Cell &cell,
                                 const correlation::analysis::StructureAnalyzer *neighbors,
                                 TetrahedralParams params) {
  const auto q_tet = calculatePerAtom(cell, neighbors, params);
  const auto &atoms = cell.atoms();
  const size_t num_atoms = atoms.size();

  correlation::analysis::Histogram hist;
  hist.title = "Chau-Hardwick Tetrahedral Order Parameter";
  hist.x_label = "Tetrahedrality q_tet";
  hist.y_label = "Probability";
  hist.x_unit = "order";
  hist.y_unit = "probability";
  hist.description = "Chau-Hardwick 4-fold angular tetrahedral order parameter";
  hist.file_suffix = "_q_tet";

  const size_t num_bins = std::max<size_t>(params.num_bins, 1);
  const real_t bin_width = (params.max_domain - params.min_domain) / static_cast<real_t>(num_bins);
  hist.bins.resize(num_bins);
  for (size_t bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
    hist.bins.at(bin_idx) =
        params.min_domain + (static_cast<real_t>(bin_idx) + static_cast<real_t>(0.5)) * bin_width;
  }

  if (num_atoms == 0) {
    return hist;
  }

  std::vector<std::string> element_symbols;
  for (const auto &elem : cell.elements()) {
    element_symbols.push_back(elem.symbol);
  }

  std::map<std::string, std::vector<real_t>> partials;
  for (const auto &sym : element_symbols) {
    partials[sym].assign(num_bins, static_cast<real_t>(0.0));
  }
  partials["Total"].assign(num_bins, static_cast<real_t>(0.0));

  for (size_t i = 0; i < num_atoms; ++i) {
    const real_t val = q_tet[i];
    const real_t clamped_val = std::clamp(val, params.min_domain, params.max_domain);
    auto bin_idx = static_cast<size_t>((clamped_val - params.min_domain) / bin_width);
    if (bin_idx >= num_bins) {
      bin_idx = num_bins - 1;
    }
    const std::string &sym = atoms[i].element().symbol;
    partials[sym][bin_idx] += static_cast<real_t>(1.0);
    partials["Total"][bin_idx] += static_cast<real_t>(1.0);
  }

  const real_t inv_total = static_cast<real_t>(1.0) / static_cast<real_t>(num_atoms);
  for (auto &[sym, counts] : partials) {
    for (auto &count : counts) {
      count *= inv_total;
    }
    hist.partials[sym] = counts;
  }

  return hist;
}

} // namespace correlation::calculators
