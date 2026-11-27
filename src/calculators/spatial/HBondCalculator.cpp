/**
 * @file HBondCalculator.cpp
 * @brief Implementation of the Hydrogen Bond Analysis calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/spatial/HBondCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "calculators/spatial/DistanceCalculator.hpp"
#include "math/Constants.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <tbb/enumerable_thread_specific.h>
#include <vector>

namespace correlation::calculators {

namespace {

struct HBondCriteria {
  real_t r_cut;
  real_t r_cut_sq;
  real_t cos_cut_sq;
};

struct DonorBond {
  correlation::math::Vector3<real_t> v_dh;
  real_t d_dh_sq{0.0};
};

// Static registration of the calculator in the factory
const bool REGISTERED = CalculatorFactory::registerTypeSafe<HBondCalculator>("HBondCalculator");

bool isElectronegative(const std::string &symbol) {
  return symbol == "O" || symbol == "N" || symbol == "F" || symbol == "S";
}

[[nodiscard]] int countDirectedHBonds(const std::vector<DonorBond> &donor_bonds,
                                      const correlation::math::Vector3<real_t> &v_da,
                                      real_t d_da_sq, real_t cos_cut_sq) noexcept {
  int matches = 0;
  for (const auto &bond : donor_bonds) {
    const real_t dot_val = bond.v_dh * v_da;
    if (dot_val > static_cast<real_t>(0.0) &&
        (dot_val * dot_val) > (cos_cut_sq * bond.d_dh_sq * d_da_sq)) {
      ++matches;
    }
  }
  return matches;
}

void findHydrogenBondsSmallBox(const correlation::core::Cell &cell,
                               const std::vector<size_t> &en_indices,
                               const std::vector<std::vector<DonorBond>> &donor_bonds,
                               const HBondCriteria &criteria, std::vector<int> &hbond_counts) {
  const auto &atoms = cell.atoms();
  for (const size_t donor_idx : en_indices) {
    if (donor_bonds[donor_idx].empty()) {
      continue;
    }
    const auto &pos_d = atoms[donor_idx].position();
    for (const size_t acc_idx : en_indices) {
      if (donor_idx == acc_idx) {
        continue;
      }
      const auto v_da = cell.minimumImage(atoms[acc_idx].position() - pos_d);
      const real_t d_da_sq = correlation::math::normSq(v_da);
      if (d_da_sq < criteria.r_cut_sq && d_da_sq >= static_cast<real_t>(1e-12)) {
        const int matches =
            countDirectedHBonds(donor_bonds[donor_idx], v_da, d_da_sq, criteria.cos_cut_sq);
        hbond_counts[donor_idx] += matches;
        hbond_counts[acc_idx] += matches;
      }
    }
  }
}

void findHydrogenBondsCellList(const correlation::core::Cell &cell,
                               const std::vector<uint8_t> &is_en,
                               const std::vector<std::vector<DonorBond>> &donor_bonds,
                               const HBondCriteria &criteria, std::vector<int> &hbond_counts) {
  const size_t num_atoms = cell.atomCount();
  tbb::enumerable_thread_specific<std::vector<int>> ets_counts(
      [num_atoms] { return std::vector<int>(num_atoms, 0); });

  DistanceCalculator::forEachPairWithin(
      cell, criteria.r_cut, true,
      [&](size_t atom_a, size_t atom_b, real_t d_sq,
          const correlation::math::Vector3<real_t> &r_ij) {
        if (atom_a >= num_atoms || atom_b >= num_atoms || is_en[atom_a] == 0U ||
            is_en[atom_b] == 0U) {
          return;
        }
        if (d_sq >= criteria.r_cut_sq || d_sq < static_cast<real_t>(1e-12)) {
          return;
        }
        const int ab_matches =
            countDirectedHBonds(donor_bonds[atom_a], r_ij, d_sq, criteria.cos_cut_sq);
        const int ba_matches =
            countDirectedHBonds(donor_bonds[atom_b], -r_ij, d_sq, criteria.cos_cut_sq);
        const int total = ab_matches + ba_matches;
        if (total > 0) {
          auto &local = ets_counts.local();
          local[atom_a] += total;
          local[atom_b] += total;
        }
      });

  for (const auto &local : ets_counts) {
    for (size_t i = 0; i < num_atoms; ++i) {
      hbond_counts[i] += local[i];
    }
  }
}

} // namespace

void HBondCalculator::calculateFrame(
    correlation::analysis::DistributionFunctions &dists,
    const correlation::analysis::AnalysisSettings & /*settings*/) const {
  dists.addHistogram("HBond", calculate(dists.cell(), dists.neighbors()));
}

correlation::analysis::Histogram
HBondCalculator::calculate(const correlation::core::Cell &cell,
                           const correlation::analysis::StructureAnalyzer *neighbors) {

  if (neighbors == nullptr) {
    return {};
  }

  const auto &atoms = cell.atoms();
  const size_t num_atoms = cell.atomCount();
  const auto &neighbor_graph = neighbors->neighborGraph();

  // Geometric criteria (R_DA < 3.5 A, angle H-D...A < 30 deg -> cos^2(alpha) > cos^2(30 deg))
  constexpr auto K_R_CUT = static_cast<real_t>(3.5);
  const real_t cos_cut =
      std::cos(static_cast<real_t>(30.0) * static_cast<real_t>(correlation::math::DEG_TO_RAD));
  const HBondCriteria criteria{
      .r_cut = K_R_CUT,
      .r_cut_sq = K_R_CUT * K_R_CUT,
      .cos_cut_sq = cos_cut * cos_cut,
  };

  std::vector<int> hbond_counts(num_atoms, 0);
  std::vector<uint8_t> is_en(num_atoms, 0U);
  std::vector<size_t> en_indices;
  en_indices.reserve(num_atoms / 4);
  std::vector<std::vector<DonorBond>> donor_bonds(num_atoms);

  for (size_t i = 0; i < num_atoms; ++i) {
    if (!isElectronegative(atoms[i].element().symbol)) {
      continue;
    }
    is_en[i] = 1U;
    en_indices.push_back(i);
    const auto &pos_d = atoms[i].position();
    for (const auto &nb_graph : neighbor_graph.getNeighbors(i)) {
      if (atoms[nb_graph.index].element().symbol == "H") {
        const auto v_dh = cell.minimumImage(atoms[nb_graph.index].position() - pos_d);
        const real_t d_dh_sq = correlation::math::normSq(v_dh);
        if (d_dh_sq >= static_cast<real_t>(1e-12)) {
          donor_bonds[i].push_back(DonorBond{.v_dh = v_dh, .d_dh_sq = d_dh_sq});
        }
      }
    }
  }

  const auto widths = cell.perpendicularWidths();
  const real_t min_width = std::min({widths[0], widths[1], widths[2]});
  if (min_width >= static_cast<real_t>(2.0) * criteria.r_cut) {
    findHydrogenBondsCellList(cell, is_en, donor_bonds, criteria, hbond_counts);
  } else {
    findHydrogenBondsSmallBox(cell, en_indices, donor_bonds, criteria, hbond_counts);
  }

  // Use pre-filtered en_indices for the final distribution loop.
  std::map<int, real_t> distribution;
  int max_hb = 0;
  const int num_en_atoms = static_cast<int>(en_indices.size());
  for (const size_t atom_idx : en_indices) {
    distribution[hbond_counts[atom_idx]] += static_cast<real_t>(1.0);
    max_hb = std::max(max_hb, hbond_counts[atom_idx]);
  }

  correlation::analysis::Histogram hist;
  hist.x_label = "N_HB";
  hist.title = "Hydrogen Bond Distribution";
  hist.y_label = "P(N_HB)";
  hist.x_unit = "bonds";
  hist.y_unit = "fraction";
  hist.description = "Distribution of Hydrogen Bonds per electronegative atom.";
  hist.file_suffix = "_HBond";

  for (int idx = 0; idx <= max_hb; ++idx) {
    hist.bins.push_back(static_cast<real_t>(idx));
    const real_t freq = (num_en_atoms > 0) ? (distribution[idx] / static_cast<real_t>(num_en_atoms))
                                           : static_cast<real_t>(0.0);
    hist.partials["Total"].push_back(freq);
  }

  return hist;
}

} // namespace correlation::calculators
