/**
 * @file HBondCalculator.cpp
 * @brief Implementation of the Hydrogen Bond Analysis calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/HBondCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Constants.hpp"
#include <algorithm>
#include <vector>

namespace correlation::calculators {

namespace {

struct HBondCriteria {
  double R_cut_sq;
  double Alpha_cut;
};

// Static registration of the calculator in the factory
// NOLINTNEXTLINE(cert-err58-cpp)
const bool registered = CalculatorFactory::registerTypeSafe<HBondCalculator>("HBondCalculator");

bool isElectronegative(const std::string &symbol) {
  return symbol == "O" || symbol == "N" || symbol == "F" || symbol == "S";
}

std::vector<size_t> findBondedHydrogens(size_t donor_idx, const std::vector<correlation::core::Atom> &atoms,
                                        const correlation::core::NeighborGraph &neighbor_graph) {
  std::vector<size_t> hydrogens;
  for (const auto &nb_graph : neighbor_graph.getNeighbors(donor_idx)) {
    if (atoms[nb_graph.index].element().symbol == "H") {
      hydrogens.push_back(nb_graph.index);
    }
  }
  return hydrogens;
}

void checkAcceptorsForHydrogen(const correlation::core::Cell &cell, size_t donor_idx, size_t h_idx,
                               const std::vector<size_t> &en_indices, const HBondCriteria &criteria,
                               std::vector<int> &hbond_counts) {
  const auto &atoms = cell.atoms();
  const auto &pos_h = atoms[h_idx].position();
  const auto &pos_d = atoms[donor_idx].position();
  correlation::math::Vector3<double> v_dh = cell.minimumImage(pos_h - pos_d);
  double d_dh_sq = correlation::math::norm_sq(v_dh);
  if (d_dh_sq < 1e-12) {
    return;
  }
  double norm_dh = std::sqrt(d_dh_sq);

  for (size_t atom_idx : en_indices) {
    if (donor_idx == atom_idx) {
      continue;
    }
    const auto &atom_j = atoms[atom_idx];
    correlation::math::Vector3<double> v_da = cell.minimumImage(atom_j.position() - pos_d);
    double d_da_sq = correlation::math::norm_sq(v_da);

    if (d_da_sq < criteria.R_cut_sq && d_da_sq >= 1e-12) {
      double dot_val = v_dh * v_da;
      double cos_alpha = dot_val / (norm_dh * std::sqrt(d_da_sq));
      double alpha = std::acos(std::max(-1.0, std::min(1.0, cos_alpha))) * correlation::math::rad_to_deg;

      if (alpha < criteria.Alpha_cut) {
        hbond_counts[donor_idx]++;
        hbond_counts[atom_idx]++;
      }
    }
  }
}

void findHydrogenBonds(const correlation::core::Cell &cell, const correlation::core::NeighborGraph &neighbor_graph,
                       const std::vector<size_t> &en_indices, const HBondCriteria &criteria,
                       std::vector<int> &hbond_counts) {
  const auto &atoms = cell.atoms();
  for (size_t donor_idx : en_indices) {
    std::vector<size_t> hydrogens = findBondedHydrogens(donor_idx, atoms, neighbor_graph);
    if (hydrogens.empty()) {
      continue;
    }

    for (size_t h_idx : hydrogens) {
      checkAcceptorsForHydrogen(cell, donor_idx, h_idx, en_indices, criteria, hbond_counts);
    }
  }
}

} // namespace

void HBondCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                     const correlation::analysis::AnalysisSettings & /*settings*/) const {
  dists.addHistogram("HBond", calculate(dists.cell(), dists.neighbors()));
}

correlation::analysis::Histogram HBondCalculator::calculate(const correlation::core::Cell &cell,
                                                            const correlation::analysis::StructureAnalyzer *neighbors) {

  if (neighbors == nullptr) {
    return {};
  }

  const auto &atoms = cell.atoms();
  const size_t num_atoms = cell.atomCount();
  const auto &neighbor_graph = neighbors->neighborGraph();

  // Geometric criteria
  const HBondCriteria criteria{
      .R_cut_sq = 3.5 * 3.5, // Donor-Acceptor distance squared
      .Alpha_cut = 30.0      // H-D...A angle
  };

  std::vector<int> hbond_counts(num_atoms, 0);

  // Identify Donors (D) and Acceptors (A) - commonly O, N, F, S.
  // Pre-filter to the electronegative subset so the inner acceptor
  // search iterates N_en atoms rather than N (typically N_en << N).
  // NOTE: for a full O(N_en) solution, StructureAnalyzer would need to
  // expose a non-bonded neighbour list at R_cut = 3.5 Å per atom.
  std::vector<size_t> en_indices;
  en_indices.reserve(num_atoms / 4);
  for (size_t i = 0; i < num_atoms; ++i) {
    if (isElectronegative(atoms[i].element().symbol)) {
      en_indices.push_back(i);
    }
  }

  findHydrogenBonds(cell, neighbor_graph, en_indices, criteria, hbond_counts);

  // Use pre-filtered en_indices for the final distribution loop.
  std::map<int, double> distribution;
  int max_hb = 0;
  int num_en_atoms = static_cast<int>(en_indices.size());
  for (size_t atom_idx : en_indices) {
    distribution[hbond_counts[atom_idx]]++;
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
    hist.bins.push_back(static_cast<double>(idx));
    double freq = (num_en_atoms > 0) ? (distribution[idx] / num_en_atoms) : 0.0;
    hist.partials["Total"].push_back(freq);
  }

  return hist;
}

} // namespace correlation::calculators
