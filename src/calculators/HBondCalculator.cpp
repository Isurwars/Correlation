/**
 * @file HBondCalculator.cpp
 * @brief Implementation of the Hydrogen Bond Analysis calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "calculators/HBondCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Constants.hpp"
#include <vector>

namespace correlation::calculators {

namespace {
bool registered = CalculatorFactory::instance().registerCalculator(
    std::make_unique<HBondCalculator>());
}

void HBondCalculator::calculateFrame(
    correlation::analysis::DistributionFunctions &df,
    const correlation::analysis::AnalysisSettings &settings) const {
  df.addHistogram("HBond", calculate(df.cell(), df.neighbors()));
}

correlation::analysis::Histogram HBondCalculator::calculate(
    const correlation::core::Cell &cell,
    const correlation::analysis::StructureAnalyzer *neighbors) {

  if (!neighbors)
    return {};

  const auto &atoms = cell.atoms();
  const size_t num_atoms = cell.atomCount();
  const auto &neighbor_graph = neighbors->neighborGraph();

  // Geometric criteria
  const double R_cut = 3.5; // Donor-Acceptor distance
  const double R_cut_sq = R_cut * R_cut;
  const double Alpha_cut = 30.0; // H-D...A angle

  size_t total_hbonds = 0;
  std::vector<int> hbond_counts(num_atoms, 0);

  // Identify Donors (D) and Acceptors (A) - commonly O, N, F, S
  auto is_electronegative = [](const std::string &symbol) {
    return symbol == "O" || symbol == "N" || symbol == "F" || symbol == "S";
  };

  for (size_t i = 0; i < num_atoms; ++i) {
    const auto &atom_i = atoms[i];
    if (!is_electronegative(atom_i.element().symbol))
      continue;

    // Atom i is a potential Donor or Acceptor.
    // Let's find all hydrogens bonded to i (i is Donor).
    std::vector<size_t> hydrogens;
    for (const auto &nb : neighbor_graph.getNeighbors(i)) {
      if (atoms[nb.index].element().symbol == "H") {
        hydrogens.push_back(nb.index);
      }
    }

    if (hydrogens.empty())
      continue;

    // For each Hydrogen, find potential Acceptors (j)
    for (size_t h_idx : hydrogens) {
      const auto &pos_h = atoms[h_idx].position();
      const auto &pos_d = atom_i.position();
      correlation::math::Vector3<double> v_dh =
          cell.minimumImage(pos_h - pos_d);

      // Guard against degenerate case where H and D are coincident
      if (correlation::math::norm_sq(v_dh) < 1e-12)
        continue;

      for (size_t j = 0; j < num_atoms; ++j) {
        if (i == j)
          continue;
        const auto &atom_j = atoms[j];
        if (!is_electronegative(atom_j.element().symbol))
          continue;

        // Atom j is a potential Acceptor.
        correlation::math::Vector3<double> v_da =
            cell.minimumImage(atom_j.position() - pos_d);
        double d_da_sq = correlation::math::norm_sq(v_da);

        if (d_da_sq < R_cut_sq) {
          // Check angle H-D...A
          double dot_val = v_dh * v_da;
          double cos_alpha = dot_val / (correlation::math::norm(v_dh) * std::sqrt(d_da_sq));
          double alpha = std::acos(std::max(-1.0, std::min(1.0, cos_alpha))) *
                         correlation::math::rad_to_deg;

          if (alpha < Alpha_cut) {
            total_hbonds++;
            hbond_counts[i]++;
            hbond_counts[j]++;
          }
        }
      }
    }
  }

  // Results: Distribution of H-bond counts per molecule (if we could identify
  // them) For now, just a distribution of H-bonds per electronegative atom.
  std::map<int, double> distribution;
  int max_hb = 0;
  int num_en_atoms = 0;
  for (size_t i = 0; i < num_atoms; ++i) {
    if (is_electronegative(atoms[i].element().symbol)) {
      distribution[hbond_counts[i]]++;
      num_en_atoms++;
      if (hbond_counts[i] > max_hb)
        max_hb = hbond_counts[i];
    }
  }

  correlation::analysis::Histogram hist;
  hist.x_label = "N_HB";
  hist.title = "Hydrogen Bond Distribution";
  hist.y_label = "P(N_HB)";
  hist.x_unit = "bonds";
  hist.y_unit = "fraction";
  hist.description = "Distribution of Hydrogen Bonds per electronegative atom.";
  hist.file_suffix = "_HBond";

  for (int n = 0; n <= max_hb; ++n) {
    hist.bins.push_back(static_cast<double>(n));
    double freq = (num_en_atoms > 0) ? (distribution[n] / num_en_atoms) : 0.0;
    hist.partials["Total"].push_back(freq);
  }

  return hist;
}

} // namespace correlation::calculators
