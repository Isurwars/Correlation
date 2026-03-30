// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/DebyeS_KCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/SIMDUtils.hpp"
#include <cmath>
#include <stdexcept>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <map>
#include <vector>
#include <string>

namespace {
bool registered = CalculatorFactory::instance().registerCalculator(
    std::make_unique<DebyeS_KCalculator>());
} // namespace

void DebyeS_KCalculator::calculateFrame(DistributionFunctions &df,
                                        const AnalysisSettings &settings) const {
  if (settings.q_bin_width <= 0 || settings.q_max <= 0) {
    throw std::invalid_argument("Q-space parameters must be positive.");
  }

  const Cell &cell = df.cell();
  const auto &atoms = cell.atoms();
  const size_t num_atoms = atoms.size();
  
  if (num_atoms == 0) return;

  const double q_max = settings.q_max;
  const double q_bin_width = settings.q_bin_width;
  const size_t num_q_bins = static_cast<size_t>(std::floor(q_max / q_bin_width));

  if (num_q_bins == 0) {
    throw std::invalid_argument("Q_max is too small for the given Q_bin_width.");
  }

  Histogram s_k_hist;
  s_k_hist.bin_label = "Q (Å⁻¹)";
  s_k_hist.bins.resize(num_q_bins);
  for (size_t i = 0; i < num_q_bins; ++i) {
    s_k_hist.bins[i] = (i + 0.5) * q_bin_width;
  }

  // 1. Group atom indices by element symbol
  std::map<std::string, std::vector<size_t>> indices_by_type;
  for (size_t i = 0; i < num_atoms; ++i) {
    indices_by_type[atoms[i].element().symbol].push_back(i);
  }

  // 2. Precompute all pairwise Euclidean distances for each partial pair
  struct PairData {
    std::string key;
    std::vector<double> distances;
    double weight; // Ashcroft-Langreth weight: c_i * c_j (or 2*c_i*c_j)
    double composition_factor; // sqrt(c_i * c_j)
    double N_A;
    double N_B;
    bool is_identical;
  };
  std::vector<PairData> pair_data_list;
  
  const auto& ashcroft_weights = df.getAshcroftWeights();

  for (auto it1 = indices_by_type.begin(); it1 != indices_by_type.end(); ++it1) {
    const auto& sym1 = it1->first;
    const auto& ids1 = it1->second;
    const double N_A = static_cast<double>(ids1.size());
    
    for (auto it2 = it1; it2 != indices_by_type.end(); ++it2) {
      const auto& sym2 = it2->first;
      const auto& ids2 = it2->second;
      const double N_B = static_cast<double>(ids2.size());
      
      bool is_identical = (sym1 == sym2);
      // Construct key consistent with DistributionFunctions::getPartialKey
      std::string key = (sym1 < sym2) ? (sym1 + "-" + sym2) : (sym2 + "-" + sym1);
      
      double weight = 0.0;
      auto wit = ashcroft_weights.find(key);
      if (wit != ashcroft_weights.end()) {
        weight = wit->second;
      }

      // Ashcroft-Langreth composition factor: sqrt(c_i * c_j)
      double comp_factor = is_identical ? std::sqrt(weight) : std::sqrt(weight / 2.0);

      std::vector<double> dists;
      if (is_identical) {
        dists.reserve(ids1.size() * (ids1.size() - 1) / 2);
        for (size_t i = 0; i < ids1.size(); ++i) {
          for (size_t j = i + 1; j < ids1.size(); ++j) {
            double dx = atoms[ids1[i]].position().x() - atoms[ids1[j]].position().x();
            double dy = atoms[ids1[i]].position().y() - atoms[ids1[j]].position().y();
            double dz = atoms[ids1[i]].position().z() - atoms[ids1[j]].position().z();
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (r > 1e-9) dists.push_back(r);
          }
        }
      } else {
        dists.reserve(ids1.size() * ids2.size());
        for (size_t i = 0; i < ids1.size(); ++i) {
          for (size_t j = 0; j < ids2.size(); ++j) {
            double dx = atoms[ids1[i]].position().x() - atoms[ids2[j]].position().x();
            double dy = atoms[ids1[i]].position().y() - atoms[ids2[j]].position().y();
            double dz = atoms[ids1[i]].position().z() - atoms[ids2[j]].position().z();
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (r > 1e-9) dists.push_back(r);
          }
        }
      }
      
      pair_data_list.push_back({key, std::move(dists), weight, comp_factor, N_A, N_B, is_identical});
      s_k_hist.partials[key].assign(num_q_bins, 0.0);
    }
  }

  std::vector<double> total_s_k(num_q_bins, 0.0);
  tbb::enumerable_thread_specific<std::vector<double>> scratch_ets;

  tbb::parallel_for(
    tbb::blocked_range<size_t>(0, num_q_bins),
    [&](const tbb::blocked_range<size_t> &range) {
      auto &scratch = scratch_ets.local();
      
      for (size_t i = range.begin(); i != range.end(); ++i) {
        const double Q = s_k_hist.bins[i];
        
        for (const auto& pair_data : pair_data_list) {
          if (scratch.size() < pair_data.distances.size()) {
            scratch.resize(pair_data.distances.size());
          }
          
          double sum = 0.0;
          if (!pair_data.distances.empty()) {
             sum = correlation::math::debye_sum(Q, pair_data.distances.data(), scratch.data(), pair_data.distances.size());
          }
          
          double partial_sk = 0.0;
          if (pair_data.is_identical) {
            // S_aa(Q) = 1 + (2 / N_a) * sum
            partial_sk = 1.0 + (2.0 / pair_data.N_A) * sum;
          } else {
            // S_ab(Q) = (1 / sqrt(N_a * N_b)) * sum
            partial_sk = (1.0 / std::sqrt(pair_data.N_A * pair_data.N_B)) * sum;
          }
          
          s_k_hist.partials[pair_data.key][i] = partial_sk;
          // Total S(Q) = sum_A c_A * S_AA(Q) + sum_{A<B} 2 * sqrt(c_A*c_B) * S_AB(Q)
          double weight_factor = pair_data.is_identical ? 1.0 : 2.0;
          total_s_k[i] += partial_sk * pair_data.composition_factor * weight_factor;
        }
      }
    });

  s_k_hist.partials["Total"] = std::move(total_s_k);
  df.addHistogram("debye_S_q", std::move(s_k_hist));
}
