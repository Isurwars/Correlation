// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/SQCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "PhysicalData.hpp"
#include "SIMDUtils.hpp"
#include <cmath>
#include <stdexcept>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace {
bool registered = CalculatorFactory::instance().registerCalculator(
    std::make_unique<SQCalculator>());
} // namespace

void SQCalculator::calculateFrame(DistributionFunctions &df,
                                  const AnalysisSettings &settings) const {
  if (df.getAllHistograms().find("g(r)") == df.getAllHistograms().end()) {
    return; // g(r) hasn't been calculated yet
  }
  df.addHistogram("S(Q)",
                  calculate(df.getHistogram("g(r)"), df.cell(),
                            df.getAshcroftWeights(), settings.q_max,
                            settings.q_bin_width, settings.r_int_max));
}

Histogram
SQCalculator::calculate(const Histogram &g_r_hist, const Cell &cell,
                        const std::map<std::string, double> &ashcroft_weights,
                        double q_max, double q_bin_width,
                        double r_integration_max) {
  if (q_bin_width <= 0 || q_max <= 0) {
    throw std::invalid_argument("Q-space parameters must be positive.");
  }
  if (r_integration_max <= 0) {
    throw std::invalid_argument("Integration cutoff must be positive.");
  }

  const auto &r_bins = g_r_hist.bins;
  if (r_bins.size() < 2) {
    throw std::logic_error("Insufficient r-bins for integration.");
  }

  const double max_r_in_hist = r_bins.back();
  const double r_max_val = std::min(r_integration_max, max_r_in_hist);

  const double dr = r_bins[1] - r_bins[0];
  const double total_rho = cell.atomCount() / cell.volume();

  Histogram s_q_hist;
  s_q_hist.bin_label = "Q (Å⁻¹)";
  const size_t num_q_bins =
      static_cast<size_t>(std::floor(q_max / q_bin_width));

  if (num_q_bins == 0) {
    throw std::invalid_argument(
        "Q_max is too small for the given Q_bin_width.");
  }

  s_q_hist.bins.resize(num_q_bins);
  for (size_t i = 0; i < num_q_bins; ++i) {
    s_q_hist.bins[i] = (i + 0.5) * q_bin_width;
  }

  size_t j_max = r_bins.size();
  for (size_t j = 0; j < r_bins.size(); ++j) {
    if (r_bins[j] > r_max_val) {
      j_max = j;
      break;
    }
  }

  std::vector<double> total_s_q_num(num_q_bins, 0.0);

  // Collect (key, partial_ref, weight, composition_factor, is_identical)
  // tuples once so that parallel Q-bins can read them without touching the map.
  struct PartialInfo {
    const std::string *key;
    const std::vector<double> *g_r_partial;
    double weight;
    double composition_sqrt_factor;
    bool is_identical;
  };
  std::vector<PartialInfo> partials;
  partials.reserve(g_r_hist.partials.size());

  for (const auto &[key, g_r_partial] : g_r_hist.partials) {
    if (key == "Total")
      continue;
    double weight = ashcroft_weights.at(key);
    bool is_identical = false;
    size_t dash_pos = key.find('-');
    if (dash_pos != std::string::npos) {
      std::string sym1 = key.substr(0, dash_pos);
      std::string sym2 = key.substr(dash_pos + 1);
      is_identical = (sym1 == sym2);
    }
    double csf = is_identical ? std::sqrt(weight) : std::sqrt(weight / 2.0);
    partials.push_back({&key, &g_r_partial, weight, csf, is_identical});

    // Allocate output partial now (serial, before parallel section)
    s_q_hist.partials[key].assign(num_q_bins, 0.0);
  }

  // Precompute integrand_term[p][j] = r*(g-1)*window*dr  (read-only in parallel)
  const size_t num_partials = partials.size();
  std::vector<std::vector<double>> integrand_terms(num_partials,
                                                   std::vector<double>(j_max));
  for (size_t p = 0; p < num_partials; ++p) {
    const auto &g_r_partial = *partials[p].g_r_partial;
    for (size_t j = 0; j < j_max; ++j) {
      const double r = r_bins[j];
      double window = 1.0;
      if (r > r_max_val * 0.8) {
        const double x = (r - 0.8 * r_max_val) / (0.2 * r_max_val);
        window = std::cos(constants::pi * x / 2.0);
      }
      integrand_terms[p][j] = r * (g_r_partial[j] - 1.0) * window * dr;
    }
  }

  // Thread-local scratch buffer for sin(Q*r) precompute (Phase 1 of kernel)
  tbb::enumerable_thread_specific<std::vector<double>> sinqr_ets(
      [&] { return std::vector<double>(j_max, 0.0); });

  // Parallel over Q-bins: each bin is fully independent
  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, num_q_bins),
      [&](const tbb::blocked_range<size_t> &range) {
        auto &sinqr = sinqr_ets.local(); // per-thread scratch

        for (size_t i = range.begin(); i != range.end(); ++i) {
          const double Q = s_q_hist.bins[i];

          for (size_t p = 0; p < num_partials; ++p) {
            const PartialInfo &pi = partials[p];
            const double integral = simd_utils::sinc_integral(
                Q, integrand_terms[p].data(), r_bins.data(),
                sinqr.data(), j_max);

            const double delta_ij = pi.is_identical ? 1.0 : 0.0;
            const double sq_val =
                delta_ij + (4.0 * constants::pi * total_rho *
                            pi.composition_sqrt_factor / Q) *
                               integral;

            s_q_hist.partials[*pi.key][i] = sq_val;
            total_s_q_num[i] += sq_val * pi.weight;
          }
        }
      });

  s_q_hist.partials["Total"] = std::move(total_s_q_num);

  return s_q_hist;
}
