// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/XRDCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Constants.hpp"
#include "math/PhysicalData.hpp"
#include "math/SIMDUtils.hpp"
#include <cmath>
#include <stdexcept>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace {
bool registered = CalculatorFactory::instance().registerCalculator(
    std::make_unique<XRDCalculator>());
} // namespace

void XRDCalculator::calculateFrame(DistributionFunctions &df,
                                   const AnalysisSettings &settings) const {
  if (df.getAllHistograms().find("g(r)") == df.getAllHistograms().end()) {
    return; // g(r) hasn't been calculated yet
  }
  df.addHistogram("XRD", calculate(df.getHistogram("g(r)"), df.cell(),
                                   df.getAshcroftWeights(), 1.5406, 5.0, 90.0,
                                   settings.q_bin_width));
}

Histogram
XRDCalculator::calculate(const Histogram &g_r_hist, const Cell &cell,
                         const std::map<std::string, double> &ashcroft_weights,
                         double lambda, double theta_min, double theta_max,
                         double bin_width) {
  if (bin_width <= 0) {
    throw std::invalid_argument("Bin width must be positive.");
  }

  const auto &r_bins = g_r_hist.bins;
  const double dr = r_bins[1] - r_bins[0];
  const double total_rho = cell.atomCount() / cell.volume();
  const double max_r = r_bins.back();

  size_t num_bins =
      static_cast<size_t>((theta_max - theta_min) / bin_width) + 1;
  Histogram xrd_hist;
  xrd_hist.bin_label = "2Theta (°)";
  xrd_hist.bins.resize(num_bins);
  std::vector<double> intensities(num_bins, 0.0);

  auto get_f_Q = [](const std::string &symbol, double Q) -> double {
    const auto &coeffs =
        correlation::math::physics::getAtomicFormFactors(symbol);
    double s = Q / (4.0 * correlation::math::constants::pi);
    double s2 = s * s;
    double f = coeffs[8];
    for (size_t i = 0; i < 4; ++i) {
      f += coeffs[2 * i] * std::exp(-coeffs[2 * i + 1] * s2);
    }
    return f;
  };

  std::map<std::string, std::vector<double>> partial_integrands;
  for (const auto &[key, g_partial] : g_r_hist.partials) {
    if (key == "Total")
      continue;
    partial_integrands[key].resize(g_partial.size());
    for (size_t k = 0; k < g_partial.size(); ++k) {
      partial_integrands[key][k] = r_bins[k] * (g_partial[k] - 1.0) * dr;
    }
  }

  std::map<std::string, double> concentrations;
  for (const auto &elem : cell.elements()) {
    double count = 0;
    for (const auto &atom : cell.atoms()) {
      if (atom.element().symbol == elem.symbol)
        count++;
    }
    concentrations[elem.symbol] = count / cell.atomCount();
  }

  // Collect partials as a flat vector for parallel access (avoids map iteration
  // inside the parallel region)
  struct PartialXRD {
    const std::string *key;
    const std::vector<double> *integrand;
    double weight;
    std::string sym1, sym2;
  };
  std::vector<PartialXRD> xrd_partials;
  xrd_partials.reserve(partial_integrands.size());
  for (const auto &[key, integ] : partial_integrands) {
    double weight = ashcroft_weights.at(key);
    size_t dash_pos = key.find('-');
    std::string s1 = key.substr(0, dash_pos);
    std::string s2 = key.substr(dash_pos + 1);
    xrd_partials.push_back({&key, &integ, weight, s1, s2});
  }
  const size_t num_xrd_partials = xrd_partials.size();

  // Thread-local sin(Q*r) scratch buffer
  const size_t r_count = r_bins.size();
  tbb::enumerable_thread_specific<std::vector<double>> sinqr_ets(
      [&] { return std::vector<double>(r_count, 0.0); });

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, num_bins),
      [&](const tbb::blocked_range<size_t> &range) {
        auto &sinqr = sinqr_ets.local();

        for (size_t i = range.begin(); i != range.end(); ++i) {
          double two_theta = theta_min + i * bin_width;
          xrd_hist.bins[i] = two_theta;

          double theta_rad =
              (two_theta / 2.0) * correlation::math::constants::deg2rad;
          double Q = 4.0 * correlation::math::constants::pi *
                     std::sin(theta_rad) / lambda;

          if (Q < 1e-6) {
            intensities[i] = 0.0;
            continue;
          }

          double I_Q = 0.0;

          for (const auto &[sym, c] : concentrations) {
            double f = get_f_Q(sym, Q);
            I_Q += c * f * f;
          }

          for (size_t p = 0; p < num_xrd_partials; ++p) {
            const PartialXRD &px = xrd_partials[p];
            // Clamp integrand to r_bins size (partial_integrands was built
            // from g_partial which may be shorter than r_bins)
            const size_t pcount = std::min(px.integrand->size(), r_count);

            double integral = correlation::math::simd::sinc_integral(
                Q, px.integrand->data(), r_bins.data(), sinqr.data(), pcount);

            double f1 = get_f_Q(px.sym1, Q);
            double f2 = get_f_Q(px.sym2, Q);

            I_Q += px.weight * f1 * f2 *
                   (4.0 * correlation::math::constants::pi * total_rho / Q) *
                   integral;
          }

          intensities[i] = I_Q;
        }
      });

  xrd_hist.partials["Total"] = std::move(intensities);

  return xrd_hist;
}
