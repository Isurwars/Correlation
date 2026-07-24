/**
 * @file XRDCalculator.cpp
 * @brief Implementation of X-ray diffraction pattern calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/XRDCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Constants.hpp"
#include "math/Precision.hpp"
#include "math/SIMDUtils.hpp"
#include "physics/PhysicalData.hpp"

#include <cmath>
#include <stdexcept>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
const bool registered = CalculatorFactory::registerTypeSafe<XRDCalculator>("XRDCalculator");
} // namespace

void XRDCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                   const correlation::analysis::AnalysisSettings &settings) const {
  if (!dists.getAllHistograms().contains("g_r")) {
    return; // g(r) hasn't been calculated yet
  }
  dists.addHistogram("XRD",
                     calculate(dists.getHistogram("g_r"), dists.cell(), dists.getAshcroftWeights(), Wavelength{1.5406},
                               MinTheta{5.0}, MaxTheta{90.0}, BinWidth{settings.q_bin_width}));
}

correlation::analysis::Histogram XRDCalculator::calculate(const correlation::analysis::Histogram &g_r_hist,
                                                          const correlation::core::Cell &cell,
                                                          const std::map<std::string, real_t> &ashcroft_weights,
                                                          Wavelength lambda, MinTheta theta_min, MaxTheta theta_max,
                                                          BinWidth bin_width) {
  real_t const lambda_val = lambda.value;
  real_t const theta_min_val = theta_min.value;
  real_t const theta_max_val = theta_max.value;
  real_t const bin_width_val = bin_width.value;

  if (bin_width_val <= 0) {
    throw std::invalid_argument("Bin width must be positive.");
  }
  if (g_r_hist.bins.size() < 2) {
    throw std::invalid_argument("g(r) histogram must contain at least 2 bins.");
  }
  if (lambda_val <= 0.0) {
    throw std::invalid_argument("Wavelength lambda must be strictly positive.");
  }

  const auto &r_bins = g_r_hist.bins;
  const real_t delta_r = r_bins[1] - r_bins[0];
  const real_t total_rho = static_cast<real_t>(cell.atomCount()) / cell.volume();
  const real_t max_r = r_bins.back();

  size_t num_bins = static_cast<size_t>((theta_max_val - theta_min_val) / bin_width_val) + 1;
  correlation::analysis::Histogram xrd_hist;
  xrd_hist.x_label = "2θ";
  xrd_hist.title = "XRD Pattern";
  xrd_hist.y_label = "Intensity";
  xrd_hist.x_unit = "°";
  xrd_hist.y_unit = "Intensity";
  xrd_hist.description = "X-Ray Diffraction Pattern";
  xrd_hist.file_suffix = "_XRD";
  xrd_hist.bins.resize(num_bins);
  std::vector<real_t> intensities(num_bins, 0.0);

  auto partial_integrands = calculatePartialIntegrands(g_r_hist, ashcroft_weights, delta_r);
  auto concentrations = calculateConcentrations(cell);

  // Collect partials as a flat vector for parallel access (avoids map iteration
  // inside the parallel region)
  struct PartialXRD {
    const std::string *key;
    const std::vector<real_t> *integrand;
    real_t weight;
    std::string sym1, sym2;
  };
  std::vector<PartialXRD> xrd_partials;
  xrd_partials.reserve(partial_integrands.size());
  for (const auto &[key, integ] : partial_integrands) {
    auto weight = static_cast<real_t>(ashcroft_weights.at(key));
    size_t dash_pos = key.find('-');
    std::string sym1 = key.substr(0, dash_pos);
    std::string sym2 = key.substr(dash_pos + 1);
    xrd_partials.push_back({
        .key = &key,
        .integrand = &integ,
        .weight = weight,
        .sym1 = sym1,
        .sym2 = sym2,
    });
  }
  const size_t num_xrd_partials = xrd_partials.size();

  // Thread-local sin(Q*r) scratch buffer
  const size_t r_count = r_bins.size();
  tbb::enumerable_thread_specific<std::vector<real_t>> sinqr_ets([&] { return std::vector<real_t>(r_count, 0.0); });

  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_bins), [&](const tbb::blocked_range<size_t> &range) {
    auto &sinqr = sinqr_ets.local();

    for (size_t i = range.begin(); i != range.end(); ++i) {
      auto const two_theta = theta_min_val + static_cast<real_t>(i) * bin_width_val;
      xrd_hist.bins[i] = two_theta;

      auto const theta_rad = static_cast<real_t>((two_theta / 2.0) * correlation::math::deg_to_rad);
      auto const q_value = static_cast<real_t>(correlation::math::four_pi * std::sin(theta_rad) / lambda_val);

      if (q_value < 1e-6) {
        intensities[i] = 0.0;
        continue;
      }

      correlation::KahanAccumulator<real_t> intensity_Q;

      for (const auto &[sym, concentration] : concentrations) {
        auto const form_factor = getAtomicFormFactor(sym, q_value);
        intensity_Q.add(static_cast<real_t>(concentration) * form_factor * form_factor);
      }

      // Precompute sinqr once per theta step (angle bin)
      for (size_t j = 0; j < r_count; ++j) {
        sinqr[j] = std::sin(q_value * r_bins[j]);
      }

      for (size_t partial_idx = 0; partial_idx < num_xrd_partials; ++partial_idx) {
        const PartialXRD &partial_xrd = xrd_partials[partial_idx];
        const size_t pcount = std::min(partial_xrd.integrand->size(), r_count);
        auto const integral =
            static_cast<real_t>(correlation::math::simd_dot(partial_xrd.integrand->data(), sinqr.data(), pcount));
        real_t const form_factor_1 = getAtomicFormFactor(partial_xrd.sym1, q_value);
        real_t const form_factor_2 = getAtomicFormFactor(partial_xrd.sym2, q_value);

        intensity_Q.add(static_cast<real_t>(form_factor_1 * form_factor_2 *
                                            (correlation::math::four_pi * total_rho / q_value) * integral));
      }

      intensities[i] = intensity_Q.value();
    }
  });

  xrd_hist.partials["Total"] = std::move(intensities);

  return xrd_hist;
}

real_t XRDCalculator::getAtomicFormFactor(const std::string &symbol, real_t q_value) {
  const auto &coeffs = correlation::physics::getAtomicFormFactors(symbol);
  auto s_value = static_cast<real_t>(q_value / correlation::math::four_pi);
  real_t s_squared = s_value * s_value;
  auto form_factor = static_cast<real_t>(coeffs.at(8));
  for (size_t i = 0; i < 4; ++i) {
    form_factor += static_cast<real_t>(coeffs.at(2 * i) * std::exp(-coeffs.at(2 * i + 1) * s_squared));
  }
  return form_factor;
}

std::map<std::string, real_t> XRDCalculator::calculateConcentrations(const correlation::core::Cell &cell) {
  std::map<std::string, real_t> concentrations;
  real_t const total_atoms = static_cast<real_t>(cell.atomCount());
  if (total_atoms == 0) {
    return concentrations;
  }
  for (const auto &elem : cell.elements()) {
    real_t count = 0;
    for (const auto &atom : cell.atoms()) {
      if (atom.element().symbol == elem.symbol) {
        count++;
      }
    }
    concentrations[elem.symbol] = count / total_atoms;
  }
  return concentrations;
}

std::map<std::string, std::vector<real_t>>
XRDCalculator::calculatePartialIntegrands(const correlation::analysis::Histogram &g_r_hist,
                                          const std::map<std::string, real_t> &ashcroft_weights, real_t delta_r) {
  std::map<std::string, std::vector<real_t>> partial_integrands;
  const auto &r_bins = g_r_hist.bins;
  for (const auto &[key, g_partial] : g_r_hist.partials) {
    if (key == "Total") {
      continue;
    }
    real_t const weight = ashcroft_weights.at(key);
    partial_integrands[key].resize(g_partial.size());
    for (size_t k = 0; k < g_partial.size(); ++k) {
      partial_integrands[key][k] = r_bins[k] * (g_partial[k] - weight) * delta_r;
    }
  }
  return partial_integrands;
}

} // namespace correlation::calculators
