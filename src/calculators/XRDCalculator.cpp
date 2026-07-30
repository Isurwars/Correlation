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

namespace correlation {
namespace calculators {

namespace {
// Static registration of the calculator in the factory
const bool registered = CalculatorFactory::registerTypeSafe<XRDCalculator>("XRDCalculator");

struct PartialInfoSq {
  std::string key;
  const std::vector<real_t> *data;
  std::string sym1;
  std::string sym2;
  bool is_identical;
  real_t c_i;
  real_t c_j;
};

struct QGrid {
  const std::vector<real_t> *bins{nullptr};
  real_t min{0.0};
  real_t max{0.0};
};

// Helper for linear interpolation of partial S_ij(Q)
real_t sampleSq(const QGrid &q_grid, const std::vector<real_t> &sq_data, real_t q_val) {
  if (q_val <= q_grid.min) {
    return sq_data.front();
  }
  if (q_val >= q_grid.max) {
    return sq_data.back();
  }
  const auto &bins = *q_grid.bins;
  auto iterator = std::lower_bound(bins.begin(), bins.end(), q_val);
  size_t idx = std::distance(bins.begin(), iterator);
  if (idx == 0) {
    return sq_data[0];
  }
  real_t q_0 = bins[idx - 1];
  real_t q_1 = bins[idx];
  real_t d_t = (q_val - q_0) / (q_1 - q_0);
  return sq_data[idx - 1] + d_t * (sq_data[idx] - sq_data[idx - 1]);
}

std::vector<PartialInfoSq> buildPartialSqList(const std::map<std::string, std::vector<real_t>> &partials,
                                              const std::map<std::string, real_t> &concentrations) {
  std::vector<PartialInfoSq> partial_sq_list;
  for (const auto &[key, sq_data] : partials) {
    if (key == "Total") {
      continue;
    }
    size_t dash_pos = key.find('-');
    std::string sym1 = key.substr(0, dash_pos);
    std::string sym2 = key.substr(dash_pos + 1);
    real_t c_1 = concentrations.contains(sym1) ? concentrations.at(sym1) : static_cast<real_t>(0.0);
    real_t c_2 = concentrations.contains(sym2) ? concentrations.at(sym2) : static_cast<real_t>(0.0);
    partial_sq_list.push_back({
        .key = key,
        .data = &sq_data,
        .sym1 = sym1,
        .sym2 = sym2,
        .is_identical = (sym1 == sym2),
        .c_i = c_1,
        .c_j = c_2,
    });
  }
  return partial_sq_list;
}

real_t calculateIntensityAtQ(real_t q_value, const std::map<std::string, real_t> &concentrations,
                             const std::vector<PartialInfoSq> &partial_sq_list, const QGrid &q_grid) {
  if (q_value < 1e-6) {
    return 0.0;
  }

  correlation::KahanAccumulator<real_t> intensity_Q;

  // 1. Self-scattering term: sum_i c_i * f_i(Q)^2
  for (const auto &[sym, concentration] : concentrations) {
    auto const form_factor = XRDCalculator::getAtomicFormFactor(sym, q_value);
    intensity_Q.add(static_cast<real_t>(concentration) * form_factor * form_factor);
  }

  // 2. Inter-atomic interference term from partial S_ij(Q):
  // sum_{i <= j} (2 - delta_ij) * sqrt(c_i * c_j) * f_i(Q) * f_j(Q) * (S_ij(Q) - delta_ij)
  for (const auto &p_sq : partial_sq_list) {
    real_t sq_val = sampleSq(q_grid, *p_sq.data, q_value);
    real_t f_1 = XRDCalculator::getAtomicFormFactor(p_sq.sym1, q_value);
    real_t f_2 = XRDCalculator::getAtomicFormFactor(p_sq.sym2, q_value);
    real_t factor = p_sq.is_identical ? static_cast<real_t>(1.0) : static_cast<real_t>(2.0);
    real_t delta_ij = p_sq.is_identical ? static_cast<real_t>(1.0) : static_cast<real_t>(0.0);
    real_t weight = factor * std::sqrt(p_sq.c_i * p_sq.c_j);

    intensity_Q.add(weight * f_1 * f_2 * (sq_val - delta_ij));
  }

  return intensity_Q.value();
}
} // namespace

void XRDCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                   const correlation::analysis::AnalysisSettings &settings) const {
  // Primary path: Use S(Q) if available
  if (dists.getAllHistograms().contains("S_q")) {
    dists.addHistogram("XRD", calculateFromSq(dists.getHistogram("S_q"), dists.cell(), dists.getAshcroftWeights(),
                                              Wavelength{1.5406}, MinTheta{10.0}, MaxTheta{140.0}, BinWidth{0.05}));
    return;
  }

  // Fallback path: Use g(r)
  if (!dists.getAllHistograms().contains("g_r")) {
    const auto *rdf_calc = CalculatorFactory::instance().getCalculator("RDF");
    if (rdf_calc != nullptr) {
      rdf_calc->calculateFrame(dists, settings);
    }
  }
  if (dists.getAllHistograms().contains("g_r")) {
    dists.addHistogram("XRD", calculate(dists.getHistogram("g_r"), dists.cell(), dists.getAshcroftWeights(),
                                        Wavelength{1.5406}, MinTheta{10.0}, MaxTheta{140.0}, BinWidth{0.05}));
  }
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

correlation::analysis::Histogram
XRDCalculator::calculateFromSq(const correlation::analysis::Histogram &s_q_hist, const correlation::core::Cell &cell,
                               const std::map<std::string, real_t> & /*ashcroft_weights*/, Wavelength lambda,
                               MinTheta theta_min, MaxTheta theta_max, BinWidth bin_width) {
  real_t const lambda_val = lambda.value;
  real_t const theta_min_val = theta_min.value;
  real_t const theta_max_val = theta_max.value;
  real_t const bin_width_val = bin_width.value;

  if (bin_width_val <= 0) {
    throw std::invalid_argument("Bin width must be positive.");
  }
  if (s_q_hist.bins.empty()) {
    throw std::invalid_argument("S(Q) histogram bins must not be empty.");
  }
  if (lambda_val <= 0.0) {
    throw std::invalid_argument("Wavelength lambda must be strictly positive.");
  }

  size_t num_bins = static_cast<size_t>((theta_max_val - theta_min_val) / bin_width_val) + 1;
  correlation::analysis::Histogram xrd_hist;
  xrd_hist.x_label = "2θ";
  xrd_hist.title = "XRD Pattern";
  xrd_hist.y_label = "Intensity";
  xrd_hist.x_unit = "°";
  xrd_hist.y_unit = "Intensity";
  xrd_hist.description = "X-Ray Diffraction Pattern (from S(Q))";
  xrd_hist.file_suffix = "_XRD";
  xrd_hist.bins.resize(num_bins);
  std::vector<real_t> intensities(num_bins, 0.0);

  auto concentrations = calculateConcentrations(cell);
  const auto &q_bins = s_q_hist.bins;
  real_t const q_min = q_bins.front();
  real_t const q_max = q_bins.back();

  QGrid q_grid{.bins = &q_bins, .min = q_min, .max = q_max};
  auto partial_sq_list = buildPartialSqList(s_q_hist.partials, concentrations);

  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_bins), [&](const tbb::blocked_range<size_t> &range) {
    for (size_t i = range.begin(); i != range.end(); ++i) {
      auto const two_theta = theta_min_val + static_cast<real_t>(i) * bin_width_val;
      xrd_hist.bins[i] = two_theta;

      auto const theta_rad = static_cast<real_t>((two_theta / 2.0) * correlation::math::deg_to_rad);
      auto const q_value = static_cast<real_t>(correlation::math::four_pi * std::sin(theta_rad) / lambda_val);

      intensities[i] = calculateIntensityAtQ(q_value, concentrations, partial_sq_list, q_grid);
    }
  });

  xrd_hist.partials["Total"] = std::move(intensities);
  return xrd_hist;
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

} // namespace calculators
} // namespace correlation
