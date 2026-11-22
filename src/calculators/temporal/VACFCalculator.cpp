/**
 * @file VACFCalculator.cpp
 * @brief Implementation of velocity autocorrelation calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/temporal/VACFCalculator.hpp"
#include "analysis/DynamicsAnalyzer.hpp"
#include "calculators/CalculatorFactory.hpp"

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
const bool REGISTERED = CalculatorFactory::registerTypeSafe<VACFCalculator>("VACFCalculator");
} // namespace

void VACFCalculator::calculateTrajectory(
    correlation::analysis::DistributionFunctions &dists, const correlation::core::Trajectory &traj,
    const correlation::analysis::AnalysisSettings &settings) const {
  if (traj.getFrameCount() < 2) {
    return;
  }
  auto results = calculate(traj, {settings.vacf_max_lag}, {0}, {static_cast<size_t>(-1)});
  for (auto &[name, histogram] : results) {
    dists.addHistogram(name, std::move(histogram));
  }
}

std::map<std::string, correlation::analysis::Histogram>
VACFCalculator::calculate(const correlation::core::Trajectory &traj,
                          correlation::analysis::MaxFrames max_correlation_frames,
                          correlation::analysis::StartFrame start_frame,
                          correlation::analysis::EndFrame end_frame) {
  std::map<std::string, correlation::analysis::Histogram> results;
  if (traj.getFrameCount() < 2) {
    return results;
  }

  const auto partial_vacf = correlation::analysis::DynamicsAnalyzer::calculatePartialVACF(
      traj, max_correlation_frames, start_frame, end_frame);
  if (partial_vacf.total_unweighted.size() < 2) {
    return results;
  }

  size_t const num_frames = partial_vacf.total_unweighted.size();
  real_t const time_step = traj.getTimeStep();

  correlation::analysis::Histogram vacf_hist;
  vacf_hist.x_label = "t";
  vacf_hist.title = "Velocity Autocorrelation";
  vacf_hist.y_label = "C(t)";
  vacf_hist.x_unit = "fs";
  vacf_hist.y_unit = "Å² fs⁻²";
  vacf_hist.description = "Velocity Autocorrelation Function";
  vacf_hist.file_suffix = "_VACF";
  vacf_hist.bins.resize(num_frames);
  for (size_t i = 0; i < num_frames; ++i) {
    vacf_hist.bins[i] = static_cast<real_t>(i) * time_step;
  }

  vacf_hist.partials["Total"] = partial_vacf.total_unweighted;
  for (const auto &[elem_sym, elem_curve] : partial_vacf.element_vacfs) {
    vacf_hist.partials[elem_sym] = elem_curve;
  }

  correlation::analysis::Histogram norm_vacf_hist;
  norm_vacf_hist.x_label = "t";
  norm_vacf_hist.title = "Normalized Velocity Autocorrelation";
  norm_vacf_hist.y_label = "C(t) / C(0)";
  norm_vacf_hist.x_unit = "fs";
  norm_vacf_hist.y_unit = "normalized";
  norm_vacf_hist.description = "Normalized Velocity Autocorrelation Function";
  norm_vacf_hist.file_suffix = "_VACF_norm";
  norm_vacf_hist.bins = vacf_hist.bins;

  const auto normalize_curve = [](const std::vector<real_t> &curve) {
    std::vector<real_t> norm = curve;
    if (!norm.empty() && norm[0] != static_cast<real_t>(0.0)) {
      const real_t inv_c0 = static_cast<real_t>(1.0) / norm[0];
      for (real_t &val : norm) {
        val *= inv_c0;
      }
    }
    return norm;
  };

  norm_vacf_hist.partials["Total"] = normalize_curve(partial_vacf.total_unweighted);
  for (const auto &[elem_sym, elem_curve] : partial_vacf.element_vacfs) {
    norm_vacf_hist.partials[elem_sym] = normalize_curve(elem_curve);
  }

  results["VACF"] = std::move(vacf_hist);
  results["vacf_normalized"] = std::move(norm_vacf_hist);
  return results;
}

} // namespace correlation::calculators
