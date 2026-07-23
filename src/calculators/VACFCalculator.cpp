/**
 * @file VACFCalculator.cpp
 * @brief Implementation of velocity autocorrelation calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/VACFCalculator.hpp"
#include "analysis/DynamicsAnalyzer.hpp"
#include "calculators/CalculatorFactory.hpp"

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
const bool registered = CalculatorFactory::registerTypeSafe<VACFCalculator>("VACFCalculator");
} // namespace

void VACFCalculator::calculateTrajectory(correlation::analysis::DistributionFunctions &dists,
                                         const correlation::core::Trajectory &traj,
                                         const correlation::analysis::AnalysisSettings & /*settings*/) const {
  auto results = calculate(traj, {-1}, {0}, {static_cast<size_t>(-1)});
  for (auto &[name, histogram] : results) {
    dists.addHistogram(name, std::move(histogram));
  }
}

std::map<std::string, correlation::analysis::Histogram>
VACFCalculator::calculate(const correlation::core::Trajectory &traj,
                          correlation::analysis::MaxFrames max_correlation_frames,
                          correlation::analysis::StartFrame start_frame, correlation::analysis::EndFrame end_frame) {
  std::map<std::string, correlation::analysis::Histogram> results;
  // Use getFrameCount() to avoid materialising a memory-mapped trajectory.
  if (traj.getFrameCount() == 0) {
    return results;
  }

  std::vector<real_t> const raw_vacf =
      correlation::analysis::DynamicsAnalyzer::calculateVACF(traj, max_correlation_frames, start_frame, end_frame);
  if (raw_vacf.empty()) {
    return results;
  }

  size_t const num_frames = raw_vacf.size();
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

  vacf_hist.partials["Total"] = raw_vacf;
  results["VACF"] = std::move(vacf_hist);

  std::vector<real_t> const norm_vacf = correlation::analysis::DynamicsAnalyzer::calculateNormalizedVACF(
      traj, max_correlation_frames, start_frame, end_frame);
  if (!norm_vacf.empty()) {
    correlation::analysis::Histogram norm_vacf_hist;
    norm_vacf_hist.x_label = "t";
    norm_vacf_hist.title = "Normalized Velocity Autocorrelation";
    norm_vacf_hist.y_label = "C(t) / C(0)";
    norm_vacf_hist.x_unit = "fs";
    norm_vacf_hist.y_unit = "normalized";
    norm_vacf_hist.description = "Normalized Velocity Autocorrelation Function";
    norm_vacf_hist.file_suffix = "_VACF_norm";
    norm_vacf_hist.bins = results["VACF"].bins;
    norm_vacf_hist.partials["Total"] = norm_vacf;
    results["Normalized VACF"] = std::move(norm_vacf_hist);
  }

  return results;
}

} // namespace correlation::calculators
