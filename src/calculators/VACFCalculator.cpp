/**
 * @file VACFCalculator.cpp
 * @brief Implementation of velocity autocorrelation calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "calculators/VACFCalculator.hpp"
#include "DynamicsAnalyzer.hpp"
#include "calculators/CalculatorFactory.hpp"

namespace {
bool registered = CalculatorFactory::instance().registerCalculator(
    std::make_unique<VACFCalculator>());
} // namespace

void VACFCalculator::calculateTrajectory(
    DistributionFunctions &df, const Trajectory &traj,
    const AnalysisSettings &settings) const {
  auto results = calculate(traj, -1, 0, static_cast<size_t>(-1));
  for (auto &[name, histogram] : results) {
    df.addHistogram(name, std::move(histogram));
  }
}

std::map<std::string, Histogram>
VACFCalculator::calculate(const Trajectory &traj, int max_correlation_frames,
                          size_t start_frame, size_t end_frame) {
  std::map<std::string, Histogram> results;
  const auto &velocities = traj.getVelocities();
  if (velocities.empty()) {
    return results;
  }

  std::vector<double> raw_vacf = DynamicsAnalyzer::calculateVACF(
      traj, max_correlation_frames, start_frame, end_frame);
  if (raw_vacf.empty())
    return results;

  size_t num_frames = raw_vacf.size();
  double dt = traj.getTimeStep();

  Histogram vacf_hist;
  vacf_hist.x_label = "t";
  vacf_hist.title = "Velocity Autocorrelation";
  vacf_hist.y_label = "C(t)";
  vacf_hist.x_unit = "fs";
  vacf_hist.y_unit = "Å²/fs²";
  vacf_hist.description = "Velocity Autocorrelation Function";
  vacf_hist.file_suffix = "_VACF";
  vacf_hist.bins.resize(num_frames);
  for (size_t i = 0; i < num_frames; ++i) {
    vacf_hist.bins[i] = i * dt;
  }

  vacf_hist.partials["Total"] = raw_vacf;
  results["VACF"] = std::move(vacf_hist);

  std::vector<double> norm_vacf = DynamicsAnalyzer::calculateNormalizedVACF(
      traj, max_correlation_frames, start_frame, end_frame);
  if (!norm_vacf.empty()) {
    Histogram norm_vacf_hist;
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
