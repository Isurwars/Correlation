/**
 * @file MSDCalculator.cpp
 * @brief Implementation of Mean Squared Displacement calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/temporal/MSDCalculator.hpp"
#include "analysis/DynamicsAnalyzer.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Precision.hpp"

#include <utility>
#include <vector>

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
const bool REGISTERED = CalculatorFactory::registerTypeSafe<MSDCalculator>("MSDCalculator");

[[nodiscard]] std::vector<real_t> computeRunningDiffusion(const std::vector<real_t> &msd_curve,
                                                          real_t time_step) {
  const size_t num_frames = msd_curve.size();
  std::vector<real_t> d_eff(num_frames, static_cast<real_t>(0.0));
  const auto dt_d = static_cast<double>(time_step);
  if (dt_d <= 0.0) {
    return d_eff;
  }
  for (size_t i = 1; i < num_frames; ++i) {
    const double time_fs = static_cast<double>(i) * dt_d;
    d_eff[i] = static_cast<real_t>(static_cast<double>(msd_curve[i]) / (6.0 * time_fs));
  }
  return d_eff;
}
} // namespace

void MSDCalculator::calculateTrajectory(
    correlation::analysis::DistributionFunctions &dists, const correlation::core::Trajectory &traj,
    const correlation::analysis::AnalysisSettings & /*settings*/) const {
  if (traj.getFrameCount() < 2) {
    return;
  }
  auto results = calculate(traj, {-1}, {0}, {static_cast<size_t>(-1)});
  for (auto &[name, histogram] : results) {
    dists.addHistogram(name, std::move(histogram));
  }
}

std::map<std::string, correlation::analysis::Histogram>
MSDCalculator::calculate(const correlation::core::Trajectory &traj,
                         correlation::analysis::MaxFrames max_correlation_frames,
                         correlation::analysis::StartFrame start_frame,
                         correlation::analysis::EndFrame end_frame) {
  std::map<std::string, correlation::analysis::Histogram> results;

  if (traj.getFrameCount() < 2) {
    return results;
  }

  const auto partial_msd = correlation::analysis::DynamicsAnalyzer::calculatePartialMSD(
      traj, max_correlation_frames, start_frame, end_frame);

  if (partial_msd.total.size() < 2) {
    return results;
  }

  const size_t num_frames = partial_msd.total.size();
  const real_t time_step = traj.getTimeStep();

  correlation::analysis::Histogram msd_hist;
  msd_hist.x_label = "t";
  msd_hist.title = "Mean Squared Displacement";
  msd_hist.y_label = "MSD";
  msd_hist.x_unit = "fs";
  msd_hist.y_unit = "Å²";
  msd_hist.description = "Mean Squared Displacement";
  msd_hist.file_suffix = "_MSD";
  msd_hist.bins.resize(num_frames);
  for (size_t i = 0; i < num_frames; ++i) {
    msd_hist.bins[i] = static_cast<real_t>(static_cast<double>(i) * static_cast<double>(time_step));
  }
  msd_hist.partials["Total"] = partial_msd.total;
  for (const auto &[elem_sym, elem_curve] : partial_msd.element_msds) {
    msd_hist.partials[elem_sym] = elem_curve;
  }

  correlation::analysis::Histogram deff_hist;
  deff_hist.x_label = "t";
  deff_hist.title = "Running Diffusion Coefficient";
  deff_hist.y_label = "D(t)";
  deff_hist.x_unit = "fs";
  deff_hist.y_unit = "Å² fs⁻¹";
  deff_hist.description = "Running Diffusion Coefficient D(t) = MSD(t) / (6t)";
  deff_hist.file_suffix = "_D_eff";
  deff_hist.bins = msd_hist.bins;
  deff_hist.partials["Total"] = computeRunningDiffusion(partial_msd.total, time_step);
  for (const auto &[elem_sym, elem_curve] : partial_msd.element_msds) {
    deff_hist.partials[elem_sym] = computeRunningDiffusion(elem_curve, time_step);
  }

  results["MSD"] = std::move(msd_hist);
  results["D_eff"] = std::move(deff_hist);
  return results;
}

} // namespace correlation::calculators
