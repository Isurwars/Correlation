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

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
const bool REGISTERED = CalculatorFactory::registerTypeSafe<MSDCalculator>("MSDCalculator");
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

  std::vector<real_t> raw_msd = correlation::analysis::DynamicsAnalyzer::calculateMSD(
      traj, max_correlation_frames, start_frame, end_frame);

  if (raw_msd.size() < 2) {
    return results;
  }

  const size_t num_frames = raw_msd.size();
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
    msd_hist.bins[i] = static_cast<real_t>(i) * time_step;
  }
  msd_hist.partials["Total"] = raw_msd;
  results["MSD"] = std::move(msd_hist);

  correlation::analysis::Histogram deff_hist;
  deff_hist.x_label = "t";
  deff_hist.title = "Running Diffusion Coefficient";
  deff_hist.y_label = "D(t)";
  deff_hist.x_unit = "fs";
  deff_hist.y_unit = "Å² fs⁻¹";
  deff_hist.description = "Running Diffusion Coefficient D(t) = MSD(t) / (6t)";
  deff_hist.file_suffix = "_MSD";
  deff_hist.bins = results["MSD"].bins; // Same time axis
  std::vector<real_t> d_eff(num_frames, 0.0);
  for (size_t i = 1; i < num_frames; ++i) {
    const real_t time = static_cast<real_t>(i) * time_step;
    if (time > 0.0) {
      d_eff[i] = raw_msd[i] / (static_cast<real_t>(6.0) * time);
    }
  }
  deff_hist.partials["Total"] = d_eff;
  results["D_eff"] = std::move(deff_hist);

  return results;
}

} // namespace correlation::calculators
