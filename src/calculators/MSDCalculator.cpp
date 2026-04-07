/**
 * @file MSDCalculator.cpp
 * @brief Implementation of Mean Squared Displacement calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "calculators/MSDCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "DynamicsAnalyzer.hpp"

namespace {
bool registered = CalculatorFactory::instance().registerCalculator(
    std::make_unique<MSDCalculator>());
} // namespace

void MSDCalculator::calculateTrajectory(DistributionFunctions &df,
                                        const Trajectory &traj,
                                        const AnalysisSettings &settings) const {
  auto results = calculate(traj, -1, 0, static_cast<size_t>(-1));
  for (auto &[name, histogram] : results) {
    df.addHistogram(name, std::move(histogram));
  }
}

std::map<std::string, Histogram>
MSDCalculator::calculate(const Trajectory &traj, int max_correlation_frames,
                         size_t start_frame, size_t end_frame) {
  std::map<std::string, Histogram> results;

  // MSD requires positional data only (frames), not velocity data
  if (traj.getFrames().empty()) {
    return results;
  }

  std::vector<double> raw_msd = DynamicsAnalyzer::calculateMSD(
      traj, max_correlation_frames, start_frame, end_frame);

  if (raw_msd.empty()) {
    return results;
  }

  const size_t num_frames = raw_msd.size();
  const double dt = traj.getTimeStep();

  // --- MSD histogram: bins = time (fs), partials["Total"] = MSD (Å²) ---
  Histogram msd_hist;
  msd_hist.bin_label = "Time";
  msd_hist.bins.resize(num_frames);
  for (size_t i = 0; i < num_frames; ++i) {
    msd_hist.bins[i] = static_cast<double>(i) * dt;
  }
  msd_hist.partials["Total"] = raw_msd;
  results["MSD"] = std::move(msd_hist);

  // --- Running D(t) = MSD(t) / (6 * t): skip lag=0 (division by zero) ---
  // We store D_eff starting from lag=1; lag=0 gets D=0 by convention.
  Histogram deff_hist;
  deff_hist.bin_label = "Time";
  deff_hist.bins = results["MSD"].bins; // Same time axis
  std::vector<double> d_eff(num_frames, 0.0);
  for (size_t i = 1; i < num_frames; ++i) {
    const double t = static_cast<double>(i) * dt;
    if (t > 0.0) {
      d_eff[i] = raw_msd[i] / (6.0 * t);
    }
  }
  deff_hist.partials["Total"] = d_eff;
  results["D_eff"] = std::move(deff_hist);

  return results;
}
