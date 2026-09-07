/**
 * @file CorrelationEngine.cpp
 * @brief Implementation of analysis pipeline orchestration engine.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "analysis/CorrelationEngine.hpp"
#include "analysis/DynamicsAnalyzer.hpp"
#include "analysis/TrajectoryAnalyzer.hpp"
#include "calculators/CalculatorFactory.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

namespace correlation::analysis {

namespace {

[[nodiscard]] std::string validateHistogramSettings(const AnalysisSettings &settings) {
  if (settings.r_max <= 0.0) {
    return "Error: r_max must be strictly positive.";
  }
  if (settings.r_bin_width <= 0.0) {
    return "Error: r_bin_width must be strictly positive.";
  }
  if (settings.r_bin_width >= settings.r_max) {
    return "Error: r_bin_width must be strictly less than r_max.";
  }
  if (settings.q_max <= 0.0) {
    return "Error: q_max must be strictly positive.";
  }
  if (settings.q_bin_width <= 0.0) {
    return "Error: q_bin_width must be strictly positive.";
  }
  if (settings.q_bin_width >= settings.q_max) {
    return "Error: q_bin_width must be strictly less than q_max.";
  }
  if (settings.r_int_max <= 0.0) {
    return "Error: r_int_max must be strictly positive.";
  }
  return "";
}

[[nodiscard]] std::string validateAngularSettings(const AnalysisSettings &settings) {
  if (settings.angle_bin_width <= 0.0) {
    return "Error: angle_bin_width must be strictly positive.";
  }
  if (settings.angle_bin_width > 180.0) {
    return "Error: angle_bin_width must be at most 180.0 degrees.";
  }
  if (settings.dihedral_bin_width <= 0.0) {
    return "Error: dihedral_bin_width must be strictly positive.";
  }
  if (settings.dihedral_bin_width > 360.0) {
    return "Error: dihedral_bin_width must be at most 360.0 degrees.";
  }
  return "";
}

[[nodiscard]] std::string validateSimulationSettings(const CorrelationEngineConfig &config) {
  if (config.time_step <= 0.0) {
    return "Error: time_step must be strictly positive.";
  }
  if (config.settings.max_ring_size < 3) {
    return "Error: max_ring_size must be an integer >= 3.";
  }
  if (config.settings.hyperuniformity_samples == 0) {
    return "Error: hyper_samples must be strictly positive.";
  }
  if (config.max_frame < -1) {
    return "Error: max_frame cannot be less than -1.";
  }
  if (config.max_frame >= 0 && config.min_frame > config.max_frame) {
    return "Error: min_frame cannot be greater than max_frame.";
  }
  return "";
}

[[nodiscard]] std::string validateSmoothingSettings(const AnalysisSettings &settings) {
  if (settings.smoothing && settings.smoothing_sigma <= 0.0) {
    return "Error: smoothing_sigma must be strictly positive when smoothing is enabled.";
  }
  if (settings.smoothing_sigma < 0.0) {
    return "Error: smoothing_sigma cannot be negative.";
  }
  if (settings.lef_cutoff <= 0.0) {
    return "Error: lef_cutoff must be strictly positive.";
  }
  if (settings.lef_sigma <= 0.0) {
    return "Error: lef_sigma must be strictly positive.";
  }
  return "";
}

[[nodiscard]] std::string validateBondCutoffs(const BondCutoffMatrix &cutoffs) {
  for (const auto &row : cutoffs) {
    for (const auto &range : row) {
      if (range.min_sq < 0.0 || range.max_sq < 0.0) {
        return "Error: Bond cutoff bounds cannot be negative.";
      }
      if (range.min_sq > range.max_sq) {
        return "Error: Minimum bond cutoff cannot be greater than maximum bond cutoff.";
      }
    }
  }
  return "";
}

void configureTrajectory(correlation::core::Trajectory &trajectory, const CorrelationEngineConfig &config,
                         size_t &start_f) {
  if (!config.bond_cutoffs.empty()) {
    trajectory.setBondCutoffs(config.bond_cutoffs);
  } else if (trajectory.getBondCutoffs().empty()) {
    trajectory.precomputeBondCutoffs();
  }

  start_f = config.min_frame >= 0 ? static_cast<size_t>(config.min_frame) : 0;
  if (start_f >= trajectory.getFrameCount()) {
    start_f = 0;
  }
  trajectory.setTimeStep(config.time_step);
}

bool requiresVelocities(const AnalysisSettings &settings) {
  const auto &factory_calcs = ::correlation::calculators::CalculatorFactory::instance().getCalculators();
  for (const auto &calc : factory_calcs) {
    bool const is_active = settings.isActive(calc->getName()) || settings.isActive(calc->getShortName());
    if (calc->isTrajectoryCalculator() && calc->isConfigured() && is_active) {
      if (calc->getName() == "VACF" || calc->getShortName() == "VACF" || calc->getName() == "vDoS" ||
          calc->getShortName() == "vDoS") {
        return true;
      }
    }
  }
  return false;
}

} // namespace

std::string CorrelationEngine::validateConfig(const CorrelationEngineConfig &config) {
  if (std::string const err = validateHistogramSettings(config.settings); !err.empty()) {
    return err;
  }
  if (std::string const err = validateAngularSettings(config.settings); !err.empty()) {
    return err;
  }
  if (std::string const err = validateSimulationSettings(config); !err.empty()) {
    return err;
  }
  if (std::string const err = validateSmoothingSettings(config.settings); !err.empty()) {
    return err;
  }
  return validateBondCutoffs(config.bond_cutoffs);
}

void CorrelationEngine::runTrajectoryCalculators(DistributionFunctions &distribution_functions,
                                                 correlation::core::Trajectory &trajectory,
                                                 const AnalysisSettings &settings) {
  if (requiresVelocities(settings)) {
    if (trajectory.getFrameCount() == 0) {
      throw std::runtime_error("VACF requires a trajectory with frames");
    }
    trajectory.calculateVelocities();
  }

  const auto &factory_calcs = ::correlation::calculators::CalculatorFactory::instance().getCalculators();
  for (const auto &calc : factory_calcs) {
    if (calc->isTrajectoryCalculator() && calc->isConfigured() &&
        (settings.isActive(calc->getName()) || settings.isActive(calc->getShortName()))) {
      calc->calculateTrajectory(distribution_functions, trajectory, settings);
    }
  }
}

void CorrelationEngine::calculateDynamicProperties(DistributionFunctions &distribution_functions) {
  const auto &hists = distribution_functions.getAllHistograms();

  auto it_msd = hists.find("MSD");
  if (it_msd != hists.end()) {
    const auto &hist = it_msd->second;
    auto it_total = hist.partials.find("Total");
    if (it_total != hist.partials.end()) {
      real_t const d_msd = DynamicsAnalyzer::computeDiffusionCoefficientMSD(hist.bins, it_total->second);
      distribution_functions.setDiffusionCoefficientMSD(d_msd);
    }
  }

  auto it_vacf = hists.find("VACF");
  if (it_vacf != hists.end()) {
    const auto &hist = it_vacf->second;
    auto it_total = hist.partials.find("Total");
    if (it_total != hist.partials.end()) {
      real_t const d_vacf = DynamicsAnalyzer::computeDiffusionCoefficientVACF(hist.bins, it_total->second);
      distribution_functions.setDiffusionCoefficientVACF(d_vacf);
    }
  }

  auto it_norm = hists.find("Normalized VACF");
  if (it_norm != hists.end()) {
    const auto &hist = it_norm->second;
    auto it_total = hist.partials.find("Total");
    if (it_total != hist.partials.end()) {
      real_t const tau = DynamicsAnalyzer::computeRelaxationTime(hist.bins, it_total->second);
      distribution_functions.setRelaxationTime(tau);
      real_t deborah = 0.0;
      if (!hist.bins.empty() && hist.bins.back() > 0.0) {
        deborah = tau / hist.bins.back();
      }
      distribution_functions.setDeborahNumber(deborah);
    }
  }

  if (distribution_functions.getDiffusionCoefficientMSD() > 0.0) {
    std::cout << "Self-diffusion coefficient (from MSD): " << distribution_functions.getDiffusionCoefficientMSD()
              << " Å²/fs\n";
  }
  if (distribution_functions.getDiffusionCoefficientVACF() > 0.0) {
    std::cout << "Self-diffusion coefficient (from VACF): " << distribution_functions.getDiffusionCoefficientVACF()
              << " Å²/fs\n";
  }
  if (distribution_functions.getRelaxationTime() > 0.0) {
    std::cout << "Relaxation time (from VACF): " << distribution_functions.getRelaxationTime() << " fs\n";
    std::cout << "Deborah number: " << distribution_functions.getDeborahNumber() << '\n';
  }
}

std::expected<std::unique_ptr<DistributionFunctions>, std::string>
CorrelationEngine::runAnalysis(correlation::core::Trajectory &trajectory, const CorrelationEngineConfig &config,
                               std::function<void(float, const std::string &)> progress_callback) {
  if (trajectory.getFrameCount() == 0) {
    return std::unexpected("Analysis aborted: No trajectory loaded.");
  }

  std::string const validation_error = validateConfig(config);
  if (!validation_error.empty()) {
    return std::unexpected(validation_error);
  }

  try {
    size_t start_f = 0;
    configureTrajectory(trajectory, config, start_f);

    const auto &active_cutoffs = !config.bond_cutoffs.empty() ? config.bond_cutoffs : trajectory.getBondCutoffs();

    if (progress_callback) {
      progress_callback(0.0F, "Starting analysis...");
    }

    auto cb_structure = [&progress_callback](float /*p*/, const std::string &msg) {
      if (progress_callback) {
        progress_callback(0.0F, msg);
      }
    };

    auto cb_dist = [&progress_callback](float progress, const std::string &msg) {
      if (progress_callback) {
        progress_callback(progress, msg);
      }
    };

    const TrajectoryAnalyzer trajectory_analyzer(trajectory, config.settings.r_max, active_cutoffs, StartFrame{start_f},
                                                 EndFrame{static_cast<size_t>(config.max_frame)}, true, cb_structure);

    auto dist_funcs =
        DistributionFunctions::computeMean(trajectory, trajectory_analyzer, start_f, config.settings, cb_dist);

    if (dist_funcs) {
      runTrajectoryCalculators(*dist_funcs, trajectory, config.settings);
      calculateDynamicProperties(*dist_funcs);
    }

    return dist_funcs;
  } catch (const std::exception &e) {
    return std::unexpected(std::string("Error during analysis: ") + e.what());
  }
}

} // namespace correlation::analysis
