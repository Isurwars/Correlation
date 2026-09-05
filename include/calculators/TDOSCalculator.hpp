/**
 * @file TDOSCalculator.hpp
 * @brief Total Density of States (TDOS) calculator from MLIP predictions.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "mlip/MLIPInterface.hpp"

#include <atomic>
#include <string>

namespace correlation::calculators {

/**
 * @struct TDOSParams
 * @brief Parameter configuration for Total Density of States evaluation.
 */
struct TDOSParams {
  real_t e_min = static_cast<real_t>(-15.0);              /**< Lower energy boundary in eV (relative to E_F). */
  real_t e_max = static_cast<real_t>(5.0);                /**< Upper energy boundary in eV (relative to E_F). */
  const correlation::mlip::MLIPInterface *model{nullptr}; /**< Pointer to MLIP evaluation model engine. */
};

/**
 * @class TDOSCalculator
 * @brief Computes Total Density of States (TDOS) from MLIP per-atom Local Density of States predictions.
 */
class TDOSCalculator : public BaseCalculator {
public:
  TDOSCalculator() = default;
  explicit TDOSCalculator(const correlation::mlip::MLIPInterface *model) noexcept : model_(model) {}

  [[nodiscard]] std::string getName() const override { return "Total Density of States (TDOS)"; }
  [[nodiscard]] std::string getShortName() const override { return "TDOS"; }
  [[nodiscard]] std::string getGroup() const override { return "Machine Learning"; }
  [[nodiscard]] std::string getDescription() const override {
    return "Computes the Total Density of States (TDOS) from MLIP per-atom Local Density of States predictions.";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return true; }

  /**
   * @brief Checks if a valid MLIP evaluation model is attached.
   * @return True if model is non-null; false otherwise.
   */
  [[nodiscard]] bool isConfigured() const override { return model_ != nullptr; }

  /**
   * @brief Sets or updates the attached MLIP model engine.
   * @param model Pointer to the MLIPInterface engine.
   */
  void setModel(const correlation::mlip::MLIPInterface *model) noexcept { model_ = model; }

  /**
   * @brief Gets the attached MLIP model engine.
   * @return Pointer to the MLIPInterface engine, or nullptr if unset.
   */
  [[nodiscard]] const correlation::mlip::MLIPInterface *getModel() const noexcept { return model_; }

  /**
   * @brief Calculates per-frame TDOS and inserts into DistributionFunctions.
   * @param[in,out] dists DistributionFunctions target container.
   * @param[in] settings Analysis configuration settings.
   */
  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Calculates trajectory-averaged TDOS and inserts into DistributionFunctions.
   * @param[in,out] dists DistributionFunctions target container.
   * @param[in] traj Atomic trajectory to evaluate.
   * @param[in] settings Analysis configuration settings.
   */
  void calculateTrajectory(correlation::analysis::DistributionFunctions &dists,
                           const correlation::core::Trajectory &traj,
                           const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Evaluates TDOS for a single unit cell.
   * @param[in] cell Atomic configuration unit cell.
   * @param[in] params Evaluation parameters and model pointer.
   * @return Generated Histogram with total and species-resolved partials.
   */
  [[nodiscard]] static correlation::analysis::Histogram calculate(const correlation::core::Cell &cell,
                                                                  const TDOSParams &params = {});

  /**
   * @brief Evaluates trajectory-averaged TDOS across multiple frames.
   * @param[in] traj Trajectory containing atomic frames.
   * @param[in] params Evaluation parameters and model pointer.
   * @param[in] cancel_flag Optional atomic flag to abort computation prematurely.
   * @return Generated Histogram averaged over processed frames.
   */
  [[nodiscard]] static correlation::analysis::Histogram
  calculateTrajectory(const correlation::core::Trajectory &traj, const TDOSParams &params = {},
                      const std::atomic<bool> *cancel_flag = nullptr);

private:
  const correlation::mlip::MLIPInterface *model_{nullptr};
};

} // namespace correlation::calculators
