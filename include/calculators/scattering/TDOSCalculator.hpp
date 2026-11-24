/**
 * @file TDOSCalculator.hpp
 * @brief Total Density of States (TDOS) and site-resolved LDOS calculator from MLIP predictions.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/AnalysisTypes.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "calculators/BaseCalculator.hpp"
#include "mlip/MLIPInterface.hpp"

#include <atomic>
#include <cstddef>
#include <vector>

namespace correlation::calculators {

/**
 * @struct TDOSParams
 * @brief Parameter configuration for Total and Local Density of States evaluation.
 */
struct TDOSParams {
  real_t e_min = static_cast<real_t>(-30.0); /**< Lower energy boundary in eV (relative to E_F). */
  real_t e_max = static_cast<real_t>(15.0);  /**< Upper energy boundary in eV (relative to E_F). */
  const correlation::mlip::MLIPInterface *model{
      nullptr}; /**< Pointer to MLIP evaluation model engine. */
  correlation::analysis::LDOSBroadeningKernel broadening_kernel{
      correlation::analysis::LDOSBroadeningKernel::None}; /**< Optional spectral broadening kernel.
                                                           */
  real_t broadening_sigma{static_cast<real_t>(0.0)};      /**< Broadening width sigma/gamma (eV). */
  real_t electronic_temp_k{static_cast<real_t>(
      300.0)}; /**< Electronic temperature for Fermi-Dirac N(E_F) smearing (K). */
  real_t gap_threshold{
      static_cast<real_t>(1e-3)}; /**< DOS threshold (states/eV/atom) for band-gap detection. */
  size_t rebin_num_bins{0}; /**< Optional target bin count for re-binning (0 = keep native). */
  bool compute_idos{true};  /**< Compute companion Integrated DOS (IDOS) histogram. */
  bool retain_per_atom_ldos{false}; /**< Retain full [N_atoms * N_bins] per-atom LDOS buffer. */
};

/**
 * @struct LDOSResult
 * @brief Comprehensive electronic DOS output containing TDOS, IDOS, diagnostics, and optional
 * site-resolved LDOS.
 */
struct LDOSResult {
  correlation::analysis::Histogram tdos; /**< Total and species-resolved TDOS. */
  correlation::analysis::Histogram idos; /**< Cumulative Integrated DOS (IDOS). */
  correlation::analysis::ElectronicDiagnostics
      diagnostics;                   /**< Fermi DOS, gap, per-atom N_i(E_F). */
  std::vector<real_t> per_atom_ldos; /**< Optional contiguous [N_atoms * ldos_bins] spectra. */
  size_t ldos_bins{0};               /**< Number of energy bins per atom in per_atom_ldos. */
};

/**
 * @class TDOSCalculator
 * @brief Computes Total Density of States (TDOS), Integrated DOS (IDOS), and site-resolved
 * electronic diagnostics from MLIP per-atom Local Density of States predictions.
 */
class TDOSCalculator : public BaseCalculator {
public:
  TDOSCalculator() = default;
  explicit TDOSCalculator(const correlation::mlip::MLIPInterface *model) noexcept : model_(model) {}

  [[nodiscard]] std::string_view getName() const override {
    return "Total Density of States (TDOS)";
  }
  [[nodiscard]] std::string_view getShortName() const override { return "TDOS"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Machine Learning"; }
  [[nodiscard]] std::string_view getDescription() const override {
    return "Computes the Total Density of States (TDOS), Integrated DOS (IDOS), and Fermi-level "
           "electronic diagnostics from MLIP per-atom Local Density of States predictions.";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return true; }

  /**
   * @brief Checks if a valid MLIP evaluation model is attached or discoverable.
   * @return True if model is attached or default ORB-v3 weights exist.
   */
  [[nodiscard]] bool isConfigured() const override;

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
   * @brief Calculates per-frame TDOS, IDOS, and electronic diagnostics into DistributionFunctions.
   * @param[in,out] dists DistributionFunctions target container.
   * @param[in] settings Analysis configuration settings.
   */
  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Calculates trajectory-averaged TDOS, IDOS, and electronic diagnostics into
   * DistributionFunctions.
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
  [[nodiscard]] static correlation::analysis::Histogram
  calculate(const correlation::core::Cell &cell, const TDOSParams &params = {});

  /**
   * @brief Evaluates full LDOS/TDOS, IDOS, and electronic diagnostics for a single unit cell.
   * @param[in] cell Atomic configuration unit cell.
   * @param[in] params Evaluation parameters and model pointer.
   * @return Populated LDOSResult containing TDOS, IDOS, diagnostics, and optional per-atom LDOS.
   */
  [[nodiscard]] static LDOSResult calculateDetailed(const correlation::core::Cell &cell,
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

  /**
   * @brief Evaluates trajectory-averaged TDOS, IDOS, and electronic diagnostics across frames.
   * @param[in] traj Trajectory containing atomic frames.
   * @param[in] params Evaluation parameters and model pointer.
   * @param[in] cancel_flag Optional atomic flag to abort computation prematurely.
   * @return Populated LDOSResult averaged over processed frames.
   */
  [[nodiscard]] static LDOSResult
  calculateTrajectoryDetailed(const correlation::core::Trajectory &traj,
                              const TDOSParams &params = {},
                              const std::atomic<bool> *cancel_flag = nullptr);

private:
  const correlation::mlip::MLIPInterface *model_{nullptr};
};

} // namespace correlation::calculators
