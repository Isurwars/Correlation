/**
 * @file EliashbergCalculator.hpp
 * @brief Hybrid MLIP-VDOS Eliashberg spectral function α²F(ω) and Allen-Dynes T_c estimator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/AnalysisTypes.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "calculators/BaseCalculator.hpp"
#include "core/Trajectory.hpp"
#include "mlip/MLIPInterface.hpp"

#include <atomic>
#include <cstddef>
#include <map>
#include <memory>
#include <span>
#include <string>

namespace correlation::calculators {

/**
 * @struct EliashbergParams
 * @brief Configuration parameters for Hybrid MLIP-VDOS Eliashberg α²F(ω) estimation.
 */
struct EliashbergParams {
  const correlation::mlip::MLIPInterface *model{
      nullptr};                              /**< Pointer to MLIP engine providing per-atom LDOS. */
  real_t mu_star{static_cast<real_t>(0.10)}; /**< Coulomb pseudopotential μ* for Allen-Dynes T_c. */
  correlation::analysis::CouplingEstimationMode coupling_mode{
      correlation::analysis::CouplingEstimationMode::TrajectoryFluctuation}; /**< Coupling mode. */
  real_t displacement_step{
      static_cast<real_t>(0.02)}; /**< Finite-displacement perturbation amplitude (Å). */
  size_t max_coupling_frames{32}; /**< Maximum trajectory frames sampled for coupling regression. */
  size_t samples_per_species{4};  /**< Representative atoms per species for finite-displacement. */
  real_t electronic_temp_k{
      static_cast<real_t>(300.0)}; /**< Electronic temperature for N_i(E_F) evaluation (K). */
  real_t omega_min_mev{
      static_cast<real_t>(0.5)}; /**< Low-frequency cutoff (meV) for 1/ω integral stability. */
};

/**
 * @struct SpeciesCouplingWeights
 * @brief Species-resolved electron-phonon deformation sensitivity proxy from MLIP LDOS.
 */
struct SpeciesCouplingWeights {
  std::map<std::string, real_t> eta_by_species; /**< Deformation sensitivity η_α (eV/Å² proxy). */
  std::map<std::string, real_t> fermi_dos_by_species; /**< Species N_α(E_F) (states/eV/atom). */
  real_t total_fermi_dos{static_cast<real_t>(0.0)};   /**< Total N(E_F) (states/eV/atom). */
};

/**
 * @struct EliashbergResult
 * @brief Output bundle containing α²F(ω) histogram, superconducting diagnostics, and coupling data.
 */
struct EliashbergResult {
  correlation::analysis::Histogram spectral_function; /**< α²F(ω) vs ħω (meV), total & species. */
  correlation::analysis::SuperconductingDiagnostics diagnostics; /**< λ, ω_log, ω_2, T_c. */
  SpeciesCouplingWeights coupling_weights; /**< Species-resolved MLIP deformation weights. */
};

/**
 * @class IElectronPhononCouplingEstimator
 * @brief Pure abstract strategy interface for estimating species electron-phonon coupling proxies.
 */
class IElectronPhononCouplingEstimator {
public:
  virtual ~IElectronPhononCouplingEstimator() = default;

  /**
   * @brief Estimates species-resolved electron-phonon deformation sensitivity from MLIP LDOS.
   * @param[in] traj Atomic trajectory.
   * @param[in] model MLIP evaluation engine producing per-atom LDOS.
   * @param[in] params Eliashberg estimation parameters.
   * @param[in] cancel_flag Optional atomic cancellation flag.
   * @return SpeciesCouplingWeights containing η_α and N_α(E_F) per species.
   */
  [[nodiscard]] virtual SpeciesCouplingWeights
  estimate(const correlation::core::Trajectory &traj, const correlation::mlip::MLIPInterface &model,
           const EliashbergParams &params,
           const std::atomic<bool> *cancel_flag = nullptr) const = 0;
};

/**
 * @class TrajectoryFluctuationCouplingEstimator
 * @brief Estimates species coupling weights by regressing per-atom δN_i(E_F, t) against thermal MD
 * displacements along sampled trajectory frames.
 */
class TrajectoryFluctuationCouplingEstimator final : public IElectronPhononCouplingEstimator {
public:
  [[nodiscard]] SpeciesCouplingWeights
  estimate(const correlation::core::Trajectory &traj, const correlation::mlip::MLIPInterface &model,
           const EliashbergParams &params,
           const std::atomic<bool> *cancel_flag = nullptr) const override;
};

/**
 * @class FiniteDisplacementCouplingEstimator
 * @brief Estimates species coupling weights via deterministic central-difference ±δ atomic
 * displacements on representative atoms per species.
 */
class FiniteDisplacementCouplingEstimator final : public IElectronPhononCouplingEstimator {
public:
  [[nodiscard]] SpeciesCouplingWeights
  estimate(const correlation::core::Trajectory &traj, const correlation::mlip::MLIPInterface &model,
           const EliashbergParams &params,
           const std::atomic<bool> *cancel_flag = nullptr) const override;
};

/**
 * @class EliashbergCalculator
 * @brief Estimates the Eliashberg spectral function α²F(ω), electron-phonon coupling constant λ,
 * logarithmic average frequency ω_log, and Allen-Dynes T_c by coupling 2PT vibrational VDOS with
 * MLIP LDOS deformation sensitivity.
 *
 * @note This calculator provides a fast MLIP-VDOS screening estimate rather than full ab-initio
 * DFPT/Wannier k-space electron-phonon matrix elements.
 */
class EliashbergCalculator : public BaseCalculator {
public:
  EliashbergCalculator() = default;
  explicit EliashbergCalculator(const correlation::mlip::MLIPInterface *model) noexcept
      : model_(model) {}

  [[nodiscard]] std::string_view getName() const override { return "Eliashberg Spectral Function"; }
  [[nodiscard]] std::string_view getShortName() const override { return "a2F"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Machine Learning"; }
  [[nodiscard]] std::string_view getDescription() const override {
    return "Estimates the Eliashberg spectral function α²F(ω), λ, ω_log, and Allen-Dynes T_c from "
           "hybrid MLIP LDOS coupling sensitivity and vibrational VDOS.";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return false; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return true; }

  /**
   * @brief Checks if an MLIP model is attached or discoverable for LDOS coupling estimation.
   * @return True if model is attached or default ORB-v3 weights exist.
   */
  [[nodiscard]] bool isConfigured() const override;

  /**
   * @brief Sets or updates the attached MLIP model engine.
   * @param[in] model Pointer to MLIPInterface engine.
   */
  void setModel(const correlation::mlip::MLIPInterface *model) noexcept { model_ = model; }

  /**
   * @brief Gets the attached MLIP model engine.
   * @return Pointer to MLIPInterface engine, or nullptr if unset.
   */
  [[nodiscard]] const correlation::mlip::MLIPInterface *getModel() const noexcept { return model_; }

  /**
   * @brief Computes Eliashberg α²F(ω) and superconducting diagnostics for a trajectory.
   * @param[in,out] dists DistributionFunctions container (uses existing VDOS if present).
   * @param[in] traj Atomic trajectory with velocities.
   * @param[in] settings Analysis configuration settings.
   */
  void calculateTrajectory(correlation::analysis::DistributionFunctions &dists,
                           const correlation::core::Trajectory &traj,
                           const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Computes superconducting diagnostics (λ, ω_log, ω_2, McMillan & Allen-Dynes T_c)
   * directly from an arbitrary Eliashberg spectral function α²F(ω) on an energy grid in meV.
   * @param[in] omega_mev Monotonically increasing phonon energy grid ħω in meV.
   * @param[in] a2f Dimensionless Eliashberg spectral function α²F(ω).
   * @param[in] mu_star Coulomb pseudopotential μ* (typically 0.10–0.15).
   * @param[in] omega_min_mev Low-frequency cutoff in meV to guard 1/ω integration.
   * @return Populated SuperconductingDiagnostics struct.
   */
  [[nodiscard]] static correlation::analysis::SuperconductingDiagnostics
  computeFromSpectrum(std::span<const real_t> omega_mev, std::span<const real_t> a2f,
                      real_t mu_star = static_cast<real_t>(0.10),
                      real_t omega_min_mev = static_cast<real_t>(0.5)) noexcept;

  /**
   * @brief Computes hybrid MLIP-VDOS Eliashberg α²F(ω) and diagnostics from a trajectory and
   * vibrational VDOS histogram.
   * @param[in] traj Atomic trajectory.
   * @param[in] vdos_vib Vibrational-only VDOS histogram g_vib(ν) in THz.
   * @param[in] params Eliashberg estimation parameters (including MLIP model pointer).
   * @param[in] estimator Optional custom coupling strategy (nullptr uses params.coupling_mode).
   * @param[in] cancel_flag Optional atomic cancellation flag.
   * @return Populated EliashbergResult.
   */
  [[nodiscard]] static EliashbergResult
  calculate(const correlation::core::Trajectory &traj,
            const correlation::analysis::Histogram &vdos_vib, const EliashbergParams &params = {},
            const IElectronPhononCouplingEstimator *estimator = nullptr,
            const std::atomic<bool> *cancel_flag = nullptr);

  /**
   * @brief Factory helper creating the concrete coupling estimator for a given mode.
   * @param[in] mode Coupling estimation strategy enum.
   * @return Unique pointer to concrete IElectronPhononCouplingEstimator.
   */
  [[nodiscard]] static std::unique_ptr<IElectronPhononCouplingEstimator>
  createCouplingEstimator(correlation::analysis::CouplingEstimationMode mode);

private:
  const correlation::mlip::MLIPInterface *model_{nullptr};
};

} // namespace correlation::calculators
