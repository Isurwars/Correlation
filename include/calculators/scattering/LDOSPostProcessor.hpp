/**
 * @file LDOSPostProcessor.hpp
 * @brief Stateless electronic LDOS/TDOS post-processing utilities (broadening, re-binning,
 * Fermi-level N(E_F), band gap detection, and Integrated DOS).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/AnalysisTypes.hpp"
#include "math/Precision.hpp"

#include <span>
#include <vector>

namespace correlation::calculators {

/**
 * @struct BandGapInfo
 * @brief Electronic band-edge and gap summary around E_F = 0 eV.
 */
struct BandGapInfo {
  real_t band_gap{static_cast<real_t>(0.0)};   ///< Energy gap E_g = max(0, CBM - VBM) in eV.
  real_t vbm_energy{static_cast<real_t>(0.0)}; ///< Valence band maximum relative to E_F (eV).
  real_t cbm_energy{static_cast<real_t>(0.0)}; ///< Conduction band minimum relative to E_F (eV).
  bool is_metallic{false};                     ///< True when DOS at E_F exceeds gap threshold.
};

/**
 * @class LDOSPostProcessor
 * @brief Pure mathematical post-processing routines for electronic density-of-states spectra.
 */
class LDOSPostProcessor {
public:
  /**
   * @brief Convolves a DOS spectrum with an area-preserving Gaussian or Lorentzian kernel.
   * @param[in] energies Monotonically increasing energy grid (eV).
   * @param[in] dos Input density of states values (states/eV/atom).
   * @param[in] kernel Broadening kernel type (None, Gaussian, Lorentzian).
   * @param[in] sigma Broadening width (standard deviation for Gaussian, HWHM for Lorentzian) in eV.
   * @return Broadened spectrum with preserved spectral integral.
   */
  [[nodiscard]] static std::vector<real_t>
  broadenSpectrum(std::span<const real_t> energies, std::span<const real_t> dos,
                  correlation::analysis::LDOSBroadeningKernel kernel, real_t sigma);

  /**
   * @brief Resamples a DOS spectrum onto a target energy grid via piecewise-linear interpolation.
   * @param[in] src_energies Monotonically increasing source energy grid (eV).
   * @param[in] src_dos Source DOS values matching @p src_energies.
   * @param[in] dst_energies Monotonically increasing target energy grid (eV).
   * @return Resampled DOS values evaluated on @p dst_energies.
   */
  [[nodiscard]] static std::vector<real_t> rebinSpectrum(std::span<const real_t> src_energies,
                                                         std::span<const real_t> src_dos,
                                                         std::span<const real_t> dst_energies);

  /**
   * @brief Evaluates the density of states at the Fermi level N(E_F = 0 eV).
   *
   * When @p electronic_temp_k <= 0, performs exact linear interpolation at E = 0.
   * When @p electronic_temp_k > 0, evaluates the Fermi-Dirac thermal smearing integral
   * N(E_F) = -∫ (∂f/∂E) N(E) dE using k_B T = 8.617333262e-5 * T (eV).
   *
   * @param[in] energies Monotonically increasing energy grid relative to E_F (eV).
   * @param[in] dos Density of states values (states/eV/atom).
   * @param[in] electronic_temp_k Electronic smearing temperature in Kelvin (0 = exact E = 0).
   * @return Evaluated N(E_F) in states/eV/atom.
   */
  [[nodiscard]] static real_t
  evaluateFermiDos(std::span<const real_t> energies, std::span<const real_t> dos,
                   real_t electronic_temp_k = static_cast<real_t>(0.0)) noexcept;

  /**
   * @brief Detects valence/conduction band edges and electronic band gap around E_F = 0 eV.
   * @param[in] energies Monotonically increasing energy grid relative to E_F (eV).
   * @param[in] dos Total density of states values (states/eV/atom).
   * @param[in] threshold Minimum DOS value considered occupied/conductive (states/eV/atom).
   * @return BandGapInfo containing VBM, CBM, E_g, and metallicity flag.
   */
  [[nodiscard]] static BandGapInfo
  detectBandGap(std::span<const real_t> energies, std::span<const real_t> dos,
                real_t threshold = static_cast<real_t>(1e-3)) noexcept;

  /**
   * @brief Computes the cumulative trapezoidal Integrated Density of States (IDOS).
   * @param[in] energies Monotonically increasing energy grid (eV).
   * @param[in] dos Density of states values (states/eV/atom).
   * @return Cumulative integral I(E) = ∫_{E_min}^{E} N(E') dE' (states/atom).
   */
  [[nodiscard]] static std::vector<real_t> computeIntegratedDos(std::span<const real_t> energies,
                                                                std::span<const real_t> dos);
};

} // namespace correlation::calculators
