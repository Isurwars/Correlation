/**
 * @file XRDCalculator.hpp
 * @brief X-ray diffraction (XRD) pattern calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"

#include <map>
#include <string>
namespace correlation::calculators {

/**
 * @brief Strong type wrapper for X-ray incident radiation wavelength in Angstroms.
 */
struct Wavelength {
  real_t value; ///< Wavelength value in Angstroms (e.g. 1.5406 for Cu K-alpha).
};

/**
 * @brief Strong type wrapper for minimum 2-theta angle in degrees.
 */
struct MinTheta {
  real_t value; ///< Minimum 2-theta scattering angle in degrees.
};

/**
 * @brief Strong type wrapper for maximum 2-theta angle in degrees.
 */
struct MaxTheta {
  real_t value; ///< Maximum 2-theta scattering angle in degrees.
};

/**
 * @brief Strong type wrapper for angular bin width in degrees.
 */
struct BinWidth {
  real_t value; ///< Angular bin resolution in degrees.
};

/**
 * @class XRDCalculator
 * @brief Computes the X-Ray Diffraction (XRD) pattern from g(r).
 */
class XRDCalculator : public BaseCalculator {
public:
  [[nodiscard]] std::string getName() const override { return "XRD"; }
  [[nodiscard]] std::string getShortName() const override { return "XRD"; }
  [[nodiscard]] std::string getGroup() const override { return "Scattering"; }
  [[nodiscard]] std::string getDescription() const override { return "Computes the X-Ray Diffraction Pattern (XRD)."; }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief High-performance computation of the X-Ray Diffraction (XRD) pattern.
   *
   * @param g_r_hist Input Radial Distribution Function g(r).
   * @param cell The periodic cell.
   * @param ashcroft_weights Composition-dependent scattering weights.
   * @param lambda X-ray wavelength (Angstrom).
   * @param theta_min Minimum 2-theta angle (degrees).
   * @param theta_max Maximum 2-theta angle (degrees).
   * @param bin_width Angular resolution (degrees).
   * @return A histogram representing intensity vs 2-theta.
   */
  static correlation::analysis::Histogram calculate(const correlation::analysis::Histogram &g_r_hist,
                                                    const correlation::core::Cell &cell,
                                                    const std::map<std::string, real_t> &ashcroft_weights,
                                                    Wavelength lambda, MinTheta theta_min, MaxTheta theta_max,
                                                    BinWidth bin_width);

private:
  static real_t getAtomicFormFactor(const std::string &symbol, real_t q_value);

  static std::map<std::string, real_t> calculateConcentrations(const correlation::core::Cell &cell);

  static std::map<std::string, std::vector<real_t>>
  calculatePartialIntegrands(const correlation::analysis::Histogram &g_r_hist,
                             const std::map<std::string, real_t> &ashcroft_weights, real_t delta_r);
};

} // namespace correlation::calculators
