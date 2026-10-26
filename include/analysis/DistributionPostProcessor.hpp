/**
 * @file DistributionPostProcessor.hpp
 * @brief Mathematical post-processing and smoothing service for distribution functions.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/AnalysisTypes.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "math/Smoothing.hpp"

#include <map>
#include <string>
#include <string_view>

namespace correlation::analysis {

/**
 * @brief Service responsible for mathematical smoothing and post-processing of histogram data.
 */
class DistributionPostProcessor {
public:
  /**
   * @brief Applies kernel smoothing to a single Histogram in-place.
   * @param[in,out] hist The histogram to be smoothed.
   * @param[in] name The name or key of the histogram (used to select angular vs. spatial domain
   * scale).
   * @param[in] sigma Smoothing width parameter.
   * @param[in] kernel Smoothing kernel type (e.g., Gaussian, Epanechnikov).
   */
  static void
  smooth(Histogram &hist, std::string_view name, real_t sigma,
         correlation::math::KernelType kernel = correlation::math::KernelType::Gaussian);

  /**
   * @brief Applies kernel smoothing to all histograms within a histogram map.
   * @param[in,out] histograms Map of name-to-Histogram entries.
   * @param[in] sigma Smoothing width parameter.
   * @param[in] kernel Smoothing kernel type.
   */
  static void
  smoothAll(std::map<std::string, Histogram> &histograms, real_t sigma,
            correlation::math::KernelType kernel = correlation::math::KernelType::Gaussian);

  /**
   * @brief Smooths a specified histogram within a DistributionFunctions container.
   * @param[in,out] dists DistributionFunctions instance containing histograms.
   * @param[in] name Name of the target histogram.
   * @param[in] sigma Smoothing width parameter.
   * @param[in] kernel Smoothing kernel type.
   */
  static void
  smooth(DistributionFunctions &dists, const std::string &name, real_t sigma,
         correlation::math::KernelType kernel = correlation::math::KernelType::Gaussian);

  /**
   * @brief Smooths all histograms within a DistributionFunctions container.
   * @param[in,out] dists DistributionFunctions instance containing histograms.
   * @param[in] sigma Smoothing width parameter.
   * @param[in] kernel Smoothing kernel type.
   */
  static void
  smoothAll(DistributionFunctions &dists, real_t sigma,
            correlation::math::KernelType kernel = correlation::math::KernelType::Gaussian);
};

} // namespace correlation::analysis
