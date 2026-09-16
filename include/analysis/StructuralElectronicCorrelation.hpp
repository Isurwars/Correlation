/**
 * @file StructuralElectronicCorrelation.hpp
 * @brief Structural-electronic correlation pipelines (CNA and Steinhardt motif-projected TDOS).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "calculators/TDOSCalculator.hpp"
#include "core/Trajectory.hpp"
#include "math/Precision.hpp"

#include <atomic>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace correlation::analysis {

/**
 * @struct MotifProjectedTDOS
 * @brief Total Density of States partitioned across structural motifs.
 */
struct MotifProjectedTDOS {
  std::vector<real_t> energies;                          /**< Energy grid values [num_bins] in eV. */
  std::map<std::string, std::vector<real_t>> motif_tdos; /**< Partial TDOS spectra by motif name. */
  std::vector<real_t> total_tdos;                        /**< Total aggregated TDOS spectrum. */
  size_t frame_count{0};                                 /**< Evaluated frame count. */

  /**
   * @brief Converts the motif-projected TDOS into a standard Histogram.
   * @param title Title for the resulting histogram.
   * @return Populated Histogram containing total and motif partials.
   */
  [[nodiscard]] Histogram toHistogram(std::string_view title = "Motif-Projected TDOS") const;
};

/**
 * @class StructuralElectronicCorrelation
 * @brief Trajectory-level pipeline correlating MLIP electronic structures with local structural motifs.
 */
class StructuralElectronicCorrelation {
public:
  /**
   * @brief Correlates trajectory LDoS with Common Neighbor Analysis (CNA) classifications.
   *
   * @param[in,out] dists DistributionFunctions target container.
   * @param[in] traj Atomic trajectory to evaluate.
   * @param[in] params TDOS evaluation parameters and MLIP model pointer.
   * @param[in] cancel_flag Optional atomic flag to abort computation prematurely.
   * @return Generated MotifProjectedTDOS structure.
   */
  static MotifProjectedTDOS correlateCNA(DistributionFunctions &dists, const core::Trajectory &traj,
                                         const calculators::TDOSParams &params,
                                         const std::atomic<bool> *cancel_flag = nullptr);

  /**
   * @brief Correlates trajectory LDoS with Steinhardt bond-orientational order parameters.
   *
   * @param[in,out] dists DistributionFunctions target container.
   * @param[in] traj Atomic trajectory to evaluate.
   * @param[in] params TDOS evaluation parameters and MLIP model pointer.
   * @param[in] cancel_flag Optional atomic flag to abort computation prematurely.
   * @return Generated MotifProjectedTDOS structure.
   */
  static MotifProjectedTDOS correlateSteinhardt(DistributionFunctions &dists, const core::Trajectory &traj,
                                                const calculators::TDOSParams &params,
                                                const std::atomic<bool> *cancel_flag = nullptr);
};

} // namespace correlation::analysis
