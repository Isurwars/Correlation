// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "math/Precision.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace correlation::analysis {

/**
 * @brief Identifiers for 1D distributions tracked on a per-frame basis.
 */
enum class DistributionMetric : uint8_t {
  RDF = 0, ///< Radial Distribution Function g(r)
  PAD = 1, ///< Planar Angle Distribution θ
  SQ = 2   ///< Structure Factor S(Q)
};

/**
 * @brief Holds 1D curve data for a specific distribution across an entire trajectory.
 */
struct MetricCurveData {
  std::vector<real_t> bins;                      ///< X-axis bin positions (shared across frames)
  std::vector<real_t> mean_curve;                ///< Y-axis ensemble mean curve (from [min_frame, max_frame])
  std::vector<std::vector<real_t>> frame_curves; ///< Y-axis per-frame curves (indexed by frame_idx)
  std::string name;                              ///< Short metric identifier ("g(r)", "PAD", "S(Q)")
  std::string title;                             ///< Full display title
  std::string x_unit;                            ///< Physical unit for X axis ("Å", "deg", "Å⁻¹")
  std::string y_unit;                            ///< Physical unit for Y axis
};

/**
 * @brief Thread-safe repository for per-frame 1D distribution curves and their ensemble mean.
 *
 * Designed to provide ultra-low latency scrub playback in 3D viewers. Holds compact 1D curves
 * (~2.5 KB per frame) without retaining heavy 3D structures or partial matrices.
 */
class PerFrameDistributionStore {
public:
  /**
   * @brief Constructs a store pre-allocated for a specific total frame count.
   * @param[in] total_frames Total number of frames across the full trajectory.
   */
  explicit PerFrameDistributionStore(std::size_t total_frames);

  /**
   * @brief Sets the shared bin coordinates for a metric.
   * @param[in] metric Target distribution identifier.
   * @param[in] bins Bin values along the x-axis.
   */
  void setBins(DistributionMetric metric, std::vector<real_t> bins);

  /**
   * @brief Sets the ensemble mean curve for a metric.
   * @param[in] metric Target distribution identifier.
   * @param[in] mean_curve Normalized mean curve values.
   */
  void setMeanCurve(DistributionMetric metric, std::vector<real_t> mean_curve);

  /**
   * @brief Sets the 1D distribution curve for a specific trajectory frame.
   * @param[in] metric Target distribution identifier.
   * @param[in] frame_idx Zero-based trajectory frame index.
   * @param[in] curve Normalized 1D curve values for the frame.
   */
  void setFrameCurve(DistributionMetric metric, std::size_t frame_idx, std::vector<real_t> curve);

  /**
   * @brief Returns the total trajectory frame count allocated in this store.
   * @return Number of frames.
   */
  [[nodiscard]] std::size_t getTotalFrames() const noexcept;

  /**
   * @brief Checks if bins and curves are available for a given metric.
   * @param[in] metric Target distribution identifier.
   * @return True if valid bins and data exist.
   */
  [[nodiscard]] bool hasData(DistributionMetric metric) const noexcept;

  /**
   * @brief Checks if curve data is populated for a specific frame index.
   * @param[in] metric Target distribution identifier.
   * @param[in] frame_idx Zero-based trajectory frame index.
   * @return True if curve data is present.
   */
  [[nodiscard]] bool hasFrame(DistributionMetric metric, std::size_t frame_idx) const noexcept;

  /**
   * @brief Retrieves the shared bin coordinates for a metric.
   * @param[in] metric Target distribution identifier.
   * @return Const reference to bins vector.
   */
  [[nodiscard]] const std::vector<real_t> &getBins(DistributionMetric metric) const;

  /**
   * @brief Retrieves the ensemble mean curve for a metric.
   * @param[in] metric Target distribution identifier.
   * @return Const reference to mean curve vector.
   */
  [[nodiscard]] const std::vector<real_t> &getMeanCurve(DistributionMetric metric) const;

  /**
   * @brief Retrieves the curve for a specific frame index.
   * @param[in] metric Target distribution identifier.
   * @param[in] frame_idx Zero-based trajectory frame index.
   * @return Const reference to the frame's curve vector (or empty vector if unpopulated).
   */
  [[nodiscard]] const std::vector<real_t> &getFrameCurve(DistributionMetric metric,
                                                        std::size_t frame_idx) const;

  /**
   * @brief Retrieves const reference to the full metric curve dataset.
   * @param[in] metric Target distribution identifier.
   * @return Const reference to MetricCurveData.
   */
  [[nodiscard]] const MetricCurveData &getMetricData(DistributionMetric metric) const;

  /**
   * @brief Retrieves mutable reference to the full metric curve dataset.
   * @param[in] metric Target distribution identifier.
   * @return Mutable reference to MetricCurveData.
   */
  [[nodiscard]] MetricCurveData &getMetricData(DistributionMetric metric);

private:
  static constexpr std::size_t kMetricCount = 3;

  [[nodiscard]] static std::size_t toIndex(DistributionMetric metric) noexcept;

  std::size_t total_frames_{0};
  std::array<MetricCurveData, kMetricCount> metrics_;
  static const std::vector<real_t> empty_curve_;
};

} // namespace correlation::analysis

