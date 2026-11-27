// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "math/Precision.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <mutex>
#include <shared_mutex>
#include <string>
#include <vector>

namespace correlation::analysis {

/**
 * @brief Identifiers for 1D distributions tracked on a per-frame basis.
 */
enum class DistributionMetric : uint8_t {
  RDF = 0, ///< Radial Distribution Function g(r)
  PAD = 1, ///< Planar Angle Distribution θ
  SQ = 2,  ///< Structure Factor S(Q)
  XRD = 3  ///< X-Ray Diffraction Pattern I(2θ)
};

/**
 * @brief Holds 1D curve data for a specific distribution across an entire trajectory.
 */
struct MetricCurveData {
  std::vector<real_t> bins;       ///< X-axis bin positions (shared across frames)
  std::vector<real_t> mean_curve; ///< Y-axis ensemble mean curve (from [min_frame, max_frame])
  std::vector<std::vector<real_t>> frame_curves; ///< Y-axis per-frame curves (indexed by frame_idx)
  std::string name;   ///< Short metric identifier ("g(r)", "PAD", "S(Q)")
  std::string title;  ///< Full display title
  std::string x_unit; ///< Physical unit for X axis ("Å", "deg", "Å⁻¹")
  std::string y_unit; ///< Physical unit for Y axis
};

/**
 * @brief Thread-safe repository for per-frame 1D distribution curves and their ensemble mean.
 *
 * Designed to provide ultra-low latency scrub playback in 3D viewers. Holds compact 1D curves
 * (~2.5 KB per frame) without retaining heavy 3D structures or partial matrices.
 */
class PerFrameDistributionStore {
public:
  static constexpr std::size_t METRIC_COUNT = 4;
  using FrameCurveArray = std::array<std::vector<real_t>, METRIC_COUNT>;
  using FrameCurveGenerator = std::function<FrameCurveArray(std::size_t)>;

  struct MeanFrameWindow {
    std::size_t start_frame{0};
    std::size_t sample_count{0};
    std::size_t stride{1};
  };

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

  /**
   * @brief Records the strided [start_frame, start_frame + sample_count * stride) mean window.
   * @param[in] window Sampled frame window parameters accumulated into the mean.
   */
  void configureMeanFrames(MeanFrameWindow window);

  /**
   * @brief Explicitly marks whether a single frame contributed to the ensemble mean.
   * @param[in] frame_idx Zero-based trajectory frame index.
   * @param[in] in_mean True if included in ensemble mean.
   */
  void markFrameInMean(std::size_t frame_idx, bool in_mean = true);

  /**
   * @brief Queries whether a trajectory frame contributed to the ensemble mean.
   * @param[in] frame_idx Zero-based trajectory frame index.
   * @return True if frame was included in ensemble mean (or if no explicit mask was configured).
   */
  [[nodiscard]] bool isFrameInMean(std::size_t frame_idx) const noexcept;

  /**
   * @brief Installs an on-demand 1D curve generator for frames outside the initial mean range.
   * @param[in] generator Callback computing [RDF, PAD, SQ, XRD] curves for a given frame index.
   */
  void setFrameGenerator(FrameCurveGenerator generator);

  /**
   * @brief Detaches any active lazy frame generator.
   */
  void clearFrameGenerator() noexcept;

  /**
   * @brief Checks whether any metric with a populated mean curve is missing a curve for frame_idx.
   * @param[in] frame_idx Zero-based trajectory frame index.
   * @return True if at least one active metric curve is still unpopulated for frame_idx.
   */
  [[nodiscard]] bool needsFrameGeneration(std::size_t frame_idx) const noexcept;

  /**
   * @brief Lazily computes and caches missing 1D curves for frame_idx if a generator is attached.
   * @param[in] frame_idx Zero-based trajectory frame index.
   */
  void ensureFrame(std::size_t frame_idx) const;

private:
  [[nodiscard]] static std::size_t toIndex(DistributionMetric metric) noexcept;
  [[nodiscard]] bool needsFrameGenerationUnlocked(std::size_t frame_idx) const noexcept;

  mutable std::shared_mutex mutex_;
  mutable std::mutex generation_mutex_;
  std::size_t total_frames_{0};
  mutable std::array<MetricCurveData, METRIC_COUNT> metrics_;
  bool has_explicit_mean_mask_{false};
  std::vector<uint8_t> in_mean_mask_;
  FrameCurveGenerator frame_generator_;
  static const std::vector<real_t> EMPTY_CURVE;
};

} // namespace correlation::analysis
