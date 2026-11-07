// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/PerFrameDistributionStore.hpp"

#include <stdexcept>
#include <utility>

namespace correlation::analysis {

const std::vector<real_t> PerFrameDistributionStore::empty_curve_{};

PerFrameDistributionStore::PerFrameDistributionStore(std::size_t total_frames)
    : total_frames_(total_frames) {
  // Initialize descriptors for the 3 supported 1D metrics
  metrics_[toIndex(DistributionMetric::RDF)] = {
      .bins = {},
      .mean_curve = {},
      .frame_curves = std::vector<std::vector<real_t>>(total_frames_),
      .name = "g(r)",
      .title = "Radial Distribution Function",
      .x_unit = "Å",
      .y_unit = ""};

  metrics_[toIndex(DistributionMetric::PAD)] = {
      .bins = {},
      .mean_curve = {},
      .frame_curves = std::vector<std::vector<real_t>>(total_frames_),
      .name = "PAD",
      .title = "Planar Angle Distribution",
      .x_unit = "deg",
      .y_unit = "deg⁻¹"};

  metrics_[toIndex(DistributionMetric::SQ)] = {
      .bins = {},
      .mean_curve = {},
      .frame_curves = std::vector<std::vector<real_t>>(total_frames_),
      .name = "S(Q)",
      .title = "Structure Factor",
      .x_unit = "Å⁻¹",
      .y_unit = ""};
}

std::size_t PerFrameDistributionStore::toIndex(DistributionMetric metric) noexcept {
  return static_cast<std::size_t>(metric);
}

void PerFrameDistributionStore::setBins(DistributionMetric metric, std::vector<real_t> bins) {
  metrics_[toIndex(metric)].bins = std::move(bins);
}

void PerFrameDistributionStore::setMeanCurve(DistributionMetric metric,
                                             std::vector<real_t> mean_curve) {
  metrics_[toIndex(metric)].mean_curve = std::move(mean_curve);
}

void PerFrameDistributionStore::setFrameCurve(DistributionMetric metric, std::size_t frame_idx,
                                              std::vector<real_t> curve) {
  auto &frame_curves = metrics_[toIndex(metric)].frame_curves;
  if (frame_idx < frame_curves.size()) {
    frame_curves[frame_idx] = std::move(curve);
  }
}

std::size_t PerFrameDistributionStore::getTotalFrames() const noexcept {
  return total_frames_;
}

bool PerFrameDistributionStore::hasData(DistributionMetric metric) const noexcept {
  const auto &data = metrics_[toIndex(metric)];
  return !data.bins.empty() && (!data.mean_curve.empty() || !data.frame_curves.empty());
}

bool PerFrameDistributionStore::hasFrame(DistributionMetric metric,
                                         std::size_t frame_idx) const noexcept {
  const auto &frame_curves = metrics_[toIndex(metric)].frame_curves;
  if (frame_idx >= frame_curves.size()) {
    return false;
  }
  return !frame_curves[frame_idx].empty();
}

const std::vector<real_t> &
PerFrameDistributionStore::getBins(DistributionMetric metric) const {
  return metrics_[toIndex(metric)].bins;
}

const std::vector<real_t> &
PerFrameDistributionStore::getMeanCurve(DistributionMetric metric) const {
  return metrics_[toIndex(metric)].mean_curve;
}

const std::vector<real_t> &
PerFrameDistributionStore::getFrameCurve(DistributionMetric metric,
                                         std::size_t frame_idx) const {
  const auto &frame_curves = metrics_[toIndex(metric)].frame_curves;
  if (frame_idx >= frame_curves.size() || frame_curves[frame_idx].empty()) {
    return empty_curve_;
  }
  return frame_curves[frame_idx];
}

const MetricCurveData &
PerFrameDistributionStore::getMetricData(DistributionMetric metric) const {
  return metrics_[toIndex(metric)];
}

MetricCurveData &PerFrameDistributionStore::getMetricData(DistributionMetric metric) {
  return metrics_[toIndex(metric)];
}

} // namespace correlation::analysis

