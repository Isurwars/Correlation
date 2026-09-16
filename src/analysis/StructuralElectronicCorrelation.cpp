/**
 * @file StructuralElectronicCorrelation.cpp
 * @brief Implementation of structural-electronic correlation pipelines.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "analysis/StructuralElectronicCorrelation.hpp"
#include "calculators/SteinhardtCalculator.hpp"
#include "math/Constants.hpp"
#include "mlip/GraphDescriptors.hpp"
#include "mlip/PeriodicGraphBuilder.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <string_view>
#include <vector>

namespace correlation::analysis {

namespace {

std::vector<real_t> buildEnergyGrid(real_t e_min, real_t e_max, size_t num_bins) {
  if (num_bins == 0) {
    return {};
  }
  std::vector<real_t> energies(num_bins, 0.0);
  const real_t delta = (e_max - e_min) / static_cast<real_t>(num_bins);
  for (size_t i = 0; i < num_bins; ++i) {
    energies[i] = e_min + (static_cast<real_t>(i) + static_cast<real_t>(0.5)) * delta;
  }
  return energies;
}

std::string_view cnaLabelToString(int label) noexcept {
  switch (static_cast<correlation::mlip::CNALabel>(label)) {
  case correlation::mlip::CNALabel::FCC:
    return "FCC";
  case correlation::mlip::CNALabel::HCP:
    return "HCP";
  case correlation::mlip::CNALabel::BCC:
    return "BCC";
  case correlation::mlip::CNALabel::ICO:
    return "ICO";
  case correlation::mlip::CNALabel::Other:
  default:
    return "Other";
  }
}

std::string_view classifySteinhardt(real_t q6_val) noexcept {
  if (q6_val < static_cast<real_t>(0.35)) {
    return "Disordered";
  }
  if (q6_val >= static_cast<real_t>(0.5)) {
    return "FCC-like";
  }
  if (q6_val >= static_cast<real_t>(0.4)) {
    return "BCC-like";
  }
  return "Ordered";
}

void accumulateVector(const std::vector<real_t> &src, std::vector<real_t> &dst) {
  const size_t num_elements = src.size();
  if (dst.size() < num_elements) {
    dst.resize(num_elements, 0.0);
  }
  for (size_t i = 0; i < num_elements; ++i) {
    dst[i] += src[i];
  }
}

void normalizeResult(MotifProjectedTDOS &result, size_t total_atoms) {
  if (total_atoms == 0) {
    return;
  }
  const auto inv_atoms = static_cast<real_t>(1.0) / static_cast<real_t>(total_atoms);
  for (auto &val : result.total_tdos) {
    val *= inv_atoms;
  }
  for (auto &pair : result.motif_tdos) {
    for (auto &val : pair.second) {
      val *= inv_atoms;
    }
  }
}

std::vector<std::vector<size_t>>
buildAtomEdgeLists(const correlation::mlip::PeriodicGraphData &graph) {
  std::vector<std::vector<size_t>> atom_edges(graph.atom_count);
  for (size_t edge_idx = 0; edge_idx < graph.edge_count; ++edge_idx) {
    const auto src = static_cast<size_t>(graph.edge_index_flat[edge_idx]);
    if (src < graph.atom_count) {
      atom_edges[src].push_back(edge_idx);
    }
  }
  return atom_edges;
}

real_t computeAtomQ6(size_t atom_idx, const correlation::mlip::PeriodicGraphData &graph,
                     const std::vector<std::vector<size_t>> &atom_edges) {
  const auto &edge_indices = atom_edges[atom_idx];
  if (edge_indices.size() < 2) {
    return static_cast<real_t>(0.0);
  }

  std::array<std::complex<real_t>, 13> q6m{};
  for (const size_t edge_idx : edge_indices) {
    const real_t delta_x = graph.edge_vectors_flat[edge_idx * 3];
    const real_t delta_y = graph.edge_vectors_flat[edge_idx * 3 + 1];
    const real_t delta_z = graph.edge_vectors_flat[edge_idx * 3 + 2];
    const real_t dist = graph.edge_distances.empty()
                            ? std::sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z)
                            : graph.edge_distances[edge_idx];
    if (dist <= static_cast<real_t>(0.0)) {
      continue;
    }

    const real_t theta =
        std::acos(std::clamp(delta_z / dist, static_cast<real_t>(-1.0), static_cast<real_t>(1.0)));
    const real_t phi = std::atan2(delta_y, delta_x);
    for (size_t array_idx = 0; array_idx < 13; ++array_idx) {
      const int m_idx = static_cast<int>(array_idx) - 6;
      q6m.at(array_idx) += calculators::SteinhardtCalculator::sphericalHarmonic(
          6, m_idx, calculators::SteinhardtCalculator::SphericalAngles{.theta = theta, .phi = phi});
    }
  }

  const auto inv_bonds = static_cast<real_t>(1.0) / static_cast<real_t>(edge_indices.size());
  real_t q6_sq = 0.0;
  for (size_t array_idx = 0; array_idx < 13; ++array_idx) {
    const auto q_val = q6m.at(array_idx) * inv_bonds;
    q6_sq += std::norm(q_val);
  }
  return std::sqrt((correlation::math::four_pi / static_cast<real_t>(13.0)) * q6_sq);
}

} // anonymous namespace

Histogram MotifProjectedTDOS::toHistogram(std::string_view title) const {
  Histogram hist;
  hist.title = std::string(title);
  hist.x_label = "Energy";
  hist.y_label = "Density of States";
  hist.x_unit = "eV";
  hist.y_unit = "states/eV/atom";
  hist.bins = energies;
  hist.partials = motif_tdos;
  hist.partials["total"] = total_tdos;
  hist.compute_count = static_cast<int>(frame_count);
  return hist;
}

MotifProjectedTDOS StructuralElectronicCorrelation::correlateCNA(
    DistributionFunctions &dists, const core::Trajectory &traj,
    const calculators::TDOSParams &params, const std::atomic<bool> *cancel_flag) {
  MotifProjectedTDOS result;
  const size_t num_frames = traj.getFrameCount();
  if (params.model == nullptr || num_frames == 0) {
    return result;
  }

  size_t total_atoms = 0;
  for (size_t frame_idx = 0; frame_idx < num_frames; ++frame_idx) {
    if (cancel_flag != nullptr && cancel_flag->load(std::memory_order_relaxed)) {
      return {};
    }

    const auto cell = traj.getFrame(frame_idx);
    const auto mlip_out = params.model->evaluate(cell);
    if (mlip_out.ldos_bins == 0 || mlip_out.ldos.empty()) {
      continue;
    }

    if (result.energies.empty()) {
      result.energies = buildEnergyGrid(params.e_min, params.e_max, mlip_out.ldos_bins);
      result.total_tdos.resize(mlip_out.ldos_bins, 0.0);
    }

    const auto graph = correlation::mlip::PeriodicGraphBuilder::buildGraph(cell);
    const auto cna_labels = correlation::mlip::GraphDescriptors::computeCNADescriptor(graph);

    const size_t n_atoms = cell.atomCount();
    for (size_t i = 0; i < n_atoms && i < mlip_out.ldos.size(); ++i) {
      const auto motif_name =
          std::string(cnaLabelToString(i < cna_labels.size() ? cna_labels[i] : 0));
      accumulateVector(mlip_out.ldos[i], result.motif_tdos[motif_name]);
      accumulateVector(mlip_out.ldos[i], result.total_tdos);
      ++total_atoms;
    }
    ++result.frame_count;
  }

  normalizeResult(result, total_atoms);
  dists.addHistogram("MotifProjectedTDOS_CNA", result.toHistogram("Motif-Projected TDOS (CNA)"));
  return result;
}

MotifProjectedTDOS StructuralElectronicCorrelation::correlateSteinhardt(
    DistributionFunctions &dists, const core::Trajectory &traj,
    const calculators::TDOSParams &params, const std::atomic<bool> *cancel_flag) {
  MotifProjectedTDOS result;
  const size_t num_frames = traj.getFrameCount();
  if (params.model == nullptr || num_frames == 0) {
    return result;
  }

  size_t total_atoms = 0;
  for (size_t frame_idx = 0; frame_idx < num_frames; ++frame_idx) {
    if (cancel_flag != nullptr && cancel_flag->load(std::memory_order_relaxed)) {
      return {};
    }

    const auto cell = traj.getFrame(frame_idx);
    const auto mlip_out = params.model->evaluate(cell);
    if (mlip_out.ldos_bins == 0 || mlip_out.ldos.empty()) {
      continue;
    }

    if (result.energies.empty()) {
      result.energies = buildEnergyGrid(params.e_min, params.e_max, mlip_out.ldos_bins);
      result.total_tdos.resize(mlip_out.ldos_bins, 0.0);
    }

    const auto graph = correlation::mlip::PeriodicGraphBuilder::buildGraph(cell);
    const auto atom_edges = buildAtomEdgeLists(graph);

    const size_t n_atoms = cell.atomCount();
    for (size_t i = 0; i < n_atoms && i < mlip_out.ldos.size(); ++i) {
      const real_t q6_val = computeAtomQ6(i, graph, atom_edges);
      const auto motif_name = std::string(classifySteinhardt(q6_val));
      accumulateVector(mlip_out.ldos[i], result.motif_tdos[motif_name]);
      accumulateVector(mlip_out.ldos[i], result.total_tdos);
      ++total_atoms;
    }
    ++result.frame_count;
  }

  normalizeResult(result, total_atoms);
  dists.addHistogram("MotifProjectedTDOS_Steinhardt",
                     result.toHistogram("Motif-Projected TDOS (Steinhardt)"));
  return result;
}

} // namespace correlation::analysis
