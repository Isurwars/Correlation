/**
 * @file RDCalculator.cpp
 * @brief Implementation of the ring distribution calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/spatial/RDCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "calculators/spatial/MotifFinder.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace correlation::calculators {

namespace {

// Static registration of the calculator in the factory
const bool REGISTERED = CalculatorFactory::registerTypeSafe<RDCalculator>("RDCalculator");

size_t countNetworkFormerAtoms(const correlation::core::Cell &cell, size_t fallback_count,
                               std::string_view former_element) {
  if (former_element.empty() || cell.atomCount() == 0) {
    if (fallback_count > 0) {
      return fallback_count;
    }
    return 1;
  }
  size_t count = 0;
  for (const auto &atom : cell.atoms()) {
    if (atom.element().symbol == former_element) {
      count++;
    }
  }
  if (count > 0) {
    return count;
  }
  if (fallback_count > 0) {
    return fallback_count;
  }
  return 1;
}

std::vector<std::vector<correlation::core::AtomID>>
extractRingCycles(const correlation::core::NeighborGraph &graph,
                  const correlation::core::Cell &cell, const RDCalculator::RDParams &params) {
  if (params.projection_mode == correlation::analysis::RingProjectionMode::BridgedProjection &&
      !params.network_former.empty() && !params.bridging_element.empty()) {
    auto const bridged =
        MotifFinder::buildBridgedGraph(graph, cell, params.network_former, params.bridging_element);
    return MotifFinder::extractAllCycles(bridged, params.max_ring_size, params.ring_type);
  }

  if (params.projection_mode == correlation::analysis::RingProjectionMode::AlternatingTracing &&
      !params.network_former.empty() && !params.bridging_element.empty()) {
    size_t const max_atom_ring =
        params.report_polyhedra_size ? (2 * params.max_ring_size) : params.max_ring_size;
    auto const all_cycles = MotifFinder::extractAllCycles(graph, max_atom_ring, params.ring_type);
    return MotifFinder::filterAlternatingCycles(all_cycles, cell, params.network_former,
                                                params.bridging_element);
  }

  return MotifFinder::extractAllCycles(graph, params.max_ring_size, params.ring_type);
}

struct RingAccumulationConfig {
  size_t num_bins{0};
  size_t min_ring_size{0};
  size_t total_former_atoms{0};
};

struct RingAccumulator {
  std::vector<real_t> raw_counts;
  std::vector<real_t> node_frequency;
  std::vector<real_t> rings_per_former;
  std::vector<real_t> fraction;
  size_t total_counts{0};
};

RingAccumulator
accumulateRingStats(const std::vector<std::vector<correlation::core::AtomID>> &cycles,
                    const correlation::core::Cell &cell, const RDCalculator::RDParams &params,
                    const RingAccumulationConfig &config) {
  RingAccumulator acc{
      .raw_counts = std::vector<real_t>(config.num_bins, 0.0),
      .node_frequency = std::vector<real_t>(config.num_bins, 0.0),
      .rings_per_former = std::vector<real_t>(config.num_bins, 0.0),
      .fraction = std::vector<real_t>(config.num_bins, 0.0),
  };

  const auto &atoms = cell.atoms();
  const size_t num_atoms = atoms.size();

  for (const auto &cycle : cycles) {
    size_t const eff_size =
        (params.projection_mode == correlation::analysis::RingProjectionMode::AlternatingTracing &&
         params.report_polyhedra_size)
            ? (cycle.size() / 2)
            : cycle.size();

    if (eff_size < config.min_ring_size || eff_size > params.max_ring_size) {
      continue;
    }

    size_t const bin = eff_size - config.min_ring_size;
    acc.raw_counts[bin] += 1.0;
    acc.total_counts++;

    size_t former_in_ring = 0;
    for (auto const atom_id : cycle) {
      if (params.network_former.empty() ||
          (atom_id < num_atoms && atoms[atom_id].element().symbol == params.network_former)) {
        former_in_ring++;
      }
    }
    if (former_in_ring == 0) {
      former_in_ring = eff_size;
    }
    acc.node_frequency[bin] += static_cast<real_t>(former_in_ring);
  }

  auto const n_former_f = static_cast<real_t>(config.total_former_atoms);
  auto const total_counts_f = static_cast<real_t>(acc.total_counts);

  for (size_t i = 0; i < config.num_bins; ++i) {
    if (acc.total_counts > 0) {
      acc.fraction[i] = acc.raw_counts[i] / total_counts_f;
    }
    acc.rings_per_former[i] = acc.raw_counts[i] / n_former_f;
    acc.node_frequency[i] /= n_former_f;
  }

  return acc;
}

} // namespace

void RDCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                  const correlation::analysis::AnalysisSettings &settings) const {
  if (dists.neighbors() == nullptr) {
    return;
  }
  RDParams const params{
      .max_ring_size = settings.max_ring_size,
      .ring_type = settings.ring_type,
      .projection_mode = settings.ring_projection_mode,
      .network_former = settings.ring_network_former,
      .bridging_element = settings.ring_bridging_element,
      .report_polyhedra_size = true,
  };
  dists.addHistogram("RD", calculate(dists.neighbors()->neighborGraph(), dists.cell(), params));
}

correlation::analysis::Histogram
RDCalculator::calculate(const correlation::core::NeighborGraph &graph, size_t max_ring_size,
                        correlation::analysis::RingType ring_type) {
  correlation::core::Cell const dummy_cell;
  RDParams const params{
      .max_ring_size = max_ring_size,
      .ring_type = ring_type,
      .projection_mode = correlation::analysis::RingProjectionMode::Direct,
      .network_former = "",
      .bridging_element = "",
      .report_polyhedra_size = true,
  };
  return calculate(graph, dummy_cell, params);
}

correlation::analysis::Histogram
RDCalculator::calculate(const correlation::core::NeighborGraph &graph,
                        const correlation::core::Cell &cell, const RDParams &params) {
  if (params.max_ring_size < 3) {
    throw std::invalid_argument("Max ring size must be at least 3");
  }

  size_t const min_ring_size =
      (params.projection_mode == correlation::analysis::RingProjectionMode::AlternatingTracing &&
       params.report_polyhedra_size)
          ? 2
          : 3;

  size_t const num_bins =
      (params.max_ring_size >= min_ring_size) ? (params.max_ring_size - min_ring_size + 1) : 0;

  bool const is_polyhedra =
      (params.projection_mode == correlation::analysis::RingProjectionMode::BridgedProjection ||
       (params.projection_mode == correlation::analysis::RingProjectionMode::AlternatingTracing &&
        params.report_polyhedra_size));

  correlation::analysis::Histogram f_motif;
  f_motif.x_label = "Ring Size";
  f_motif.title = (params.ring_type == correlation::analysis::RingType::Franzblau)
                      ? "Franzblau Primitive Ring Distribution"
                      : "Ring Distribution";
  f_motif.y_label = "Frequency";
  f_motif.x_unit = is_polyhedra ? "polyhedra" : "atoms";
  f_motif.y_unit = "fraction";
  f_motif.description = (params.ring_type == correlation::analysis::RingType::Franzblau)
                            ? "Franzblau Shortest-Path Ring Distribution"
                            : "Ring Distribution";
  f_motif.file_suffix = "_RD";

  f_motif.bins.resize(num_bins);
  for (size_t i = 0; i < num_bins; ++i) {
    f_motif.bins[i] = static_cast<real_t>(i + min_ring_size);
  }

  auto const cycles = extractRingCycles(graph, cell, params);
  size_t const total_former =
      countNetworkFormerAtoms(cell, graph.nodeCount(), params.network_former);
  auto const acc = accumulateRingStats(cycles, cell, params,
                                       RingAccumulationConfig{
                                           .num_bins = num_bins,
                                           .min_ring_size = min_ring_size,
                                           .total_former_atoms = total_former,
                                       });

  f_motif.partials["Rings"] = acc.fraction;
  f_motif.partials["Total"] = acc.fraction;
  f_motif.partials["Fraction"] = acc.fraction;
  f_motif.partials["RawCounts"] = acc.raw_counts;
  f_motif.partials["RingsPerFormer"] = acc.rings_per_former;
  f_motif.partials["NodeFrequency"] = acc.node_frequency;

  return f_motif;
}

} // namespace correlation::calculators
