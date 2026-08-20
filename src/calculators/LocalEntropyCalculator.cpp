/**
 * @file LocalEntropyCalculator.cpp
 * @brief Implementation of the Local Entropic Fingerprint order parameter.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/LocalEntropyCalculator.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Constants.hpp"
#include "math/LinearAlgebra.hpp"
#include "math/Precision.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
const bool registered = CalculatorFactory::registerTypeSafe<LocalEntropyCalculator>("LocalEntropyCalculator");

struct SearchGridConfig {
  int K_x = 1;
  int K_y = 1;
  int K_z = 1;
  int max_dx = 1;
  int max_dy = 1;
  int max_dz = 1;
};

[[nodiscard]] inline real_t wrapUnit(real_t val) noexcept {
  real_t const res = std::fmod(val, static_cast<real_t>(1.0));
  return (res < 0.0) ? res + static_cast<real_t>(1.0) : res;
}

[[nodiscard]] inline correlation::math::Vector3<real_t>
wrapFractional(const correlation::math::Vector3<real_t> &frac) noexcept {
  return {wrapUnit(frac.x()), wrapUnit(frac.y()), wrapUnit(frac.z())};
}

struct PeriodicShift {
  int wrap = 0;
  int shift = 0;
};

[[nodiscard]] inline PeriodicShift computePeriodicShift(int index, int grid_k) noexcept {
  int const shift = (index >= 0) ? (index / grid_k) : ((index - grid_k + 1) / grid_k);
  int const wrap = index - shift * grid_k;
  return {.wrap = wrap, .shift = shift};
}

/**
 * @brief Spatial cell-grid structure for dedicated geometric neighbor distance queries.
 */
struct MetricNeighborGrid {
  SearchGridConfig grid;
  std::vector<int> atom_bin;
  std::vector<size_t> bin_offsets;
  std::vector<size_t> bin_indices;
  std::vector<correlation::math::Vector3<real_t>> wrapped_positions;

  struct QueryContext {
    size_t atom_idx = 0;
    correlation::math::Vector3<real_t> pos_i;
    real_t cutoff_sq = 0.0;
    bool ignore_self_periodic = false;
  };

  void build(const correlation::core::Cell &cell, real_t cutoff) {
    const size_t num_atoms = cell.atomCount();
    atom_bin.resize(num_atoms);
    wrapped_positions.resize(num_atoms);
    if (num_atoms == 0 || cutoff <= 0.0) {
      return;
    }

    const auto &lattice = cell.latticeVectors();
    const auto &inv_lattice = cell.inverseLatticeVectors();

    real_t const vol = cell.volume();
    const correlation::math::Vector3<real_t> v_a = {lattice(0, 0), lattice(0, 1), lattice(0, 2)};
    const correlation::math::Vector3<real_t> v_b = {lattice(1, 0), lattice(1, 1), lattice(1, 2)};
    const correlation::math::Vector3<real_t> v_c = {lattice(2, 0), lattice(2, 1), lattice(2, 2)};

    const correlation::math::Vector3<real_t> cross_bc = correlation::math::cross(v_b, v_c);
    const correlation::math::Vector3<real_t> cross_ca = correlation::math::cross(v_c, v_a);
    const correlation::math::Vector3<real_t> cross_ab = correlation::math::cross(v_a, v_b);

    real_t const h_x = vol / correlation::math::norm(cross_bc);
    real_t const h_y = vol / correlation::math::norm(cross_ca);
    real_t const h_z = vol / correlation::math::norm(cross_ab);

    grid.K_x = std::max(1, static_cast<int>(std::floor(h_x / cutoff)));
    grid.K_y = std::max(1, static_cast<int>(std::floor(h_y / cutoff)));
    grid.K_z = std::max(1, static_cast<int>(std::floor(h_z / cutoff)));

    grid.max_dx = static_cast<int>(std::ceil(cutoff / (h_x / static_cast<real_t>(grid.K_x))));
    grid.max_dy = static_cast<int>(std::ceil(cutoff / (h_y / static_cast<real_t>(grid.K_y))));
    grid.max_dz = static_cast<int>(std::ceil(cutoff / (h_z / static_cast<real_t>(grid.K_z))));

    const size_t total_bins =
        static_cast<size_t>(grid.K_x) * static_cast<size_t>(grid.K_y) * static_cast<size_t>(grid.K_z);
    std::vector<size_t> bin_counts(total_bins, 0);

    const auto &atoms = cell.atoms();
    for (size_t i = 0; i < num_atoms; ++i) {
      correlation::math::Vector3<real_t> const frac = wrapFractional(inv_lattice * atoms[i].position());
      wrapped_positions[i] = lattice * frac;

      int const b_x = std::clamp(static_cast<int>(frac.x() * static_cast<real_t>(grid.K_x)), 0, grid.K_x - 1);
      int const b_y = std::clamp(static_cast<int>(frac.y() * static_cast<real_t>(grid.K_y)), 0, grid.K_y - 1);
      int const b_z = std::clamp(static_cast<int>(frac.z() * static_cast<real_t>(grid.K_z)), 0, grid.K_z - 1);

      int const bin_idx = b_x * (grid.K_y * grid.K_z) + b_y * grid.K_z + b_z;
      atom_bin[i] = bin_idx;
      bin_counts[bin_idx]++;
    }

    bin_offsets.resize(total_bins + 1, 0);
    for (size_t b_idx = 0; b_idx < total_bins; ++b_idx) {
      bin_offsets[b_idx + 1] = bin_offsets[b_idx] + bin_counts[b_idx];
    }

    bin_indices.resize(num_atoms);
    std::vector<size_t> current_offsets = bin_offsets;
    for (size_t atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
      int const bin_idx = atom_bin[atom_idx];
      bin_indices[current_offsets[bin_idx]++] = atom_idx;
    }
  }

  void scanBinDistances(const QueryContext &ctx, int bin_idx, const correlation::math::Vector3<real_t> &disp,
                        bool zero_disp, std::vector<real_t> &distances_out) const {
    size_t const start = bin_offsets[bin_idx];
    size_t const end = bin_offsets[bin_idx + 1];

    for (size_t idx = start; idx < end; ++idx) {
      size_t const j_atom = bin_indices[idx];
      if (ctx.atom_idx == j_atom && (zero_disp || ctx.ignore_self_periodic)) {
        continue;
      }

      correlation::math::Vector3<real_t> const wrapped_pos_j = wrapped_positions[j_atom] + disp;
      correlation::math::Vector3<real_t> const r_vec = wrapped_pos_j - ctx.pos_i;
      real_t const d_sq = r_vec.x() * r_vec.x() + r_vec.y() * r_vec.y() + r_vec.z() * r_vec.z();

      if (d_sq <= ctx.cutoff_sq && d_sq > static_cast<real_t>(1e-12)) {
        distances_out.push_back(std::sqrt(d_sq));
      }
    }
  }

  void queryDistances(size_t atom_idx, const correlation::core::Cell &cell, real_t cutoff,
                      std::vector<real_t> &distances_out, bool ignore_self_periodic) const {
    distances_out.clear();
    const auto &lattice = cell.latticeVectors();
    QueryContext const ctx = {
        .atom_idx = atom_idx,
        .pos_i = wrapped_positions[atom_idx],
        .cutoff_sq = cutoff * cutoff,
        .ignore_self_periodic = ignore_self_periodic,
    };

    int const center_x = atom_bin[atom_idx] / (grid.K_y * grid.K_z);
    int const center_y = (atom_bin[atom_idx] / grid.K_z) % grid.K_y;
    int const center_z = atom_bin[atom_idx] % grid.K_z;

    for (int dx = -grid.max_dx; dx <= grid.max_dx; ++dx) {
      auto const [wrap_x, shift_x] = computePeriodicShift(center_x + dx, grid.K_x);

      for (int dy = -grid.max_dy; dy <= grid.max_dy; ++dy) {
        auto const [wrap_y, shift_y] = computePeriodicShift(center_y + dy, grid.K_y);

        for (int dz = -grid.max_dz; dz <= grid.max_dz; ++dz) {
          auto const [wrap_z, shift_z] = computePeriodicShift(center_z + dz, grid.K_z);

          int const n_bin = wrap_x * (grid.K_y * grid.K_z) + wrap_y * grid.K_z + wrap_z;
          correlation::math::Vector3<real_t> const disp =
              lattice * correlation::math::Vector3<real_t>(static_cast<real_t>(shift_x), static_cast<real_t>(shift_y),
                                                           static_cast<real_t>(shift_z));
          bool const zero_disp = (shift_x == 0 && shift_y == 0 && shift_z == 0);

          scanBinDistances(ctx, n_bin, disp, zero_disp, distances_out);
        }
      }
    }
  }
};

real_t computeSingleAtomEntropy(const std::vector<real_t> &atom_distances, const correlation::core::Cell &cell,
                                const LocalEntropyParams &params) {
  real_t const cutoff = params.cutoff;
  real_t const sigma = params.sigma;

  real_t const vol = cell.volume();
  if (vol <= 0.0) {
    return 0.0;
  }
  real_t const density = static_cast<real_t>(cell.atomCount()) / vol;

  // We integrate from r = dr to Rc with step size dr
  real_t const dr_val = 0.02;
  auto const num_steps = static_cast<size_t>(std::ceil(cutoff / dr_val));

  real_t integral = 0.0;

  // Gaussian prefactor: 1 / sqrt(2 * pi * sigma^2)
  real_t const gaussian_prefactor =
      static_cast<real_t>(1.0) / (sigma * std::sqrt(static_cast<real_t>(2.0) * correlation::math::pi));

  const real_t two_sigma_sq = static_cast<real_t>(2.0) * sigma * sigma;

  for (size_t step = 1; step <= num_steps; ++step) {
    real_t const r_val = static_cast<real_t>(step) * dr_val;
    if (r_val > cutoff) {
      break;
    }

    real_t g_val = 0.0;
    for (const real_t r_ij : atom_distances) {
      // Gaussian kernel with periodic boundary corrections at r=0
      real_t const diff = r_val - r_ij;
      real_t const sum_r = r_val + r_ij;
      auto const term1 = static_cast<real_t>(std::exp(-(diff * diff) / two_sigma_sq));
      auto const term2 = static_cast<real_t>(std::exp(-(sum_r * sum_r) / two_sigma_sq));
      g_val += term1 + term2;
    }

    g_val *= gaussian_prefactor;

    // Normalization: 4 * pi * density * r^2
    real_t const norm = static_cast<real_t>(4.0) * correlation::math::pi * density * r_val * r_val;
    if (norm > 0.0) {
      g_val /= norm;
    } else {
      g_val = 0.0;
    }

    real_t integrand = 0.0;
    if (g_val > 1e-10) {
      integrand = (g_val * std::log(g_val) - g_val + static_cast<real_t>(1.0)) * r_val * r_val;
    } else {
      integrand = r_val * r_val;
    }

    // Trapezoidal rule integration
    real_t weight = dr_val;
    if (step == num_steps) {
      weight = static_cast<real_t>(0.5) * dr_val;
    }
    integral += integrand * weight;
  }

  return static_cast<real_t>(-2.0) * correlation::math::pi * density * integral;
}

struct BinningConfig {
  real_t min_val;
  real_t max_val;
  real_t d_val;
};

void initHistogramMap(std::map<std::string, std::vector<real_t>> &partials,
                      const std::vector<std::string> &element_symbols, size_t bins) {
  for (const auto &sym : element_symbols) {
    partials[sym].assign(bins, 0.0);
  }
  partials["Total"].assign(bins, 0.0);
}

void addValueToHistogram(std::map<std::string, std::vector<real_t>> &partials, const std::string &symbol, real_t val,
                         BinningConfig config) {
  if (val >= config.min_val && val < config.max_val) {
    auto const bin_idx = static_cast<size_t>((val - config.min_val) / config.d_val);
    partials[symbol][bin_idx] += 1.0;
    partials["Total"][bin_idx] += 1.0;
  }
}

void accumulateHistogramMap(std::map<std::string, std::vector<real_t>> &dest,
                            const std::map<std::string, std::vector<real_t>> &src) {
  for (const auto &[key, vec] : src) {
    auto &dest_vec = dest[key];
    for (size_t bin_idx = 0; bin_idx < vec.size(); ++bin_idx) {
      dest_vec[bin_idx] += vec[bin_idx];
    }
  }
}

void normalizeHistogramMap(std::map<std::string, std::vector<real_t>> &partials, real_t factor) {
  if (factor <= 0.0) {
    return;
  }
  for (auto &[key, vec] : partials) {
    for (auto &val : vec) {
      val /= factor;
    }
  }
}

void copyPartialsToHistogram(correlation::analysis::Histogram &hist,
                             const std::map<std::string, std::vector<real_t>> &partials) {
  for (const auto &[key, vec] : partials) {
    hist.partials[key] = vec;
  }
}

} // namespace

void LocalEntropyCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                            const correlation::analysis::AnalysisSettings &settings) const {
  auto hist = calculate(dists.cell(), dists.neighbors(),
                        LocalEntropyParams{
                            .cutoff = settings.lef_cutoff,
                            .sigma = settings.lef_sigma,
                        });
  dists.addHistogram("LEF", std::move(hist));
}

correlation::analysis::Histogram
LocalEntropyCalculator::calculate(const correlation::core::Cell &cell,
                                  const correlation::analysis::StructureAnalyzer *neighbors,
                                  LocalEntropyParams params) {
  const auto &atoms = cell.atoms();
  size_t const num_atoms = atoms.size();
  if (num_atoms == 0) {
    return {};
  }

  bool const ignore_self_periodic = (neighbors != nullptr) ? neighbors->getIgnorePeriodicSelfInteractions() : false;

  MetricNeighborGrid grid;
  grid.build(cell, params.cutoff);

  struct ThreadLocalScratch {
    std::vector<real_t> distances;
  };

  tbb::enumerable_thread_specific<ThreadLocalScratch> thread_scratch;

  // Compute local entropy for all atoms
  std::vector<real_t> entropies(num_atoms, 0.0);
  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_atoms), [&](const tbb::blocked_range<size_t> &range) {
    auto &scratch = thread_scratch.local();
    for (size_t i = range.begin(); i != range.end(); ++i) {
      grid.queryDistances(i, cell, params.cutoff, scratch.distances, ignore_self_periodic);
      entropies[i] = computeSingleAtomEntropy(scratch.distances, cell, params);
    }
  });

  // Setup histogram configuration
  size_t const bins = 150;
  real_t const min_val = -15.0;
  real_t const max_val = 0.0;
  real_t const d_val = (max_val - min_val) / static_cast<real_t>(bins);

  correlation::analysis::Histogram hist;
  hist.title = "Local Entropic Fingerprint";
  hist.x_label = "Local Entropy (s_i / k_B)";
  hist.y_label = "Probability";
  hist.x_unit = "k_B";
  hist.y_unit = "probability";
  hist.description = "Local Entropic Fingerprint order parameter distribution";
  hist.file_suffix = "_lef";
  hist.bins.resize(bins);
  for (size_t i = 0; i < bins; ++i) {
    hist.bins[i] = min_val + (static_cast<real_t>(i) + static_cast<real_t>(0.5)) * d_val;
  }

  // Pre-size thread-local maps for element symbols
  std::vector<std::string> element_symbols;
  for (const auto &elem : cell.elements()) {
    element_symbols.push_back(elem.symbol);
  }

  struct ThreadLocalHist {
    std::map<std::string, std::vector<real_t>> partials;
    real_t num_atoms_f = 0.0;
  };

  tbb::enumerable_thread_specific<ThreadLocalHist> ets([&]() {
    ThreadLocalHist local;
    initHistogramMap(local.partials, element_symbols, bins);
    return local;
  });

  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_atoms), [&](const tbb::blocked_range<size_t> &range) {
    auto &local = ets.local();
    for (size_t i = range.begin(); i != range.end(); ++i) {
      local.num_atoms_f += 1.0;
      const std::string &symbol = atoms[i].element().symbol;
      addValueToHistogram(local.partials, symbol, entropies[i],
                          {
                              .min_val = min_val,
                              .max_val = max_val,
                              .d_val = d_val,
                          });
    }
  });

  // Reduce thread-local histograms
  std::map<std::string, std::vector<real_t>> partials;
  initHistogramMap(partials, element_symbols, bins);
  real_t total_atoms_f = 0.0;

  for (const auto &local_hist : ets) {
    accumulateHistogramMap(partials, local_hist.partials);
    total_atoms_f += local_hist.num_atoms_f;
  }

  // Normalize histograms to get probability distribution
  normalizeHistogramMap(partials, total_atoms_f);
  copyPartialsToHistogram(hist, partials);

  return hist;
}

} // namespace correlation::calculators
