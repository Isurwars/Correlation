/**
 * @file StructureFactorCalculator.cpp
 * @brief Implementation of structure factor S(Q) calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/StructureFactorCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Constants.hpp"
#include "math/SIMDUtils.hpp"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
// NOLINTNEXTLINE(cert-err58-cpp)
const bool registered = CalculatorFactory::registerTypeSafe<StructureFactorCalculator>("StructureFactorCalculator");

struct ReciprocalVector {
  real_t x;
  real_t y;
  real_t z;
};

struct ReciprocalBasis {
  ReciprocalVector b1;
  ReciprocalVector b2;
  ReciprocalVector b3;
  int hmax, kmax, lmax;
};

struct QVector {
  int h, k, l;
  real_t qmag;
};

struct TypeBlock {
  std::string symbol;
  size_t offset{}; // start index in flat position arrays
  size_t count{};
};

struct PartialInfo {
  std::string key;
  size_t typeA_idx; // offset into flat per-type phase arrays
  size_t typeB_idx;
  size_t N_A;
  size_t N_B;
  bool is_identical;
  real_t weight;
};

struct CoordinateArrays {
  const real_t *x;
  const real_t *y;
  const real_t *z;
};

struct PhaseArrays {
  std::vector<real_t> *cos;
  std::vector<real_t> *sin;
};

struct PrecomputedPhases {
  const real_t *E1_cos;
  const real_t *E1_sin;
  const real_t *E2_cos;
  const real_t *E2_sin;
  const real_t *E3_cos;
  const real_t *E3_sin;
};

struct ThreadAccumulators {
  std::vector<real_t> *total_sum;
  std::vector<size_t> *total_count;
  std::vector<std::vector<real_t>> *partial_sums;
  std::vector<std::vector<size_t>> *partial_counts;
  std::vector<real_t> *c_total_sum;
  std::vector<std::vector<real_t>> *c_partial_sums;
};

struct QBinning {
  real_t width;
  size_t num_bins;
};

ReciprocalBasis computeReciprocalBasis(const correlation::core::Cell &cell, real_t q_max) {
  const auto &inv = cell.inverseLatticeVectors();
  const auto two_pi_val = static_cast<real_t>(correlation::math::two_pi);
  const real_t bx_x = two_pi_val * inv(0, 0);
  const real_t bx_y = two_pi_val * inv(0, 1);
  const real_t bx_z = two_pi_val * inv(0, 2);
  const real_t by_x = two_pi_val * inv(1, 0);
  const real_t by_y = two_pi_val * inv(1, 1);
  const real_t by_z = two_pi_val * inv(1, 2);
  const real_t bz_x = two_pi_val * inv(2, 0);
  const real_t bz_y = two_pi_val * inv(2, 1);
  const real_t bz_z = two_pi_val * inv(2, 2);

  const real_t b1_norm = std::sqrt(bx_x * bx_x + bx_y * bx_y + bx_z * bx_z);
  const real_t b2_norm = std::sqrt(by_x * by_x + by_y * by_y + by_z * by_z);
  const real_t b3_norm = std::sqrt(bz_x * bz_x + bz_y * bz_y + bz_z * bz_z);

  const int hmax = (b1_norm > static_cast<real_t>(1e-10)) ? static_cast<int>(std::ceil(q_max / b1_norm)) : 0;
  const int kmax = (b2_norm > static_cast<real_t>(1e-10)) ? static_cast<int>(std::ceil(q_max / b2_norm)) : 0;
  const int lmax = (b3_norm > static_cast<real_t>(1e-10)) ? static_cast<int>(std::ceil(q_max / b3_norm)) : 0;

  return {{bx_x, bx_y, bx_z}, {by_x, by_y, by_z}, {bz_x, bz_y, bz_z}, hmax, kmax, lmax};
}

std::vector<QVector> generateQVectors(const ReciprocalBasis &basis, real_t q_max) {
  std::vector<QVector> q_vectors;
  q_vectors.reserve(8 * basis.hmax * basis.kmax * basis.lmax + 1);

  const real_t q_max_sq = q_max * q_max;
  for (int h_idx = -basis.hmax; h_idx <= basis.hmax; ++h_idx) {
    for (int k_idx = -basis.kmax; k_idx <= basis.kmax; ++k_idx) {
      for (int l_idx = -basis.lmax; l_idx <= basis.lmax; ++l_idx) {
        if (h_idx == 0 && k_idx == 0 && l_idx == 0) {
          continue;
        }
        const real_t q_x = static_cast<real_t>(h_idx) * basis.b1.x + static_cast<real_t>(k_idx) * basis.b2.x +
                           static_cast<real_t>(l_idx) * basis.b3.x;
        const real_t q_y = static_cast<real_t>(h_idx) * basis.b1.y + static_cast<real_t>(k_idx) * basis.b2.y +
                           static_cast<real_t>(l_idx) * basis.b3.y;
        const real_t q_z = static_cast<real_t>(h_idx) * basis.b1.z + static_cast<real_t>(k_idx) * basis.b2.z +
                           static_cast<real_t>(l_idx) * basis.b3.z;
        const real_t qmag_sq = q_x * q_x + q_y * q_y + q_z * q_z;
        if (qmag_sq <= q_max_sq) {
          q_vectors.push_back({.h = h_idx, .k = k_idx, .l = l_idx, .qmag = std::sqrt(qmag_sq)});
        }
      }
    }
  }
  return q_vectors;
}

void buildTypeBlocks(const std::vector<correlation::core::Atom> &atoms, std::vector<TypeBlock> &type_blocks,
                     std::vector<real_t> &x_s, std::vector<real_t> &y_s, std::vector<real_t> &z_s) {
  std::map<std::string, std::vector<size_t>> indices_by_type;
  for (size_t j = 0; j < atoms.size(); ++j) {
    indices_by_type[atoms[j].element().symbol].push_back(j);
  }

  type_blocks.reserve(indices_by_type.size());
  const size_t num_atoms = atoms.size();
  x_s.resize(num_atoms);
  y_s.resize(num_atoms);
  z_s.resize(num_atoms);

  size_t flat_offset = 0;
  for (const auto &[sym, idxs] : indices_by_type) {
    TypeBlock t_b;
    t_b.symbol = sym;
    t_b.offset = flat_offset;
    t_b.count = idxs.size();
    for (size_t j = 0; j < idxs.size(); ++j) {
      const auto &pos = atoms[idxs[j]].position();
      x_s[flat_offset + j] = pos.x();
      y_s[flat_offset + j] = pos.y();
      z_s[flat_offset + j] = pos.z();
    }
    flat_offset += idxs.size();
    type_blocks.push_back(t_b);
  }
}

std::vector<PartialInfo> buildPartialsInfo(const std::vector<TypeBlock> &type_blocks,
                                           const std::map<std::string, real_t> &ashcroft_weights) {
  std::vector<PartialInfo> partials_info;
  for (size_t ti = 0; ti < type_blocks.size(); ++ti) {
    for (size_t tj = ti; tj < type_blocks.size(); ++tj) {
      const auto &tbA = type_blocks[ti];
      const auto &tbB = type_blocks[tj];
      bool const is_identical = (ti == tj);
      std::string const key =
          (tbA.symbol < tbB.symbol) ? (tbA.symbol + "-" + tbB.symbol) : (tbB.symbol + "-" + tbA.symbol);
      real_t weight = static_cast<real_t>(0.0);
      auto wit = ashcroft_weights.find(key);
      if (wit != ashcroft_weights.end()) {
        weight = wit->second;
      }
      partials_info.push_back({.key = key,
                               .typeA_idx = ti,
                               .typeB_idx = tj,
                               .N_A = tbA.count,
                               .N_B = tbB.count,
                               .is_identical = is_identical,
                               .weight = weight});
    }
  }
  return partials_info;
}

void precomputePhases(int max_idx, const ReciprocalVector &rec_vec, size_t num_atoms, const CoordinateArrays &coords,
                      PhaseArrays phases) {
  phases.cos->resize((2 * max_idx + 1) * num_atoms);
  phases.sin->resize((2 * max_idx + 1) * num_atoms);

  real_t *cos_zero = &(*phases.cos)[max_idx * num_atoms];
  real_t *sin_zero = &(*phases.sin)[max_idx * num_atoms];
  std::fill_n(cos_zero, num_atoms, static_cast<real_t>(1.0));
  std::fill_n(sin_zero, num_atoms, static_cast<real_t>(0.0));

  if (max_idx == 0) {
    return;
  }

  std::vector<real_t> Cosine(num_atoms);
  std::vector<real_t> Sine(num_atoms);
  for (size_t j = 0; j < num_atoms; ++j) {
    real_t const phase = rec_vec.x * coords.x[j] + rec_vec.y * coords.y[j] + rec_vec.z * coords.z[j];
    Cosine[j] = std::cos(phase);
    Sine[j] = std::sin(phase);
  }

  real_t *cos_one = &(*phases.cos)[(1 + max_idx) * num_atoms];
  real_t *sin_one = &(*phases.sin)[(1 + max_idx) * num_atoms];
  std::copy(Cosine.begin(), Cosine.end(), cos_one);
  std::copy(Sine.begin(), Sine.end(), sin_one);

  real_t *cos_neg_one = &(*phases.cos)[(-1 + max_idx) * num_atoms];
  real_t *sin_neg_one = &(*phases.sin)[(-1 + max_idx) * num_atoms];
  std::copy(Cosine.begin(), Cosine.end(), cos_neg_one);
  for (size_t j = 0; j < num_atoms; ++j) {
    sin_neg_one[j] = -Sine[j];
  }

  for (int h_idx = 2; h_idx <= max_idx; ++h_idx) {
    const real_t *cos_prev = &(*phases.cos)[(h_idx - 1 + max_idx) * num_atoms];
    const real_t *sin_prev = &(*phases.sin)[(h_idx - 1 + max_idx) * num_atoms];
    real_t *cos_curr = &(*phases.cos)[(h_idx + max_idx) * num_atoms];
    real_t *sin_curr = &(*phases.sin)[(h_idx + max_idx) * num_atoms];
    real_t *cos_neg = &(*phases.cos)[(-h_idx + max_idx) * num_atoms];
    real_t *sin_neg = &(*phases.sin)[(-h_idx + max_idx) * num_atoms];

    for (size_t j = 0; j < num_atoms; ++j) {
      real_t const c_val = cos_prev[j] * Cosine[j] - sin_prev[j] * Sine[j];
      real_t const s_val = sin_prev[j] * Cosine[j] + cos_prev[j] * Sine[j];
      cos_curr[j] = c_val;
      sin_curr[j] = s_val;
      cos_neg[j] = c_val;
      sin_neg[j] = -s_val;
    }
  }
}

inline void processSingleQVector(const QVector &q_vec, QBinning binning, size_t num_atoms, const ReciprocalBasis &basis,
                                 const std::vector<TypeBlock> &type_blocks,
                                 const std::vector<PartialInfo> &partials_info, const PrecomputedPhases &phases,
                                 ThreadAccumulators accum) {
  const auto bin = static_cast<size_t>(q_vec.qmag / binning.width);
  if (bin >= binning.num_bins) {
    return;
  }

  // Compute cos/sin sums for each type block by multiplying the
  // per-index contributions (already factored by Miller h,k,l).
  std::vector<real_t> type_cos(type_blocks.size(), static_cast<real_t>(0.0));
  std::vector<real_t> type_sin(type_blocks.size(), static_cast<real_t>(0.0));

  for (size_t ti = 0; ti < type_blocks.size(); ++ti) {
    const size_t off = type_blocks[ti].offset;
    const size_t cnt = type_blocks[ti].count;
    real_t c_sum = static_cast<real_t>(0.0);
    real_t s_sum = static_cast<real_t>(0.0);
    correlation::math::miller_phase_sum(phases.E1_cos + (q_vec.h + basis.hmax) * num_atoms + off,
                                        phases.E1_sin + (q_vec.h + basis.hmax) * num_atoms + off,
                                        phases.E2_cos + (q_vec.k + basis.kmax) * num_atoms + off,
                                        phases.E2_sin + (q_vec.k + basis.kmax) * num_atoms + off,
                                        phases.E3_cos + (q_vec.l + basis.lmax) * num_atoms + off,
                                        phases.E3_sin + (q_vec.l + basis.lmax) * num_atoms + off, cnt, c_sum, s_sum);
    type_cos[ti] = c_sum;
    type_sin[ti] = s_sum;
  }

  // Accumulate total S(Q) via the full sum |rho(q)|^2 / N.
  real_t full_cos = static_cast<real_t>(0.0);
  real_t full_sin = static_cast<real_t>(0.0);
  for (size_t ti = 0; ti < type_blocks.size(); ++ti) {
    full_cos += type_cos[ti];
    full_sin += type_sin[ti];
  }
  const real_t sq_total = (full_cos * full_cos + full_sin * full_sin) / static_cast<real_t>(num_atoms);
  {
    real_t const y_total = sq_total - (*accum.c_total_sum)[bin];
    real_t const t_total = (*accum.total_sum)[bin] + y_total;
    (*accum.c_total_sum)[bin] = (t_total - (*accum.total_sum)[bin]) - y_total;
    (*accum.total_sum)[bin] = t_total;
  }
  (*accum.total_count)[bin] += 1;

  // Per-partial: S_AB(Q) = Re[rho_A(q)* · rho_B(q)] / sqrt(N_A * N_B)
  const size_t num_partials = partials_info.size();
  for (size_t pi = 0; pi < num_partials; ++pi) {
    const auto &pinfo = partials_info[pi];
    const size_t tiA = pinfo.typeA_idx;
    const size_t tiB = pinfo.typeB_idx;
    // Re[rho_A* rho_B] = cosA*cosB + sinA*sinB
    real_t const cross = type_cos[tiA] * type_cos[tiB] + type_sin[tiA] * type_sin[tiB];
    real_t const denom = std::sqrt(static_cast<real_t>(pinfo.N_A) * static_cast<real_t>(pinfo.N_B));
    real_t const sq_partial = (denom > static_cast<real_t>(0.0)) ? cross / denom : static_cast<real_t>(0.0);
    real_t const y_partial = sq_partial - (*accum.c_partial_sums)[pi][bin];
    real_t const t_partial = (*accum.partial_sums)[pi][bin] + y_partial;
    (*accum.c_partial_sums)[pi][bin] = (t_partial - (*accum.partial_sums)[pi][bin]) - y_partial;
    (*accum.partial_sums)[pi][bin] = t_partial;
    (*accum.partial_counts)[pi][bin] += 1;
  }
}

using LocalAccumulatorTuple =
    std::tuple<std::vector<real_t>, std::vector<size_t>, std::vector<std::vector<real_t>>,
               std::vector<std::vector<size_t>>, std::vector<real_t>, std::vector<std::vector<real_t>>>;

correlation::analysis::Histogram createInitialHistogram(const QBinning &binning,
                                                        const std::vector<PartialInfo> &partials_info) {
  correlation::analysis::Histogram s_q_hist;
  s_q_hist.bins.resize(binning.num_bins);
  s_q_hist.x_label = "Q";
  s_q_hist.title = "S(Q) — Structure Factor";
  s_q_hist.y_label = "S(Q)";
  s_q_hist.x_unit = "Å⁻¹";
  s_q_hist.y_unit = "arbitrary units";
  s_q_hist.description = "Structure Factor S(Q)";
  s_q_hist.file_suffix = "_S";
  for (size_t i = 0; i < binning.num_bins; ++i) {
    s_q_hist.bins[i] = (static_cast<real_t>(i) + static_cast<real_t>(0.5)) * binning.width;
  }

  for (const auto &partial : partials_info) {
    s_q_hist.partials[partial.key].assign(binning.num_bins, static_cast<real_t>(0.0));
  }
  return s_q_hist;
}

void combineAccumulators(tbb::enumerable_thread_specific<LocalAccumulatorTuple> &local_ets, size_t num_q_bins,
                         size_t num_partials, std::vector<real_t> &sq_total_sum, std::vector<size_t> &sq_total_count,
                         std::vector<std::vector<real_t>> &partial_sums,
                         std::vector<std::vector<size_t>> &partial_counts) {
  std::vector<real_t> c_sq_total_sum(num_q_bins, static_cast<real_t>(0.0));
  std::vector<std::vector<real_t>> c_partial_sums(num_partials,
                                                  std::vector<real_t>(num_q_bins, static_cast<real_t>(0.0)));

  local_ets.combine_each([&](const auto &local) {
    const auto &[lt_sum, lt_count, lp_sums, lp_counts, lc_sum, lcp_sums] = local;
    for (size_t i = 0; i < num_q_bins; ++i) {
      real_t const y_val = lt_sum[i] - c_sq_total_sum[i];
      real_t const t_val = sq_total_sum[i] + y_val;
      c_sq_total_sum[i] = (t_val - sq_total_sum[i]) - y_val;
      sq_total_sum[i] = t_val;

      sq_total_count[i] += lt_count[i];
    }
    for (size_t pi = 0; pi < num_partials; ++pi) {
      for (size_t i = 0; i < num_q_bins; ++i) {
        real_t const y_val = lp_sums[pi][i] - c_partial_sums[pi][i];
        real_t const t_val = partial_sums[pi][i] + y_val;
        c_partial_sums[pi][i] = (t_val - partial_sums[pi][i]) - y_val;
        partial_sums[pi][i] = t_val;

        partial_counts[pi][i] += lp_counts[pi][i];
      }
    }
  });
}

void normalizeAndStoreResults(correlation::analysis::Histogram &s_q_hist, const std::vector<real_t> &sq_total_sum,
                              const std::vector<size_t> &sq_total_count,
                              const std::vector<std::vector<real_t>> &partial_sums,
                              const std::vector<std::vector<size_t>> &partial_counts,
                              const std::vector<PartialInfo> &partials_info, size_t num_q_bins) {
  const size_t num_partials = partials_info.size();
  std::vector<real_t> total_sq(num_q_bins, static_cast<real_t>(0.0));
  for (size_t i = 0; i < num_q_bins; ++i) {
    if (sq_total_count[i] > 0) {
      total_sq[i] = sq_total_sum[i] / static_cast<real_t>(sq_total_count[i]);
    }
  }

  for (size_t pi = 0; pi < num_partials; ++pi) {
    auto &out = s_q_hist.partials[partials_info[pi].key];
    for (size_t i = 0; i < num_q_bins; ++i) {
      if (partial_counts[pi][i] > 0) {
        out[i] = partial_sums[pi][i] / static_cast<real_t>(partial_counts[pi][i]);
      }
    }
  }

  s_q_hist.partials["Total"] = std::move(total_sq);
}

} // namespace

void StructureFactorCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                               const correlation::analysis::AnalysisSettings &settings) const {
  if (settings.q_bin_width <= 0 || settings.q_max <= 0) {
    throw std::invalid_argument("Q-space parameters must be positive.");
  }

  const correlation::core::Cell &cell = dists.cell();
  const auto &atoms = cell.atoms();
  const size_t num_atoms = atoms.size();
  if (num_atoms == 0) {
    return;
  }

  const real_t q_max = static_cast<real_t>(settings.q_max);
  const real_t q_bin_width = static_cast<real_t>(settings.q_bin_width);
  const auto num_q_bins = static_cast<size_t>(std::floor(q_max / q_bin_width));
  if (num_q_bins == 0) {
    throw std::invalid_argument("Q_max is too small for the given Q_bin_width.");
  }

  ReciprocalBasis const basis = computeReciprocalBasis(cell, q_max);
  std::vector<QVector> const q_vectors = generateQVectors(basis, q_max);

  if (q_vectors.empty()) {
    return;
  }

  std::vector<TypeBlock> type_blocks;
  std::vector<real_t> x_s;
  std::vector<real_t> y_s;
  std::vector<real_t> z_s;
  buildTypeBlocks(atoms, type_blocks, x_s, y_s, z_s);

  const auto &ashcroft_weights = dists.getAshcroftWeights();
  std::vector<PartialInfo> const partials_info = buildPartialsInfo(type_blocks, ashcroft_weights);

  std::vector<real_t> E1_cos;
  std::vector<real_t> E1_sin;
  std::vector<real_t> E2_cos;
  std::vector<real_t> E2_sin;
  std::vector<real_t> E3_cos;
  std::vector<real_t> E3_sin;

  precomputePhases(basis.hmax, basis.b1, num_atoms, {.x = x_s.data(), .y = y_s.data(), .z = z_s.data()},
                   {.cos = &E1_cos, .sin = &E1_sin});
  precomputePhases(basis.kmax, basis.b2, num_atoms, {.x = x_s.data(), .y = y_s.data(), .z = z_s.data()},
                   {.cos = &E2_cos, .sin = &E2_sin});
  precomputePhases(basis.lmax, basis.b3, num_atoms, {.x = x_s.data(), .y = y_s.data(), .z = z_s.data()},
                   {.cos = &E3_cos, .sin = &E3_sin});

  QBinning const binning{.width = q_bin_width, .num_bins = num_q_bins};
  correlation::analysis::Histogram s_q_hist = createInitialHistogram(binning, partials_info);

  std::vector<real_t> sq_total_sum(num_q_bins, static_cast<real_t>(0.0));
  std::vector<size_t> sq_total_count(num_q_bins, 0);

  const size_t num_partials = partials_info.size();
  std::vector<std::vector<real_t>> partial_sums(num_partials,
                                                std::vector<real_t>(num_q_bins, static_cast<real_t>(0.0)));
  std::vector<std::vector<size_t>> partial_counts(num_partials, std::vector<size_t>(num_q_bins, 0));

  tbb::enumerable_thread_specific<LocalAccumulatorTuple> local_ets([&] {
    return std::make_tuple(
        std::vector<real_t>(num_q_bins, static_cast<real_t>(0.0)), std::vector<size_t>(num_q_bins, 0),
        std::vector<std::vector<real_t>>(num_partials, std::vector<real_t>(num_q_bins, static_cast<real_t>(0.0))),
        std::vector<std::vector<size_t>>(num_partials, std::vector<size_t>(num_q_bins, 0)),
        std::vector<real_t>(num_q_bins, static_cast<real_t>(0.0)),
        std::vector<std::vector<real_t>>(num_partials, std::vector<real_t>(num_q_bins, static_cast<real_t>(0.0))));
  });

  tbb::parallel_for(tbb::blocked_range<size_t>(0, q_vectors.size()), [&](const tbb::blocked_range<size_t> &range) {
    auto &[local_total_sum, local_total_count, local_partial_sums, local_partial_counts, local_c_total_sum,
           local_c_partial_sums] = local_ets.local();

    for (size_t qi = range.begin(); qi != range.end(); ++qi) {
      processSingleQVector(q_vectors[qi], binning, num_atoms, basis, type_blocks, partials_info,
                           {.E1_cos = E1_cos.data(),
                            .E1_sin = E1_sin.data(),
                            .E2_cos = E2_cos.data(),
                            .E2_sin = E2_sin.data(),
                            .E3_cos = E3_cos.data(),
                            .E3_sin = E3_sin.data()},
                           {.total_sum = &local_total_sum,
                            .total_count = &local_total_count,
                            .partial_sums = &local_partial_sums,
                            .partial_counts = &local_partial_counts,
                            .c_total_sum = &local_c_total_sum,
                            .c_partial_sums = &local_c_partial_sums});
    }
  });

  combineAccumulators(local_ets, num_q_bins, num_partials, sq_total_sum, sq_total_count, partial_sums, partial_counts);
  normalizeAndStoreResults(s_q_hist, sq_total_sum, sq_total_count, partial_sums, partial_counts, partials_info,
                           num_q_bins);

  dists.addHistogram("S_q", std::move(s_q_hist));
}

} // namespace correlation::calculators
