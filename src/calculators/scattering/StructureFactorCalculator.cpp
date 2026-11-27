/**
 * @file StructureFactorCalculator.cpp
 * @brief Implementation of structure factor S(Q) calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/scattering/StructureFactorCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Constants.hpp"
#include "math/SIMDUtils.hpp"

#if defined(CORRELATION_USE_CUDA) || defined(CORRELATION_USE_HIP)
#include "calculators/gpu/GPUSQCalculator.hpp"
#endif

#include <algorithm>
#include <map>
#include <stdexcept>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
const bool REGISTERED =
    CalculatorFactory::registerTypeSafe<StructureFactorCalculator>("StructureFactorCalculator");

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
  const real_t *e1_cos;
  const real_t *e1_sin;
  const real_t *e2_cos;
  const real_t *e2_sin;
  const real_t *e3_cos;
  const real_t *e3_sin;
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
  const auto two_pi_val = static_cast<real_t>(correlation::math::TWO_PI);
  const real_t bx_x = two_pi_val * inv(0, 0);
  const real_t bx_y = two_pi_val * inv(0, 1);
  const real_t bx_z = two_pi_val * inv(0, 2);
  const real_t by_x = two_pi_val * inv(1, 0);
  const real_t by_y = two_pi_val * inv(1, 1);
  const real_t by_z = two_pi_val * inv(1, 2);
  const real_t bz_x = two_pi_val * inv(2, 0);
  const real_t bz_y = two_pi_val * inv(2, 1);
  const real_t bz_z = two_pi_val * inv(2, 2);

  const auto &lattice = cell.latticeVectors();
  const real_t a_norm = correlation::math::norm(lattice[0]);
  const real_t b_norm = correlation::math::norm(lattice[1]);
  const real_t c_norm = correlation::math::norm(lattice[2]);

  // Exact Miller-index bounding box for arbitrary triclinic/hexagonal cells:
  // |h| = |a . q| / (2*pi) <= |a| * q_max / (2*pi)
  const int hmax = (a_norm > static_cast<real_t>(1e-10))
                       ? static_cast<int>(std::ceil(q_max * a_norm / two_pi_val))
                       : 0;
  const int kmax = (b_norm > static_cast<real_t>(1e-10))
                       ? static_cast<int>(std::ceil(q_max * b_norm / two_pi_val))
                       : 0;
  const int lmax = (c_norm > static_cast<real_t>(1e-10))
                       ? static_cast<int>(std::ceil(q_max * c_norm / two_pi_val))
                       : 0;

  return ReciprocalBasis{
      .b1 = {.x = bx_x, .y = bx_y, .z = bx_z},
      .b2 = {.x = by_x, .y = by_y, .z = by_z},
      .b3 = {.x = bz_x, .y = bz_y, .z = bz_z},
      .hmax = hmax,
      .kmax = kmax,
      .lmax = lmax,
  };
}

constexpr size_t K_MAX_DIRECT_QVECTORS = 25000;
constexpr size_t K_SMALL_CELL_ATOM_THRESHOLD = 1000;
constexpr size_t K_STRATIFIED_QVECTORS_PER_BIN = 512;

[[nodiscard]] double computeDeterministicCellPhase(const correlation::core::Cell &cell) noexcept {
  const auto &atoms = cell.atoms();
  if (atoms.empty()) {
    return 0.0;
  }
  const auto &pos0 = atoms.front().position();
  const auto &pos_last = atoms.back().position();
  const double mix = static_cast<double>(pos0.x()) * 12.9898 +
                     static_cast<double>(pos0.y()) * 78.233 +
                     static_cast<double>(pos0.z()) * 37.719 +
                     static_cast<double>(pos_last.x()) * 45.164;
  const double frac = std::abs(std::sin(mix) * 43758.5453);
  return (frac - std::floor(frac)) * correlation::math::TWO_PI;
}

[[nodiscard]] bool isFriedelHalfSpace(int h_idx, int k_idx, int l_idx) noexcept {
  return (h_idx > 0) || (h_idx == 0 && k_idx > 0) || (h_idx == 0 && k_idx == 0 && l_idx > 0);
}

[[nodiscard]] uint64_t packMillerKey(int h_idx, int k_idx, int l_idx) noexcept {
  constexpr int64_t BIAS = 1048576;
  const auto h_u = static_cast<uint64_t>(static_cast<int64_t>(h_idx) + BIAS) & 0x1FFFFFULL;
  const auto k_u = static_cast<uint64_t>(static_cast<int64_t>(k_idx) + BIAS) & 0x1FFFFFULL;
  const auto l_u = static_cast<uint64_t>(static_cast<int64_t>(l_idx) + BIAS) & 0x1FFFFFULL;
  return (h_u << 42U) | (k_u << 21U) | l_u;
}

std::vector<QVector> generateQVectors(const ReciprocalBasis &basis, real_t q_max) {
  std::vector<QVector> q_vectors;
  const size_t box_est = static_cast<size_t>(2 * basis.hmax + 1) *
                         static_cast<size_t>(2 * basis.kmax + 1) *
                         static_cast<size_t>(2 * basis.lmax + 1) / 2;
  q_vectors.reserve(std::min<size_t>(box_est, K_MAX_DIRECT_QVECTORS * 2 + 4096));

  const real_t q_max_sq = q_max * q_max;
  for (int h_idx = 0; h_idx <= basis.hmax; ++h_idx) {
    for (int k_idx = -basis.kmax; k_idx <= basis.kmax; ++k_idx) {
      for (int l_idx = -basis.lmax; l_idx <= basis.lmax; ++l_idx) {
        if (!isFriedelHalfSpace(h_idx, k_idx, l_idx)) {
          continue;
        }
        const real_t q_x = static_cast<real_t>(h_idx) * basis.b1.x +
                           static_cast<real_t>(k_idx) * basis.b2.x +
                           static_cast<real_t>(l_idx) * basis.b3.x;
        const real_t q_y = static_cast<real_t>(h_idx) * basis.b1.y +
                           static_cast<real_t>(k_idx) * basis.b2.y +
                           static_cast<real_t>(l_idx) * basis.b3.y;
        const real_t q_z = static_cast<real_t>(h_idx) * basis.b1.z +
                           static_cast<real_t>(k_idx) * basis.b2.z +
                           static_cast<real_t>(l_idx) * basis.b3.z;
        const real_t qmag_sq = q_x * q_x + q_y * q_y + q_z * q_z;
        if (qmag_sq <= q_max_sq) {
          q_vectors.push_back({
              .h = h_idx,
              .k = k_idx,
              .l = l_idx,
              .qmag = std::sqrt(qmag_sq),
          });
        }
      }
    }
  }
  return q_vectors;
}

struct FibonacciShellParams {
  double phi_offset;
  double q_lo;
  double q_hi;
  real_t q_exact_max;
  double q_max_d;
};

[[nodiscard]] bool trySampleFibonacciMiller(size_t sample_idx, const FibonacciShellParams &params,
                                            const correlation::math::Matrix3<real_t> &lattice,
                                            const ReciprocalBasis &full_basis,
                                            QVector &out_vec) noexcept {
  constexpr double GOLDEN_ANGLE = 2.39996322972865332; // pi * (3 - sqrt(5))
  const double inv_two_pi = 1.0 / correlation::math::TWO_PI;
  const double frac = (static_cast<double>(sample_idx) + 0.5) /
                      static_cast<double>(K_STRATIFIED_QVECTORS_PER_BIN);
  const double cos_theta = frac;
  const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
  const double phi = params.phi_offset + static_cast<double>(sample_idx) * GOLDEN_ANGLE;
  const double q_target = params.q_lo + frac * (params.q_hi - params.q_lo);

  const double qx_t = q_target * sin_theta * std::cos(phi);
  const double qy_t = q_target * sin_theta * std::sin(phi);
  const double qz_t = q_target * cos_theta;

  int h_idx = static_cast<int>(std::llround(
      (static_cast<double>(lattice[0].x()) * qx_t + static_cast<double>(lattice[0].y()) * qy_t +
       static_cast<double>(lattice[0].z()) * qz_t) *
      inv_two_pi));
  int k_idx = static_cast<int>(std::llround(
      (static_cast<double>(lattice[1].x()) * qx_t + static_cast<double>(lattice[1].y()) * qy_t +
       static_cast<double>(lattice[1].z()) * qz_t) *
      inv_two_pi));
  int l_idx = static_cast<int>(std::llround(
      (static_cast<double>(lattice[2].x()) * qx_t + static_cast<double>(lattice[2].y()) * qy_t +
       static_cast<double>(lattice[2].z()) * qz_t) *
      inv_two_pi));

  if (h_idx == 0 && k_idx == 0 && l_idx == 0) {
    return false;
  }
  if (!isFriedelHalfSpace(h_idx, k_idx, l_idx)) {
    h_idx = -h_idx;
    k_idx = -k_idx;
    l_idx = -l_idx;
  }
  if (std::abs(h_idx) > full_basis.hmax || std::abs(k_idx) > full_basis.kmax ||
      std::abs(l_idx) > full_basis.lmax) {
    return false;
  }

  const real_t q_x = static_cast<real_t>(h_idx) * full_basis.b1.x +
                     static_cast<real_t>(k_idx) * full_basis.b2.x +
                     static_cast<real_t>(l_idx) * full_basis.b3.x;
  const real_t q_y = static_cast<real_t>(h_idx) * full_basis.b1.y +
                     static_cast<real_t>(k_idx) * full_basis.b2.y +
                     static_cast<real_t>(l_idx) * full_basis.b3.y;
  const real_t q_z = static_cast<real_t>(h_idx) * full_basis.b1.z +
                     static_cast<real_t>(k_idx) * full_basis.b2.z +
                     static_cast<real_t>(l_idx) * full_basis.b3.z;
  const real_t qmag = std::sqrt(q_x * q_x + q_y * q_y + q_z * q_z);
  if (qmag <= params.q_exact_max || static_cast<double>(qmag) > params.q_max_d) {
    return false;
  }

  out_vec = {.h = h_idx, .k = k_idx, .l = l_idx, .qmag = qmag};
  return true;
}

std::vector<QVector> generateStratifiedQVectors(const correlation::core::Cell &cell,
                                                const ReciprocalBasis &full_basis,
                                                real_t q_exact_max, QBinning binning) {
  const ReciprocalBasis low_basis = computeReciprocalBasis(cell, q_exact_max);
  std::vector<QVector> q_vectors = generateQVectors(low_basis, q_exact_max);

  const auto start_bin = static_cast<size_t>(std::floor(q_exact_max / binning.width));
  if (start_bin >= binning.num_bins) {
    return q_vectors;
  }

  const double q_max_d = static_cast<double>(binning.num_bins) * static_cast<double>(binning.width);
  const double cell_phase = computeDeterministicCellPhase(cell);
  const auto &lattice = cell.latticeVectors();
  q_vectors.reserve(q_vectors.size() +
                    (binning.num_bins - start_bin) * K_STRATIFIED_QVECTORS_PER_BIN);

  std::vector<uint64_t> seen_keys;
  seen_keys.reserve(K_STRATIFIED_QVECTORS_PER_BIN);

  for (size_t bin_idx = start_bin; bin_idx < binning.num_bins; ++bin_idx) {
    seen_keys.clear();
    const double q_lo = static_cast<double>(bin_idx) * static_cast<double>(binning.width);
    const FibonacciShellParams shell_params{
        .phi_offset = cell_phase + static_cast<double>(bin_idx) * 1.6180339887498948,
        .q_lo = q_lo,
        .q_hi = q_lo + static_cast<double>(binning.width),
        .q_exact_max = q_exact_max,
        .q_max_d = q_max_d,
    };

    for (size_t s_idx = 0; s_idx < K_STRATIFIED_QVECTORS_PER_BIN; ++s_idx) {
      QVector cand{};
      if (!trySampleFibonacciMiller(s_idx, shell_params, lattice, full_basis, cand)) {
        continue;
      }
      const uint64_t key = packMillerKey(cand.h, cand.k, cand.l);
      if (std::ranges::find(seen_keys, key) != seen_keys.end()) {
        continue;
      }
      seen_keys.push_back(key);
      q_vectors.push_back(cand);
    }
  }
  return q_vectors;
}

void buildTypeBlocks(const std::vector<correlation::core::Atom> &atoms,
                     std::vector<TypeBlock> &type_blocks, std::vector<real_t> &x_s,
                     std::vector<real_t> &y_s, std::vector<real_t> &z_s) {
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
      const auto &tb_a = type_blocks[ti];
      const auto &tb_b = type_blocks[tj];
      bool const is_identical = (ti == tj);
      std::string const key = (tb_a.symbol < tb_b.symbol) ? (tb_a.symbol + "-" + tb_b.symbol)
                                                          : (tb_b.symbol + "-" + tb_a.symbol);
      auto weight = static_cast<real_t>(0.0);
      const auto wit = ashcroft_weights.find(key);
      if (wit != ashcroft_weights.end()) {
        weight = wit->second;
      }
      partials_info.push_back({
          .key = key,
          .typeA_idx = ti,
          .typeB_idx = tj,
          .N_A = tb_a.count,
          .N_B = tb_b.count,
          .is_identical = is_identical,
          .weight = weight,
      });
    }
  }
  return partials_info;
}

void precomputePhases(int max_idx, const ReciprocalVector &rec_vec, size_t num_atoms,
                      const CoordinateArrays &coords, PhaseArrays phases) {
  phases.cos->resize((2 * max_idx + 1) * num_atoms);
  phases.sin->resize((2 * max_idx + 1) * num_atoms);

  real_t *cos_zero = &(*phases.cos)[max_idx * num_atoms];
  real_t *sin_zero = &(*phases.sin)[max_idx * num_atoms];
  std::fill_n(cos_zero, num_atoms, static_cast<real_t>(1.0));
  std::fill_n(sin_zero, num_atoms, static_cast<real_t>(0.0));

  if (max_idx == 0) {
    return;
  }

  std::vector<real_t> cosine(num_atoms);
  std::vector<real_t> sine(num_atoms);
  for (size_t j = 0; j < num_atoms; ++j) {
    real_t const phase =
        rec_vec.x * coords.x[j] + rec_vec.y * coords.y[j] + rec_vec.z * coords.z[j];
    cosine[j] = std::cos(phase);
    sine[j] = std::sin(phase);
  }

  real_t *cos_one = &(*phases.cos)[(1 + max_idx) * num_atoms];
  real_t *sin_one = &(*phases.sin)[(1 + max_idx) * num_atoms];
  std::ranges::copy(cosine, cos_one);
  std::ranges::copy(sine, sin_one);

  real_t *cos_neg_one = &(*phases.cos)[(-1 + max_idx) * num_atoms];
  real_t *sin_neg_one = &(*phases.sin)[(-1 + max_idx) * num_atoms];
  std::ranges::copy(cosine, cos_neg_one);
  for (size_t j = 0; j < num_atoms; ++j) {
    sin_neg_one[j] = -sine[j];
  }

  // Run angle-addition recurrence in double precision to prevent O(hmax * eps) drift in float mode
  std::vector<double> prev_cos_d(num_atoms);
  std::vector<double> prev_sin_d(num_atoms);
  std::vector<double> base_cos_d(num_atoms);
  std::vector<double> base_sin_d(num_atoms);
  for (size_t j = 0; j < num_atoms; ++j) {
    base_cos_d[j] = static_cast<double>(cosine[j]);
    base_sin_d[j] = static_cast<double>(sine[j]);
    prev_cos_d[j] = base_cos_d[j];
    prev_sin_d[j] = base_sin_d[j];
  }

  for (int h_idx = 2; h_idx <= max_idx; ++h_idx) {
    real_t *cos_curr = &(*phases.cos)[(h_idx + max_idx) * num_atoms];
    real_t *sin_curr = &(*phases.sin)[(h_idx + max_idx) * num_atoms];
    real_t *cos_neg = &(*phases.cos)[(-h_idx + max_idx) * num_atoms];
    real_t *sin_neg = &(*phases.sin)[(-h_idx + max_idx) * num_atoms];

    for (size_t j = 0; j < num_atoms; ++j) {
      const double c_val = prev_cos_d[j] * base_cos_d[j] - prev_sin_d[j] * base_sin_d[j];
      const double s_val = prev_sin_d[j] * base_cos_d[j] + prev_cos_d[j] * base_sin_d[j];
      prev_cos_d[j] = c_val;
      prev_sin_d[j] = s_val;
      const auto c_out = static_cast<real_t>(c_val);
      const auto s_out = static_cast<real_t>(s_val);
      cos_curr[j] = c_out;
      sin_curr[j] = s_out;
      cos_neg[j] = c_out;
      sin_neg[j] = -s_out;
    }
  }
}

inline void processSingleQVector(const QVector &q_vec, QBinning binning, size_t num_atoms,
                                 const ReciprocalBasis &basis,
                                 const std::vector<TypeBlock> &type_blocks,
                                 const std::vector<PartialInfo> &partials_info,
                                 const PrecomputedPhases &phases, ThreadAccumulators accum,
                                 std::vector<real_t> &type_cos, std::vector<real_t> &type_sin) {
  const auto bin = static_cast<size_t>(q_vec.qmag / binning.width);
  if (bin >= binning.num_bins) {
    return;
  }

  // Compute cos/sin sums for each type block into thread-local reusable scratch buffers.
  type_cos.assign(type_blocks.size(), static_cast<real_t>(0.0));
  type_sin.assign(type_blocks.size(), static_cast<real_t>(0.0));

  for (size_t ti = 0; ti < type_blocks.size(); ++ti) {
    const size_t off = type_blocks[ti].offset;
    const size_t cnt = type_blocks[ti].count;
    const auto result =
        correlation::math::millerPhaseSum(correlation::math::MillerPhaseSumParams<real_t>{
            .cos1 = phases.e1_cos + (q_vec.h + basis.hmax) * num_atoms + off,
            .sin1 = phases.e1_sin + (q_vec.h + basis.hmax) * num_atoms + off,
            .cos2 = phases.e2_cos + (q_vec.k + basis.kmax) * num_atoms + off,
            .sin2 = phases.e2_sin + (q_vec.k + basis.kmax) * num_atoms + off,
            .cos3 = phases.e3_cos + (q_vec.l + basis.lmax) * num_atoms + off,
            .sin3 = phases.e3_sin + (q_vec.l + basis.lmax) * num_atoms + off,
            .count = cnt,
        });
    type_cos[ti] = result.cos_sum;
    type_sin[ti] = result.sin_sum;
  }

  // Accumulate total S(Q) via the full sum |rho(q)|^2 / N.
  auto full_cos = static_cast<real_t>(0.0);
  auto full_sin = static_cast<real_t>(0.0);
  for (size_t ti = 0; ti < type_blocks.size(); ++ti) {
    full_cos += type_cos[ti];
    full_sin += type_sin[ti];
  }
  const real_t sq_total =
      (full_cos * full_cos + full_sin * full_sin) / static_cast<real_t>(num_atoms);
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
    const size_t ti_a = pinfo.typeA_idx;
    const size_t ti_b = pinfo.typeB_idx;
    // Re[rho_A* rho_B] = cosA*cosB + sinA*sinB
    real_t const cross = type_cos[ti_a] * type_cos[ti_b] + type_sin[ti_a] * type_sin[ti_b];
    real_t const denom = std::sqrt(static_cast<real_t>(pinfo.N_A) * static_cast<real_t>(pinfo.N_B));
    real_t const sq_partial =
        (denom > static_cast<real_t>(0.0)) ? cross / denom : static_cast<real_t>(0.0);
    real_t const y_partial = sq_partial - (*accum.c_partial_sums)[pi][bin];
    real_t const t_partial = (*accum.partial_sums)[pi][bin] + y_partial;
    (*accum.c_partial_sums)[pi][bin] = (t_partial - (*accum.partial_sums)[pi][bin]) - y_partial;
    (*accum.partial_sums)[pi][bin] = t_partial;
    (*accum.partial_counts)[pi][bin] += 1;
  }
}

using LocalAccumulatorTuple =
    std::tuple<std::vector<real_t>, std::vector<size_t>, std::vector<std::vector<real_t>>,
               std::vector<std::vector<size_t>>, std::vector<real_t>,
               std::vector<std::vector<real_t>>>;

correlation::analysis::Histogram
createInitialHistogram(const QBinning &binning, const std::vector<PartialInfo> &partials_info) {
  correlation::analysis::Histogram s_q_hist;
  s_q_hist.bins.resize(binning.num_bins);
  s_q_hist.x_label = "Q";
  s_q_hist.title = "S(Q) — Structure Factor";
  s_q_hist.y_label = "S(Q)";
  s_q_hist.x_unit = "Å⁻¹";
  s_q_hist.y_unit = "";
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

void combineAccumulators(tbb::enumerable_thread_specific<LocalAccumulatorTuple> &local_ets,
                         size_t num_q_bins, size_t num_partials, std::vector<real_t> &sq_total_sum,
                         std::vector<size_t> &sq_total_count,
                         std::vector<std::vector<real_t>> &partial_sums,
                         std::vector<std::vector<size_t>> &partial_counts) {
  std::vector<real_t> c_sq_total_sum(num_q_bins, static_cast<real_t>(0.0));
  std::vector<std::vector<real_t>> c_partial_sums(
      num_partials, std::vector<real_t>(num_q_bins, static_cast<real_t>(0.0)));

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

void normalizeAndStoreResults(correlation::analysis::Histogram &s_q_hist,
                              const std::vector<real_t> &sq_total_sum,
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

[[nodiscard]] real_t computeDirectQCutoff(real_t cell_volume, real_t q_max,
                                          real_t q_bin_width) noexcept {
  const double vol = std::abs(static_cast<double>(cell_volume));
  const auto q_max_d = static_cast<double>(q_max);
  if (vol <= 1e-12 || q_max_d <= 0.0) {
    return q_max;
  }
  const double twelve_pi_sq = 12.0 * correlation::math::PI * correlation::math::PI;
  const double est_half_vectors = (vol * q_max_d * q_max_d * q_max_d) / twelve_pi_sq;
  if (est_half_vectors <= static_cast<double>(K_MAX_DIRECT_QVECTORS)) {
    return q_max;
  }
  const double q_switch =
      std::cbrt((twelve_pi_sq * static_cast<double>(K_MAX_DIRECT_QVECTORS)) / vol);
  return std::clamp(static_cast<real_t>(q_switch), q_bin_width, q_max);
}

} // namespace

void StructureFactorCalculator::calculateFrame(
    correlation::analysis::DistributionFunctions &dists,
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

#if defined(CORRELATION_USE_CUDA) || defined(CORRELATION_USE_HIP)
  static const GPUSQCalculator GPU_SQ_CALC;
  if (GPU_SQ_CALC.hasGPU() && cell.elements().size() == 1) {
    GPU_SQ_CALC.calculateFrame(dists, settings);
    return;
  }
#endif

  const auto q_max = static_cast<real_t>(settings.q_max);
  const auto q_bin_width = static_cast<real_t>(settings.q_bin_width);
  const auto num_q_bins = static_cast<size_t>(std::floor(q_max / q_bin_width));
  if (num_q_bins == 0) {
    throw std::invalid_argument("Q_max is too small for the given Q_bin_width.");
  }

  const real_t cell_volume = cell.volume();
  const bool force_full_direct =
      settings.sq_exact_direct || (num_atoms < K_SMALL_CELL_ATOM_THRESHOLD);
  const real_t q_direct_max =
      force_full_direct ? q_max : computeDirectQCutoff(cell_volume, q_max, q_bin_width);
  const bool use_stratified = (q_direct_max < q_max);

  QBinning const binning{
      .width = q_bin_width,
      .num_bins = num_q_bins,
  };
  ReciprocalBasis const basis = computeReciprocalBasis(cell, q_max);
  std::vector<QVector> const q_vectors =
      use_stratified ? generateStratifiedQVectors(cell, basis, q_direct_max, binning)
                     : generateQVectors(basis, q_max);

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

  correlation::analysis::Histogram s_q_hist = createInitialHistogram(binning, partials_info);

  std::vector<real_t> sq_total_sum(num_q_bins, static_cast<real_t>(0.0));
  std::vector<size_t> sq_total_count(num_q_bins, 0);

  const size_t num_partials = partials_info.size();
  std::vector<std::vector<real_t>> partial_sums(
      num_partials, std::vector<real_t>(num_q_bins, static_cast<real_t>(0.0)));
  std::vector<std::vector<size_t>> partial_counts(num_partials, std::vector<size_t>(num_q_bins, 0));

  if (!q_vectors.empty()) {
    std::vector<real_t> e1_cos;
    std::vector<real_t> e1_sin;
    std::vector<real_t> e2_cos;
    std::vector<real_t> e2_sin;
    std::vector<real_t> e3_cos;
    std::vector<real_t> e3_sin;

    precomputePhases(basis.hmax, basis.b1, num_atoms,
                     {.x = x_s.data(), .y = y_s.data(), .z = z_s.data()},
                     {.cos = &e1_cos, .sin = &e1_sin});
    precomputePhases(basis.kmax, basis.b2, num_atoms,
                     {.x = x_s.data(), .y = y_s.data(), .z = z_s.data()},
                     {.cos = &e2_cos, .sin = &e2_sin});
    precomputePhases(basis.lmax, basis.b3, num_atoms,
                     {.x = x_s.data(), .y = y_s.data(), .z = z_s.data()},
                     {.cos = &e3_cos, .sin = &e3_sin});

    tbb::enumerable_thread_specific<LocalAccumulatorTuple> local_ets([&] {
      return std::make_tuple(
          std::vector<real_t>(num_q_bins, static_cast<real_t>(0.0)),
          std::vector<size_t>(num_q_bins, 0),
          std::vector<std::vector<real_t>>(
              num_partials, std::vector<real_t>(num_q_bins, static_cast<real_t>(0.0))),
          std::vector<std::vector<size_t>>(num_partials, std::vector<size_t>(num_q_bins, 0)),
          std::vector<real_t>(num_q_bins, static_cast<real_t>(0.0)),
          std::vector<std::vector<real_t>>(
              num_partials, std::vector<real_t>(num_q_bins, static_cast<real_t>(0.0))));
    });

    tbb::parallel_for(tbb::blocked_range<size_t>(0, q_vectors.size()),
                      [&](const tbb::blocked_range<size_t> &range) {
                        auto &[local_total_sum, local_total_count, local_partial_sums,
                               local_partial_counts, local_c_total_sum, local_c_partial_sums] =
                            local_ets.local();
                        std::vector<real_t> type_cos(type_blocks.size(), static_cast<real_t>(0.0));
                        std::vector<real_t> type_sin(type_blocks.size(), static_cast<real_t>(0.0));

                        for (size_t qi = range.begin(); qi != range.end(); ++qi) {
                          processSingleQVector(q_vectors[qi], binning, num_atoms, basis,
                                               type_blocks, partials_info,
                                               {
                                                   .e1_cos = e1_cos.data(),
                                                   .e1_sin = e1_sin.data(),
                                                   .e2_cos = e2_cos.data(),
                                                   .e2_sin = e2_sin.data(),
                                                   .e3_cos = e3_cos.data(),
                                                   .e3_sin = e3_sin.data(),
                                               },
                                               {
                                                   .total_sum = &local_total_sum,
                                                   .total_count = &local_total_count,
                                                   .partial_sums = &local_partial_sums,
                                                   .partial_counts = &local_partial_counts,
                                                   .c_total_sum = &local_c_total_sum,
                                                   .c_partial_sums = &local_c_partial_sums,
                                               },
                                               type_cos, type_sin);
                        }
                      });

    combineAccumulators(local_ets, num_q_bins, num_partials, sq_total_sum, sq_total_count,
                        partial_sums, partial_counts);
  }

  normalizeAndStoreResults(s_q_hist, sq_total_sum, sq_total_count, partial_sums, partial_counts,
                           partials_info, num_q_bins);

  dists.addHistogram("S_q", std::move(s_q_hist));
}

} // namespace correlation::calculators
