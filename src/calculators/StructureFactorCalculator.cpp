/**
 * @file StructureFactorCalculator.cpp
 * @brief Implementation of structure factor S(Q) calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "calculators/StructureFactorCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Constants.hpp"
#include "math/SIMDUtils.hpp"

#include <cmath>
#include <map>
#include <stdexcept>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <vector>

namespace {
bool registered = CalculatorFactory::instance().registerCalculator(
    std::make_unique<StructureFactorCalculator>());
} // namespace

void StructureFactorCalculator::calculateFrame(
    DistributionFunctions &df, const AnalysisSettings &settings) const {
  if (settings.q_bin_width <= 0 || settings.q_max <= 0) {
    throw std::invalid_argument("Q-space parameters must be positive.");
  }

  const Cell &cell = df.cell();
  const auto &atoms = cell.atoms();
  const size_t N = atoms.size();
  if (N == 0)
    return;

  const double q_max = settings.q_max;
  const double q_bin_width = settings.q_bin_width;
  const size_t num_q_bins =
      static_cast<size_t>(std::floor(q_max / q_bin_width));
  if (num_q_bins == 0) {
    throw std::invalid_argument(
        "Q_max is too small for the given Q_bin_width.");
  }

  // -----------------------------------------------------------------------
  // Build the reciprocal lattice basis vectors.
  // The Cell's inverseLatticeVectors() gives M^{-1} such that M * M^{-1} = I
  // where M is the lattice matrix. The reciprocal vectors are b_i = 2*pi *
  // (M^{-T})_i. We extract each reciprocal basis vector from the transpose of
  // the inverse.
  // -----------------------------------------------------------------------
  const auto &inv = cell.inverseLatticeVectors();
  const double bx_x = correlation::math::two_pi * inv(0, 0);
  const double bx_y = correlation::math::two_pi * inv(0, 1);
  const double bx_z = correlation::math::two_pi * inv(0, 2);
  const double by_x = correlation::math::two_pi * inv(1, 0);
  const double by_y = correlation::math::two_pi * inv(1, 1);
  const double by_z = correlation::math::two_pi * inv(1, 2);
  const double bz_x = correlation::math::two_pi * inv(2, 0);
  const double bz_y = correlation::math::two_pi * inv(2, 1);
  const double bz_z = correlation::math::two_pi * inv(2, 2);

  // -----------------------------------------------------------------------
  // Generate all (h, k, l) reciprocal lattice vectors with |q| <= q_max.
  // -----------------------------------------------------------------------
  const double b1_norm = std::sqrt(bx_x * bx_x + bx_y * bx_y + bx_z * bx_z);
  const double b2_norm = std::sqrt(by_x * by_x + by_y * by_y + by_z * by_z);
  const double b3_norm = std::sqrt(bz_x * bz_x + bz_y * bz_y + bz_z * bz_z);

  const int hmax =
      (b1_norm > 1e-10) ? static_cast<int>(std::ceil(q_max / b1_norm)) : 0;
  const int kmax =
      (b2_norm > 1e-10) ? static_cast<int>(std::ceil(q_max / b2_norm)) : 0;
  const int lmax =
      (b3_norm > 1e-10) ? static_cast<int>(std::ceil(q_max / b3_norm)) : 0;

  struct QVector {
    int h, k, l;
    double qmag;
  };
  std::vector<QVector> q_vectors;
  q_vectors.reserve(8 * hmax * kmax * lmax + 1);

  const double q_max_sq = q_max * q_max;
  for (int h = -hmax; h <= hmax; ++h) {
    for (int k = -kmax; k <= kmax; ++k) {
      for (int l = -lmax; l <= lmax; ++l) {
        if (h == 0 && k == 0 && l == 0)
          continue;
        const double qx = h * bx_x + k * by_x + l * bz_x;
        const double qy = h * bx_y + k * by_y + l * bz_y;
        const double qz = h * bx_z + k * by_z + l * bz_z;
        const double qmag_sq = qx * qx + qy * qy + qz * qz;
        if (qmag_sq <= q_max_sq) {
          q_vectors.push_back({h, k, l, std::sqrt(qmag_sq)});
        }
      }
    }
  }

  if (q_vectors.empty()) {
    return;
  }

  // -----------------------------------------------------------------------
  // Group atom indices by element type for partial S(Q) calculation.
  // -----------------------------------------------------------------------
  std::map<std::string, std::vector<size_t>> indices_by_type;
  for (size_t j = 0; j < N; ++j) {
    indices_by_type[atoms[j].element().symbol].push_back(j);
  }

  // Collect element pairs for partial S(Q) output.
  struct PartialInfo {
    std::string key;
    size_t typeA_idx; // offset into flat per-type phase arrays
    size_t typeB_idx;
    size_t N_A;
    size_t N_B;
    bool is_identical;
    double weight;
  };
  std::vector<PartialInfo> partials_info;
  const auto &ashcroft_weights = df.getAshcroftWeights();

  // Build a flat array of per-type start offsets into the global xs/ys/zs.
  // We'll store positions contiguously: all of type0, then type1, etc.
  struct TypeBlock {
    std::string symbol;
    size_t offset; // start index in flat position arrays
    size_t count;
  };
  std::vector<TypeBlock> type_blocks;
  type_blocks.reserve(indices_by_type.size());
  std::vector<double> xs(N), ys(N), zs(N);

  size_t flat_offset = 0;
  for (const auto &[sym, idxs] : indices_by_type) {
    TypeBlock tb;
    tb.symbol = sym;
    tb.offset = flat_offset;
    tb.count = idxs.size();
    for (size_t j = 0; j < idxs.size(); ++j) {
      const auto &pos = atoms[idxs[j]].position();
      xs[flat_offset + j] = pos.x();
      ys[flat_offset + j] = pos.y();
      zs[flat_offset + j] = pos.z();
    }
    flat_offset += idxs.size();
    type_blocks.push_back(tb);
  }

  // Build partials_info for every (i >= j) pair of type blocks.
  for (size_t ti = 0; ti < type_blocks.size(); ++ti) {
    for (size_t tj = ti; tj < type_blocks.size(); ++tj) {
      const auto &tbA = type_blocks[ti];
      const auto &tbB = type_blocks[tj];
      bool is_identical = (ti == tj);
      std::string key = (tbA.symbol < tbB.symbol)
                            ? (tbA.symbol + "-" + tbB.symbol)
                            : (tbB.symbol + "-" + tbA.symbol);
      double w = 0.0;
      auto wit = ashcroft_weights.find(key);
      if (wit != ashcroft_weights.end())
        w = wit->second;
      partials_info.push_back(
          {key, ti, tj, tbA.count, tbB.count, is_identical, w});
    }
  }

  // -----------------------------------------------------------------------
  // Precompute phase factors for Miller index decomposition.
  // -----------------------------------------------------------------------
  auto precompute_phases = [&](int max_idx, double bx, double by, double bz,
                                std::vector<double> &E_cos,
                                std::vector<double> &E_sin) {
    E_cos.resize((2 * max_idx + 1) * N);
    E_sin.resize((2 * max_idx + 1) * N);
    for (int h = -max_idx; h <= max_idx; ++h) {
      double *cos_h = &E_cos[(h + max_idx) * N];
      double *sin_h = &E_sin[(h + max_idx) * N];
      for (size_t j = 0; j < N; ++j) {
        double phase = h * (bx * xs[j] + by * ys[j] + bz * zs[j]);
        cos_h[j] = std::cos(phase);
        sin_h[j] = std::sin(phase);
      }
    }
  };

  std::vector<double> E1_cos, E1_sin, E2_cos, E2_sin, E3_cos, E3_sin;
  precompute_phases(hmax, bx_x, bx_y, bx_z, E1_cos, E1_sin);
  precompute_phases(kmax, by_x, by_y, by_z, E2_cos, E2_sin);
  precompute_phases(lmax, bz_x, bz_y, bz_z, E3_cos, E3_sin);

  // -----------------------------------------------------------------------
  // Allocate output histogram with partials.
  // -----------------------------------------------------------------------
  Histogram s_q_hist;
  s_q_hist.bin_label = "Q (Å⁻¹)";
  s_q_hist.bins.resize(num_q_bins);
  for (size_t i = 0; i < num_q_bins; ++i)
    s_q_hist.bins[i] = (i + 0.5) * q_bin_width;

  for (const auto &pi : partials_info)
    s_q_hist.partials[pi.key].assign(num_q_bins, 0.0);

  std::vector<double> sq_total_sum(num_q_bins, 0.0);
  std::vector<size_t> sq_total_count(num_q_bins, 0);

  // Per-partial accumulators: sum and count.
  const size_t np = partials_info.size();
  std::vector<std::vector<double>> partial_sums(np,
                                                std::vector<double>(num_q_bins, 0.0));
  std::vector<std::vector<size_t>> partial_counts(
      np, std::vector<size_t>(num_q_bins, 0));

  // Thread-local storage: one pair<sum, count> per partial + one for total.
  using LocalPartials =
      std::pair<std::vector<std::vector<double>>, std::vector<std::vector<size_t>>>;
  tbb::enumerable_thread_specific<
      std::tuple<std::vector<double>, std::vector<size_t>, // total sum/count
                 std::vector<std::vector<double>>,          // partial sums
                 std::vector<std::vector<size_t>>>>         // partial counts
      local_ets([&] {
        return std::make_tuple(
            std::vector<double>(num_q_bins, 0.0),
            std::vector<size_t>(num_q_bins, 0),
            std::vector<std::vector<double>>(np,
                                             std::vector<double>(num_q_bins, 0.0)),
            std::vector<std::vector<size_t>>(
                np, std::vector<size_t>(num_q_bins, 0)));
      });

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, q_vectors.size()),
      [&](const tbb::blocked_range<size_t> &range) {
        auto &[local_total_sum, local_total_count, local_partial_sums,
               local_partial_counts] = local_ets.local();

        for (size_t qi = range.begin(); qi != range.end(); ++qi) {
          const auto &q = q_vectors[qi];
          const size_t bin = static_cast<size_t>(q.qmag / q_bin_width);
          if (bin >= num_q_bins)
            continue;

          // Compute cos/sin sums for each type block by multiplying the
          // per-index contributions (already factored by Miller h,k,l).
          std::vector<double> type_cos(type_blocks.size(), 0.0);
          std::vector<double> type_sin(type_blocks.size(), 0.0);

          for (size_t ti = 0; ti < type_blocks.size(); ++ti) {
            const size_t off = type_blocks[ti].offset;
            const size_t cnt = type_blocks[ti].count;
            double c_sum = 0.0, s_sum = 0.0;
            correlation::math::miller_phase_sum(
                &E1_cos[(q.h + hmax) * N] + off,
                &E1_sin[(q.h + hmax) * N] + off,
                &E2_cos[(q.k + kmax) * N] + off,
                &E2_sin[(q.k + kmax) * N] + off,
                &E3_cos[(q.l + lmax) * N] + off,
                &E3_sin[(q.l + lmax) * N] + off, cnt, c_sum, s_sum);
            type_cos[ti] = c_sum;
            type_sin[ti] = s_sum;
          }

          // Accumulate total S(Q) via the full sum |rho(q)|^2 / N.
          double full_cos = 0.0, full_sin = 0.0;
          for (size_t ti = 0; ti < type_blocks.size(); ++ti) {
            full_cos += type_cos[ti];
            full_sin += type_sin[ti];
          }
          const double sq_total =
              (full_cos * full_cos + full_sin * full_sin) /
              static_cast<double>(N);
          local_total_sum[bin] += sq_total;
          local_total_count[bin] += 1;

          // Per-partial: S_AB(Q) = Re[rho_A(q)* · rho_B(q)] / sqrt(N_A * N_B)
          for (size_t pi = 0; pi < np; ++pi) {
            const auto &pinfo = partials_info[pi];
            const size_t tiA = pinfo.typeA_idx;
            const size_t tiB = pinfo.typeB_idx;
            // Re[rho_A* rho_B] = cosA*cosB + sinA*sinB
            double cross = type_cos[tiA] * type_cos[tiB] +
                           type_sin[tiA] * type_sin[tiB];
            double denom = std::sqrt(static_cast<double>(pinfo.N_A) *
                                     static_cast<double>(pinfo.N_B));
            double sq_partial = (denom > 0.0) ? cross / denom : 0.0;
            local_partial_sums[pi][bin] += sq_partial;
            local_partial_counts[pi][bin] += 1;
          }
        }
      });

  // Reduce thread-local results.
  local_ets.combine_each([&](const auto &local) {
    const auto &[lt_sum, lt_count, lp_sums, lp_counts] = local;
    for (size_t i = 0; i < num_q_bins; ++i) {
      sq_total_sum[i] += lt_sum[i];
      sq_total_count[i] += lt_count[i];
    }
    for (size_t pi = 0; pi < np; ++pi) {
      for (size_t i = 0; i < num_q_bins; ++i) {
        partial_sums[pi][i] += lp_sums[pi][i];
        partial_counts[pi][i] += lp_counts[pi][i];
      }
    }
  });

  // Average and store.
  std::vector<double> total_sq(num_q_bins, 0.0);
  for (size_t i = 0; i < num_q_bins; ++i) {
    if (sq_total_count[i] > 0)
      total_sq[i] = sq_total_sum[i] / static_cast<double>(sq_total_count[i]);
  }

  for (size_t pi = 0; pi < np; ++pi) {
    auto &out = s_q_hist.partials[partials_info[pi].key];
    for (size_t i = 0; i < num_q_bins; ++i) {
      if (partial_counts[pi][i] > 0)
        out[i] = partial_sums[pi][i] /
                 static_cast<double>(partial_counts[pi][i]);
    }
  }

  s_q_hist.partials["Total"] = std::move(total_sq);
  df.addHistogram("S_q", std::move(s_q_hist));
}
