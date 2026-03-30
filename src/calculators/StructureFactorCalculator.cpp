// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

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
  // inv stores M^{-1} column-major. Row i of inv is the i-th reciprocal vector
  // direction (without the 2*pi). So bx = 2*pi * row0(inv), etc.
  // Row 0 of inv = (inv(0,0), inv(0,1), inv(0,2))
  // inv(r,c) = inv column c, row r
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
  // Estimate max Miller index needed.
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
          continue; // skip Gamma point
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
  // Precompute phase factors for Miller index decomposition
  // arrays memory size: (2*hmax+1) * N doubles
  // -----------------------------------------------------------------------
  std::vector<double> xs(N), ys(N), zs(N);
  for (size_t j = 0; j < N; ++j) {
    const auto &pos = atoms[j].position();
    xs[j] = pos.x();
    ys[j] = pos.y();
    zs[j] = pos.z();
  }

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
  // For each q-vector, compute S(q) = |rho(q)|^2 / N, then bin by |q|.
  // Parallel over q-vectors.
  // -----------------------------------------------------------------------
  Histogram s_q_hist;
  s_q_hist.bin_label = "Q (Å⁻¹)";
  s_q_hist.bins.resize(num_q_bins);
  for (size_t i = 0; i < num_q_bins; ++i) {
    s_q_hist.bins[i] = (i + 0.5) * q_bin_width;
  }

  std::vector<double> sq_sum(num_q_bins, 0.0);
  std::vector<size_t> sq_count(num_q_bins, 0);

  tbb::enumerable_thread_specific<
      std::pair<std::vector<double>, std::vector<size_t>>>
      local_ets([&] {
        return std::make_pair(std::vector<double>(num_q_bins, 0.0),
                              std::vector<size_t>(num_q_bins, 0));
      });

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, q_vectors.size()),
      [&](const tbb::blocked_range<size_t> &range) {
        auto &[local_sum, local_count] = local_ets.local();
        for (size_t qi = range.begin(); qi != range.end(); ++qi) {
          const auto &q = q_vectors[qi];
          double cos_sum = 0.0, sin_sum = 0.0;

          correlation::math::miller_phase_sum(
              &E1_cos[(q.h + hmax) * N], &E1_sin[(q.h + hmax) * N],
              &E2_cos[(q.k + kmax) * N], &E2_sin[(q.k + kmax) * N],
              &E3_cos[(q.l + lmax) * N], &E3_sin[(q.l + lmax) * N], N, cos_sum,
              sin_sum);

          const double sq =
              (cos_sum * cos_sum + sin_sum * sin_sum) / static_cast<double>(N);
          const size_t bin = static_cast<size_t>(q.qmag / q_bin_width);
          if (bin < num_q_bins) {
            local_sum[bin] += sq;
            local_count[bin] += 1;
          }
        }
      });

  local_ets.combine_each([&](const auto &local) {
    const auto &[local_sum, local_count] = local;
    for (size_t i = 0; i < num_q_bins; ++i) {
      sq_sum[i] += local_sum[i];
      sq_count[i] += local_count[i];
    }
  });

  std::vector<double> total_sq(num_q_bins, 0.0);
  for (size_t i = 0; i < num_q_bins; ++i) {
    if (sq_count[i] > 0) {
      total_sq[i] = sq_sum[i] / static_cast<double>(sq_count[i]);
    }
  }

  s_q_hist.partials["Total"] = std::move(total_sq);
  df.addHistogram("S_q", std::move(s_q_hist));
}
