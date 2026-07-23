/**
 * @file SIMDDistance.hpp
 * @brief SIMD-accelerated distance and position block calculation routines.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/simd/SIMDTypes.hpp"
#include "math/simd/detail/AVX2Kernels.hpp"
#include "math/simd/detail/AVX512Kernels.hpp" // IWYU pragma: export
#include "math/simd/detail/ScalarKernels.hpp" // IWYU pragma: export

namespace correlation::math {

/**
 * @brief Computes squared distances from a reference point to a SoA block of positions.
 * @tparam T Floating-point coordinate type (float or double).
 */
template <typename T>
inline void compute_dsq_block(T ref_x, T ref_y, T ref_z, const PositionBlockT<T> &block,
                              T *CORRELATION_RESTRICT out_dsq) noexcept {
#ifdef CORRELATION_SIMD_AVX512
  detail::avx512::compute_dsq_block(ref_x, ref_y, ref_z, block, out_dsq);
#elif defined(CORRELATION_SIMD_AVX2)
  detail::avx2::compute_dsq_block(ref_x, ref_y, ref_z, block, out_dsq);
#else
  detail::scalar::compute_dsq_block(ref_x, ref_y, ref_z, block, out_dsq);
#endif
}

/**
 * @brief Populates vectors with standard atom block positions for SIMD.
 * @tparam AtomRange A range of atom objects (e.g., std::vector<Atom>).
 * @tparam T Coordinate scalar type (defaults to real_t).
 */
template <typename AtomRange, typename T = real_t>
inline std::size_t fill_position_block(const FillPositionBlockParams<AtomRange, T> &params) noexcept {
  if (params.atoms == nullptr || params.x_s == nullptr || params.y_s == nullptr || params.z_s == nullptr) {
    return 0;
  }
  const std::size_t count = params.end_idx - params.begin_idx;
  params.x_s->resize(count);
  params.y_s->resize(count);
  params.z_s->resize(count);
  for (std::size_t idx = 0; idx < count; ++idx) {
    const auto &pos = (*params.atoms)[params.begin_idx + idx].position();
    (*params.x_s)[idx] = static_cast<T>(pos.x());
    (*params.y_s)[idx] = static_cast<T>(pos.y());
    (*params.z_s)[idx] = static_cast<T>(pos.z());
  }
  return count;
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
template <typename AtomRange, typename T = real_t>
inline std::size_t fill_position_block(const AtomRange &atoms, std::size_t begin_idx, std::size_t end_idx,
                                       std::vector<T> &x_s, std::vector<T> &y_s, std::vector<T> &z_s) noexcept {
  return fill_position_block(FillPositionBlockParams<AtomRange, T>{
      .atoms = &atoms, .begin_idx = begin_idx, .end_idx = end_idx, .x_s = &x_s, .y_s = &y_s, .z_s = &z_s});
}

} // namespace correlation::math
