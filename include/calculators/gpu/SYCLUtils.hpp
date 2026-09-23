/**
 * @file SYCLUtils.hpp
 * @brief Multi-vendor GPU utilities and queue management using SYCL/oneAPI.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#ifdef CORRELATION_USE_SYCL
#include <sycl/sycl.hpp>
#endif

namespace correlation::calculators::sycl_gpu {

/**
 * @brief Checks if a compatible SYCL GPU device (NVIDIA, AMD, Intel) is available at runtime.
 */
[[nodiscard]] bool hasSyclGpuDevice() noexcept;

#ifdef CORRELATION_USE_SYCL
/**
 * @brief Returns a reference to the shared SYCL queue instance.
 */
[[nodiscard]] sycl::queue &get_sycl_queue();
#endif

} // namespace correlation::calculators::sycl_gpu
