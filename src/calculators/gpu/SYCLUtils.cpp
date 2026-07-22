/**
 * @file SYCLUtils.cpp
 * @brief Implementation of SYCL device selection and queue management.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/gpu/SYCLUtils.hpp"

namespace correlation::calculators::sycl_gpu {

bool has_sycl_gpu_device() noexcept {
#if defined(CORRELATION_USE_SYCL)
  try {
    const auto devices = sycl::device::get_devices(sycl::info::device_type::gpu);
    return !devices.empty();
  } catch (...) {
    return false;
  }
#else
  return false;
#endif
}

#if defined(CORRELATION_USE_SYCL)
sycl::queue &get_sycl_queue() {
  static sycl::queue queue(sycl::gpu_selector_v, sycl::property::queue::in_order());
  return queue;
}
#endif

} // namespace correlation::calculators::sycl_gpu
