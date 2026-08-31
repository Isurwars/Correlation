/**
 * @file DeviceBuffer.hpp
 * @brief RAII wrapper for GPU device memory allocations (CUDA/HIP).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/GPUErrorCheck.hpp"
#include "core/GPUPortability.hpp"

#include <cstddef>

namespace correlation::core::gpu {

/**
 * @brief RAII wrapper for a typed GPU device allocation.
 *
 * Allocates device memory on construction and releases it on destruction.
 * Move-only — cannot be copied. Replaces raw hipMalloc/hipFree pairs to
 * guarantee exception-safe cleanup.
 *
 * @tparam T Element type stored in the device buffer.
 */
template <typename T> class DeviceBuffer {
public:
  DeviceBuffer() = default;

  /**
   * @brief Allocates device memory for @p count elements of type T.
   * @param[in] count Number of elements to allocate.
   * @throws GPUError If the device allocation fails.
   */
  explicit DeviceBuffer(std::size_t count) : count_(count) {
    if (count_ > 0) {
      hipCheck(hipMalloc(&ptr_, count_ * sizeof(T)));
    }
  }

  ~DeviceBuffer() { release(); }

  // Move-only semantics
  DeviceBuffer(DeviceBuffer &&other) noexcept : ptr_(other.ptr_), count_(other.count_) {
    other.ptr_ = nullptr;
    other.count_ = 0;
  }

  DeviceBuffer &operator=(DeviceBuffer &&other) noexcept {
    if (this != &other) {
      release();
      ptr_ = other.ptr_;
      count_ = other.count_;
      other.ptr_ = nullptr;
      other.count_ = 0;
    }
    return *this;
  }

  DeviceBuffer(const DeviceBuffer &) = delete;
  DeviceBuffer &operator=(const DeviceBuffer &) = delete;

  /**
   * @brief Copies host data to this device buffer.
   * @param[in] host_data Pointer to host source data.
   * @param[in] num_elements Number of elements to copy (must be ≤ count_).
   * @throws GPUError If the memory copy fails.
   */
  void copyFromHost(const T *host_data, std::size_t num_elements) {
    hipCheck(hipMemcpy(ptr_, host_data, num_elements * sizeof(T), hipMemcpyHostToDevice));
  }

  /**
   * @brief Copies device data back to a host buffer.
   * @param[out] host_data Pointer to host destination buffer.
   * @param[in] num_elements Number of elements to copy (must be ≤ count_).
   * @throws GPUError If the memory copy fails.
   */
  void copyToHost(T *host_data, std::size_t num_elements) const {
    hipCheck(hipMemcpy(host_data, ptr_, num_elements * sizeof(T), hipMemcpyDeviceToHost));
  }

  /**
   * @brief Sets the entire buffer to a value via host-to-device copy of a single element.
   * @param[in] value The value to write (copied to device for the first element).
   */
  void setScalar(const T &value) { hipCheck(hipMemcpy(ptr_, &value, sizeof(T), hipMemcpyHostToDevice)); }

  /**
   * @brief Returns the raw device pointer.
   * @return Device pointer (may be nullptr if default-constructed or moved-from).
   */
  [[nodiscard]] T *get() noexcept { return ptr_; }

  /**
   * @brief Returns the raw device pointer (const overload).
   * @return Const device pointer.
   */
  [[nodiscard]] const T *get() const noexcept { return ptr_; }

  /**
   * @brief Returns the number of allocated elements.
   * @return Element count.
   */
  [[nodiscard]] std::size_t size() const noexcept { return count_; }

private:
  void release() noexcept {
    if (ptr_ != nullptr) {
      hipFree(ptr_);
      ptr_ = nullptr;
      count_ = 0;
    }
  }

  T *ptr_{nullptr};
  std::size_t count_{0};
};

} // namespace correlation::core::gpu
