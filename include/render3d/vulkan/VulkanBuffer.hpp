/**
 * @file VulkanBuffer.hpp
 * @brief RAII encapsulation of a Vulkan buffer and its associated device memory.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#ifdef CORRELATION_USE_VULKAN
#include <vulkan/vulkan.h>
#endif

#include <cstddef>
#include <span>
#include <utility>

namespace correlation::render3d::vulkan {

#ifdef CORRELATION_USE_VULKAN

/**
 * @struct BufferAllocationDesc
 * @brief Parameters defining size, usage, and memory properties for a Vulkan buffer.
 */
struct BufferAllocationDesc {
  VkDeviceSize size{0};
  VkBufferUsageFlags usage{0};
  VkMemoryPropertyFlags properties{0};
};

/**
 * @class VulkanBuffer
 * @brief Manages lifecycle, allocation, mapping, and resizing of a VkBuffer.
 */
class VulkanBuffer {
public:
  VulkanBuffer() noexcept = default;
  ~VulkanBuffer() noexcept { destroy(); }

  VulkanBuffer(const VulkanBuffer &) = delete;
  VulkanBuffer &operator=(const VulkanBuffer &) = delete;

  VulkanBuffer(VulkanBuffer &&other) noexcept { moveFrom(other); }

  VulkanBuffer &operator=(VulkanBuffer &&other) noexcept {
    if (this != &other) {
      destroy();
      moveFrom(other);
    }
    return *this;
  }

  /**
   * @brief Allocates and binds a Vulkan buffer with specified descriptor.
   */
  void create(VkDevice device, VkPhysicalDevice physical_device, const BufferAllocationDesc &desc);

  /**
   * @brief Frees underlying buffer and device memory.
   */
  void destroy() noexcept;

  /**
   * @brief Maps buffer memory into host virtual address space.
   * @return Pointer to mapped memory.
   */
  [[nodiscard]] void *map();

  /**
   * @brief Unmaps host virtual address space.
   */
  void unmap() noexcept;

  /**
   * @brief Uploads contiguous byte data into host-visible mapped buffer.
   * @param data Byte span to copy.
   */
  void upload(std::span<const std::byte> data);

  /**
   * @brief Ensures buffer capacity meets required size, resizing to next power-of-two if needed.
   * @return True if buffer was reallocated.
   */
  bool ensureCapacity(VkDevice device, VkPhysicalDevice physical_device,
                      const BufferAllocationDesc &desc);

  [[nodiscard]] VkBuffer buffer() const noexcept { return buffer_; }
  [[nodiscard]] VkDeviceMemory memory() const noexcept { return memory_; }
  [[nodiscard]] VkDeviceSize size() const noexcept { return size_; }
  [[nodiscard]] void *mappedData() const noexcept { return mapped_data_; }
  [[nodiscard]] bool isValid() const noexcept { return buffer_ != VK_NULL_HANDLE; }

private:
  void moveFrom(VulkanBuffer &other) noexcept {
    device_ = std::exchange(other.device_, VK_NULL_HANDLE);
    buffer_ = std::exchange(other.buffer_, VK_NULL_HANDLE);
    memory_ = std::exchange(other.memory_, VK_NULL_HANDLE);
    size_ = std::exchange(other.size_, 0);
    mapped_data_ = std::exchange(other.mapped_data_, nullptr);
  }

  VkDevice device_{VK_NULL_HANDLE};
  VkBuffer buffer_{VK_NULL_HANDLE};
  VkDeviceMemory memory_{VK_NULL_HANDLE};
  VkDeviceSize size_{0};
  void *mapped_data_{nullptr};
};

#endif // CORRELATION_USE_VULKAN

} // namespace correlation::render3d::vulkan
