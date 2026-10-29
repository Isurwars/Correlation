/**
 * @file VulkanDevice.hpp
 * @brief RAII management of Vulkan instance, physical device, logical device, and command pool.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#ifdef CORRELATION_USE_VULKAN
#include <vulkan/vulkan.h>
#endif

#include <cstdint>
#include <string>

namespace correlation::render3d::vulkan {

#ifdef CORRELATION_USE_VULKAN

/**
 * @class VulkanDevice
 * @brief Manages Vulkan instance, GPU device selection, logical queues, and command allocation.
 */
class VulkanDevice {
public:
  VulkanDevice() noexcept = default;
  ~VulkanDevice() noexcept { destroy(); }

  VulkanDevice(const VulkanDevice &) = delete;
  VulkanDevice &operator=(const VulkanDevice &) = delete;

  VulkanDevice(VulkanDevice &&other) noexcept { moveFrom(other); }
  VulkanDevice &operator=(VulkanDevice &&other) noexcept {
    if (this != &other) {
      destroy();
      moveFrom(other);
    }
    return *this;
  }

  /**
   * @brief Discovers GPUs and initializes Vulkan instance, device, and command pool.
   * @return True if a compatible physical device was found and initialized.
   */
  [[nodiscard]] bool initialize();

  /**
   * @brief Frees all allocated Vulkan device and instance resources in safe reverse order.
   */
  void destroy() noexcept;

  /**
   * @brief Allocates and begins a transient primary command buffer for transfer/setup.
   * @return Begun command buffer handle.
   */
  [[nodiscard]] VkCommandBuffer beginSingleTimeCommands();

  /**
   * @brief Ends, submits, and waits for a transient command buffer to complete.
   * @param command_buffer Command buffer to submit and free.
   */
  void endSingleTimeCommands(VkCommandBuffer command_buffer);

  [[nodiscard]] VkInstance instance() const noexcept { return instance_; }
  [[nodiscard]] VkPhysicalDevice physicalDevice() const noexcept { return physical_device_; }
  [[nodiscard]] VkDevice device() const noexcept { return device_; }
  [[nodiscard]] VkQueue graphicsQueue() const noexcept { return graphics_queue_; }
  [[nodiscard]] uint32_t queueFamilyIndex() const noexcept { return queue_family_index_; }
  [[nodiscard]] VkCommandPool commandPool() const noexcept { return command_pool_; }
  [[nodiscard]] bool isValid() const noexcept { return device_ != VK_NULL_HANDLE; }
  [[nodiscard]] std::string physicalDeviceName() const;

private:
  void moveFrom(VulkanDevice &other) noexcept;
  [[nodiscard]] bool createInstance();
  [[nodiscard]] bool pickPhysicalDevice();
  [[nodiscard]] bool createLogicalDevice();
  [[nodiscard]] bool createCommandPool();

  VkInstance instance_{VK_NULL_HANDLE};
  VkPhysicalDevice physical_device_{VK_NULL_HANDLE};
  VkDevice device_{VK_NULL_HANDLE};
  VkQueue graphics_queue_{VK_NULL_HANDLE};
  uint32_t queue_family_index_{0};
  VkCommandPool command_pool_{VK_NULL_HANDLE};
};

#endif // CORRELATION_USE_VULKAN

} // namespace correlation::render3d::vulkan
