/**
 * @file VulkanContext.cpp
 * @brief Implementation of VulkanContext runtime device probing.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/vulkan/VulkanContext.hpp"

#ifdef CORRELATION_USE_VULKAN
#include <cstdint>
#include <vector>
#include <vulkan/vulkan.h>
#endif

namespace correlation::render3d {

#ifdef CORRELATION_USE_VULKAN

bool VulkanContext::probeVulkanDevice() noexcept {
  static int cached_result = -1;
  if (cached_result != -1) {
    return cached_result == 1;
  }

  const VkApplicationInfo app_info{
      .sType = VK_STRUCTURE_TYPE_APPLICATION_INFO,
      .pNext = nullptr,
      .pApplicationName = "CorrelationProbe",
      .applicationVersion = VK_MAKE_VERSION(1, 0, 0),
      .pEngineName = "Probe",
      .engineVersion = VK_MAKE_VERSION(1, 0, 0),
      .apiVersion = VK_API_VERSION_1_2,
  };

  const VkInstanceCreateInfo create_info{
      .sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .pApplicationInfo = &app_info,
      .enabledLayerCount = 0,
      .ppEnabledLayerNames = nullptr,
      .enabledExtensionCount = 0,
      .ppEnabledExtensionNames = nullptr,
  };

  VkInstance instance = VK_NULL_HANDLE;
  if (vkCreateInstance(&create_info, nullptr, &instance) != VK_SUCCESS) {
    cached_result = 0;
    return false;
  }

  uint32_t device_count = 0;
  bool found_compatible = false;
  if (vkEnumeratePhysicalDevices(instance, &device_count, nullptr) == VK_SUCCESS &&
      device_count > 0) {
    std::vector<VkPhysicalDevice> devices(device_count);
    if (vkEnumeratePhysicalDevices(instance, &device_count, devices.data()) == VK_SUCCESS) {
      for (VkPhysicalDevice dev : devices) {
        uint32_t queue_count = 0;
        vkGetPhysicalDeviceQueueFamilyProperties(dev, &queue_count, nullptr);
        std::vector<VkQueueFamilyProperties> queue_families(queue_count);
        vkGetPhysicalDeviceQueueFamilyProperties(dev, &queue_count, queue_families.data());

        for (const auto &family : queue_families) {
          if ((family.queueFlags & VK_QUEUE_GRAPHICS_BIT) != 0) {
            found_compatible = true;
            break;
          }
        }
        if (found_compatible) {
          break;
        }
      }
    }
  }

  vkDestroyInstance(instance, nullptr);
  cached_result = found_compatible ? 1 : 0;
  return found_compatible;
}

#endif // CORRELATION_USE_VULKAN

} // namespace correlation::render3d
