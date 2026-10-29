/**
 * @file VulkanCommon.hpp
 * @brief Common types, macros, and RAII helpers for the Vulkan rendering subsystem.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#ifdef CORRELATION_USE_VULKAN
#include <vulkan/vulkan.h>
#endif

#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>

namespace correlation::render3d::vulkan {

#ifdef CORRELATION_USE_VULKAN

/**
 * @brief Evaluates a VkResult and throws std::runtime_error upon failure.
 * @param result Vulkan status result code.
 * @param message Context description for the exception message.
 */
inline void checkVkResult(VkResult result, std::string_view message) {
  if (result != VK_SUCCESS) {
    throw std::runtime_error(std::string(message) +
                             " (VkResult: " + std::to_string(static_cast<int>(result)) + ")");
  }
}

/**
 * @brief Locates a suitable memory type index on a physical device.
 * @param physical_device Physical device handle.
 * @param type_filter Bitmask of permitted memory types.
 * @param properties Required memory property flags.
 * @return Memory type index matching the filter and properties.
 */
inline uint32_t findMemoryType(VkPhysicalDevice physical_device, uint32_t type_filter,
                               VkMemoryPropertyFlags properties) {
  VkPhysicalDeviceMemoryProperties mem_properties{};
  vkGetPhysicalDeviceMemoryProperties(physical_device, &mem_properties);

  const auto memory_types =
      std::span(mem_properties.memoryTypes).subspan(0, mem_properties.memoryTypeCount);
  for (uint32_t i = 0; i < mem_properties.memoryTypeCount; ++i) {
    const bool is_type_set = (type_filter & (1U << i)) != 0;
    const bool has_properties = (memory_types[i].propertyFlags & properties) == properties;
    if (is_type_set && has_properties) {
      return i;
    }
  }

  throw std::runtime_error("Failed to find suitable Vulkan device memory type");
}

#endif // CORRELATION_USE_VULKAN

} // namespace correlation::render3d::vulkan
