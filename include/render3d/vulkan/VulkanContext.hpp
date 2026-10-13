/**
 * @file VulkanContext.hpp
 * @brief Headless Vulkan runtime device probing and capability detection.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include <string_view>

namespace correlation::render3d {

/**
 * @class VulkanContext
 * @brief Inspects host Vulkan drivers and physical devices for headless rendering.
 */
class VulkanContext {
public:
  /**
   * @brief Checks if a functional Vulkan driver and physical device are present.
   * @return True if headless Vulkan offscreen rendering is supported.
   */
  [[nodiscard]] static bool isSupported() noexcept {
#if defined(CORRELATION_USE_VULKAN)
    return probeVulkanDevice();
#else
    return false;
#endif
  }

  /**
   * @brief Descriptive status of the Vulkan subsystem.
   * @return Status description string.
   */
  [[nodiscard]] static std::string_view statusMessage() noexcept {
#if defined(CORRELATION_USE_VULKAN)
    return probeVulkanDevice() ? "Vulkan headless GPU available"
                               : "Vulkan drivers found but no compatible physical device";
#else
    return "Vulkan disabled at build time (CPU rasterizer active)";
#endif
  }

private:
#if defined(CORRELATION_USE_VULKAN)
  [[nodiscard]] static bool probeVulkanDevice() noexcept;
#endif
};

} // namespace correlation::render3d
