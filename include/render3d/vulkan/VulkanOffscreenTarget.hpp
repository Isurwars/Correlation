/**
 * @file VulkanOffscreenTarget.hpp
 * @brief Headless offscreen render target, color/depth framebuffers, and staging readback.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#ifdef CORRELATION_USE_VULKAN
#include <vulkan/vulkan.h>
#endif

#include "render3d/vulkan/VulkanBuffer.hpp"
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>

namespace correlation::render3d::vulkan {

#ifdef CORRELATION_USE_VULKAN

/**
 * @class VulkanOffscreenTarget
 * @brief Manages offscreen framebuffer, depth/color attachments, and double-buffered staging
 * readout.
 */
class VulkanOffscreenTarget {
public:
  static constexpr std::size_t STAGING_BUFFER_COUNT = 2;

  VulkanOffscreenTarget() noexcept = default;
  ~VulkanOffscreenTarget() noexcept { destroy(); }

  VulkanOffscreenTarget(const VulkanOffscreenTarget &) = delete;
  VulkanOffscreenTarget &operator=(const VulkanOffscreenTarget &) = delete;

  VulkanOffscreenTarget(VulkanOffscreenTarget &&other) noexcept { moveFrom(other); }
  VulkanOffscreenTarget &operator=(VulkanOffscreenTarget &&other) noexcept {
    if (this != &other) {
      destroy();
      moveFrom(other);
    }
    return *this;
  }

  /**
   * @brief Creates the offscreen target with specified device handles and resolution.
   * @param device Vulkan logical device handle.
   * @param physical_device Vulkan physical device handle.
   * @param width Viewport width in pixels.
   * @param height Viewport height in pixels.
   */
  void create(VkDevice device, VkPhysicalDevice physical_device, uint32_t width, uint32_t height);

  /**
   * @brief Resizes the target attachments and recreates framebuffers if dimensions change.
   * @param width New viewport width in pixels.
   * @param height New viewport height in pixels.
   */
  void resize(uint32_t width, uint32_t height);

  /**
   * @brief Frees all Vulkan render pass, framebuffer, image, and staging buffer resources.
   */
  void destroy() noexcept;

  /**
   * @brief Begins the offscreen render pass, setting viewports and clearing attachments.
   * @param cmd Vulkan command buffer in recording state.
   * @param clear_color Clear color for the RGBA color attachment.
   */
  void beginRenderPass(VkCommandBuffer cmd,
                       VkClearColorValue clear_color = {{0.1F, 0.1F, 0.12F, 1.0F}});

  /**
   * @brief Ends the active offscreen render pass.
   * @param cmd Vulkan command buffer in recording state.
   */
  static void endRenderPass(VkCommandBuffer cmd);

  /**
   * @brief Records a transfer command copying the rendered color image to a host-visible staging
   * buffer.
   * @param cmd Vulkan command buffer in recording state.
   * @param staging_index Ring index of the staging buffer to receive pixel data.
   */
  void copyToStaging(VkCommandBuffer cmd, uint32_t staging_index);

  /**
   * @brief Obtains a mapped byte span of the requested staging buffer's pixel contents.
   * @param staging_index Ring index of the staging buffer to read.
   * @return Contiguous span of RGBA8 pixels, or empty span if invalid.
   */
  [[nodiscard]] std::span<const uint8_t> readoutStagingSpan(uint32_t staging_index) const noexcept;

  [[nodiscard]] VkRenderPass renderPass() const noexcept { return render_pass_; }
  [[nodiscard]] VkFramebuffer framebuffer() const noexcept { return framebuffer_; }
  [[nodiscard]] VkExtent2D extent() const noexcept { return extent_; }
  [[nodiscard]] uint32_t width() const noexcept { return extent_.width; }
  [[nodiscard]] uint32_t height() const noexcept { return extent_.height; }
  [[nodiscard]] VkFormat colorFormat() const noexcept { return color_format_; }
  [[nodiscard]] VkFormat depthFormat() const noexcept { return depth_format_; }
  [[nodiscard]] VkImage colorImage() const noexcept { return color_image_; }
  [[nodiscard]] VkImage depthImage() const noexcept { return depth_image_; }
  [[nodiscard]] bool isValid() const noexcept { return framebuffer_ != VK_NULL_HANDLE; }

private:
  void moveFrom(VulkanOffscreenTarget &other) noexcept;
  void createRenderPass();
  void createColorResources();
  void createDepthResources();
  void createFramebuffer();
  void createStagingBuffers();
  void destroyFrameResources() noexcept;

  [[nodiscard]] VkFormat findDepthFormat() const;
  void createImage(uint32_t width, uint32_t height, VkFormat format, VkImageTiling tiling,
                   VkImageUsageFlags usage, VkMemoryPropertyFlags properties, VkImage &image,
                   VkDeviceMemory &memory);
  [[nodiscard]] VkImageView createImageView(VkImage image, VkFormat format,
                                            VkImageAspectFlags aspect_flags) const;

  VkDevice device_{VK_NULL_HANDLE};
  VkPhysicalDevice physical_device_{VK_NULL_HANDLE};
  VkExtent2D extent_{.width = 0, .height = 0};

  VkFormat color_format_{VK_FORMAT_R8G8B8A8_UNORM};
  VkFormat depth_format_{VK_FORMAT_UNDEFINED};

  VkRenderPass render_pass_{VK_NULL_HANDLE};
  VkImage color_image_{VK_NULL_HANDLE};
  VkDeviceMemory color_memory_{VK_NULL_HANDLE};
  VkImageView color_image_view_{VK_NULL_HANDLE};

  VkImage depth_image_{VK_NULL_HANDLE};
  VkDeviceMemory depth_memory_{VK_NULL_HANDLE};
  VkImageView depth_image_view_{VK_NULL_HANDLE};

  VkFramebuffer framebuffer_{VK_NULL_HANDLE};

  std::array<VulkanBuffer, STAGING_BUFFER_COUNT> staging_buffers_{};
};

#endif // CORRELATION_USE_VULKAN

} // namespace correlation::render3d::vulkan
