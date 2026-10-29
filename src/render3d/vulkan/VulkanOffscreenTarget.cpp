/**
 * @file VulkanOffscreenTarget.cpp
 * @brief Implementation of VulkanOffscreenTarget offscreen rendering and staging readout.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/vulkan/VulkanOffscreenTarget.hpp"

#ifdef CORRELATION_USE_VULKAN
#include "render3d/vulkan/VulkanCommon.hpp"
#include <algorithm>
#include <stdexcept>
#include <utility>

namespace correlation::render3d::vulkan {

void VulkanOffscreenTarget::create(VkDevice device, VkPhysicalDevice physical_device,
                                   uint32_t width, uint32_t height) {
  destroy();

  device_ = device;
  physical_device_ = physical_device;
  extent_ = {.width = std::max(width, 1U), .height = std::max(height, 1U)};

  depth_format_ = findDepthFormat();
  createRenderPass();
  createColorResources();
  createDepthResources();
  createFramebuffer();
  createStagingBuffers();
}

void VulkanOffscreenTarget::resize(uint32_t width, uint32_t height) {
  const uint32_t target_width = std::max(width, 1U);
  const uint32_t target_height = std::max(height, 1U);

  if (target_width == extent_.width && target_height == extent_.height && isValid()) {
    return;
  }

  if (device_ != VK_NULL_HANDLE) {
    vkDeviceWaitIdle(device_);
  }

  destroyFrameResources();
  extent_ = {.width = target_width, .height = target_height};

  createColorResources();
  createDepthResources();
  createFramebuffer();
  createStagingBuffers();
}

void VulkanOffscreenTarget::destroy() noexcept {
  if (device_ == VK_NULL_HANDLE) {
    return;
  }

  destroyFrameResources();

  for (auto &staging_buffer : staging_buffers_) {
    staging_buffer.destroy();
  }

  if (render_pass_ != VK_NULL_HANDLE) {
    vkDestroyRenderPass(device_, render_pass_, nullptr);
    render_pass_ = VK_NULL_HANDLE;
  }

  device_ = VK_NULL_HANDLE;
  physical_device_ = VK_NULL_HANDLE;
  extent_ = {.width = 0, .height = 0};
}

void VulkanOffscreenTarget::destroyFrameResources() noexcept {
  if (device_ == VK_NULL_HANDLE) {
    return;
  }

  if (framebuffer_ != VK_NULL_HANDLE) {
    vkDestroyFramebuffer(device_, framebuffer_, nullptr);
    framebuffer_ = VK_NULL_HANDLE;
  }

  if (color_image_view_ != VK_NULL_HANDLE) {
    vkDestroyImageView(device_, color_image_view_, nullptr);
    color_image_view_ = VK_NULL_HANDLE;
  }
  if (color_image_ != VK_NULL_HANDLE) {
    vkDestroyImage(device_, color_image_, nullptr);
    color_image_ = VK_NULL_HANDLE;
  }
  if (color_memory_ != VK_NULL_HANDLE) {
    vkFreeMemory(device_, color_memory_, nullptr);
    color_memory_ = VK_NULL_HANDLE;
  }

  if (depth_image_view_ != VK_NULL_HANDLE) {
    vkDestroyImageView(device_, depth_image_view_, nullptr);
    depth_image_view_ = VK_NULL_HANDLE;
  }
  if (depth_image_ != VK_NULL_HANDLE) {
    vkDestroyImage(device_, depth_image_, nullptr);
    depth_image_ = VK_NULL_HANDLE;
  }
  if (depth_memory_ != VK_NULL_HANDLE) {
    vkFreeMemory(device_, depth_memory_, nullptr);
    depth_memory_ = VK_NULL_HANDLE;
  }
}

void VulkanOffscreenTarget::moveFrom(VulkanOffscreenTarget &other) noexcept {
  device_ = std::exchange(other.device_, VK_NULL_HANDLE);
  physical_device_ = std::exchange(other.physical_device_, VK_NULL_HANDLE);
  extent_ = std::exchange(other.extent_, {.width = 0, .height = 0});
  color_format_ = other.color_format_;
  depth_format_ = other.depth_format_;
  render_pass_ = std::exchange(other.render_pass_, VK_NULL_HANDLE);
  color_image_ = std::exchange(other.color_image_, VK_NULL_HANDLE);
  color_memory_ = std::exchange(other.color_memory_, VK_NULL_HANDLE);
  color_image_view_ = std::exchange(other.color_image_view_, VK_NULL_HANDLE);
  depth_image_ = std::exchange(other.depth_image_, VK_NULL_HANDLE);
  depth_memory_ = std::exchange(other.depth_memory_, VK_NULL_HANDLE);
  depth_image_view_ = std::exchange(other.depth_image_view_, VK_NULL_HANDLE);
  framebuffer_ = std::exchange(other.framebuffer_, VK_NULL_HANDLE);
  staging_buffers_ = std::move(other.staging_buffers_);
}

void VulkanOffscreenTarget::beginRenderPass(VkCommandBuffer cmd, VkClearColorValue clear_color) {
  std::array<VkClearValue, 2> clear_values{};
  clear_values.at(0).color = clear_color;
  clear_values.at(1).depthStencil = {.depth = 1.0F, .stencil = 0};

  const VkRenderPassBeginInfo begin_info{
      .sType = VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO,
      .pNext = nullptr,
      .renderPass = render_pass_,
      .framebuffer = framebuffer_,
      .renderArea = {.offset = {.x = 0, .y = 0}, .extent = extent_},
      .clearValueCount = static_cast<uint32_t>(clear_values.size()),
      .pClearValues = clear_values.data(),
  };

  vkCmdBeginRenderPass(cmd, &begin_info, VK_SUBPASS_CONTENTS_INLINE);

  const VkViewport viewport{
      .x = 0.0F,
      .y = 0.0F,
      .width = static_cast<float>(extent_.width),
      .height = static_cast<float>(extent_.height),
      .minDepth = 0.0F,
      .maxDepth = 1.0F,
  };
  vkCmdSetViewport(cmd, 0, 1, &viewport);

  const VkRect2D scissor{
      .offset = {.x = 0, .y = 0},
      .extent = extent_,
  };
  vkCmdSetScissor(cmd, 0, 1, &scissor);
}

void VulkanOffscreenTarget::endRenderPass(VkCommandBuffer cmd) { vkCmdEndRenderPass(cmd); }

void VulkanOffscreenTarget::copyToStaging(VkCommandBuffer cmd, uint32_t staging_index) {
  const size_t idx = staging_index % STAGING_BUFFER_COUNT;
  const VkBufferImageCopy region{
      .bufferOffset = 0,
      .bufferRowLength = 0,
      .bufferImageHeight = 0,
      .imageSubresource =
          {
              .aspectMask = VK_IMAGE_ASPECT_COLOR_BIT,
              .mipLevel = 0,
              .baseArrayLayer = 0,
              .layerCount = 1,
          },
      .imageOffset = {.x = 0, .y = 0, .z = 0},
      .imageExtent = {.width = extent_.width, .height = extent_.height, .depth = 1},
  };

  vkCmdCopyImageToBuffer(cmd, color_image_, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
                         staging_buffers_.at(idx).buffer(), 1, &region);

  const VkBufferMemoryBarrier barrier{
      .sType = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER,
      .pNext = nullptr,
      .srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT,
      .dstAccessMask = VK_ACCESS_HOST_READ_BIT,
      .srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED,
      .dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED,
      .buffer = staging_buffers_.at(idx).buffer(),
      .offset = 0,
      .size = VK_WHOLE_SIZE,
  };

  vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_HOST_BIT, 0, 0,
                       nullptr, 1, &barrier, 0, nullptr);
}

std::span<const uint8_t>
VulkanOffscreenTarget::readoutStagingSpan(uint32_t staging_index) const noexcept {
  const size_t idx = staging_index % STAGING_BUFFER_COUNT;
  if (!staging_buffers_.at(idx).isValid() || staging_buffers_.at(idx).mappedData() == nullptr) {
    return {};
  }
  const size_t byte_size = static_cast<size_t>(extent_.width) * extent_.height * 4;
  return {static_cast<const uint8_t *>(staging_buffers_.at(idx).mappedData()), byte_size};
}

VkFormat VulkanOffscreenTarget::findDepthFormat() const {
  constexpr std::array<VkFormat, 3> CANDIDATES{
      VK_FORMAT_D32_SFLOAT,
      VK_FORMAT_D32_SFLOAT_S8_UINT,
      VK_FORMAT_D24_UNORM_S8_UINT,
  };

  for (const auto format : CANDIDATES) {
    VkFormatProperties props{};
    vkGetPhysicalDeviceFormatProperties(physical_device_, format, &props);
    if ((props.optimalTilingFeatures & VK_FORMAT_FEATURE_DEPTH_STENCIL_ATTACHMENT_BIT) != 0) {
      return format;
    }
  }

  throw std::runtime_error("Failed to find supported Vulkan depth attachment format");
}

void VulkanOffscreenTarget::createRenderPass() {
  const VkAttachmentDescription color_attachment{
      .flags = 0,
      .format = color_format_,
      .samples = VK_SAMPLE_COUNT_1_BIT,
      .loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR,
      .storeOp = VK_ATTACHMENT_STORE_OP_STORE,
      .stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE,
      .stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE,
      .initialLayout = VK_IMAGE_LAYOUT_UNDEFINED,
      .finalLayout = VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
  };

  const VkAttachmentDescription depth_attachment{
      .flags = 0,
      .format = depth_format_,
      .samples = VK_SAMPLE_COUNT_1_BIT,
      .loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR,
      .storeOp = VK_ATTACHMENT_STORE_OP_DONT_CARE,
      .stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE,
      .stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE,
      .initialLayout = VK_IMAGE_LAYOUT_UNDEFINED,
      .finalLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL,
  };

  const VkAttachmentReference color_ref{
      .attachment = 0,
      .layout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL,
  };

  const VkAttachmentReference depth_ref{
      .attachment = 1,
      .layout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL,
  };

  const VkSubpassDescription subpass{
      .flags = 0,
      .pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS,
      .inputAttachmentCount = 0,
      .pInputAttachments = nullptr,
      .colorAttachmentCount = 1,
      .pColorAttachments = &color_ref,
      .pResolveAttachments = nullptr,
      .pDepthStencilAttachment = &depth_ref,
      .preserveAttachmentCount = 0,
      .pPreserveAttachments = nullptr,
  };

  const std::array<VkSubpassDependency, 2> dependencies{
      VkSubpassDependency{
          .srcSubpass = VK_SUBPASS_EXTERNAL,
          .dstSubpass = 0,
          .srcStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT |
                          VK_PIPELINE_STAGE_EARLY_FRAGMENT_TESTS_BIT,
          .dstStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT |
                          VK_PIPELINE_STAGE_EARLY_FRAGMENT_TESTS_BIT,
          .srcAccessMask = 0,
          .dstAccessMask =
              VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT | VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT,
          .dependencyFlags = 0,
      },
      VkSubpassDependency{
          .srcSubpass = 0,
          .dstSubpass = VK_SUBPASS_EXTERNAL,
          .srcStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT,
          .dstStageMask = VK_PIPELINE_STAGE_TRANSFER_BIT,
          .srcAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT,
          .dstAccessMask = VK_ACCESS_TRANSFER_READ_BIT,
          .dependencyFlags = 0,
      },
  };

  const std::array<VkAttachmentDescription, 2> attachments{color_attachment, depth_attachment};

  const VkRenderPassCreateInfo render_pass_info{
      .sType = VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .attachmentCount = static_cast<uint32_t>(attachments.size()),
      .pAttachments = attachments.data(),
      .subpassCount = 1,
      .pSubpasses = &subpass,
      .dependencyCount = static_cast<uint32_t>(dependencies.size()),
      .pDependencies = dependencies.data(),
  };

  checkVkResult(vkCreateRenderPass(device_, &render_pass_info, nullptr, &render_pass_),
                "Failed to create offscreen render pass");
}

void VulkanOffscreenTarget::createColorResources() {
  createImage(extent_.width, extent_.height, color_format_, VK_IMAGE_TILING_OPTIMAL,
              VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT,
              VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, color_image_, color_memory_);
  color_image_view_ = createImageView(color_image_, color_format_, VK_IMAGE_ASPECT_COLOR_BIT);
}

void VulkanOffscreenTarget::createDepthResources() {
  createImage(extent_.width, extent_.height, depth_format_, VK_IMAGE_TILING_OPTIMAL,
              VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
              depth_image_, depth_memory_);
  depth_image_view_ = createImageView(depth_image_, depth_format_, VK_IMAGE_ASPECT_DEPTH_BIT);
}

void VulkanOffscreenTarget::createFramebuffer() {
  const std::array<VkImageView, 2> attachments{color_image_view_, depth_image_view_};

  const VkFramebufferCreateInfo framebuffer_info{
      .sType = VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .renderPass = render_pass_,
      .attachmentCount = static_cast<uint32_t>(attachments.size()),
      .pAttachments = attachments.data(),
      .width = extent_.width,
      .height = extent_.height,
      .layers = 1,
  };

  checkVkResult(vkCreateFramebuffer(device_, &framebuffer_info, nullptr, &framebuffer_),
                "Failed to create offscreen framebuffer");
}

void VulkanOffscreenTarget::createStagingBuffers() {
  const VkDeviceSize required_bytes =
      std::max<VkDeviceSize>(static_cast<VkDeviceSize>(extent_.width) * extent_.height * 4, 1024);
  const BufferAllocationDesc desc{
      .size = required_bytes,
      .usage = VK_BUFFER_USAGE_TRANSFER_DST_BIT,
      .properties = VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
  };

  for (auto &staging_buffer : staging_buffers_) {
    staging_buffer.ensureCapacity(device_, physical_device_, desc);
    if (staging_buffer.mappedData() == nullptr) {
      static_cast<void>(staging_buffer.map());
    }
  }
}

void VulkanOffscreenTarget::createImage(uint32_t width, uint32_t height, VkFormat format,
                                        VkImageTiling tiling, VkImageUsageFlags usage,
                                        VkMemoryPropertyFlags properties, VkImage &image,
                                        VkDeviceMemory &memory) {
  const VkImageCreateInfo image_info{
      .sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .imageType = VK_IMAGE_TYPE_2D,
      .format = format,
      .extent = {.width = width, .height = height, .depth = 1},
      .mipLevels = 1,
      .arrayLayers = 1,
      .samples = VK_SAMPLE_COUNT_1_BIT,
      .tiling = tiling,
      .usage = usage,
      .sharingMode = VK_SHARING_MODE_EXCLUSIVE,
      .queueFamilyIndexCount = 0,
      .pQueueFamilyIndices = nullptr,
      .initialLayout = VK_IMAGE_LAYOUT_UNDEFINED,
  };

  checkVkResult(vkCreateImage(device_, &image_info, nullptr, &image),
                "Failed to create Vulkan image");

  VkMemoryRequirements mem_requirements{};
  vkGetImageMemoryRequirements(device_, image, &mem_requirements);

  const VkMemoryAllocateInfo alloc_info{
      .sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO,
      .pNext = nullptr,
      .allocationSize = mem_requirements.size,
      .memoryTypeIndex =
          findMemoryType(physical_device_, mem_requirements.memoryTypeBits, properties),
  };

  checkVkResult(vkAllocateMemory(device_, &alloc_info, nullptr, &memory),
                "Failed to allocate Vulkan image memory");
  checkVkResult(vkBindImageMemory(device_, image, memory, 0), "Failed to bind Vulkan image memory");
}

VkImageView VulkanOffscreenTarget::createImageView(VkImage image, VkFormat format,
                                                   VkImageAspectFlags aspect_flags) const {
  const VkImageViewCreateInfo view_info{
      .sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .image = image,
      .viewType = VK_IMAGE_VIEW_TYPE_2D,
      .format = format,
      .components =
          {
              .r = VK_COMPONENT_SWIZZLE_IDENTITY,
              .g = VK_COMPONENT_SWIZZLE_IDENTITY,
              .b = VK_COMPONENT_SWIZZLE_IDENTITY,
              .a = VK_COMPONENT_SWIZZLE_IDENTITY,
          },
      .subresourceRange =
          {
              .aspectMask = aspect_flags,
              .baseMipLevel = 0,
              .levelCount = 1,
              .baseArrayLayer = 0,
              .layerCount = 1,
          },
  };

  VkImageView image_view = VK_NULL_HANDLE;
  checkVkResult(vkCreateImageView(device_, &view_info, nullptr, &image_view),
                "Failed to create Vulkan image view");
  return image_view;
}

} // namespace correlation::render3d::vulkan

#endif // CORRELATION_USE_VULKAN
