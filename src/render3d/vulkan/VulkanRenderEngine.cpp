/**
 * @file VulkanRenderEngine.cpp
 * @brief Implementation of VulkanRenderEngine headless pipeline.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/vulkan/VulkanRenderEngine.hpp"
#include "render3d/vulkan/VulkanBufferManager.hpp"
#include "render3d/vulkan/VulkanContext.hpp"
#include "render3d/vulkan/VulkanDevice.hpp"
#include "render3d/vulkan/VulkanOffscreenTarget.hpp"
#include "render3d/vulkan/VulkanPipelineManager.hpp"
#include <array>
#include <cstdint>

namespace correlation::render3d {

VulkanRenderEngine::VulkanRenderEngine() = default;

VulkanRenderEngine::~VulkanRenderEngine() { cleanup(); }

void VulkanRenderEngine::cleanup() noexcept {
#ifdef CORRELATION_USE_VULKAN
  if (device_ && device_->isValid()) {
    vkDeviceWaitIdle(device_->device());
    if (render_fence_ != VK_NULL_HANDLE) {
      vkDestroyFence(device_->device(), render_fence_, nullptr);
      render_fence_ = VK_NULL_HANDLE;
    }
    render_command_buffer_ = VK_NULL_HANDLE;
    buffer_manager_.reset();
    pipeline_manager_.reset();
    target_.reset();
    device_.reset();
  }
#endif
  initialized_ = false;
}

bool VulkanRenderEngine::initialize(const ViewportDimensions &dimensions) {
  cleanup();

  if (!VulkanContext::isSupported()) {
    return false;
  }

#ifdef CORRELATION_USE_VULKAN
  device_ = std::make_unique<vulkan::VulkanDevice>();
  if (!device_->initialize()) {
    cleanup();
    return false;
  }

  target_ = std::make_unique<vulkan::VulkanOffscreenTarget>();
  target_->create(device_->device(), device_->physicalDevice(), dimensions.width,
                  dimensions.height);
  if (!target_->isValid()) {
    cleanup();
    return false;
  }

  pipeline_manager_ = std::make_unique<vulkan::VulkanPipelineManager>();
  if (!pipeline_manager_->initialize(device_->device(), target_->renderPass())) {
    cleanup();
    return false;
  }

  buffer_manager_ = std::make_unique<vulkan::VulkanBufferManager>();
  if (!buffer_manager_->initialize(device_->device(), device_->physicalDevice(),
                                   pipeline_manager_->descriptorSetLayout())) {
    cleanup();
    return false;
  }

  const VkCommandBufferAllocateInfo alloc_info{
      .sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO,
      .pNext = nullptr,
      .commandPool = device_->commandPool(),
      .level = VK_COMMAND_BUFFER_LEVEL_PRIMARY,
      .commandBufferCount = 1,
  };

  if (vkAllocateCommandBuffers(device_->device(), &alloc_info, &render_command_buffer_) !=
      VK_SUCCESS) {
    cleanup();
    return false;
  }

  const VkFenceCreateInfo fence_info{
      .sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO,
      .pNext = nullptr,
      .flags = VK_FENCE_CREATE_SIGNALED_BIT,
  };

  if (vkCreateFence(device_->device(), &fence_info, nullptr, &render_fence_) != VK_SUCCESS) {
    cleanup();
    return false;
  }

  dimensions_ = dimensions;
  active_staging_index_ = 0;
  initialized_ = true;
  return true;
#else
  (void)dimensions;
  return false;
#endif
}

void VulkanRenderEngine::resize(const ViewportDimensions &dimensions) {
  if (dimensions.width == dimensions_.width && dimensions.height == dimensions_.height) {
    return;
  }

  dimensions_ = dimensions;

#ifdef CORRELATION_USE_VULKAN
  if (initialized_ && device_ && target_) {
    vkDeviceWaitIdle(device_->device());
    target_->resize(dimensions.width, dimensions.height);
  }
#endif

  const std::size_t pixel_bytes =
      static_cast<std::size_t>(dimensions_.width) * dimensions_.height * 4;
  fallback_pixels_.resize(pixel_bytes, 0);
}

void VulkanRenderEngine::setScene(std::shared_ptr<const Scene3D> scene) {
  scene_ = std::move(scene);
}

void VulkanRenderEngine::setCamera(const Camera3D &camera) { camera_ = camera; }

void VulkanRenderEngine::setFeatureVisibility(bool render_bonds, bool render_wireframe) {
  render_bonds_ = render_bonds;
  render_wireframe_ = render_wireframe;
}

void VulkanRenderEngine::bindViewportAndDescriptors() {
#ifdef CORRELATION_USE_VULKAN
  const VkViewport viewport{
      .x = 0.0F,
      .y = 0.0F,
      .width = static_cast<float>(dimensions_.width),
      .height = static_cast<float>(dimensions_.height),
      .minDepth = 0.0F,
      .maxDepth = 1.0F,
  };
  vkCmdSetViewport(render_command_buffer_, 0, 1, &viewport);

  const VkRect2D scissor{.offset = {.x = 0, .y = 0}, .extent = target_->extent()};
  vkCmdSetScissor(render_command_buffer_, 0, 1, &scissor);

  VkDescriptorSet desc_set = buffer_manager_->descriptorSet();
  vkCmdBindDescriptorSets(render_command_buffer_, VK_PIPELINE_BIND_POINT_GRAPHICS,
                          pipeline_manager_->pipelineLayout(), 0, 1, &desc_set, 0, nullptr);
#endif
}

void VulkanRenderEngine::drawWireframe() {
#ifdef CORRELATION_USE_VULKAN
  if (render_wireframe_ && buffer_manager_->wireframeVertexCount() > 0) {
    vkCmdBindPipeline(render_command_buffer_, VK_PIPELINE_BIND_POINT_GRAPHICS,
                      pipeline_manager_->wireframePipeline());
    VkBuffer line_buf = buffer_manager_->wireframeBuffer();
    const VkDeviceSize offset = 0;
    vkCmdBindVertexBuffers(render_command_buffer_, 0, 1, &line_buf, &offset);
    vkCmdDraw(render_command_buffer_, buffer_manager_->wireframeVertexCount(), 1, 0, 0);
  }
#endif
}

void VulkanRenderEngine::drawBonds() {
#ifdef CORRELATION_USE_VULKAN
  if (render_bonds_ && buffer_manager_->bondInstanceCount() > 0) {
    vkCmdBindPipeline(render_command_buffer_, VK_PIPELINE_BIND_POINT_GRAPHICS,
                      pipeline_manager_->cylinderPipeline());
    const std::array<VkBuffer, 2> vertex_buffers{
        buffer_manager_->cylinderQuadBuffer(),
        buffer_manager_->bondInstanceBuffer(),
    };
    const std::array<VkDeviceSize, 2> offsets{0, 0};
    vkCmdBindVertexBuffers(render_command_buffer_, 0, 2, vertex_buffers.data(), offsets.data());
    vkCmdDraw(render_command_buffer_, 4, buffer_manager_->bondInstanceCount(), 0, 0);
  }
#endif
}

void VulkanRenderEngine::drawAtoms() {
#ifdef CORRELATION_USE_VULKAN
  if (buffer_manager_->atomInstanceCount() > 0) {
    vkCmdBindPipeline(render_command_buffer_, VK_PIPELINE_BIND_POINT_GRAPHICS,
                      pipeline_manager_->spherePipeline());
    const std::array<VkBuffer, 2> vertex_buffers{
        buffer_manager_->sphereQuadBuffer(),
        buffer_manager_->atomInstanceBuffer(),
    };
    const std::array<VkDeviceSize, 2> offsets{0, 0};
    vkCmdBindVertexBuffers(render_command_buffer_, 0, 2, vertex_buffers.data(), offsets.data());
    vkCmdDraw(render_command_buffer_, 4, buffer_manager_->atomInstanceCount(), 0, 0);
  }
#endif
}

void VulkanRenderEngine::recordRenderCommands() {
#ifdef CORRELATION_USE_VULKAN
  vkResetCommandBuffer(render_command_buffer_, 0);

  const VkCommandBufferBeginInfo begin_info{
      .sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO,
      .pNext = nullptr,
      .flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT,
      .pInheritanceInfo = nullptr,
  };
  vkBeginCommandBuffer(render_command_buffer_, &begin_info);

  target_->beginRenderPass(render_command_buffer_, {{0.1F, 0.1F, 0.12F, 1.0F}});

  bindViewportAndDescriptors();
  drawWireframe();
  drawBonds();
  drawAtoms();

  vulkan::VulkanOffscreenTarget::endRenderPass(render_command_buffer_);
  target_->copyToStaging(render_command_buffer_, active_staging_index_);

  vkEndCommandBuffer(render_command_buffer_);
#endif
}

std::span<const uint8_t> VulkanRenderEngine::renderFrame() {
  if (!initialized_) {
    const std::size_t pixel_bytes =
        static_cast<std::size_t>(dimensions_.width) * dimensions_.height * 4;
    fallback_pixels_.resize(pixel_bytes, 0);
    return fallback_pixels_;
  }

#ifdef CORRELATION_USE_VULKAN
  vkWaitForFences(device_->device(), 1, &render_fence_, VK_TRUE, UINT64_MAX);
  vkResetFences(device_->device(), 1, &render_fence_);

  if (scene_) {
    buffer_manager_->uploadScene(*scene_, camera_, dimensions_);
  } else {
    const Scene3D empty_scene;
    buffer_manager_->uploadScene(empty_scene, camera_, dimensions_);
  }

  recordRenderCommands();

  const VkSubmitInfo submit_info{
      .sType = VK_STRUCTURE_TYPE_SUBMIT_INFO,
      .pNext = nullptr,
      .waitSemaphoreCount = 0,
      .pWaitSemaphores = nullptr,
      .pWaitDstStageMask = nullptr,
      .commandBufferCount = 1,
      .pCommandBuffers = &render_command_buffer_,
      .signalSemaphoreCount = 0,
      .pSignalSemaphores = nullptr,
  };

  vkQueueSubmit(device_->graphicsQueue(), 1, &submit_info, render_fence_);
  vkWaitForFences(device_->device(), 1, &render_fence_, VK_TRUE, UINT64_MAX);

  const auto pixels = target_->readoutStagingSpan(active_staging_index_);
  active_staging_index_ =
      (active_staging_index_ + 1) % vulkan::VulkanOffscreenTarget::STAGING_BUFFER_COUNT;
  return pixels;
#else
  return fallback_pixels_;
#endif
}

ViewportDimensions VulkanRenderEngine::dimensions() const noexcept { return dimensions_; }

std::string_view VulkanRenderEngine::backendName() const noexcept {
  return "Vulkan (Headless Offscreen)";
}

bool VulkanRenderEngine::isHardwareAccelerated() const noexcept { return true; }

} // namespace correlation::render3d
