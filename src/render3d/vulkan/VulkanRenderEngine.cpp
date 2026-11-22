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
#include <algorithm>
#include <array>
#include <cstdint>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

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
    compute_service_.reset();
    buffer_manager_.reset();
    pipeline_manager_.reset();
    target_.reset();
    device_.reset();
  }
#endif
  initialized_ = false;
  last_uploaded_scene_ = nullptr;
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

  compute_service_ = std::make_unique<vulkan::VulkanComputeService>();
  if (!compute_service_->initialize(device_->device(), device_->physicalDevice(),
                                    device_->graphicsQueue(), device_->queueFamilyIndex())) {
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

void VulkanRenderEngine::setGpuHaloConfig(const GpuHaloConfig &config) {
  gpu_halos_enabled_ = config.enabled;
  halo_settings_.descriptor = config.descriptor;
  halo_settings_.min_threshold = config.min_threshold;
  halo_settings_.max_threshold = config.max_threshold;
  halo_settings_.halo_radius_ratio = config.halo_radius_ratio;
  halo_settings_.halo_intensity = config.halo_intensity;
  halo_settings_.halo_color = config.halo_color;
  halo_settings_.dim_non_matching = config.dim_non_matching;
  halo_settings_.hide_non_matching = config.hide_non_matching;
}

void VulkanRenderEngine::setRenderStyle(const RenderStyleConfig &style) { render_style_ = style; }

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

void VulkanRenderEngine::drawHalos() {
#ifdef CORRELATION_USE_VULKAN
  if (buffer_manager_->atomInstanceCount() > 0 &&
      pipeline_manager_->haloPipeline() != VK_NULL_HANDLE) {
    vkCmdBindPipeline(render_command_buffer_, VK_PIPELINE_BIND_POINT_GRAPHICS,
                      pipeline_manager_->haloPipeline());
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

  constexpr float INV_255 = 1.0F / 255.0F;
  const VkClearColorValue clear_rgba{{
      static_cast<float>(render_style_.background_color.r) * INV_255,
      static_cast<float>(render_style_.background_color.g) * INV_255,
      static_cast<float>(render_style_.background_color.b) * INV_255,
      static_cast<float>(render_style_.background_color.a) * INV_255,
  }};
  target_->beginRenderPass(render_command_buffer_, clear_rgba);

  bindViewportAndDescriptors();
  drawWireframe();
  drawBonds();
  drawAtoms();
  drawHalos();

  vulkan::VulkanOffscreenTarget::endRenderPass(render_command_buffer_);
  target_->copyToStaging(render_command_buffer_, active_staging_index_, render_style_.ssao_enabled);

  vkEndCommandBuffer(render_command_buffer_);
#endif
}

namespace {

void applyDepthSsaoPass(std::span<uint8_t> rgba, std::span<const float> depth, uint32_t width,
                        uint32_t height) noexcept {
  if (rgba.size() < static_cast<std::size_t>(width) * height * 4 ||
      depth.size() < static_cast<std::size_t>(width) * height || width < 16 || height < 16) {
    return;
  }
  constexpr std::array<std::pair<int, int>, 8> K_OFFSETS{{
      {3, 0},
      {-3, 0},
      {0, 3},
      {0, -3},
      {5, 5},
      {-5, 5},
      {5, -5},
      {-5, -5},
  }};
  constexpr float K_BIAS = 0.00035F;
  constexpr float K_RANGE = 0.025F;
  constexpr float K_STRENGTH = 0.42F;

  const auto max_x = static_cast<int>(width) - 1;
  const auto max_y = static_cast<int>(height) - 1;

  tbb::parallel_for(
      tbb::blocked_range<uint32_t>(0, height), [&](const tbb::blocked_range<uint32_t> &rows) {
        for (uint32_t py = rows.begin(); py < rows.end(); ++py) {
          const std::size_t row_off = static_cast<std::size_t>(py) * width;
          for (uint32_t px = 0; px < width; ++px) {
            const float center_d = depth[row_off + px];
            if (center_d >= 0.9999F) {
              continue;
            }
            float occ_sum = 0.0F;
            for (const auto &[off_x, off_y] : K_OFFSETS) {
              const auto sample_x =
                  static_cast<uint32_t>(std::clamp(static_cast<int>(px) + off_x, 0, max_x));
              const auto sample_y =
                  static_cast<uint32_t>(std::clamp(static_cast<int>(py) + off_y, 0, max_y));
              const float delta =
                  center_d - depth[static_cast<std::size_t>(sample_y) * width + sample_x];
              if (delta > K_BIAS && delta < K_RANGE) {
                occ_sum += 1.0F - (delta / K_RANGE);
              }
            }
            const float factor = 1.0F - (occ_sum * 0.125F) * K_STRENGTH;
            const std::size_t p_idx = (row_off + px) * 4;
            rgba[p_idx + 0] = static_cast<uint8_t>(static_cast<float>(rgba[p_idx + 0]) * factor);
            rgba[p_idx + 1] = static_cast<uint8_t>(static_cast<float>(rgba[p_idx + 1]) * factor);
            rgba[p_idx + 2] = static_cast<uint8_t>(static_cast<float>(rgba[p_idx + 2]) * factor);
          }
        }
      });
}

} // namespace

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
    if (scene_.get() != last_uploaded_scene_) {
      buffer_manager_->uploadScene(*scene_, camera_, dimensions_, render_style_);
      last_uploaded_scene_ = scene_.get();
    } else {
      buffer_manager_->uploadUbo(camera_, dimensions_, render_style_);
    }
  } else {
    const Scene3D empty_scene;
    buffer_manager_->uploadScene(empty_scene, camera_, dimensions_, render_style_);
    last_uploaded_scene_ = nullptr;
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

  if (render_style_.ssao_enabled && target_->hasDepthReadout()) {
    const auto mut_pixels = target_->readoutStagingMutableSpan(active_staging_index_);
    const auto depth_span = target_->readoutDepthStagingSpan(active_staging_index_);
    applyDepthSsaoPass(mut_pixels, depth_span, dimensions_.width, dimensions_.height);
  }
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

std::string VulkanRenderEngine::deviceName() const {
#ifdef CORRELATION_USE_VULKAN
  if (device_ && device_->isValid()) {
    return device_->physicalDeviceName();
  }
#endif
  return "Vulkan (Uninitialized)";
}

bool VulkanRenderEngine::isHardwareAccelerated() const noexcept { return true; }

} // namespace correlation::render3d
