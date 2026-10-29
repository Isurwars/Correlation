/**
 * @file VulkanDevice.cpp
 * @brief Implementation of VulkanDevice physical/logical GPU initialization.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/vulkan/VulkanDevice.hpp"
#include "render3d/vulkan/VulkanCommon.hpp"
#include <span>
#include <vector>

namespace correlation::render3d::vulkan {

#ifdef CORRELATION_USE_VULKAN

namespace {

int scorePhysicalDevice(VkPhysicalDevice device, uint32_t &out_queue_family) {
  VkPhysicalDeviceProperties properties{};
  vkGetPhysicalDeviceProperties(device, &properties);

  uint32_t queue_family_count = 0;
  vkGetPhysicalDeviceQueueFamilyProperties(device, &queue_family_count, nullptr);
  std::vector<VkQueueFamilyProperties> queue_families(queue_family_count);
  vkGetPhysicalDeviceQueueFamilyProperties(device, &queue_family_count, queue_families.data());

  int graphics_family = -1;
  for (uint32_t idx = 0; idx < queue_family_count; ++idx) {
    if ((queue_families[idx].queueFlags & VK_QUEUE_GRAPHICS_BIT) != 0) {
      graphics_family = static_cast<int>(idx);
      break;
    }
  }

  if (graphics_family < 0) {
    return -1;
  }
  out_queue_family = static_cast<uint32_t>(graphics_family);

  int score = 10;
  if (properties.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU) {
    score += 1000;
  } else if (properties.deviceType == VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU) {
    score += 100;
  }
  return score;
}

} // namespace

void VulkanDevice::moveFrom(VulkanDevice &other) noexcept {
  instance_ = other.instance_;
  physical_device_ = other.physical_device_;
  device_ = other.device_;
  graphics_queue_ = other.graphics_queue_;
  queue_family_index_ = other.queue_family_index_;
  command_pool_ = other.command_pool_;

  other.instance_ = VK_NULL_HANDLE;
  other.physical_device_ = VK_NULL_HANDLE;
  other.device_ = VK_NULL_HANDLE;
  other.graphics_queue_ = VK_NULL_HANDLE;
  other.queue_family_index_ = 0;
  other.command_pool_ = VK_NULL_HANDLE;
}

bool VulkanDevice::initialize() {
  destroy();

  if (!createInstance()) {
    return false;
  }
  if (!pickPhysicalDevice()) {
    destroy();
    return false;
  }
  if (!createLogicalDevice()) {
    destroy();
    return false;
  }
  if (!createCommandPool()) {
    destroy();
    return false;
  }

  return true;
}

void VulkanDevice::destroy() noexcept {
  if (device_ != VK_NULL_HANDLE) {
    vkDeviceWaitIdle(device_);
    if (command_pool_ != VK_NULL_HANDLE) {
      vkDestroyCommandPool(device_, command_pool_, nullptr);
      command_pool_ = VK_NULL_HANDLE;
    }
    vkDestroyDevice(device_, nullptr);
    device_ = VK_NULL_HANDLE;
  }

  if (instance_ != VK_NULL_HANDLE) {
    vkDestroyInstance(instance_, nullptr);
    instance_ = VK_NULL_HANDLE;
  }

  physical_device_ = VK_NULL_HANDLE;
  graphics_queue_ = VK_NULL_HANDLE;
  queue_family_index_ = 0;
}

bool VulkanDevice::createInstance() {
  const VkApplicationInfo app_info{
      .sType = VK_STRUCTURE_TYPE_APPLICATION_INFO,
      .pNext = nullptr,
      .pApplicationName = "Correlation",
      .applicationVersion = VK_MAKE_VERSION(4, 1, 2),
      .pEngineName = "CorrelationRender3D",
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

  return vkCreateInstance(&create_info, nullptr, &instance_) == VK_SUCCESS;
}

bool VulkanDevice::pickPhysicalDevice() {
  uint32_t device_count = 0;
  if (vkEnumeratePhysicalDevices(instance_, &device_count, nullptr) != VK_SUCCESS ||
      device_count == 0) {
    return false;
  }

  std::vector<VkPhysicalDevice> devices(device_count);
  if (vkEnumeratePhysicalDevices(instance_, &device_count, devices.data()) != VK_SUCCESS) {
    return false;
  }

  int best_score = -1;
  VkPhysicalDevice best_device = VK_NULL_HANDLE;
  uint32_t best_queue_family = 0;

  for (VkPhysicalDevice dev : devices) {
    uint32_t queue_family = 0;
    const int score = scorePhysicalDevice(dev, queue_family);
    if (score > best_score) {
      best_score = score;
      best_device = dev;
      best_queue_family = queue_family;
    }
  }

  if (best_device == VK_NULL_HANDLE) {
    return false;
  }

  physical_device_ = best_device;
  queue_family_index_ = best_queue_family;
  return true;
}

bool VulkanDevice::createLogicalDevice() {
  const float queue_priority = 1.0F;
  const VkDeviceQueueCreateInfo queue_info{
      .sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .queueFamilyIndex = queue_family_index_,
      .queueCount = 1,
      .pQueuePriorities = &queue_priority,
  };

  VkPhysicalDeviceFeatures device_features{};
  device_features.fillModeNonSolid = VK_TRUE;

  const VkDeviceCreateInfo create_info{
      .sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .queueCreateInfoCount = 1,
      .pQueueCreateInfos = &queue_info,
      .enabledLayerCount = 0,
      .ppEnabledLayerNames = nullptr,
      .enabledExtensionCount = 0,
      .ppEnabledExtensionNames = nullptr,
      .pEnabledFeatures = &device_features,
  };

  if (vkCreateDevice(physical_device_, &create_info, nullptr, &device_) != VK_SUCCESS) {
    return false;
  }

  vkGetDeviceQueue(device_, queue_family_index_, 0, &graphics_queue_);
  return true;
}

bool VulkanDevice::createCommandPool() {
  const VkCommandPoolCreateInfo pool_info{
      .sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO,
      .pNext = nullptr,
      .flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT,
      .queueFamilyIndex = queue_family_index_,
  };

  return vkCreateCommandPool(device_, &pool_info, nullptr, &command_pool_) == VK_SUCCESS;
}

VkCommandBuffer VulkanDevice::beginSingleTimeCommands() {
  const VkCommandBufferAllocateInfo alloc_info{
      .sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO,
      .pNext = nullptr,
      .commandPool = command_pool_,
      .level = VK_COMMAND_BUFFER_LEVEL_PRIMARY,
      .commandBufferCount = 1,
  };

  VkCommandBuffer command_buffer = VK_NULL_HANDLE;
  checkVkResult(vkAllocateCommandBuffers(device_, &alloc_info, &command_buffer),
                "Failed to allocate transient command buffer");

  const VkCommandBufferBeginInfo begin_info{
      .sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO,
      .pNext = nullptr,
      .flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT,
      .pInheritanceInfo = nullptr,
  };

  checkVkResult(vkBeginCommandBuffer(command_buffer, &begin_info),
                "Failed to begin transient command buffer");
  return command_buffer;
}

void VulkanDevice::endSingleTimeCommands(VkCommandBuffer command_buffer) {
  checkVkResult(vkEndCommandBuffer(command_buffer), "Failed to end transient command buffer");

  const VkSubmitInfo submit_info{
      .sType = VK_STRUCTURE_TYPE_SUBMIT_INFO,
      .pNext = nullptr,
      .waitSemaphoreCount = 0,
      .pWaitSemaphores = nullptr,
      .pWaitDstStageMask = nullptr,
      .commandBufferCount = 1,
      .pCommandBuffers = &command_buffer,
      .signalSemaphoreCount = 0,
      .pSignalSemaphores = nullptr,
  };

  checkVkResult(vkQueueSubmit(graphics_queue_, 1, &submit_info, VK_NULL_HANDLE),
                "Failed to submit transient command buffer");
  checkVkResult(vkQueueWaitIdle(graphics_queue_),
                "Failed to wait for transient command queue idle");

  vkFreeCommandBuffers(device_, command_pool_, 1, &command_buffer);
}

std::string VulkanDevice::physicalDeviceName() const {
  if (physical_device_ == VK_NULL_HANDLE) {
    return "None";
  }
  VkPhysicalDeviceProperties properties{};
  vkGetPhysicalDeviceProperties(physical_device_, &properties);
  const std::span name_span(properties.deviceName);
  return {name_span.data()};
}

#endif // CORRELATION_USE_VULKAN

} // namespace correlation::render3d::vulkan
