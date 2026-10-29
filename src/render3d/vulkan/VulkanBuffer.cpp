/**
 * @file VulkanBuffer.cpp
 * @brief Implementation of VulkanBuffer RAII lifecycle and memory operations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/vulkan/VulkanBuffer.hpp"
#include "render3d/vulkan/VulkanCommon.hpp"
#include <algorithm>
#include <cstring>

namespace correlation::render3d::vulkan {

#ifdef CORRELATION_USE_VULKAN

void VulkanBuffer::create(VkDevice device, VkPhysicalDevice physical_device,
                          const BufferAllocationDesc &desc) {
  destroy();

  device_ = device;
  size_ = desc.size;

  const VkBufferCreateInfo buffer_info{
      .sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .size = size_,
      .usage = desc.usage,
      .sharingMode = VK_SHARING_MODE_EXCLUSIVE,
      .queueFamilyIndexCount = 0,
      .pQueueFamilyIndices = nullptr,
  };

  checkVkResult(vkCreateBuffer(device_, &buffer_info, nullptr, &buffer_),
                "Failed to create VkBuffer");

  VkMemoryRequirements mem_requirements{};
  vkGetBufferMemoryRequirements(device_, buffer_, &mem_requirements);

  const uint32_t mem_type_index =
      findMemoryType(physical_device, mem_requirements.memoryTypeBits, desc.properties);

  const VkMemoryAllocateInfo alloc_info{
      .sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO,
      .pNext = nullptr,
      .allocationSize = mem_requirements.size,
      .memoryTypeIndex = mem_type_index,
  };

  checkVkResult(vkAllocateMemory(device_, &alloc_info, nullptr, &memory_),
                "Failed to allocate VkDeviceMemory for buffer");
  checkVkResult(vkBindBufferMemory(device_, buffer_, memory_, 0), "Failed to bind VkBuffer memory");
}

void VulkanBuffer::destroy() noexcept {
  if (device_ == VK_NULL_HANDLE) {
    return;
  }

  if (mapped_data_ != nullptr) {
    vkUnmapMemory(device_, memory_);
    mapped_data_ = nullptr;
  }

  if (buffer_ != VK_NULL_HANDLE) {
    vkDestroyBuffer(device_, buffer_, nullptr);
    buffer_ = VK_NULL_HANDLE;
  }

  if (memory_ != VK_NULL_HANDLE) {
    vkFreeMemory(device_, memory_, nullptr);
    memory_ = VK_NULL_HANDLE;
  }

  size_ = 0;
  device_ = VK_NULL_HANDLE;
}

void *VulkanBuffer::map() {
  if (mapped_data_ == nullptr && memory_ != VK_NULL_HANDLE) {
    checkVkResult(vkMapMemory(device_, memory_, 0, size_, 0, &mapped_data_),
                  "Failed to map VulkanBuffer memory");
  }
  return mapped_data_;
}

void VulkanBuffer::unmap() noexcept {
  if (mapped_data_ != nullptr && device_ != VK_NULL_HANDLE && memory_ != VK_NULL_HANDLE) {
    vkUnmapMemory(device_, memory_);
    mapped_data_ = nullptr;
  }
}

void VulkanBuffer::upload(std::span<const std::byte> data) {
  if (data.empty()) {
    return;
  }
  if (data.size() > static_cast<size_t>(size_)) {
    throw std::runtime_error("Upload data size exceeds VulkanBuffer capacity");
  }

  void *dst = map();
  std::memcpy(dst, data.data(), data.size());
}

bool VulkanBuffer::ensureCapacity(VkDevice device, VkPhysicalDevice physical_device,
                                  const BufferAllocationDesc &desc) {
  if (desc.size <= size_) {
    return false;
  }

  VkDeviceSize new_size = std::max(size_ * 2, static_cast<VkDeviceSize>(1024));
  while (new_size < desc.size) {
    new_size *= 2;
  }

  const BufferAllocationDesc expanded_desc{
      .size = new_size,
      .usage = desc.usage,
      .properties = desc.properties,
  };

  create(device, physical_device, expanded_desc);
  return true;
}

#endif // CORRELATION_USE_VULKAN

} // namespace correlation::render3d::vulkan
