/**
 * @file VulkanBufferManager.hpp
 * @brief Manages GPU buffers, UBO descriptors, and instance geometry uploads for Vulkan rendering.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#ifdef CORRELATION_USE_VULKAN
#include <vulkan/vulkan.h>
#endif

#include "render3d/Camera3D.hpp"
#include "render3d/IRenderEngine.hpp"
#include "render3d/Scene3D.hpp"
#include "render3d/vulkan/VulkanBuffer.hpp"
#include <array>
#include <cstdint>
#include <vector>

namespace correlation::render3d::vulkan {

#ifdef CORRELATION_USE_VULKAN

/**
 * @struct SceneUboData
 * @brief Memory layout for the uniform buffer binding 0 in scene shaders.
 */
struct alignas(16) SceneUboData {
  std::array<float, 16> view{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  std::array<float, 16> proj{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  std::array<float, 4> camera_pos{0.0F, 0.0F, 0.0F, 1.0F};
  std::array<float, 4> light_dir{0.0F, 0.5F, 1.0F, 0.0F};
  std::array<float, 2> viewport_size{800.0F, 600.0F};
  float fov{45.0F};
  float padding{0.0F};
};

/**
 * @struct AtomInstanceData
 * @brief Per-instance attributes for billboard sphere impostors.
 */
struct AtomInstanceData {
  std::array<float, 3> center{0.0F, 0.0F, 0.0F};
  float radius{0.5F};
  std::array<float, 4> color{1.0F, 1.0F, 1.0F, 1.0F};
  std::array<float, 4> halo_color{0.0F, 0.0F, 0.0F, 0.0F};
  std::array<float, 2> style{1.35F, 1.0F}; // x = halo_radius_ratio, y = dim_factor
};

/**
 * @struct BondInstanceData
 * @brief Per-instance attributes for cylinder bond impostors.
 */
struct BondInstanceData {
  std::array<float, 3> start{0.0F, 0.0F, 0.0F};
  std::array<float, 3> end{0.0F, 0.0F, 0.0F};
  float radius{0.12F};
  std::array<float, 4> color{1.0F, 1.0F, 1.0F, 1.0F};
};

/**
 * @struct WireframeVertexData
 * @brief Vertex attributes for unit cell and bounding box wireframe lines.
 */
struct WireframeVertexData {
  std::array<float, 3> position{0.0F, 0.0F, 0.0F};
  std::array<float, 4> color{1.0F, 1.0F, 1.0F, 1.0F};
};

/**
 * @class VulkanBufferManager
 * @brief Coordinates GPU memory allocation, dynamic growth, and scene geometric data marshaling.
 */
class VulkanBufferManager {
public:
  VulkanBufferManager() noexcept = default;
  ~VulkanBufferManager() noexcept { destroy(); }

  VulkanBufferManager(const VulkanBufferManager &) = delete;
  VulkanBufferManager &operator=(const VulkanBufferManager &) = delete;

  VulkanBufferManager(VulkanBufferManager &&other) noexcept { moveFrom(other); }
  VulkanBufferManager &operator=(VulkanBufferManager &&other) noexcept {
    if (this != &other) {
      destroy();
      moveFrom(other);
    }
    return *this;
  }

  /**
   * @brief Allocates scene UBO, descriptor pool, descriptor set, and static quad vertex buffers.
   * @param device Vulkan logical device handle.
   * @param physical_device Vulkan physical device handle.
   * @param ubo_layout Descriptor set layout for the scene uniform buffer.
   * @return True if all resources were allocated successfully.
   */
  [[nodiscard]] bool initialize(VkDevice device, VkPhysicalDevice physical_device,
                                VkDescriptorSetLayout ubo_layout);

  /**
   * @brief Destroys all allocated GPU buffers and descriptor pools.
   */
  void destroy() noexcept;

  /**
   * @brief Marshals and uploads scene geometry, camera matrices, and lighting uniforms.
   * @param scene Molecular 3D scene representation.
   * @param camera View and projection camera state.
   * @param dimensions Active viewport dimensions.
   */
  void uploadScene(const Scene3D &scene, const Camera3D &camera,
                   const ViewportDimensions &dimensions);

  [[nodiscard]] VkDescriptorSet descriptorSet() const noexcept { return descriptor_set_; }
  [[nodiscard]] VkBuffer sphereQuadBuffer() const noexcept { return sphere_quad_buffer_.buffer(); }
  [[nodiscard]] VkBuffer cylinderQuadBuffer() const noexcept {
    return cylinder_quad_buffer_.buffer();
  }
  [[nodiscard]] VkBuffer atomInstanceBuffer() const noexcept {
    return atom_instance_buffer_.buffer();
  }
  [[nodiscard]] uint32_t atomInstanceCount() const noexcept { return atom_instance_count_; }
  [[nodiscard]] VkBuffer bondInstanceBuffer() const noexcept {
    return bond_instance_buffer_.buffer();
  }
  [[nodiscard]] uint32_t bondInstanceCount() const noexcept { return bond_instance_count_; }
  [[nodiscard]] VkBuffer wireframeBuffer() const noexcept {
    return wireframe_vertex_buffer_.buffer();
  }
  [[nodiscard]] uint32_t wireframeVertexCount() const noexcept { return wireframe_vertex_count_; }
  [[nodiscard]] bool isValid() const noexcept { return descriptor_set_ != VK_NULL_HANDLE; }

private:
  void moveFrom(VulkanBufferManager &other) noexcept;
  [[nodiscard]] bool createSceneUbo();
  [[nodiscard]] bool createDescriptorPoolAndSet(VkDescriptorSetLayout ubo_layout);
  [[nodiscard]] bool createStaticQuadBuffers();
  void uploadUbo(const Camera3D &camera, const ViewportDimensions &dimensions);
  void uploadAtomInstances(const std::vector<AtomInstance> &atoms);
  void uploadBondInstances(const std::vector<BondInstance> &bonds);
  void uploadWireframeLines(const std::vector<LineSegment> &lines);

  VkDevice device_{VK_NULL_HANDLE};
  VkPhysicalDevice physical_device_{VK_NULL_HANDLE};
  VkDescriptorPool descriptor_pool_{VK_NULL_HANDLE};
  VkDescriptorSet descriptor_set_{VK_NULL_HANDLE};

  VulkanBuffer scene_ubo_;
  VulkanBuffer sphere_quad_buffer_;
  VulkanBuffer cylinder_quad_buffer_;
  VulkanBuffer atom_instance_buffer_;
  VulkanBuffer bond_instance_buffer_;
  VulkanBuffer wireframe_vertex_buffer_;

  uint32_t atom_instance_count_{0};
  uint32_t bond_instance_count_{0};
  uint32_t wireframe_vertex_count_{0};
};

#endif // CORRELATION_USE_VULKAN

} // namespace correlation::render3d::vulkan
