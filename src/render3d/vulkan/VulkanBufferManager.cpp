/**
 * @file VulkanBufferManager.cpp
 * @brief Implementation of VulkanBufferManager GPU allocations, instances, and UBO updates.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/vulkan/VulkanBufferManager.hpp"
#include <algorithm>
#include <cmath>
#include <utility>

namespace correlation::render3d::vulkan {

#ifdef CORRELATION_USE_VULKAN

namespace {

constexpr float INV_255 = 1.0F / 255.0F;

[[nodiscard]] std::array<float, 4> toFloatRgba(const RgbaColor &color) noexcept {
  return {
      static_cast<float>(color.r) * INV_255,
      static_cast<float>(color.g) * INV_255,
      static_cast<float>(color.b) * INV_255,
      static_cast<float>(color.a) * INV_255,
  };
}

} // namespace

void VulkanBufferManager::moveFrom(VulkanBufferManager &other) noexcept {
  device_ = std::exchange(other.device_, VK_NULL_HANDLE);
  physical_device_ = std::exchange(other.physical_device_, VK_NULL_HANDLE);
  descriptor_pool_ = std::exchange(other.descriptor_pool_, VK_NULL_HANDLE);
  descriptor_set_ = std::exchange(other.descriptor_set_, VK_NULL_HANDLE);

  scene_ubo_ = std::move(other.scene_ubo_);
  sphere_quad_buffer_ = std::move(other.sphere_quad_buffer_);
  cylinder_quad_buffer_ = std::move(other.cylinder_quad_buffer_);
  atom_instance_buffer_ = std::move(other.atom_instance_buffer_);
  bond_instance_buffer_ = std::move(other.bond_instance_buffer_);
  wireframe_vertex_buffer_ = std::move(other.wireframe_vertex_buffer_);

  atom_instance_count_ = std::exchange(other.atom_instance_count_, 0);
  bond_instance_count_ = std::exchange(other.bond_instance_count_, 0);
  wireframe_vertex_count_ = std::exchange(other.wireframe_vertex_count_, 0);
}

bool VulkanBufferManager::initialize(VkDevice device, VkPhysicalDevice physical_device,
                                     VkDescriptorSetLayout ubo_layout) {
  destroy();
  device_ = device;
  physical_device_ = physical_device;

  if (!createSceneUbo()) {
    destroy();
    return false;
  }

  if (!createDescriptorPoolAndSet(ubo_layout)) {
    destroy();
    return false;
  }

  if (!createStaticQuadBuffers()) {
    destroy();
    return false;
  }

  return true;
}

void VulkanBufferManager::destroy() noexcept {
  atom_instance_count_ = 0;
  bond_instance_count_ = 0;
  wireframe_vertex_count_ = 0;

  atom_instance_buffer_.destroy();
  bond_instance_buffer_.destroy();
  wireframe_vertex_buffer_.destroy();
  sphere_quad_buffer_.destroy();
  cylinder_quad_buffer_.destroy();
  scene_ubo_.destroy();

  if (device_ != VK_NULL_HANDLE) {
    if (descriptor_pool_ != VK_NULL_HANDLE) {
      vkDestroyDescriptorPool(device_, descriptor_pool_, nullptr);
      descriptor_pool_ = VK_NULL_HANDLE;
      descriptor_set_ = VK_NULL_HANDLE;
    }
    device_ = VK_NULL_HANDLE;
  }
  physical_device_ = VK_NULL_HANDLE;
}

bool VulkanBufferManager::createSceneUbo() {
  const BufferAllocationDesc ubo_desc{
      .size = sizeof(SceneUboData),
      .usage = VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT,
      .properties = VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
  };

  try {
    scene_ubo_.create(device_, physical_device_, ubo_desc);
  } catch (...) {
    return false;
  }
  return scene_ubo_.isValid();
}

bool VulkanBufferManager::createDescriptorPoolAndSet(VkDescriptorSetLayout ubo_layout) {
  const VkDescriptorPoolSize pool_size{
      .type = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
      .descriptorCount = 1,
  };

  const VkDescriptorPoolCreateInfo pool_info{
      .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .maxSets = 1,
      .poolSizeCount = 1,
      .pPoolSizes = &pool_size,
  };

  if (vkCreateDescriptorPool(device_, &pool_info, nullptr, &descriptor_pool_) != VK_SUCCESS) {
    return false;
  }

  const VkDescriptorSetAllocateInfo alloc_info{
      .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO,
      .pNext = nullptr,
      .descriptorPool = descriptor_pool_,
      .descriptorSetCount = 1,
      .pSetLayouts = &ubo_layout,
  };

  if (vkAllocateDescriptorSets(device_, &alloc_info, &descriptor_set_) != VK_SUCCESS) {
    return false;
  }

  const VkDescriptorBufferInfo buffer_info{
      .buffer = scene_ubo_.buffer(),
      .offset = 0,
      .range = sizeof(SceneUboData),
  };

  const VkWriteDescriptorSet descriptor_write{
      .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
      .pNext = nullptr,
      .dstSet = descriptor_set_,
      .dstBinding = 0,
      .dstArrayElement = 0,
      .descriptorCount = 1,
      .descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
      .pImageInfo = nullptr,
      .pBufferInfo = &buffer_info,
      .pTexelBufferView = nullptr,
  };

  vkUpdateDescriptorSets(device_, 1, &descriptor_write, 0, nullptr);
  return true;
}

bool VulkanBufferManager::createStaticQuadBuffers() {
  const std::array<float, 8> sphere_quad_vertices{
      -1.0F, -1.0F, 1.0F, -1.0F, -1.0F, 1.0F, 1.0F, 1.0F,
  };
  const std::array<float, 8> cylinder_quad_vertices{
      -1.0F, 0.0F, 1.0F, 0.0F, -1.0F, 1.0F, 1.0F, 1.0F,
  };

  const BufferAllocationDesc quad_desc{
      .size = sizeof(float) * 8,
      .usage = VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
      .properties = VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
  };

  try {
    sphere_quad_buffer_.create(device_, physical_device_, quad_desc);
    sphere_quad_buffer_.upload(std::as_bytes(std::span{sphere_quad_vertices}));

    cylinder_quad_buffer_.create(device_, physical_device_, quad_desc);
    cylinder_quad_buffer_.upload(std::as_bytes(std::span{cylinder_quad_vertices}));
  } catch (...) {
    return false;
  }

  return sphere_quad_buffer_.isValid() && cylinder_quad_buffer_.isValid();
}

void VulkanBufferManager::uploadScene(const Scene3D &scene, const Camera3D &camera,
                                      const ViewportDimensions &dimensions) {
  uploadUbo(camera, dimensions);
  uploadAtomInstances(scene.atoms());
  uploadBondInstances(scene.bonds());
  uploadWireframeLines(scene.unitCellLines());
}

void VulkanBufferManager::uploadUbo(const Camera3D &camera, const ViewportDimensions &dimensions) {
  Matrix4<float> proj = camera.projectionMatrix();
  proj(1, 1) = -proj(1, 1);

  const math::Vector3<float> eye = camera.eyePosition();
  const SceneUboData ubo_data{
      .view = camera.viewMatrix().data(),
      .proj = proj.data(),
      .camera_pos = {eye.x(), eye.y(), eye.z(), 1.0F},
      .light_dir = {0.0F, 0.5F, 1.0F, 0.0F},
      .viewport_size =
          {
              static_cast<float>(dimensions.width),
              static_cast<float>(dimensions.height),
          },
      .fov = camera.fov(),
      .padding = 0.0F,
  };

  const auto byte_span = std::as_bytes(std::span{&ubo_data, 1});
  scene_ubo_.upload(byte_span);
}

void VulkanBufferManager::uploadAtomInstances(const std::vector<AtomInstance> &atoms) {
  atom_instance_count_ = static_cast<uint32_t>(atoms.size());
  if (atoms.empty()) {
    return;
  }

  std::vector<AtomInstanceData> instances;
  instances.reserve(atoms.size());
  for (const auto &atom : atoms) {
    instances.push_back(AtomInstanceData{
        .center = {atom.position.x(), atom.position.y(), atom.position.z()},
        .radius = atom.radius,
        .color = toFloatRgba(atom.color),
        .halo_color = toFloatRgba(atom.halo_color),
        .style = {atom.halo_radius_ratio, atom.dim_factor},
    });
  }

  const auto byte_span = std::as_bytes(std::span{instances.data(), instances.size()});
  const BufferAllocationDesc desc{
      .size = static_cast<VkDeviceSize>(byte_span.size()),
      .usage = VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
      .properties = VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
  };

  atom_instance_buffer_.ensureCapacity(device_, physical_device_, desc);
  atom_instance_buffer_.upload(byte_span);
}

void VulkanBufferManager::uploadBondInstances(const std::vector<BondInstance> &bonds) {
  bond_instance_count_ = static_cast<uint32_t>(bonds.size());
  if (bonds.empty()) {
    return;
  }

  std::vector<BondInstanceData> instances;
  instances.reserve(bonds.size());
  for (const auto &bond : bonds) {
    instances.push_back(BondInstanceData{
        .start = {bond.start.x(), bond.start.y(), bond.start.z()},
        .end = {bond.end.x(), bond.end.y(), bond.end.z()},
        .radius = bond.radius,
        .color = toFloatRgba(bond.color),
    });
  }

  const auto byte_span = std::as_bytes(std::span{instances.data(), instances.size()});
  const BufferAllocationDesc desc{
      .size = static_cast<VkDeviceSize>(byte_span.size()),
      .usage = VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
      .properties = VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
  };

  bond_instance_buffer_.ensureCapacity(device_, physical_device_, desc);
  bond_instance_buffer_.upload(byte_span);
}

void VulkanBufferManager::uploadWireframeLines(const std::vector<LineSegment> &lines) {
  wireframe_vertex_count_ = static_cast<uint32_t>(lines.size() * 2);
  if (lines.empty()) {
    return;
  }

  std::vector<WireframeVertexData> vertices;
  vertices.reserve(lines.size() * 2);
  for (const auto &line : lines) {
    const std::array<float, 4> color = toFloatRgba(line.color);
    vertices.push_back(WireframeVertexData{
        .position = {line.start.x(), line.start.y(), line.start.z()},
        .color = color,
    });
    vertices.push_back(WireframeVertexData{
        .position = {line.end.x(), line.end.y(), line.end.z()},
        .color = color,
    });
  }

  const auto byte_span = std::as_bytes(std::span{vertices.data(), vertices.size()});
  const BufferAllocationDesc desc{
      .size = static_cast<VkDeviceSize>(byte_span.size()),
      .usage = VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
      .properties = VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
  };

  wireframe_vertex_buffer_.ensureCapacity(device_, physical_device_, desc);
  wireframe_vertex_buffer_.upload(byte_span);
}

#endif // CORRELATION_USE_VULKAN

} // namespace correlation::render3d::vulkan
