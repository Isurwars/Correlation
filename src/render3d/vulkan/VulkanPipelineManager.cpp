/**
 * @file VulkanPipelineManager.cpp
 * @brief Implementation of graphics pipelines and layouts for 3D impostors and wireframes.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/vulkan/VulkanPipelineManager.hpp"
#include "render3d/vulkan/VulkanCommon.hpp"
#include "render3d/vulkan/shaders/CylinderImpostorSpv.hpp"
#include "render3d/vulkan/shaders/LineWireframeSpv.hpp"
#include "render3d/vulkan/shaders/SphereImpostorSpv.hpp"
#include <array>
#include <cstdint>
#include <utility>

namespace correlation::render3d::vulkan {

#ifdef CORRELATION_USE_VULKAN

namespace {

VkPipelineColorBlendAttachmentState createBlendAttachmentState() noexcept {
  const auto write_mask = static_cast<uint32_t>(VK_COLOR_COMPONENT_R_BIT) |
                          static_cast<uint32_t>(VK_COLOR_COMPONENT_G_BIT) |
                          static_cast<uint32_t>(VK_COLOR_COMPONENT_B_BIT) |
                          static_cast<uint32_t>(VK_COLOR_COMPONENT_A_BIT);
  return VkPipelineColorBlendAttachmentState{
      .blendEnable = VK_TRUE,
      .srcColorBlendFactor = VK_BLEND_FACTOR_SRC_ALPHA,
      .dstColorBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA,
      .colorBlendOp = VK_BLEND_OP_ADD,
      .srcAlphaBlendFactor = VK_BLEND_FACTOR_ONE,
      .dstAlphaBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA,
      .alphaBlendOp = VK_BLEND_OP_ADD,
      .colorWriteMask = write_mask,
  };
}

VkPipelineRasterizationStateCreateInfo createRasterizationState() noexcept {
  return VkPipelineRasterizationStateCreateInfo{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .depthClampEnable = VK_FALSE,
      .rasterizerDiscardEnable = VK_FALSE,
      .polygonMode = VK_POLYGON_MODE_FILL,
      .cullMode = VK_CULL_MODE_NONE,
      .frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE,
      .depthBiasEnable = VK_FALSE,
      .depthBiasConstantFactor = 0.0F,
      .depthBiasClamp = 0.0F,
      .depthBiasSlopeFactor = 0.0F,
      .lineWidth = 1.0F,
  };
}

VkPipelineMultisampleStateCreateInfo createMultisampleState() noexcept {
  return VkPipelineMultisampleStateCreateInfo{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .rasterizationSamples = VK_SAMPLE_COUNT_1_BIT,
      .sampleShadingEnable = VK_FALSE,
      .minSampleShading = 1.0F,
      .pSampleMask = nullptr,
      .alphaToCoverageEnable = VK_FALSE,
      .alphaToOneEnable = VK_FALSE,
  };
}

VkPipelineDepthStencilStateCreateInfo createDepthStencilState() noexcept {
  return VkPipelineDepthStencilStateCreateInfo{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .depthTestEnable = VK_TRUE,
      .depthWriteEnable = VK_TRUE,
      .depthCompareOp = VK_COMPARE_OP_LESS_OR_EQUAL,
      .depthBoundsTestEnable = VK_FALSE,
      .stencilTestEnable = VK_FALSE,
      .front = {},
      .back = {},
      .minDepthBounds = 0.0F,
      .maxDepthBounds = 1.0F,
  };
}

} // namespace

void VulkanPipelineManager::moveFrom(VulkanPipelineManager &other) noexcept {
  device_ = std::exchange(other.device_, VK_NULL_HANDLE);
  descriptor_set_layout_ = std::exchange(other.descriptor_set_layout_, VK_NULL_HANDLE);
  pipeline_layout_ = std::exchange(other.pipeline_layout_, VK_NULL_HANDLE);
  sphere_pipeline_ = std::exchange(other.sphere_pipeline_, VK_NULL_HANDLE);
  cylinder_pipeline_ = std::exchange(other.cylinder_pipeline_, VK_NULL_HANDLE);
  wireframe_pipeline_ = std::exchange(other.wireframe_pipeline_, VK_NULL_HANDLE);
}

bool VulkanPipelineManager::initialize(VkDevice device, VkRenderPass render_pass) {
  destroy();
  device_ = device;

  if (!createDescriptorSetLayout() || !createPipelineLayout()) {
    destroy();
    return false;
  }

  if (!createSpherePipeline(render_pass) || !createCylinderPipeline(render_pass) ||
      !createWireframePipeline(render_pass)) {
    destroy();
    return false;
  }

  return true;
}

void VulkanPipelineManager::destroy() noexcept {
  if (device_ != VK_NULL_HANDLE) {
    if (sphere_pipeline_ != VK_NULL_HANDLE) {
      vkDestroyPipeline(device_, sphere_pipeline_, nullptr);
      sphere_pipeline_ = VK_NULL_HANDLE;
    }
    if (cylinder_pipeline_ != VK_NULL_HANDLE) {
      vkDestroyPipeline(device_, cylinder_pipeline_, nullptr);
      cylinder_pipeline_ = VK_NULL_HANDLE;
    }
    if (wireframe_pipeline_ != VK_NULL_HANDLE) {
      vkDestroyPipeline(device_, wireframe_pipeline_, nullptr);
      wireframe_pipeline_ = VK_NULL_HANDLE;
    }
    if (pipeline_layout_ != VK_NULL_HANDLE) {
      vkDestroyPipelineLayout(device_, pipeline_layout_, nullptr);
      pipeline_layout_ = VK_NULL_HANDLE;
    }
    if (descriptor_set_layout_ != VK_NULL_HANDLE) {
      vkDestroyDescriptorSetLayout(device_, descriptor_set_layout_, nullptr);
      descriptor_set_layout_ = VK_NULL_HANDLE;
    }
    device_ = VK_NULL_HANDLE;
  }
}

VkShaderModule VulkanPipelineManager::createShaderModule(std::span<const uint32_t> spv_code) const {
  const VkShaderModuleCreateInfo create_info{
      .sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .codeSize = spv_code.size_bytes(),
      .pCode = spv_code.data(),
  };

  VkShaderModule shader_module{VK_NULL_HANDLE};
  if (vkCreateShaderModule(device_, &create_info, nullptr, &shader_module) != VK_SUCCESS) {
    return VK_NULL_HANDLE;
  }
  return shader_module;
}

bool VulkanPipelineManager::createDescriptorSetLayout() {
  const VkDescriptorSetLayoutBinding ubo_binding{
      .binding = 0,
      .descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
      .descriptorCount = 1,
      .stageFlags = VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT,
      .pImmutableSamplers = nullptr,
  };

  const VkDescriptorSetLayoutCreateInfo layout_info{
      .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .bindingCount = 1,
      .pBindings = &ubo_binding,
  };

  return vkCreateDescriptorSetLayout(device_, &layout_info, nullptr, &descriptor_set_layout_) ==
         VK_SUCCESS;
}

bool VulkanPipelineManager::createPipelineLayout() {
  const VkPipelineLayoutCreateInfo pipeline_layout_info{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .setLayoutCount = 1,
      .pSetLayouts = &descriptor_set_layout_,
      .pushConstantRangeCount = 0,
      .pPushConstantRanges = nullptr,
  };

  return vkCreatePipelineLayout(device_, &pipeline_layout_info, nullptr, &pipeline_layout_) ==
         VK_SUCCESS;
}

bool VulkanPipelineManager::createSpherePipeline(VkRenderPass render_pass) {
  VkShaderModule vert_module = createShaderModule(shaders::getSphereImpostorVertSpv());
  VkShaderModule frag_module = createShaderModule(shaders::getSphereImpostorFragSpv());
  if (vert_module == VK_NULL_HANDLE || frag_module == VK_NULL_HANDLE) {
    if (vert_module != VK_NULL_HANDLE) {
      vkDestroyShaderModule(device_, vert_module, nullptr);
    }
    if (frag_module != VK_NULL_HANDLE) {
      vkDestroyShaderModule(device_, frag_module, nullptr);
    }
    return false;
  }

  const std::array<VkPipelineShaderStageCreateInfo, 2> stages{
      VkPipelineShaderStageCreateInfo{
          .sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
          .pNext = nullptr,
          .flags = 0,
          .stage = VK_SHADER_STAGE_VERTEX_BIT,
          .module = vert_module,
          .pName = "main",
          .pSpecializationInfo = nullptr,
      },
      VkPipelineShaderStageCreateInfo{
          .sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
          .pNext = nullptr,
          .flags = 0,
          .stage = VK_SHADER_STAGE_FRAGMENT_BIT,
          .module = frag_module,
          .pName = "main",
          .pSpecializationInfo = nullptr,
      },
  };

  const std::array<VkVertexInputBindingDescription, 2> bindings{
      VkVertexInputBindingDescription{
          .binding = 0,
          .stride = 8,
          .inputRate = VK_VERTEX_INPUT_RATE_VERTEX,
      },
      VkVertexInputBindingDescription{
          .binding = 1,
          .stride = 56,
          .inputRate = VK_VERTEX_INPUT_RATE_INSTANCE,
      },
  };

  const std::array<VkVertexInputAttributeDescription, 6> attributes{
      VkVertexInputAttributeDescription{
          .location = 0,
          .binding = 0,
          .format = VK_FORMAT_R32G32_SFLOAT,
          .offset = 0,
      },
      VkVertexInputAttributeDescription{
          .location = 1,
          .binding = 1,
          .format = VK_FORMAT_R32G32B32_SFLOAT,
          .offset = 0,
      },
      VkVertexInputAttributeDescription{
          .location = 2,
          .binding = 1,
          .format = VK_FORMAT_R32_SFLOAT,
          .offset = 12,
      },
      VkVertexInputAttributeDescription{
          .location = 3,
          .binding = 1,
          .format = VK_FORMAT_R32G32B32A32_SFLOAT,
          .offset = 16,
      },
      VkVertexInputAttributeDescription{
          .location = 4,
          .binding = 1,
          .format = VK_FORMAT_R32G32B32A32_SFLOAT,
          .offset = 32,
      },
      VkVertexInputAttributeDescription{
          .location = 5,
          .binding = 1,
          .format = VK_FORMAT_R32G32_SFLOAT,
          .offset = 48,
      },
  };

  const VkPipelineVertexInputStateCreateInfo vertex_input{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .vertexBindingDescriptionCount = static_cast<uint32_t>(bindings.size()),
      .pVertexBindingDescriptions = bindings.data(),
      .vertexAttributeDescriptionCount = static_cast<uint32_t>(attributes.size()),
      .pVertexAttributeDescriptions = attributes.data(),
  };

  const VkPipelineInputAssemblyStateCreateInfo input_assembly{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_STRIP,
      .primitiveRestartEnable = VK_FALSE,
  };

  const std::array<VkDynamicState, 2> dynamic_states{
      VK_DYNAMIC_STATE_VIEWPORT,
      VK_DYNAMIC_STATE_SCISSOR,
  };
  const VkPipelineDynamicStateCreateInfo dynamic_state{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .dynamicStateCount = static_cast<uint32_t>(dynamic_states.size()),
      .pDynamicStates = dynamic_states.data(),
  };

  const VkPipelineViewportStateCreateInfo viewport_state{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .viewportCount = 1,
      .pViewports = nullptr,
      .scissorCount = 1,
      .pScissors = nullptr,
  };

  const VkPipelineRasterizationStateCreateInfo rasterization = createRasterizationState();
  const VkPipelineMultisampleStateCreateInfo multisample = createMultisampleState();
  const VkPipelineDepthStencilStateCreateInfo depth_stencil = createDepthStencilState();
  const VkPipelineColorBlendAttachmentState color_blend_attachment = createBlendAttachmentState();

  const VkPipelineColorBlendStateCreateInfo color_blend{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .logicOpEnable = VK_FALSE,
      .logicOp = VK_LOGIC_OP_COPY,
      .attachmentCount = 1,
      .pAttachments = &color_blend_attachment,
      .blendConstants = {0.0F, 0.0F, 0.0F, 0.0F},
  };

  const VkGraphicsPipelineCreateInfo pipeline_info{
      .sType = VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .stageCount = static_cast<uint32_t>(stages.size()),
      .pStages = stages.data(),
      .pVertexInputState = &vertex_input,
      .pInputAssemblyState = &input_assembly,
      .pTessellationState = nullptr,
      .pViewportState = &viewport_state,
      .pRasterizationState = &rasterization,
      .pMultisampleState = &multisample,
      .pDepthStencilState = &depth_stencil,
      .pColorBlendState = &color_blend,
      .pDynamicState = &dynamic_state,
      .layout = pipeline_layout_,
      .renderPass = render_pass,
      .subpass = 0,
      .basePipelineHandle = VK_NULL_HANDLE,
      .basePipelineIndex = -1,
  };

  const VkResult result = vkCreateGraphicsPipelines(device_, VK_NULL_HANDLE, 1, &pipeline_info,
                                                    nullptr, &sphere_pipeline_);

  vkDestroyShaderModule(device_, vert_module, nullptr);
  vkDestroyShaderModule(device_, frag_module, nullptr);

  return result == VK_SUCCESS;
}

bool VulkanPipelineManager::createCylinderPipeline(VkRenderPass render_pass) {
  VkShaderModule vert_module = createShaderModule(shaders::getCylinderImpostorVertSpv());
  VkShaderModule frag_module = createShaderModule(shaders::getCylinderImpostorFragSpv());
  if (vert_module == VK_NULL_HANDLE || frag_module == VK_NULL_HANDLE) {
    if (vert_module != VK_NULL_HANDLE) {
      vkDestroyShaderModule(device_, vert_module, nullptr);
    }
    if (frag_module != VK_NULL_HANDLE) {
      vkDestroyShaderModule(device_, frag_module, nullptr);
    }
    return false;
  }

  const std::array<VkPipelineShaderStageCreateInfo, 2> stages{
      VkPipelineShaderStageCreateInfo{
          .sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
          .pNext = nullptr,
          .flags = 0,
          .stage = VK_SHADER_STAGE_VERTEX_BIT,
          .module = vert_module,
          .pName = "main",
          .pSpecializationInfo = nullptr,
      },
      VkPipelineShaderStageCreateInfo{
          .sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
          .pNext = nullptr,
          .flags = 0,
          .stage = VK_SHADER_STAGE_FRAGMENT_BIT,
          .module = frag_module,
          .pName = "main",
          .pSpecializationInfo = nullptr,
      },
  };

  const std::array<VkVertexInputBindingDescription, 2> bindings{
      VkVertexInputBindingDescription{
          .binding = 0,
          .stride = 8,
          .inputRate = VK_VERTEX_INPUT_RATE_VERTEX,
      },
      VkVertexInputBindingDescription{
          .binding = 1,
          .stride = 44,
          .inputRate = VK_VERTEX_INPUT_RATE_INSTANCE,
      },
  };

  const std::array<VkVertexInputAttributeDescription, 5> attributes{
      VkVertexInputAttributeDescription{
          .location = 0,
          .binding = 0,
          .format = VK_FORMAT_R32G32_SFLOAT,
          .offset = 0,
      },
      VkVertexInputAttributeDescription{
          .location = 1,
          .binding = 1,
          .format = VK_FORMAT_R32G32B32_SFLOAT,
          .offset = 0,
      },
      VkVertexInputAttributeDescription{
          .location = 2,
          .binding = 1,
          .format = VK_FORMAT_R32G32B32_SFLOAT,
          .offset = 12,
      },
      VkVertexInputAttributeDescription{
          .location = 3,
          .binding = 1,
          .format = VK_FORMAT_R32_SFLOAT,
          .offset = 24,
      },
      VkVertexInputAttributeDescription{
          .location = 4,
          .binding = 1,
          .format = VK_FORMAT_R32G32B32A32_SFLOAT,
          .offset = 28,
      },
  };

  const VkPipelineVertexInputStateCreateInfo vertex_input{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .vertexBindingDescriptionCount = static_cast<uint32_t>(bindings.size()),
      .pVertexBindingDescriptions = bindings.data(),
      .vertexAttributeDescriptionCount = static_cast<uint32_t>(attributes.size()),
      .pVertexAttributeDescriptions = attributes.data(),
  };

  const VkPipelineInputAssemblyStateCreateInfo input_assembly{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_STRIP,
      .primitiveRestartEnable = VK_FALSE,
  };

  const std::array<VkDynamicState, 2> dynamic_states{
      VK_DYNAMIC_STATE_VIEWPORT,
      VK_DYNAMIC_STATE_SCISSOR,
  };
  const VkPipelineDynamicStateCreateInfo dynamic_state{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .dynamicStateCount = static_cast<uint32_t>(dynamic_states.size()),
      .pDynamicStates = dynamic_states.data(),
  };

  const VkPipelineViewportStateCreateInfo viewport_state{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .viewportCount = 1,
      .pViewports = nullptr,
      .scissorCount = 1,
      .pScissors = nullptr,
  };

  const VkPipelineRasterizationStateCreateInfo rasterization = createRasterizationState();
  const VkPipelineMultisampleStateCreateInfo multisample = createMultisampleState();
  const VkPipelineDepthStencilStateCreateInfo depth_stencil = createDepthStencilState();
  const VkPipelineColorBlendAttachmentState color_blend_attachment = createBlendAttachmentState();

  const VkPipelineColorBlendStateCreateInfo color_blend{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .logicOpEnable = VK_FALSE,
      .logicOp = VK_LOGIC_OP_COPY,
      .attachmentCount = 1,
      .pAttachments = &color_blend_attachment,
      .blendConstants = {0.0F, 0.0F, 0.0F, 0.0F},
  };

  const VkGraphicsPipelineCreateInfo pipeline_info{
      .sType = VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .stageCount = static_cast<uint32_t>(stages.size()),
      .pStages = stages.data(),
      .pVertexInputState = &vertex_input,
      .pInputAssemblyState = &input_assembly,
      .pTessellationState = nullptr,
      .pViewportState = &viewport_state,
      .pRasterizationState = &rasterization,
      .pMultisampleState = &multisample,
      .pDepthStencilState = &depth_stencil,
      .pColorBlendState = &color_blend,
      .pDynamicState = &dynamic_state,
      .layout = pipeline_layout_,
      .renderPass = render_pass,
      .subpass = 0,
      .basePipelineHandle = VK_NULL_HANDLE,
      .basePipelineIndex = -1,
  };

  const VkResult result = vkCreateGraphicsPipelines(device_, VK_NULL_HANDLE, 1, &pipeline_info,
                                                    nullptr, &cylinder_pipeline_);

  vkDestroyShaderModule(device_, vert_module, nullptr);
  vkDestroyShaderModule(device_, frag_module, nullptr);

  return result == VK_SUCCESS;
}

bool VulkanPipelineManager::createWireframePipeline(VkRenderPass render_pass) {
  VkShaderModule vert_module = createShaderModule(shaders::getLineWireframeVertSpv());
  VkShaderModule frag_module = createShaderModule(shaders::getLineWireframeFragSpv());
  if (vert_module == VK_NULL_HANDLE || frag_module == VK_NULL_HANDLE) {
    if (vert_module != VK_NULL_HANDLE) {
      vkDestroyShaderModule(device_, vert_module, nullptr);
    }
    if (frag_module != VK_NULL_HANDLE) {
      vkDestroyShaderModule(device_, frag_module, nullptr);
    }
    return false;
  }

  const std::array<VkPipelineShaderStageCreateInfo, 2> stages{
      VkPipelineShaderStageCreateInfo{
          .sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
          .pNext = nullptr,
          .flags = 0,
          .stage = VK_SHADER_STAGE_VERTEX_BIT,
          .module = vert_module,
          .pName = "main",
          .pSpecializationInfo = nullptr,
      },
      VkPipelineShaderStageCreateInfo{
          .sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
          .pNext = nullptr,
          .flags = 0,
          .stage = VK_SHADER_STAGE_FRAGMENT_BIT,
          .module = frag_module,
          .pName = "main",
          .pSpecializationInfo = nullptr,
      },
  };

  const std::array<VkVertexInputBindingDescription, 1> bindings{
      VkVertexInputBindingDescription{
          .binding = 0,
          .stride = 28,
          .inputRate = VK_VERTEX_INPUT_RATE_VERTEX,
      },
  };

  const std::array<VkVertexInputAttributeDescription, 2> attributes{
      VkVertexInputAttributeDescription{
          .location = 0,
          .binding = 0,
          .format = VK_FORMAT_R32G32B32_SFLOAT,
          .offset = 0,
      },
      VkVertexInputAttributeDescription{
          .location = 1,
          .binding = 0,
          .format = VK_FORMAT_R32G32B32A32_SFLOAT,
          .offset = 12,
      },
  };

  const VkPipelineVertexInputStateCreateInfo vertex_input{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .vertexBindingDescriptionCount = static_cast<uint32_t>(bindings.size()),
      .pVertexBindingDescriptions = bindings.data(),
      .vertexAttributeDescriptionCount = static_cast<uint32_t>(attributes.size()),
      .pVertexAttributeDescriptions = attributes.data(),
  };

  const VkPipelineInputAssemblyStateCreateInfo input_assembly{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .topology = VK_PRIMITIVE_TOPOLOGY_LINE_LIST,
      .primitiveRestartEnable = VK_FALSE,
  };

  const std::array<VkDynamicState, 2> dynamic_states{
      VK_DYNAMIC_STATE_VIEWPORT,
      VK_DYNAMIC_STATE_SCISSOR,
  };
  const VkPipelineDynamicStateCreateInfo dynamic_state{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .dynamicStateCount = static_cast<uint32_t>(dynamic_states.size()),
      .pDynamicStates = dynamic_states.data(),
  };

  const VkPipelineViewportStateCreateInfo viewport_state{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .viewportCount = 1,
      .pViewports = nullptr,
      .scissorCount = 1,
      .pScissors = nullptr,
  };

  const VkPipelineRasterizationStateCreateInfo rasterization = createRasterizationState();
  const VkPipelineMultisampleStateCreateInfo multisample = createMultisampleState();
  const VkPipelineDepthStencilStateCreateInfo depth_stencil = createDepthStencilState();
  const VkPipelineColorBlendAttachmentState color_blend_attachment = createBlendAttachmentState();

  const VkPipelineColorBlendStateCreateInfo color_blend{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .logicOpEnable = VK_FALSE,
      .logicOp = VK_LOGIC_OP_COPY,
      .attachmentCount = 1,
      .pAttachments = &color_blend_attachment,
      .blendConstants = {0.0F, 0.0F, 0.0F, 0.0F},
  };

  const VkGraphicsPipelineCreateInfo pipeline_info{
      .sType = VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .stageCount = static_cast<uint32_t>(stages.size()),
      .pStages = stages.data(),
      .pVertexInputState = &vertex_input,
      .pInputAssemblyState = &input_assembly,
      .pTessellationState = nullptr,
      .pViewportState = &viewport_state,
      .pRasterizationState = &rasterization,
      .pMultisampleState = &multisample,
      .pDepthStencilState = &depth_stencil,
      .pColorBlendState = &color_blend,
      .pDynamicState = &dynamic_state,
      .layout = pipeline_layout_,
      .renderPass = render_pass,
      .subpass = 0,
      .basePipelineHandle = VK_NULL_HANDLE,
      .basePipelineIndex = -1,
  };

  const VkResult result = vkCreateGraphicsPipelines(device_, VK_NULL_HANDLE, 1, &pipeline_info,
                                                    nullptr, &wireframe_pipeline_);

  vkDestroyShaderModule(device_, vert_module, nullptr);
  vkDestroyShaderModule(device_, frag_module, nullptr);

  return result == VK_SUCCESS;
}

#endif // CORRELATION_USE_VULKAN

} // namespace correlation::render3d::vulkan
