/**
 * @file SoftwareFramebuffer.hpp
 * @brief Cache-aligned RGBA8 color buffer and 32-bit floating-point depth buffer.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "render3d/ElementColors.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace correlation::render3d {

/**
 * @class SoftwareFramebuffer
 * @brief Manages RGBA8 color pixels and Z-buffer for software rasterization.
 */
class SoftwareFramebuffer {
public:
  /**
   * @brief Constructs an empty framebuffer.
   */
  SoftwareFramebuffer() = default;

  /**
   * @brief Constructs a framebuffer with specific pixel dimensions.
   * @param width Viewport pixel width.
   * @param height Viewport pixel height.
   */
  SoftwareFramebuffer(uint32_t width, uint32_t height) { resize(width, height); }

  /**
   * @brief Resizes color and depth buffers.
   * @param width New pixel width.
   * @param height New pixel height.
   */
  void resize(uint32_t width, uint32_t height) {
    width_ = width;
    height_ = height;
    const std::size_t pixel_count = static_cast<std::size_t>(width_) * height_;
    color_buffer_.resize(pixel_count * 4, 0);
    depth_buffer_.resize(pixel_count, 1.0F);
  }

  /**
   * @brief Clears color buffer to specified RGBA color.
   * @param color Fill color.
   */
  void clearColor(RgbaColor color) noexcept {
    const std::size_t pixel_count = static_cast<std::size_t>(width_) * height_;
    for (std::size_t i = 0; i < pixel_count; ++i) {
      const std::size_t offset = i * 4;
      color_buffer_[offset] = color.r;
      color_buffer_[offset + 1] = color.g;
      color_buffer_[offset + 2] = color.b;
      color_buffer_[offset + 3] = color.a;
    }
  }

  /**
   * @brief Clears depth buffer to specified value.
   * @param depth Depth value (default 1.0F).
   */
  void clearDepth(float depth = 1.0F) noexcept {
    std::fill(depth_buffer_.begin(), depth_buffer_.end(), depth);
  }

  /**
   * @brief Clears both color and depth buffers.
   * @param color Background fill color.
   * @param depth Clear depth value.
   */
  void clear(RgbaColor color, float depth = 1.0F) noexcept {
    clearColor(color);
    clearDepth(depth);
  }

  /**
   * @brief Tests depth and writes pixel if closer than current depth.
   * @param x X pixel coordinate.
   * @param y Y pixel coordinate.
   * @param depth Fragment depth value.
   * @param color Fragment RGBA color.
   * @return True if fragment passed depth test and was written.
   */
  bool testAndSetPixel(uint32_t x, uint32_t y, float depth, RgbaColor color) noexcept {
    if (x >= width_ || y >= height_) {
      return false;
    }
    const std::size_t idx = static_cast<std::size_t>(y) * width_ + x;
    if (depth < depth_buffer_[idx]) {
      depth_buffer_[idx] = depth;
      const std::size_t color_idx = idx * 4;
      color_buffer_[color_idx] = color.r;
      color_buffer_[color_idx + 1] = color.g;
      color_buffer_[color_idx + 2] = color.b;
      color_buffer_[color_idx + 3] = color.a;
      return true;
    }
    return false;
  }

  /**
   * @brief Reads depth value at pixel.
   * @param x X pixel coordinate.
   * @param y Y pixel coordinate.
   * @return Depth buffer value.
   */
  [[nodiscard]] float getDepth(uint32_t x, uint32_t y) const noexcept {
    if (x >= width_ || y >= height_) {
      return 1.0F;
    }
    return depth_buffer_[static_cast<std::size_t>(y) * width_ + x];
  }

  [[nodiscard]] uint32_t width() const noexcept { return width_; }
  [[nodiscard]] uint32_t height() const noexcept { return height_; }

  [[nodiscard]] std::span<const uint8_t> colorBytes() const noexcept { return color_buffer_; }

  [[nodiscard]] const std::vector<uint8_t> &colorBuffer() const noexcept { return color_buffer_; }

  [[nodiscard]] const std::vector<float> &depthBuffer() const noexcept { return depth_buffer_; }

private:
  uint32_t width_{0};
  uint32_t height_{0};
  std::vector<uint8_t> color_buffer_;
  std::vector<float> depth_buffer_;
};

} // namespace correlation::render3d
