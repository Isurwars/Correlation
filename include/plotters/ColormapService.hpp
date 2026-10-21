/**
 * @file ColormapService.hpp
 * @brief High-throughput scientific colormap engine and precomputed LUT service.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "render3d/ElementColors.hpp"

#include <cstdint>
#include <optional>
#include <span>
#include <string_view>

namespace correlation::plotters {

using render3d::RgbaColor;

/**
 * @enum ColormapType
 * @brief Supported scientific colormaps for scalar mapping.
 */
enum class ColormapType : uint8_t { Viridis = 0, Plasma = 1, Coolwarm = 2, Turbo = 3 };

/**
 * @class ColormapService
 * @brief High-throughput, perceptually uniform scientific colormap engine.
 *
 * Provides O(1) zero-allocation scalar sampling and access to 256-entry
 * precomputed look-up tables (LUTs) for continuous scalar field rendering.
 */
class ColormapService {
public:
  ColormapService() = delete;

  /**
   * @brief Samples a color from the specified colormap at normalized parameter value in [0.0, 1.0].
   * Values outside [0.0, 1.0] are clamped to boundaries.
   * @param value Normalized scalar value in [0.0, 1.0].
   * @param type Colormap scheme to sample from.
   * @return 8-bit RGBA color representation.
   */
  [[nodiscard]] static RgbaColor sample(double value, ColormapType type) noexcept;

  /**
   * @brief Returns a non-owning span over the 256-entry precomputed look-up table.
   * @param type Colormap scheme.
   * @return 256-entry span of RgbaColor.
   */
  [[nodiscard]] static std::span<const RgbaColor, 256> getLut(ColormapType type) noexcept;

  /**
   * @brief Returns the human-readable display name of a colormap.
   * @param type Colormap scheme.
   * @return String view of the name (e.g. "Viridis").
   */
  [[nodiscard]] static std::string_view name(ColormapType type) noexcept;

  /**
   * @brief Parses a colormap from its display name (case-insensitive).
   * @param name String view of the name.
   * @return Optional ColormapType if recognized.
   */
  [[nodiscard]] static std::optional<ColormapType> fromString(std::string_view name) noexcept;
};

} // namespace correlation::plotters
