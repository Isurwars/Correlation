/**
 * @file LodController.hpp
 * @brief Dynamic level-of-detail (LOD) policy controller for atomic 3D rendering.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include <cstddef>

namespace correlation::render3d {

/**
 * @enum class LodLevel
 * @brief Level of detail tiers for geometry and shading complexity.
 */
enum class LodLevel {
  High,   ///< Systems < 5k atoms: Full analytical ray-cast spheres and cylinder bonds.
  Medium, ///< Systems 5k-25k atoms: Sphere impostors and line bonds.
  Low     ///< Systems > 25k atoms: Simplified point sprites, bonds suppressed during motion.
};

/**
 * @struct LodRenderPolicy
 * @brief Effective rendering policy derived from LOD state and camera motion.
 */
struct LodRenderPolicy {
  LodLevel level{LodLevel::High}; ///< Current active LOD tier.
  bool render_bonds{true};        ///< Whether bond geometry should be rendered.
  bool render_wireframe{true};    ///< Whether simulation box lines should be rendered.
  float quality_scale{1.0F};      ///< Sphere quality/resolution factor in [0.0, 1.0].
};

/**
 * @class LodController
 * @brief Evaluates active LOD policies based on atom count and camera interaction state.
 */
class LodController {
public:
  constexpr LodController() noexcept = default;

  /**
   * @brief Evaluates the active rendering policy.
   * @param atom_count Total atom count in current scene.
   * @param is_interacting True if camera is actively orbiting or panning.
   * @return Computed LodRenderPolicy.
   */
  [[nodiscard]] LodRenderPolicy evaluatePolicy(std::size_t atom_count,
                                               bool is_interacting = false) const noexcept;

  [[nodiscard]] std::size_t highThreshold() const noexcept { return high_threshold_; }
  [[nodiscard]] std::size_t mediumThreshold() const noexcept { return medium_threshold_; }

  void setHighThreshold(std::size_t threshold) noexcept { high_threshold_ = threshold; }
  void setMediumThreshold(std::size_t threshold) noexcept { medium_threshold_ = threshold; }

private:
  std::size_t high_threshold_{5000};
  std::size_t medium_threshold_{25000};
};

} // namespace correlation::render3d
