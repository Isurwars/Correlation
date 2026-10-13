/**
 * @file LodController.cpp
 * @brief Implementation of dynamic level-of-detail (LOD) policy controller.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/LodController.hpp"

namespace correlation::render3d {

LodRenderPolicy LodController::evaluatePolicy(std::size_t atom_count,
                                              bool is_interacting) const noexcept {
  if (atom_count <= high_threshold_) {
    return LodRenderPolicy{
        .level = LodLevel::High,
        .render_bonds = true,
        .render_wireframe = true,
        .quality_scale = 1.0F,
    };
  }

  if (atom_count <= medium_threshold_) {
    return LodRenderPolicy{
        .level = LodLevel::Medium,
        .render_bonds = !is_interacting,
        .render_wireframe = true,
        .quality_scale = is_interacting ? 0.6F : 0.8F,
    };
  }

  return LodRenderPolicy{
      .level = LodLevel::Low,
      .render_bonds = false,
      .render_wireframe = true,
      .quality_scale = is_interacting ? 0.35F : 0.5F,
  };
}

} // namespace correlation::render3d
