/**
 * @file Camera3D.hpp
 * @brief Orbit and arcball camera model for interactive 3D molecular visualization.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "render3d/Math3D.hpp"

namespace correlation::render3d {

/**
 * @class Camera3D
 * @brief Manages camera positioning, orbit rotations, panning, zooming, and view-projection
 * matrices.
 */
class Camera3D {
public:
  /**
   * @brief Constructs camera with default orbital parameters.
   */
  Camera3D() noexcept;

  /**
   * @brief Orbits the camera around the current target.
   * @param delta_azimuth_deg Horizontal angle change in degrees.
   * @param delta_elevation_deg Vertical angle change in degrees.
   */
  void orbit(float delta_azimuth_deg, float delta_elevation_deg) noexcept;

  /**
   * @brief Pans the camera target along the view-aligned screen plane.
   * @param delta_x Horizontal translation delta.
   * @param delta_y Vertical translation delta.
   */
  void pan(float delta_x, float delta_y) noexcept;

  /**
   * @brief Zooms camera by adjusting orbital distance.
   * @param factor Zoom factor (positive zooms in, negative zooms out).
   */
  void zoom(float factor) noexcept;

  /**
   * @brief Automatically frames the camera around a 3D bounding box.
   * @param min_pt Minimum bounding box corner.
   * @param max_pt Maximum bounding box corner.
   */
  void fitToBox(const math::Vector3<float> &min_pt, const math::Vector3<float> &max_pt) noexcept;

  /**
   * @brief Updates viewport aspect ratio.
   * @param width Viewport pixel width.
   * @param height Viewport pixel height.
   */
  void setViewport(float width, float height) noexcept;

  /**
   * @brief Sets vertical field of view.
   * @param fov_y_deg FOV in degrees.
   */
  void setFov(float fov_y_deg) noexcept;

  /**
   * @brief Sets target look-at point.
   * @param target World coordinates of the focus target.
   */
  void setTarget(const math::Vector3<float> &target) noexcept;

  /**
   * @brief Computes eye position in world space.
   * @return Eye coordinate.
   */
  [[nodiscard]] math::Vector3<float> eyePosition() const noexcept;

  [[nodiscard]] const math::Vector3<float> &target() const noexcept { return target_; }
  [[nodiscard]] float distance() const noexcept { return distance_; }
  [[nodiscard]] float azimuth() const noexcept { return azimuth_deg_; }
  [[nodiscard]] float elevation() const noexcept { return elevation_deg_; }
  [[nodiscard]] float fov() const noexcept { return fov_y_deg_; }
  [[nodiscard]] float aspectRatio() const noexcept { return aspect_ratio_; }

  /**
   * @brief Generates current view matrix.
   * @return 4x4 View matrix.
   */
  [[nodiscard]] Matrix4<float> viewMatrix() const;

  /**
   * @brief Generates perspective projection matrix.
   * @return 4x4 Projection matrix.
   */
  [[nodiscard]] Matrix4<float> projectionMatrix() const;

  /**
   * @brief Generates combined view-projection matrix (P * V).
   * @return 4x4 View-projection matrix.
   */
  [[nodiscard]] Matrix4<float> viewProjectionMatrix() const;

private:
  math::Vector3<float> target_{0.0F, 0.0F, 0.0F};
  math::Vector3<float> world_up_{0.0F, 1.0F, 0.0F};
  float distance_{20.0F};
  float azimuth_deg_{0.0F};
  float elevation_deg_{20.0F};
  float fov_y_deg_{45.0F};
  float aspect_ratio_{1.3333F};
  float z_near_{0.1F};
  float z_far_{2000.0F};
};

} // namespace correlation::render3d
