/**
 * @file TrajectoryStreamer.hpp
 * @brief Zero-copy trajectory frame navigation and playback streamer for 3D scenes.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/Trajectory.hpp"
#include "render3d/Scene3D.hpp"
#include <cstddef>
#include <memory>

namespace correlation::render3d {

/**
 * @class TrajectoryStreamer
 * @brief Streams trajectory snapshot frames into Scene3D with zero vector reallocation.
 */
class TrajectoryStreamer {
public:
  /**
   * @brief Constructs streamer referencing a trajectory dataset.
   * @param[in] trajectory Shared pointer to simulation trajectory.
   */
  explicit TrajectoryStreamer(
      std::shared_ptr<const core::Trajectory> trajectory = nullptr) noexcept;

  /**
   * @brief Sets or replaces the active trajectory source.
   * @param[in] trajectory Shared pointer to trajectory.
   */
  void setTrajectory(std::shared_ptr<const core::Trajectory> trajectory) noexcept;

  /**
   * @brief Seeks to a specific frame and updates the scene in-place.
   * @param[in] frame_index Zero-based frame index.
   * @param[in,out] scene Target scene to update.
   * @param[in] recompute_bonds If true, recalculates bond pairs for this frame.
   * @param[in] bond_cutoff_scale Multiplier for covalent radius bond cutoff.
   * @param[in] in_cell If true, wraps coordinates into [0, 1) and generates boundary-clipped bonds.
   * @return True if frame exists and scene was updated.
   */
  bool seekFrame(std::size_t frame_index, Scene3D &scene, bool recompute_bonds = false,
                 double bond_cutoff_scale = 1.15, bool in_cell = false);

  /**
   * @brief Advances to the next frame.
   * @param[in,out] scene Target scene.
   * @param[in] loop If true, loops back to frame 0 upon reaching the end.
   * @param[in] recompute_bonds Whether to re-evaluate bonds.
   * @param[in] bond_cutoff_scale Multiplier for covalent radius bond cutoff.
   * @param[in] in_cell If true, wraps coordinates into [0, 1) and generates boundary-clipped bonds.
   * @return True if advance succeeded.
   */
  bool stepNext(Scene3D &scene, bool loop = true, bool recompute_bonds = false,
                double bond_cutoff_scale = 1.15, bool in_cell = false);

  /**
   * @brief Steps backward to the previous frame.
   * @param[in,out] scene Target scene.
   * @param[in] loop If true, wraps to the last frame from frame 0.
   * @param[in] recompute_bonds Whether to re-evaluate bonds.
   * @param[in] bond_cutoff_scale Multiplier for covalent radius bond cutoff.
   * @param[in] in_cell If true, wraps coordinates into [0, 1) and generates boundary-clipped bonds.
   * @return True if step succeeded.
   */
  bool stepPrev(Scene3D &scene, bool loop = true, bool recompute_bonds = false,
                double bond_cutoff_scale = 1.15, bool in_cell = false);

  /**
   * @brief Gets current active frame index.
   * @return Zero-based frame index.
   */
  [[nodiscard]] std::size_t currentFrameIndex() const noexcept { return current_frame_; }

  /**
   * @brief Gets total frame count in the active trajectory.
   * @return Total number of frames, or 0 if trajectory is unset.
   */
  [[nodiscard]] std::size_t totalFrames() const noexcept;

  /**
   * @brief Checks if a valid trajectory dataset is attached.
   * @return True if trajectory is loaded.
   */
  [[nodiscard]] bool hasTrajectory() const noexcept { return trajectory_ != nullptr; }

  /**
   * @brief Sets target playback frame rate.
   * @param[in] fps Playback rate in frames per second.
   */
  void setPlaybackFps(float fps) noexcept { playback_fps_ = fps; }

  /**
   * @brief Gets current playback frame rate.
   * @return Playback rate in frames per second.
   */
  [[nodiscard]] float playbackFps() const noexcept { return playback_fps_; }

private:
  std::shared_ptr<const core::Trajectory> trajectory_;
  std::size_t current_frame_{0};
  float playback_fps_{30.0F};
};

} // namespace correlation::render3d
