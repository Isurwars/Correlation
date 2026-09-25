/**
 * @file TrajectoryLoader.hpp
 * @brief Service responsible for loading trajectories and extracting structural metadata.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "math/Precision.hpp"

#include <expected>
#include <functional>
#include <map>
#include <memory>
#include <string>

namespace correlation::app {

/**
 * @class TrajectoryLoader
 * @brief Manages loading of structure/trajectory files and provides atomic/cell metadata.
 */
class TrajectoryLoader {
public:
  using ProgressCallback = std::function<void(float, const std::string &)>;

  TrajectoryLoader() = default;
  ~TrajectoryLoader() = default;

  TrajectoryLoader(const TrajectoryLoader &) = delete;
  TrajectoryLoader &operator=(const TrajectoryLoader &) = delete;
  TrajectoryLoader(TrajectoryLoader &&) noexcept = default;
  TrajectoryLoader &operator=(TrajectoryLoader &&) noexcept = default;

  /**
   * @brief Loads a structure or trajectory file from the given path.
   * @param[in] path Absolute or relative file path.
   * @param[in] progress Optional progress reporting callback.
   * @return Success message string or error description string.
   */
  std::expected<std::string, std::string> loadFile(const std::string &path,
                                                   const ProgressCallback &progress = nullptr);

  /**
   * @brief Gets a pointer to the loaded trajectory.
   * @return Const pointer to Trajectory, or nullptr if none loaded.
   */
  [[nodiscard]] const correlation::core::Trajectory *trajectory() const noexcept {
    return trajectory_.get();
  }

  /**
   * @brief Gets a mutable pointer to the loaded trajectory.
   * @return Pointer to Trajectory, or nullptr if none loaded.
   */
  [[nodiscard]] correlation::core::Trajectory *trajectoryMut() noexcept {
    return trajectory_.get();
  }

  /**
   * @brief Gets a pointer to the first cell frame.
   * @return Const pointer to Cell, or nullptr if none loaded.
   */
  [[nodiscard]] const correlation::core::Cell *cell() const noexcept {
    if (trajectory_ && trajectory_->getFrameCount() > 0) {
      return &trajectory_->firstFrame();
    }
    return nullptr;
  }

  /**
   * @brief Computes atomic species counts for the active cell.
   * @return Map of element symbol to count.
   */
  [[nodiscard]] std::map<std::string, int> getAtomCounts() const;

  /**
   * @brief Total number of frames in the loaded trajectory.
   */
  [[nodiscard]] size_t getFrameCount() const noexcept {
    return trajectory_ ? trajectory_->getFrameCount() : 0;
  }

  /**
   * @brief Total number of atoms in the primary frame.
   */
  [[nodiscard]] size_t getTotalAtomCount() const noexcept {
    if (!trajectory_ || trajectory_->getFrameCount() == 0) {
      return 0;
    }
    return trajectory_->firstFrame().atomCount();
  }

  /**
   * @brief Number of skipped/removed frames during parsing.
   */
  [[nodiscard]] size_t getRemovedFrameCount() const noexcept {
    return trajectory_ ? trajectory_->getRemovedFrameCount() : 0;
  }

  /**
   * @brief Simulation timestep of the trajectory in femtoseconds.
   */
  [[nodiscard]] real_t getTimeStep() const noexcept {
    return trajectory_ ? trajectory_->getTimeStep() : 1.0;
  }

  /**
   * @brief Recommended timestep derived from atomic masses in the active cell.
   */
  [[nodiscard]] real_t getRecommendedTimeStep() const;

  /**
   * @brief Injects an existing trajectory instance.
   * @param[in] traj Trajectory instance to own.
   */
  void setTrajectory(std::unique_ptr<correlation::core::Trajectory> traj) noexcept {
    trajectory_ = std::move(traj);
  }

  /**
   * @brief Clears currently loaded trajectory and resets state.
   */
  void clear() noexcept { trajectory_.reset(); }

private:
  std::unique_ptr<correlation::core::Trajectory> trajectory_;
};

} // namespace correlation::app
