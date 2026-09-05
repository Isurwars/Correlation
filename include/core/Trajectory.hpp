/**
 * @file Trajectory.hpp
 * @brief Trajectory data structure for time-series atomistic simulations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/AnalysisTypes.hpp"
#include "core/Cell.hpp"

#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <vector>

namespace correlation::core {

class MappedFile;

/**
 * @brief This class stores a series of snapshots of a system.
 *
 * It ensures that the number of atoms and their identities are consistent
 * across all frames. It provides methods to analyze dynamic properties
 * of the system, such as velocity distributions.
 */
class Trajectory {

public:
  /** @name Constructors */
  ///@{
  Trajectory();

  /**
   * @brief Constructs a trajectory from a vector of frames.
   * @param frames Vector of Cell objects representing simulation snapshots.
   * @param time_step The time interval between consecutive frames.
   */
  Trajectory(std::vector<Cell> frames, real_t time_step);

  using FrameParser = std::function<Cell(const char *data, size_t size)>;

  /**
   * @brief Constructs a lazy-loading trajectory from a memory-mapped file.
   * @param mapped_file Shared pointer to the mapped file.
   * @param frame_offsets Byte offsets of the frames in the file.
   * @param parser Function to parse a frame from a slice of memory.
   * @param time_step The time interval between consecutive frames.
   */
  Trajectory(std::shared_ptr<MappedFile> mapped_file, std::vector<size_t> frame_offsets, FrameParser parser,
             real_t time_step);

  /**
   * @brief Appends a new frame to the trajectory.
   * @param frame The Cell object to add.
   */
  void addFrame(const Cell &frame);

  ///@}

  /** @name Accessors */
  ///@{
  /**
   * @brief Gets a mutable reference to the frames.
   * @return A vector of Cell objects representing the frames.
   */
  [[nodiscard]] std::vector<Cell> &getFrames();

  /**
   * @brief Gets a constant reference to the frames.
   * @return A vector of Cell objects representing the frames.
   */
  [[nodiscard]] const std::vector<Cell> &getFrames() const;

  /**
   * @brief Gets the total number of frames in the trajectory.
   * @return The number of frames.
   */
  [[nodiscard]] size_t getFrameCount() const noexcept;

  /**
   * @brief Gets a constant reference to the first frame of the trajectory.
   * @return A Cell object representing the first frame.
   * @throws std::runtime_error if the trajectory is empty.
   */
  [[nodiscard]] const Cell &firstFrame() const;

  /**
   * @brief Gets a copy of a frame at a specific index.
   * @param index The index of the frame to retrieve.
   * @return A Cell object representing the frame.
   * @throws std::out_of_range if the index is invalid.
   */
  [[nodiscard]] Cell getFrame(size_t index) const;

  /**
   * @brief Gets the time step between frames.
   * @return The time step in simulation units.
   */
  [[nodiscard]] real_t getTimeStep() const noexcept { return time_step_; }

  /**
   * @brief Sets the simulation time step.
   * @param step Time step in femtoseconds.
   */
  void setTimeStep(real_t step) noexcept { time_step_ = step; }

  /**
   * @brief Gets the bond cutoff range for two given element types.
   * @param type1 ID of the first element type.
   * @param type2 ID of the second element type.
   * @return Const reference to the BondCutoffRange.
   */
  [[nodiscard]] const correlation::analysis::BondCutoffRange &getBondCutoffRange(size_t type1, size_t type2) const;

  /**
   * @brief Gets a squared bond cutoff distance for two given element types.
   * @param type1 ID of the first element type.
   * @param type2 ID of the second element type.
   * @return The squared bond cutoff distance.
   */
  [[nodiscard]] real_t getBondCutoffSQ(size_t type1, size_t type2) const;

  /**
   * @brief Gets a squared minimum bond cutoff distance for two given element types.
   * @param type1 ID of the first element type.
   * @param type2 ID of the second element type.
   * @return The squared minimum bond cutoff distance.
   */
  [[nodiscard]] real_t getMinBondCutoffSQ(size_t type1, size_t type2) const;

  /**
   * @brief Gets the maximum bond cutoff distance for two given element types.
   * @param type1 ID of the first element type.
   * @param type2 ID of the second element type.
   * @return The bond cutoff distance.
   */
  [[nodiscard]] real_t getBondCutoff(size_t type1, size_t type2) const;

  /**
   * @brief Gets the minimum bond cutoff distance for two given element types.
   * @param type1 ID of the first element type.
   * @param type2 ID of the second element type.
   * @return The minimum bond cutoff distance.
   */
  [[nodiscard]] real_t getMinBondCutoff(size_t type1, size_t type2) const;

  /**
   * @brief Sets the bond cutoff matrix for neighbor searching.
   * @param cutoffs Matrix of bond cutoff ranges [type1][type2].
   * @throws std::invalid_argument If any cutoff bound is negative or min_sq > max_sq.
   */
  void setBondCutoffs(const correlation::analysis::BondCutoffMatrix &cutoffs);

  /**
   * @brief Removes consecutive duplicated frames from the trajectory.
   * Frames are considered duplicates if all atom positions are identical
   * (within a small tolerance).
   */
  void removeDuplicatedFrames();

  /**
   * @brief Pre-calculates the squared bond cutoff distances for every pair of
   * element types using the elements from the first frame.
   */
  void precomputeBondCutoffs() const;

  /**
   * @brief Calculates velocities for all atoms in the trajectory using
   * finite differences, accounting for periodic boundary conditions.
   * Assumes constant time step.
   * Populates the internal velocities_ vector.
   */
  void calculateVelocities();

  /**
   * @brief Gets the matrix of bond cutoff ranges for element pairs.
   * @return The bond cutoff matrix.
   */
  [[nodiscard]] const correlation::analysis::BondCutoffMatrix &getBondCutoffs() const noexcept { return bond_cutoffs_; }

  /**
   * @brief Backward-compatible alias returning the bond cutoff matrix.
   * @return The bond cutoff matrix.
   */
  [[nodiscard]] const correlation::analysis::BondCutoffMatrix &getBondCutoffsSQ() const noexcept { return bond_cutoffs_; }

  /**
   * @brief Returns the number of frames removed during deduplication.
   * @return Count of duplicate frames found and removed.
   */
  [[nodiscard]] size_t getRemovedFrameCount() const noexcept { return removed_frames_count_; }

  ///@}

private:
  /**
   * @brief Internal helper to validate that a new frame matches the trajectory
   * topology.
   * @param new_frame The frame to check.
   * @throws std::invalid_argument if frame is inconsistent.
   */
  void validateFrame(const Cell &new_frame) const;
  void ensureMaterialized() const;

  mutable std::unique_ptr<std::mutex> init_mutex_{std::make_unique<std::mutex>()};

  mutable std::vector<Cell> frames_; ///< Collection of simulation snapshots.
  mutable std::optional<Cell> first_frame_;
  mutable correlation::analysis::BondCutoffMatrix bond_cutoffs_; ///< Cached bond cutoff ranges.
  real_t time_step_;                                             ///< Time between snapshots.
  size_t removed_frames_count_{0};                               ///< Counter for deduplicated frames.

  std::shared_ptr<MappedFile> mapped_file_;
  std::vector<size_t> frame_offsets_;
  FrameParser parser_;
};

} // namespace correlation::core
