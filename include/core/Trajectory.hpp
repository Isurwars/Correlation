/**
 * @file Trajectory.hpp
 * @brief Trajectory data structure for time-series atomistic simulations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "core/Cell.hpp"
#include "math/LinearAlgebra.hpp"

#include <vector>

namespace correlation::core {

/**
 * @brief This class stores a series of snapshots of a system.
 *
 * It ensures that the number of atoms and their identities are consistent
 * across all frames. It provides methods to analyze dynamic properties
 * of the system, such as velocity distributions.
 */
class Trajectory {

public:
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//
  Trajectory();

  /**
   * @brief Constructs a trajectory from a vector of frames.
   * @param frames Vector of Cell objects representing simulation snapshots.
   * @param time_step The time interval between consecutive frames.
   */
  Trajectory(std::vector<Cell> frames, double time_step);

  /**
   * @brief Appends a new frame to the trajectory.
   * @param frame The Cell object to add.
   */
  void addFrame(const Cell &frame);

  //-------------------------------------------------------------------------//
  //------------------------------- Accessors -------------------------------//
  //-------------------------------------------------------------------------//
  /**
   * @brief Gets a mutable reference to the frames.
   * @return A vector of Cell objects representing the frames.
   */
  [[nodiscard]] std::vector<Cell> &getFrames() { return frames_; }

  /**
   * @brief Gets a constant reference to the frames.
   * @return A vector of Cell objects representing the frames.
   */
  [[nodiscard]] const std::vector<Cell> &getFrames() const { return frames_; }

  /**
   * @brief Gets the total number of frames in the trajectory.
   * @return The number of frames.
   */
  [[nodiscard]] size_t getFrameCount() const { return frames_.size(); }

  /**
   * @brief Gets the time step between frames.
   * @return The time step in simulation units.
   */
  [[nodiscard]] double getTimeStep() const { return time_step_; }

  /**
   * @brief Sets the time step between frames.
   * @param dt The new time step value.
   */
  void setTimeStep(double dt) { time_step_ = dt; }

  /**
   * @brief Gets a linear bond cutoff distance for two given element types.
   * @param type1 ID of the first element type.
   * @param type2 ID of the second element type.
   * @return The linear bond cutoff distance.
   */
  [[nodiscard]] double getBondCutoff(int type1, int type2) const;

  /**
   * @brief Gets a squared bond cutoff distance for two given element types.
   * @param type1 ID of the first element type.
   * @param type2 ID of the second element type.
   * @return The squared bond cutoff distance.
   */
  [[nodiscard]] double getBondCutoffSQ(int type1, int type2) const;

  /**
   * @brief Sets the squared bond cutoffs for neighbor searching.
   * @param cutoffs Matrix of squared distance cutoffs [type1][type2].
   */
  void setBondCutoffsSQ(const std::vector<std::vector<double>> &cutoffs) {
    bond_cutoffs_sq_ = cutoffs;
  }

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
   * @brief Gets the pre-calculated velocities for all atoms.
   * @return A vector of velocity vectors, indexed by `[frame][atom_index]`.
   */
  [[nodiscard]] const std::vector<std::vector<math::Vector3<double>>> &
  getVelocities() const {
    return velocities_;
  }

  /**
   * @brief Gets the matrix of squared bond cutoffs for element pairs.
   * @return The squared bond cutoff matrix.
   */
  [[nodiscard]] const std::vector<std::vector<double>> &
  getBondCutoffsSQ() const {
    return bond_cutoffs_sq_;
  }

  /**
   * @brief Returns the number of frames removed during deduplication.
   * @return Count of duplicate frames found and removed.
   */
  [[nodiscard]] size_t getRemovedFrameCount() const {
    return removed_frames_count_;
  }

private:
  //-------------------------------------------------------------------------//
  //--------------------------- Private Methods -----------------------------//
  //-------------------------------------------------------------------------//
  /**
   * @brief Internal helper to validate that a new frame matches the trajectory
   * topology.
   * @param new_frame The frame to check.
   * @throws std::invalid_argument if frame is inconsistent.
   */
  void validateFrame(const Cell &new_frame) const;

  std::vector<Cell> frames_; ///< Collection of simulation snapshots.
  mutable std::vector<std::vector<double>> bond_cutoffs_sq_; ///< Cached squared bond cutoffs.
  std::vector<std::vector<math::Vector3<double>>> velocities_; ///< Inter-frame velocities.
  double time_step_; ///< Time between snapshots.
  size_t removed_frames_count_{0}; ///< Counter for deduplicated frames.
};

} // namespace correlation::core
