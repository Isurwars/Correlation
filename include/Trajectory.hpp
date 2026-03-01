// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <vector>

#include "../include/Cell.hpp"

class Trajectory {
  /**
   * @brief This class stores a series of snapshots of a system.
   *
   * It ensures that the number of atoms and their identities are consistent
   * across all frames. It provides methods to analyze dynamic properties
   * of the system, such as velocity distributions.
   */

public:
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//
  Trajectory();
  Trajectory(std::vector<Cell> frames, double time_step);
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
  [[nodiscard]] const std::vector<std::vector<linalg::Vector3<double>>> &
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

  [[nodiscard]] size_t getRemovedFrameCount() const {
    return removed_frames_count_;
  }

private:
  //-------------------------------------------------------------------------//
  //--------------------------- Private Methods -----------------------------//
  //-------------------------------------------------------------------------//
  void validateFrame(const Cell &new_frame) const;
  std::vector<Cell> frames_;
  mutable std::vector<std::vector<double>> bond_cutoffs_sq_;
  // Stores calculate velocities Vector<Atom<Vector<Velocity>>>
  std::vector<std::vector<linalg::Vector3<double>>> velocities_;
  double time_step_;
  size_t removed_frames_count_{0};
};
