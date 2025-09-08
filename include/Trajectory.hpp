#ifndef INCLUDE_TRAJECTORY_HPP_
#define INCLUDE_TRAJECTORY_HPP_
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include <vector>

#include "../include/Cell.hpp"

class Trajectory {
  /* --------------------------------------------------------------------------
   * This class stores a series of snapshots of a system, ensuring that the
   * number of atoms and their identities are consistent across all frames. It
   * provides methods to analyze the dynamic properties of the system, such as
   * velocity distributions.
   * --------------------------------------------------------------------------
   */

private:
  void validateFrame(const Cell &new_frame) const;
  std::vector<Cell> frames_;
  double time_step_; // Time step between frames

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
  const std::vector<Cell> &getFrames() const { return frames_; }
  size_t getFrameCount() const { return frames_.size(); }
  double getTimeStep() const { return time_step_; }
};

#endif // INCLUDE_TRAJECTORY_HPP_
