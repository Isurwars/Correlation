// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

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
  std::vector<Cell> &getFrames() { return frames_; }
  const std::vector<Cell> &getFrames() const { return frames_; }
  size_t getFrameCount() const { return frames_.size(); }
  double getTimeStep() const { return time_step_; }
};
