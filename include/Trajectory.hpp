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
  
  /**
   * @brief get a bond cut off for two given elements.
   *@param type1
   *@param type2
   **/
  double getBondCutoff(int type1, int type2);

  void setBondCutoffs(const std::vector<std::vector<double>> &cutoffs) {
    bond_cutoffs_sq_ = cutoffs;
  }
  
  /**
   * @brief Pre-calculates the squared bond cutoff distances for every pair of
   * element types using the elements from the first frame.
   */
  void precomputeBondCutoffs();
  
  const std::vector<std::vector<double>>& getBondCutoffs() const { return bond_cutoffs_sq_; }

private:
  void validateFrame(const Cell &new_frame) const;
  std::vector<Cell> frames_;
  std::vector<std::vector<double>> bond_cutoffs_sq_;
  double time_step_;
};
