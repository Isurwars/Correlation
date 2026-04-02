// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#pragma once

#include <functional>
#include <string>
#include <vector>

#include "Cell.hpp"
#include "Trajectory.hpp"
#include "BaseReader.hpp"

namespace FileReader {

/**
 * @brief Reads one or more frames from a LAMMPS custom dump file.
 *
 * Supports the common `ITEM: ATOMS id type x y z` format as well as the
 * element-name variant (`element x y z`) and scaled / wrapped coordinates
 * (`xs ys zs`). When the file contains more than one timestep the reader
 * returns every frame, making it suitable for LAMMPS MD trajectories.
 */
class LammpsDumpReader : public BaseReader {
public:
  std::string getName() const override { return "LAMMPS Dump"; }
  std::vector<std::string> getExtensions() const override {
    return {"dump", "lammpstrj"};
  }
  bool isTrajectory() const override { return true; }

  Cell readStructure(const std::string &filename,
                     std::function<void(float, const std::string &)>
                         progress_callback = nullptr) override;

  Trajectory readTrajectory(const std::string &filename,
                            std::function<void(float, const std::string &)>
                                progress_callback = nullptr) override;

  /**
   * @brief Reads all frames from a LAMMPS dump file.
   *
   * @param file_name        Path to the dump file.
   * @param progress_callback Optional callback invoked with (fraction, message).
   * @return A vector of Cell objects, one per timestep found in the file.
   */
  static std::vector<Cell>
  read(const std::string &file_name,
       std::function<void(float, const std::string &)> progress_callback =
           nullptr);
};

} // namespace FileReader
