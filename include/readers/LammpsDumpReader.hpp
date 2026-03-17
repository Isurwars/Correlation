// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#pragma once

#include <string>

#include "BaseReader.hpp"

namespace FileReader {

/**
 * @brief Reads a LAMMPS atomic dump trajectory format.
 */
class LammpsDumpReader : public BaseReader {
public:
  std::string getName() const override { return "LAMMPS Dump"; }
  std::vector<std::string> getExtensions() const override { return {"dump", "lammpstrj"}; }
  bool isTrajectory() const override { return false; } // Current impl only reads one frame? Or should it be traj?

  Cell readStructure(const std::string &filename,
                     std::function<void(float, const std::string &)>
                         progress_callback = nullptr) override;

  Trajectory readTrajectory(const std::string &filename,
                            std::function<void(float, const std::string &)>
                                progress_callback = nullptr) override;

  static Cell read(const std::string &file_name);
};

} // namespace FileReader
