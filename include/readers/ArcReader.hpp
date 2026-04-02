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
 * @brief Reads an ARC file and returns a vector of Cell objects.
 */
class ArcReader : public BaseReader {
public:
  std::string getName() const override { return "Accelrys ARC"; }
  std::vector<std::string> getExtensions() const override { return {"arc"}; }
  bool isTrajectory() const override { return true; }

  Cell readStructure(const std::string &filename,
                     std::function<void(float, const std::string &)>
                         progress_callback = nullptr) override;

  Trajectory readTrajectory(const std::string &filename,
                            std::function<void(float, const std::string &)>
                                progress_callback = nullptr) override;

  static std::vector<Cell> read(const std::string &file_name,
                                std::function<void(float, const std::string &)>
                                    progress_callback = nullptr);
};

} // namespace FileReader
