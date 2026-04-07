/**
 * @file OutmolReader.hpp
 * @brief Reader for CASTEP .outmol trajectory files.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */
#pragma once

#include <functional>
#include <string>
#include <vector>

#include "Cell.hpp"
#include "Trajectory.hpp"

#include "BaseReader.hpp"

namespace FileReader {

/**
 * @brief Reads a DMol3 .outmol output format.
 */
class OutmolReader : public BaseReader {
public:
  std::string getName() const override { return "Outmol"; }
  std::vector<std::string> getExtensions() const override { return {"outmol"}; }
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
