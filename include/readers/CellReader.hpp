/**
 * @file CellReader.hpp
 * @brief Reader for CASTEP CELL structure files.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */
#pragma once

#include <string>

#include "BaseReader.hpp"

namespace correlation::readers {

/**
 * @brief Reads a CASTEP .cell physical structure format.
 */
class CellReader : public BaseReader {
public:
  std::string getName() const override { return "CASTEP CELL"; }
  std::vector<std::string> getExtensions() const override { return {"cell"}; }
  bool isTrajectory() const override { return false; }

  Cell readStructure(const std::string &filename,
                     std::function<void(float, const std::string &)>
                         progress_callback = nullptr) override;

  Trajectory readTrajectory(const std::string &filename,
                            std::function<void(float, const std::string &)>
                                progress_callback = nullptr) override;

  static Cell read(const std::string &file_name);
};

} // namespace correlation::readers
