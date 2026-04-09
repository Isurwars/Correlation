/**
 * @file CifReader.hpp
 * @brief Reader for Crystallographic Information File (CIF) format.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */
#pragma once

#include <string>

#include "BaseReader.hpp"

namespace correlation::readers {

/**
 * @brief Reads a Crystallographic Information File (.cif).
 */
class CifReader : public BaseReader {
public:
  std::string getName() const override { return "CIF"; }
  std::vector<std::string> getExtensions() const override { return {"cif"}; }
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
