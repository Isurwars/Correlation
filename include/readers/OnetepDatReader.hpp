/**
 * @file OnetepDatReader.hpp
 * @brief Reader for ONETEP .dat structure files.
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
 * @brief Reads a ONETEP .dat structure format.
 */
class OnetepDatReader : public BaseReader {
public:
  std::string getName() const override { return "ONETEP DAT"; }
  std::vector<std::string> getExtensions() const override { return {"dat"}; }
  bool isTrajectory() const override { return false; }

  Cell readStructure(const std::string &filename,
                     std::function<void(float, const std::string &)>
                         progress_callback = nullptr) override;

  Trajectory readTrajectory(const std::string &filename,
                            std::function<void(float, const std::string &)>
                                progress_callback = nullptr) override;

  static Cell read(const std::string &file_name);
};

} // namespace FileReader
