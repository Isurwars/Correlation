// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#pragma once

#include <functional>
#include <string>
#include <vector>

#include "../Cell.hpp"
#include "../Trajectory.hpp"

#include "BaseReader.hpp"

namespace FileReader {

/**
 * @brief Reads a Materials Studio .car file and returns a Cell object.
 */
class CarReader : public BaseReader {
public:
  std::string getName() const override { return "Accelrys CAR"; }
  std::vector<std::string> getExtensions() const override { return {"car"}; }
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
