/**
 * @file CarReader.hpp
 * @brief Reader for BIOVIA/Accelrys CAR structure files.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */
#pragma once

#include "BaseReader.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include <functional>
#include <string>
#include <vector>

namespace correlation::readers {

/**
 * @brief Reads a Materials Studio .car file and returns a
 * correlation::core::Cell object.
 */
class CarReader : public BaseReader {
public:
  std::string getName() const override { return "Accelrys CAR"; }
  std::vector<std::string> getExtensions() const override { return {"car"}; }
  bool isTrajectory() const override { return false; }

  correlation::core::Cell readStructure(
      const std::string &filename,
      std::function<void(float, const std::string &)> progress_callback =
          nullptr) override;

  correlation::core::Trajectory readTrajectory(
      const std::string &filename,
      std::function<void(float, const std::string &)> progress_callback =
          nullptr) override;

  /**
   * @brief Low-level parser for the specific file format.
   * @param file_name Path to the source file.
   * @return The parsed simulation cell.
   */
  static correlation::core::Cell read(const std::string &file_name);
};

} // namespace correlation::readers
