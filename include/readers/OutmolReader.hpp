/**
 * @file OutmolReader.hpp
 * @brief Reader for CASTEP .outmol trajectory files.
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
 * @brief Reads a DMol3 .outmol output format.
 */
class OutmolReader : public BaseReader {
public:
  std::string getName() const override { return "Outmol"; }
  std::vector<std::string> getExtensions() const override { return {"outmol"}; }
  bool isTrajectory() const override { return true; }

  correlation::core::Cell readStructure(
      const std::string &filename,
      std::function<void(float, const std::string &)> progress_callback =
          nullptr) override;

  correlation::core::Trajectory readTrajectory(
      const std::string &filename,
      std::function<void(float, const std::string &)> progress_callback =
          nullptr) override;

  static std::vector<correlation::core::Cell>
  read(const std::string &file_name,
       std::function<void(float, const std::string &)> progress_callback =
           nullptr);
};

} // namespace correlation::readers
