/**
 * @file XdatcarReader.hpp
 * @brief Reader for VASP XDATCAR trajectory files.
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
 * @brief Reads VASP XDATCAR trajectory files.
 *
 * XDATCAR files contain a shared header (lattice, species, counts)
 * followed by multiple frames of atomic positions in Direct (fractional)
 * coordinates, each preceded by a "Direct configuration= N" separator.
 */
class XdatcarReader : public BaseReader {
public:
  std::string getName() const override { return "VASP XDATCAR"; }
  std::vector<std::string> getExtensions() const override {
    return {"xdatcar"};
  }
  bool isTrajectory() const override { return true; }

  correlation::core::Cell readStructure(
      const std::string &filename,
      std::function<void(float, const std::string &)> progress_callback =
          nullptr) override;

  correlation::core::Trajectory readTrajectory(
      const std::string &filename,
      std::function<void(float, const std::string &)> progress_callback =
          nullptr) override;

  /**
   * @brief Low-level parser for XDATCAR files.
   * @param file_name Path to the source file.
   * @param progress_callback Optional callback for progress updates.
   * @return A vector of parsed frames.
   */
  static std::vector<correlation::core::Cell>
  read(const std::string &file_name,
       std::function<void(float, const std::string &)> progress_callback =
           nullptr);
};

} // namespace correlation::readers
