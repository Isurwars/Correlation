/**
 * @file CastepMdReader.hpp
 * @brief Reader for CASTEP molecular dynamics trajectory files.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
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
 * @brief Reads a sequence of frames from a CASTEP .md trajectory file.
 */
class CastepMdReader : public BaseReader {
public:
  std::string getName() const override { return "CASTEP MD"; }
  std::vector<std::string> getExtensions() const override { return {"md"}; }
  bool isTrajectory() const override { return true; }

  correlation::core::Cell
  readStructure(const std::string &filename,
                std::function<void(float, const std::string &)> progress_callback = nullptr) override;

  correlation::core::Trajectory
  readTrajectory(const std::string &filename,
                 std::function<void(float, const std::string &)> progress_callback = nullptr) override;

  /**
   * @brief Low-level parser for the specific file format.
   * @param file_name Path to the source file.
   * @param progress_callback Optional callback for progress updates.
   * @return A vector of parsed frames.
   */
  static std::vector<correlation::core::Cell>
  read(const std::string &file_name,
       const std::function<void(float, const std::string &)> &progress_callback = nullptr);

private:
  static void updateProgress(std::streampos current_pos, std::streampos file_size, std::streampos &last_progress_pos,
                             size_t update_interval,
                             const std::function<void(float, const std::string &)> &progress_callback);

  static void parseEnergyLine(const std::string &line, real_t &current_energy, correlation::core::Cell &tempCell,
                              bool &cell_has_atoms, std::vector<correlation::core::Cell> &frames);

  static void parseLatticeLine(std::ifstream &myfile, const std::string &line, correlation::core::Cell &tempCell);

  static void parseAtomLine(const std::string &line, real_t current_energy, correlation::core::Cell &tempCell,
                            bool &cell_has_atoms);
};

} // namespace correlation::readers
