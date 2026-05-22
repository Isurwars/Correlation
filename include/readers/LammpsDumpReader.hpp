/**
 * @file LammpsDumpReader.hpp
 * @brief Reader for LAMMPS dump trajectory files.
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
 * @brief Reads one or more frames from a LAMMPS custom dump file.
 *
 * Supports the common `ITEM: ATOMS id type x y z` format as well as the
 * element-name variant (`element x y z`) and scaled / wrapped coordinates
 * (`xs ys zs`). When the file contains more than one timestep the reader
 * returns every frame, making it suitable for LAMMPS MD trajectories.
 *
 * Uses memory-mapped I/O for lazy frame loading on large trajectories.
 */
class LammpsDumpReader : public BaseReader {
public:
  std::string getName() const override { return "LAMMPS Dump"; }
  std::vector<std::string> getExtensions() const override {
    return {"dump", "lammpstrj"};
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
   * @brief Parses a single LAMMPS dump frame from a memory region.
   *
   * The region must contain a complete frame starting with "ITEM: TIMESTEP".
   *
   * @param data Pointer to the start of the frame data.
   * @param size Number of bytes in the frame region.
   * @return A parsed Cell object.
   */
  static correlation::core::Cell parseDumpFrame(const char *data, size_t size);
};

} // namespace correlation::readers
