/**
 * @file XdatcarReader.hpp
 * @brief Reader for VASP XDATCAR trajectory files with lazy loading.
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
 * @brief Reads VASP XDATCAR trajectory files with memory-mapped lazy loading.
 *
 * XDATCAR files contain a shared header (lattice, species, counts)
 * followed by multiple frames of atomic positions in Direct (fractional)
 * coordinates, each preceded by a "Direct configuration= N" separator.
 */
class XdatcarReader : public BaseReader {
public:
  std::string getName() const override { return "VASP XDATCAR"; }
  std::vector<std::string> getExtensions() const override { return {"xdatcar"}; }
  bool isTrajectory() const override { return true; }

  correlation::core::Cell
  readStructure(const std::string &filename,
                std::function<void(float, const std::string &)> progress_callback = nullptr) override;

  correlation::core::Trajectory
  readTrajectory(const std::string &filename,
                 std::function<void(float, const std::string &)> progress_callback = nullptr) override;
};

} // namespace correlation::readers
