/**
 * @file GromacsReader.hpp
 * @brief Reader for GROMACS (.gro) files.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */
#pragma once

#include "BaseReader.hpp"

namespace correlation::readers {

/**
 * @class GromacsReader
 * @brief Handles parsing of GROMACS (.gro) structure files.
 */
class GromacsReader : public BaseReader {
public:
    std::string getName() const override { return "GROMACS Reader"; }
    std::vector<std::string> getExtensions() const override { return {".gro"}; }
    bool isTrajectory() const override { return false; }

    correlation::core::Cell readStructure(
        const std::string& filename,
        std::function<void(float, const std::string&)> progress_callback = nullptr) override;

    correlation::core::Trajectory readTrajectory(
        const std::string& filename,
        std::function<void(float, const std::string&)> progress_callback = nullptr) override;
};

} // namespace correlation::readers
