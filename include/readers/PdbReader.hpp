/**
 * @file PdbReader.hpp
 * @brief Reader for Protein Data Bank (.pdb) files.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */
#pragma once

#include "BaseReader.hpp"

namespace correlation::readers {

/**
 * @class PdbReader
 * @brief Handles parsing of PDB (.pdb) structure files.
 */
class PdbReader : public BaseReader {
public:
    std::string getName() const override { return "PDB Reader"; }
    std::vector<std::string> getExtensions() const override { return {".pdb", ".ent"}; }
    bool isTrajectory() const override { return true; }

    correlation::core::Cell readStructure(
        const std::string& filename,
        std::function<void(float, const std::string&)> progress_callback = nullptr) override;

    correlation::core::Trajectory readTrajectory(
        const std::string& filename,
        std::function<void(float, const std::string&)> progress_callback = nullptr) override;
};

} // namespace correlation::readers
