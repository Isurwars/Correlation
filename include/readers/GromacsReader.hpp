/**
 * @file GromacsReader.hpp
 * @brief Reader for GROMACS (.gro) files with multi-frame lazy loading.
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
 *
 * Supports multi-frame `.gro` files via memory-mapped lazy loading.
 * Each frame consists of: title line, atom count, N atom lines, box line.
 */
class GromacsReader : public BaseReader {
public:
    std::string getName() const override { return "GROMACS Reader"; }
    std::vector<std::string> getExtensions() const override { return {".gro"}; }
    bool isTrajectory() const override { return true; }

    correlation::core::Cell readStructure(
        const std::string& filename,
        std::function<void(float, const std::string&)> progress_callback = nullptr) override;

    correlation::core::Trajectory readTrajectory(
        const std::string& filename,
        std::function<void(float, const std::string&)> progress_callback = nullptr) override;

    /**
     * @brief Parses a single GROMACS .gro frame from a memory region.
     *
     * @param data Pointer to the start of the frame data.
     * @param size Number of bytes in the frame region.
     * @return A parsed Cell object.
     */
    static correlation::core::Cell parseGroFrame(const char *data, size_t size);
};

} // namespace correlation::readers
