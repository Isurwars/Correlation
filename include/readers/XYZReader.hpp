/**
 * @file XYZReader.hpp
 * @brief Reader for Extended XYZ (.xyz, .exyz) structure files.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */
#pragma once

#include "BaseReader.hpp"

namespace correlation::readers {

/**
 * @class XYZReader
 * @brief Reads structure data from Extended XYZ format files.
 *
 * Supports the standard XYZ format (atom count, comment line, coordinate
 * lines) as well as the Extended XYZ convention where lattice vectors are
 * encoded in the comment line via `Lattice="..."`.
 *
 * Only structural data (species, positions, lattice) is parsed.
 * Custom per-atom properties and arbitrary key-value metadata are ignored.
 */
class XYZReader : public BaseReader {
public:
    std::string getName() const override { return "XYZ Reader"; }
    std::vector<std::string> getExtensions() const override { return {".xyz", ".exyz"}; }
    bool isTrajectory() const override { return true; }

    correlation::core::Cell readStructure(
        const std::string& filename,
        std::function<void(float, const std::string&)> progress_callback = nullptr) override;

    correlation::core::Trajectory readTrajectory(
        const std::string& filename,
        std::function<void(float, const std::string&)> progress_callback = nullptr) override;

private:
    /**
     * @brief Parses lattice vectors from an Extended XYZ comment line.
     *
     * Looks for `Lattice="v1x v1y v1z v2x v2y v2z v3x v3y v3z"` and
     * returns the 3×3 lattice matrix. If no Lattice key is found, returns
     * std::nullopt.
     *
     * @param comment The comment line string.
     * @return An optional array of 9 doubles (row-major: a1,a2,a3, b1,b2,b3, c1,c2,c3).
     */
    static std::optional<std::array<double, 9>> parseLattice(const std::string& comment);
    static correlation::core::Cell parseXYZFrame(const char *data, size_t size);
};

} // namespace correlation::readers
