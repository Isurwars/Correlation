/**
 * @file XYZReader.hpp
 * @brief Reader for Extended XYZ (.xyz, .exyz) structure files.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
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

  correlation::core::Cell
  readStructure(const std::string &filename,
                std::function<void(float, const std::string &)> progress_callback = nullptr) override;

  correlation::core::Trajectory
  readTrajectory(const std::string &filename,
                 std::function<void(float, const std::string &)> progress_callback = nullptr) override;

private:
  struct CommentData {
    std::optional<std::array<double, 9>> lattice;
    std::optional<double> energy;
    int species_col = 0;
    int pos_x_col = 1;
    int pos_y_col = 2;
    int pos_z_col = 3;
  };

  static void parseLattice(const std::string &comment, CommentData &data);
  static void parseEnergy(const std::string &comment, CommentData &data);
  static void parseProperties(const std::string &comment, CommentData &data);
  static void parsePropertiesParts(const std::vector<std::string> &parts, CommentData &data);

  /**
   * @brief Parses an Extended XYZ comment line for properties, lattice and energy.
   *
   * @param comment The comment line string.
   * @return CommentData struct containing parsed info.
   */
  static CommentData parseCommentLine(const std::string &comment);
  static correlation::core::Cell parseXYZFrame(const char *data, size_t size);
};

} // namespace correlation::readers
