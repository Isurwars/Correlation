/**
 * @file MaceReader.hpp
 * @brief Reader for MACE (Machine Learning Interatomic Potential) Extended XYZ trajectory files.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */
#pragma once

#include "readers/BaseReader.hpp"

#include <optional>

namespace correlation::readers {

/**
 * @class MaceReader
 * @brief Reads atomic structures and trajectories from MACE Extended XYZ outputs.
 *
 * MACE outputs follow the Extended XYZ convention, writing simulation frames with
 * `Lattice="..."`, energy fields (`mace_energy`, `energy`, `free_energy`), and
 * per-atom properties in `Properties=...` (such as `species:S:1:pos:R:3:mace_forces:R:3`).
 *
 * When atomic forces (`mace_forces`, `forces`, `MACE_forces`) are present, they are
 * captured and stored in the atom's velocity vector for downstream structural and
 * dynamical analysis.
 */
class MaceReader : public BaseReader {
public:
  [[nodiscard]] std::string getName() const override { return "MACE Reader"; }
  [[nodiscard]] std::vector<std::string> getExtensions() const override { return {".mace", ".xyz", ".extxyz"}; }
  [[nodiscard]] bool isTrajectory() const override { return true; }

  correlation::core::Cell
  readStructure(const std::string &filename,
                std::function<void(float, const std::string &)> progress_callback = nullptr) override;

  correlation::core::Trajectory
  readTrajectory(const std::string &filename,
                 std::function<void(float, const std::string &)> progress_callback = nullptr) override;

  /**
   * @struct CommentData
   * @brief Intermediate storage for parsed MACE Extended XYZ comment line metadata and column layout.
   */
  struct CommentData {
    std::optional<std::array<real_t, 9>> lattice; ///< Optional 3x3 lattice vector matrix (Angstroms).
    std::optional<real_t> energy;                 ///< Optional total potential energy.
    int species_col = 0;                          ///< 0-based column index for chemical species.
    int pos_x_col = 1;                            ///< 0-based column index for X coordinate.
    int pos_y_col = 2;                            ///< 0-based column index for Y coordinate.
    int pos_z_col = 3;                            ///< 0-based column index for Z coordinate.
    int force_x_col = -1;                         ///< 0-based column index for X force (or -1 if absent).
    int force_y_col = -1;                         ///< 0-based column index for Y force (or -1 if absent).
    int force_z_col = -1;                         ///< 0-based column index for Z force (or -1 if absent).
  };

private:
  static void parseLattice(const std::string &comment, CommentData &data);
  static void parseEnergy(const std::string &comment, CommentData &data);
  static void parseProperties(const std::string &comment, CommentData &data);
  static void parsePropertiesParts(const std::vector<std::string> &parts, CommentData &data);
  static CommentData parseCommentLine(const std::string &comment);
  static correlation::core::Cell parseMaceFrame(const char *data, size_t size);
};

} // namespace correlation::readers
