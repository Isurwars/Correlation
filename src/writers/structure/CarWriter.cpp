/**
 * @file CarWriter.cpp
 * @brief Implementation of the Accelrys CAR structure file writer.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "writers/structure/CarWriter.hpp"
#include "writers/structure/StructureWriterFactory.hpp"

#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <string>

namespace correlation::writers::structure {

// Static self-registration in factory
const bool CAR_REGISTERED = StructureWriterFactory::registerTypeSafe<CarWriter>("CarWriter");

namespace {

void validateCell(const correlation::core::Cell &cell, const std::string &filepath) {
  if (cell.volume() <= 1e-9) {
    throw std::runtime_error("Cannot write CAR to '" + filepath +
                             "': cell volume is degenerate or non-periodic.");
  }
  if (cell.atoms().empty()) {
    throw std::runtime_error("Cannot write CAR to '" + filepath + "': cell contains no atoms.");
  }
}

void writePbcHeader(std::ostream &out, const correlation::core::Cell &cell) {
  out << "!BIOSYM archive 3\n";
  out << "PBC=ON\n";

  const auto &params = cell.latticeParameters();
  out << std::fixed << std::setprecision(4);
  out << "PBC " << std::setw(9) << params[0] << " " << std::setw(9) << params[1] << " "
      << std::setw(9) << params[2] << " " << std::setw(9) << params[3] << " " << std::setw(9)
      << params[4] << " " << std::setw(9) << params[5] << "\n";
}

void writeAtomRecords(std::ostream &out, const correlation::core::Cell &cell) {
  std::map<std::string, std::size_t> element_counts;
  out << std::fixed << std::setprecision(4);

  for (std::size_t idx = 0; idx < cell.atoms().size(); ++idx) {
    const auto &atom = cell.atoms()[idx];
    const std::string &sym = atom.element().symbol;
    if (sym.empty()) {
      throw std::runtime_error("Cannot write CAR: atom index " + std::to_string(idx) +
                               " has an empty chemical symbol.");
    }

    element_counts[sym]++;
    std::string const label = sym + std::to_string(element_counts[sym]);
    const auto &pos = atom.position();

    out << std::left << std::setw(6) << label << " " << std::right << std::setw(9) << pos[0] << " "
        << std::setw(9) << pos[1] << " " << std::setw(9) << pos[2] << " " << std::left
        << "XXXX 1    xx    " << std::setw(4) << sym << " 0.000\n";
  }

  out << "end\nend\n";
}

} // namespace

void CarWriter::write(const std::string &filepath, const correlation::core::Cell &cell,
                      const StructureWriterOptions & /*options*/) const {
  validateCell(cell, filepath);

  std::ofstream out(filepath);
  if (!out.is_open()) {
    throw std::runtime_error("Unable to open destination file for writing: " + filepath);
  }

  writePbcHeader(out, cell);
  writeAtomRecords(out, cell);
}

} // namespace correlation::writers::structure
