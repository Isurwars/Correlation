/**
 * @file PdbReader.cpp
 * @brief Implementation of the PDB (.pdb) file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/PdbReader.hpp"
#include "readers/ReaderFactory.hpp"

#include <fstream>
#include <optional>
#include <stdexcept>

namespace correlation::readers {

namespace {

// Automatic registration
const bool registered = ReaderFactory::registerTypeSafe<PdbReader>("PdbReader");

struct PdbCrystParams {
  real_t a = 0.0;
  real_t b = 0.0;
  real_t c = 0.0;
};

struct PdbAtomData {
  std::string symbol;
  real_t x = 0.0;
  real_t y = 0.0;
  real_t z = 0.0;
};

std::optional<PdbCrystParams> parsePdbCrystLine(const std::string &line) {
  try {
    real_t const param_a = static_cast<real_t>(std::stod(line.substr(6, 9)));
    real_t const param_b = static_cast<real_t>(std::stod(line.substr(15, 9)));
    real_t const param_c = static_cast<real_t>(std::stod(line.substr(24, 9)));
    return PdbCrystParams{
        .a = param_a,
        .b = param_b,
        .c = param_c,
    };
  } catch (...) {
    return std::nullopt;
  }
}

std::optional<PdbAtomData> parsePdbAtomLine(const std::string &line) {
  try {
    std::string symbol;
    if (line.length() >= 78) {
      symbol = line.substr(76, 2);
      symbol.erase(0, symbol.find_first_not_of(' '));
      symbol.erase(symbol.find_last_not_of(' ') + 1);
    }

    if (symbol.empty()) {
      // Fallback to atom name columns 12-16
      std::string atom_name = line.substr(12, 4);
      atom_name.erase(0, atom_name.find_first_not_of(' '));
      atom_name.erase(atom_name.find_last_not_of(' ') + 1);

      // Common PDB fixes (e.g. 1H -> H)
      std::string parsed_name = atom_name;
      if (parsed_name.length() > 1 && (std::isdigit(parsed_name[0]) != 0)) {
        parsed_name = parsed_name.substr(1);
      }

      for (char const chr : parsed_name) {
        if (std::isalpha(chr) != 0) {
          symbol += chr;
        } else {
          break;
        }
      }
    }

    real_t const frac_x = static_cast<real_t>(std::stod(line.substr(30, 8)));
    real_t const frac_y = static_cast<real_t>(std::stod(line.substr(38, 8)));
    real_t const frac_z = static_cast<real_t>(std::stod(line.substr(46, 8)));
    return PdbAtomData{
        .symbol = symbol,
        .x = frac_x,
        .y = frac_y,
        .z = frac_z,
    };
  } catch (...) {
    return std::nullopt;
  }
}

} // namespace

correlation::core::Cell PdbReader::readStructure(const std::string &filename,
                                                 std::function<void(float, const std::string &)> progress_callback) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  correlation::core::Cell cell;
  std::string line;
  bool has_box = false;

  while (std::getline(file, line)) {
    if (line.starts_with("CRYST1")) {
      if (auto params = parsePdbCrystLine(line)) {
        cell.setLatticeParameters({params->a, params->b, params->c, 90.0, 90.0, 90.0});
        has_box = true;
      }
    } else if (line.starts_with("ATOM") || line.starts_with("HETATM")) {
      if (auto atom = parsePdbAtomLine(line)) {
        cell.addAtom(atom->symbol, correlation::math::Vector3<real_t>(atom->x, atom->y, atom->z));
      }
    } else if (line.starts_with("ENDMDL") || line.starts_with("END")) {
      break;
    }
  }

  if (progress_callback) {
    progress_callback(1.0F, "PDB structure loaded.");
  }
  return cell;
}

correlation::core::Trajectory
PdbReader::readTrajectory(const std::string &filename,
                          std::function<void(float, const std::string &)> progress_callback) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  std::vector<correlation::core::Cell> frames;
  std::string line;
  correlation::core::Cell current_cell;
  bool in_model = false;
  bool has_box = false;
  real_t param_a = 0.0;
  real_t param_b = 0.0;
  real_t param_c = 0.0;

  while (std::getline(file, line)) {
    if (line.starts_with("CRYST1")) {
      if (auto params = parsePdbCrystLine(line)) {
        param_a = params->a;
        param_b = params->b;
        param_c = params->c;
        has_box = true;
      }
    } else if (line.starts_with("MODEL")) {
      current_cell = correlation::core::Cell();
      if (has_box) {
        current_cell.setLatticeParameters({param_a, param_b, param_c, 90.0, 90.0, 90.0});
      }
      in_model = true;
    } else if (line.starts_with("ATOM") || line.starts_with("HETATM")) {
      if (auto atom = parsePdbAtomLine(line)) {
        current_cell.addAtom(atom->symbol, correlation::math::Vector3<real_t>(atom->x, atom->y, atom->z));
      }
    } else if (line.starts_with("ENDMDL")) {
      frames.push_back(std::move(current_cell));
      current_cell = correlation::core::Cell();
      in_model = false;
    }
  }

  // If no MODEL tags were found, but we have atoms, it's a single frame PDB
  if (frames.empty()) {
    file.clear();
    file.seekg(0);
    frames.push_back(readStructure(filename, progress_callback));
  }

  if (progress_callback) {
    progress_callback(1.0F, "PDB trajectory loaded.");
  }
  return {frames, 1.0};
}

} // namespace correlation::readers
