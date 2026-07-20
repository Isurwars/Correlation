/**
 * @file CifReader.cpp
 * @brief Implementation of the CIF file format reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */
#include "readers/CifReader.hpp"

#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "math/LinearAlgebra.hpp"
#include "readers/ReaderFactory.hpp"
#include <math.h>

#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace correlation::readers {

// Automatic registration
// NOLINTNEXTLINE(cert-err58-cpp, bugprone-throwing-static-initialization)
static const bool registered = ReaderFactory::instance().registerReader(std::make_unique<CifReader>());

correlation::core::Cell
CifReader::readStructure(const std::string &filename,
                         std::function<void(float, const std::string &)> /*progress_callback*/) {
  return read(filename);
}

correlation::core::Trajectory
CifReader::readTrajectory(const std::string & /*filename*/,
                          std::function<void(float, const std::string &)> /*progress_callback*/) {
  throw std::runtime_error("CIF files are structures, use readStructure.");
}

namespace {

// --- Helper Struct for CIF Symmetry Operations ---
struct SymmetryOp {
  correlation::math::Matrix3<real_t> rotation{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  correlation::math::Vector3<real_t> translation{0, 0, 0};

  // Applies the operation: new_pos = rotation * old_pos + translation
  [[nodiscard]] correlation::math::Vector3<real_t> apply(const correlation::math::Vector3<real_t> &pos) const {
    return rotation * pos + translation;
  }
};

// Helper to trim whitespace and remove trailing parentheses with uncertainties
std::string cleanCifValue(std::string str) {
  // Trim leading whitespace
  str.erase(0, str.find_first_not_of(" \t\n\r"));
  // Trim trailing whitespace
  str.erase(str.find_last_not_of(" \t\n\r") + 1);
  // Remove uncertainty in parentheses, e.g., "1.234(5)" -> "1.234"
  size_t const p_pos = str.find('(');
  if (p_pos != std::string::npos) {
    str.erase(p_pos);
  }
  // Remove surrounding quotes
  if (!str.empty() && ((str.front() == '\'' && str.back() == '\'') || (str.front() == '"' && str.back() == '"'))) {
    return str.substr(1, str.length() - 2);
  }
  return str;
}

// Parses a single component of a symmetry string like "-y+1/2"
void parseSymmetryComponent(std::string comp_str, int row, SymmetryOp &sym_op) {
  std::erase_if(comp_str, ::isspace);

  real_t sign = 1.0;
  size_t current_pos = 0;

  while (current_pos < comp_str.length()) {
    if (comp_str[current_pos] == '+') {
      sign = 1.0;
      current_pos++;
    } else if (comp_str[current_pos] == '-') {
      sign = -1.0;
      current_pos++;
    }

    // Check for x, y, z rotation/permutation part
    if (comp_str[current_pos] == 'x' || comp_str[current_pos] == 'X') {
      sym_op.rotation(row, 0) = sign;
      current_pos++;
    } else if (comp_str[current_pos] == 'y' || comp_str[current_pos] == 'Y') {
      sym_op.rotation(row, 1) = sign;
      current_pos++;
    } else if (comp_str[current_pos] == 'z' || comp_str[current_pos] == 'Z') {
      sym_op.rotation(row, 2) = sign;
      current_pos++;
    }
    // Check for translation part
    else if (isdigit(comp_str[current_pos]) != 0) {
      size_t next_pos = 0;
      real_t const num = static_cast<real_t>(std::stod(comp_str.substr(current_pos), &next_pos));
      current_pos += next_pos;
      if (current_pos < comp_str.length() && comp_str[current_pos] == '/') {
        current_pos++; // Skip '/'
        real_t const den = static_cast<real_t>(std::stod(comp_str.substr(current_pos), &next_pos));
        current_pos += next_pos;
        sym_op.translation[row] += sign * (num / den);
      } else {
        sym_op.translation[row] += sign * num;
      }
    } else {
      // Should not happen with valid CIF
      current_pos++;
    }
  }
}

// Parses a full symmetry operation string like 'x, y, z+1/2'
SymmetryOp parseSymmetryString(const std::string &op_str) {
  SymmetryOp sym_op;
  // Set rotation to zero initially
  sym_op.rotation = correlation::math::Matrix3<real_t>({0, 0, 0}, {0, 0, 0}, {0, 0, 0});
  std::stringstream str_stream(op_str);
  std::string component;
  int row = 0;

  while (std::getline(str_stream, component, ',') && row < 3) {
    parseSymmetryComponent(component, row, sym_op);
    row++;
  }
  return sym_op;
}

enum class TokenizerState : std::uint8_t { OUTSIDE_TOKEN, INSIDE_UNQUOTED, INSIDE_QUOTED };

void handleOutsideToken(char chr, char &quote_char, TokenizerState &state, std::string &current_token) {
  if (chr == '\'' || chr == '"') {
    quote_char = chr;
    state = TokenizerState::INSIDE_QUOTED;
  } else if (::isspace(static_cast<unsigned char>(chr)) == 0) {
    current_token += chr;
    state = TokenizerState::INSIDE_UNQUOTED;
  }
}

void handleInsideUnquoted(char chr, TokenizerState &state, std::string &current_token,
                          std::vector<std::string> &tokens) {
  if (::isspace(static_cast<unsigned char>(chr)) != 0) {
    tokens.push_back(current_token);
    current_token.clear();
    state = TokenizerState::OUTSIDE_TOKEN;
  } else if (chr != '\'' && chr != '"') {
    current_token += chr;
  }
}

void handleInsideQuoted(char chr, char &quote_char, TokenizerState &state, std::string &current_token) {
  if (chr == quote_char) {
    quote_char = '\0';
    state = TokenizerState::INSIDE_UNQUOTED;
  } else {
    current_token += chr;
  }
}

// A more robust tokenizer that handles quoted strings with spaces
std::vector<std::string> tokenizeCifLine(const std::string &line) {
  std::vector<std::string> tokens;
  std::string current_token;
  TokenizerState state = TokenizerState::OUTSIDE_TOKEN;
  char quote_char = '\0';

  for (char const chr : line) {
    switch (state) {
    case TokenizerState::OUTSIDE_TOKEN:
      handleOutsideToken(chr, quote_char, state, current_token);
      break;
    case TokenizerState::INSIDE_UNQUOTED:
      handleInsideUnquoted(chr, state, current_token, tokens);
      break;
    case TokenizerState::INSIDE_QUOTED:
      handleInsideQuoted(chr, quote_char, state, current_token);
      break;
    }
  }
  if (state != TokenizerState::OUTSIDE_TOKEN) {
    tokens.push_back(current_token);
  }
  return tokens;
}

struct AsymmetricAtom {
  std::string symbol;
  correlation::math::Vector3<real_t> frac_pos;
};

enum class ParseState : std::uint8_t { GLOBAL, LOOP_HEADER, LOOP_DATA };

void processLoopDataLine(const std::string &line, const std::vector<std::string> &loop_headers,
                         std::vector<AsymmetricAtom> &asymmetric_atoms, std::vector<SymmetryOp> &symmetry_ops) {
  // Map headers to their column index for efficient lookup
  std::map<std::string, size_t> header_map;
  for (size_t i = 0; i < loop_headers.size(); ++i) {
    header_map[loop_headers[i]] = i;
  }

  // Check which loop we're in by looking for key headers
  bool const is_atom_loop = header_map.contains("_atom_site_fract_x");
  bool const is_symm_loop =
      header_map.contains("_symmetry_equiv_pos_as_xyz") || header_map.contains("_space_group_symop_operation_xyz");

  auto tokens = tokenizeCifLine(line);
  if (tokens.empty()) {
    return;
  }

  if (is_atom_loop) {
    if (tokens.size() < 3) {
      return; // Malformed line
    }
    try {
      std::string element = cleanCifValue(tokens.at(header_map.at("_atom_site_type_symbol")));
      std::erase_if(element, ::isdigit);

      correlation::math::Vector3<real_t> const pos = {
          static_cast<real_t>(std::stod(cleanCifValue(tokens.at(header_map.at("_atom_site_fract_x"))))),
          static_cast<real_t>(std::stod(cleanCifValue(tokens.at(header_map.at("_atom_site_fract_y"))))),
          static_cast<real_t>(std::stod(cleanCifValue(tokens.at(header_map.at("_atom_site_fract_z")))))};
      asymmetric_atoms.push_back({.symbol = element, .frac_pos = pos});
    } catch (const std::out_of_range &oor) {
      throw std::runtime_error("CIF Error: Missing required atom site data "
                               "(e.g., _atom_site_fract_x).");
    }
  } else if (is_symm_loop) {
    // The symmetry operation might be the only token on the line
    std::string const op_key = (static_cast<unsigned int>(header_map.contains("_symmetry_equiv_pos_as_xyz")) != 0U)
                                   ? "_symmetry_equiv_pos_as_xyz"
                                   : "_space_group_symop_operation_xyz";
    symmetry_ops.push_back(parseSymmetryString(tokens.at(header_map.at(op_key))));
  }
}

void processCifLine(const std::string &line, ParseState &state, std::vector<std::string> &loop_headers,
                    std::map<std::string, std::string> &cif_data, std::vector<AsymmetricAtom> &asymmetric_atoms,
                    std::vector<SymmetryOp> &symmetry_ops) {
  // A new global tag or a new loop definition ends a previous loop's data section
  if (state == ParseState::LOOP_DATA && (line[0] == '_' || line.starts_with("loop_"))) {
    state = ParseState::GLOBAL;
    loop_headers.clear();
  }

  if (line.starts_with("loop_")) {
    state = ParseState::LOOP_HEADER;
    loop_headers.clear();
    return;
  }

  if (state == ParseState::GLOBAL) {
    if (line[0] == '_') {
      std::stringstream str_stream(line);
      std::string key;
      std::string value;
      str_stream >> key;
      std::getline(str_stream, value); // The rest of the line is the value
      cif_data[key] = cleanCifValue(value);
    }
  } else if (state == ParseState::LOOP_HEADER) {
    if (line[0] == '_') {
      loop_headers.push_back(line);
    } else {
      state = ParseState::LOOP_DATA;
      // Fall through to process this line as the first data line
    }
  }

  if (state == ParseState::LOOP_DATA) {
    processLoopDataLine(line, loop_headers, asymmetric_atoms, symmetry_ops);
  }
}

void parseCifFile(std::ifstream &file, std::map<std::string, std::string> &cif_data,
                  std::vector<AsymmetricAtom> &asymmetric_atoms, std::vector<SymmetryOp> &symmetry_ops) {
  std::string line;
  ParseState state = ParseState::GLOBAL;
  std::vector<std::string> loop_headers;

  while (std::getline(file, line)) {
    // Trim leading whitespace
    line.erase(0, line.find_first_not_of(" \t\n\r"));
    if (line.empty() || line[0] == '#') {
      continue;
    }
    processCifLine(line, state, loop_headers, cif_data, asymmetric_atoms, symmetry_ops);
  }
}

void setupLatticeParameters(correlation::core::Cell &cell, const std::map<std::string, std::string> &cif_data) {
  try {
    std::array<real_t, 6> const params = {static_cast<real_t>(std::stod(cif_data.at("_cell_length_a"))),
                                          static_cast<real_t>(std::stod(cif_data.at("_cell_length_b"))),
                                          static_cast<real_t>(std::stod(cif_data.at("_cell_length_c"))),
                                          static_cast<real_t>(std::stod(cif_data.at("_cell_angle_alpha"))),
                                          static_cast<real_t>(std::stod(cif_data.at("_cell_angle_beta"))),
                                          static_cast<real_t>(std::stod(cif_data.at("_cell_angle_gamma")))};
    cell.setLatticeParameters(params);
  } catch (const std::exception &e) {
    throw std::runtime_error("CIF Error: Missing or invalid cell parameters: " + std::string(e.what()));
  }
}

bool isDuplicateAtom(const correlation::math::Vector3<real_t> &pos1, const correlation::math::Vector3<real_t> &pos2,
                     real_t tolerance) {
  correlation::math::Vector3<real_t> diff = pos1 - pos2;
  // Account for periodic boundary wrapping
  diff.x() = static_cast<real_t>(std::fmod(diff.x(), 1.0));
  if (std::abs(diff.x()) > 0.5) {
    diff.x() -= static_cast<real_t>(std::copysign(1.0, diff.x()));
  }
  diff.y() = static_cast<real_t>(std::fmod(diff.y(), 1.0));
  if (std::abs(diff.y()) > 0.5) {
    diff.y() -= static_cast<real_t>(std::copysign(1.0, diff.y()));
  }
  diff.z() = static_cast<real_t>(std::fmod(diff.z(), 1.0));
  if (std::abs(diff.z()) > 0.5) {
    diff.z() -= static_cast<real_t>(std::copysign(1.0, diff.z()));
  }
  return correlation::math::norm(diff) < tolerance;
}

std::vector<AsymmetricAtom> generateSymmetryAtoms(const std::vector<AsymmetricAtom> &asymmetric_atoms,
                                                  std::vector<SymmetryOp> &symmetry_ops) {
  if (symmetry_ops.empty()) { // If no symmetry specified, 'x,y,z' is implicit
    symmetry_ops.push_back(parseSymmetryString("x,y,z"));
  }

  std::vector<AsymmetricAtom> final_atoms;
  const real_t tolerance = 1e-4;

  for (const auto &atom : asymmetric_atoms) {
    for (const auto &sym_op : symmetry_ops) {
      correlation::math::Vector3<real_t> frac_pos = sym_op.apply(atom.frac_pos);

      // Normalize fractional coordinates to be within [0-epsilon, 1-epsilon)
      auto normalizeCoord = [](real_t val) -> real_t {
        auto res = static_cast<real_t>(std::fmod(val, 1.0));
        if (res < 0) {
          res += 1.0;
        }
        return res;
      };
      frac_pos.x() = normalizeCoord(frac_pos.x());
      frac_pos.y() = normalizeCoord(frac_pos.y());
      frac_pos.z() = normalizeCoord(frac_pos.z());

      // Avoid adding duplicate atoms
      bool exists = false;
      for (const auto &final_atom : final_atoms) {
        if (final_atom.symbol == atom.symbol && isDuplicateAtom(frac_pos, final_atom.frac_pos, tolerance)) {
          exists = true;
          break;
        }
      }
      if (!exists) {
        final_atoms.push_back({.symbol = atom.symbol, .frac_pos = frac_pos});
      }
    }
  }
  return final_atoms;
}

} // namespace

correlation::core::Cell CifReader::read(const std::string &file_name) {
  std::ifstream file(file_name);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open CIF file: " + file_name);
  }

  correlation::core::Cell tempCell;
  std::map<std::string, std::string> cif_data;
  std::vector<AsymmetricAtom> asymmetric_atoms;
  std::vector<SymmetryOp> symmetry_ops;

  // 1. Parse CIF file
  parseCifFile(file, cif_data, asymmetric_atoms, symmetry_ops);

  // 2. Set Lattice Parameters
  setupLatticeParameters(tempCell, cif_data);

  // 3. Generate all atoms by applying symmetry and filtering duplicates
  std::vector<AsymmetricAtom> final_atoms = generateSymmetryAtoms(asymmetric_atoms, symmetry_ops);

  // 4. Convert to Cartesian coordinates and add to the cell
  const auto &lattice = tempCell.latticeVectors();
  for (const auto &atom : final_atoms) {
    tempCell.addAtom(atom.symbol, lattice * atom.frac_pos);
  }

  return tempCell;
}

} // namespace correlation::readers
