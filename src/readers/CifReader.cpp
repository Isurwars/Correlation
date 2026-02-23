// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "readers/CifReader.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "LinearAlgebra.hpp"

namespace {

// --- Helper Struct for CIF Symmetry Operations ---
struct SymmetryOp {
  linalg::Matrix3<double> rotation{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  linalg::Vector3<double> translation{0, 0, 0};

  // Applies the operation: new_pos = rotation * old_pos + translation
  linalg::Vector3<double> apply(const linalg::Vector3<double> &pos) const {
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
  size_t p_pos = str.find('(');
  if (p_pos != std::string::npos) {
    str.erase(p_pos);
  }
  // Remove surrounding quotes
  if (!str.empty() && ((str.front() == '\'' && str.back() == '\'') ||
                       (str.front() == '"' && str.back() == '"'))) {
    return str.substr(1, str.length() - 2);
  }
  return str;
}

// Parses a single component of a symmetry string like "-y+1/2"
void parseSymmetryComponent(std::string comp_str, int row, SymmetryOp &op) {
  comp_str.erase(std::remove_if(comp_str.begin(), comp_str.end(), ::isspace),
                 comp_str.end());

  double sign = 1.0;
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
      op.rotation(row, 0) = sign;
      current_pos++;
    } else if (comp_str[current_pos] == 'y' || comp_str[current_pos] == 'Y') {
      op.rotation(row, 1) = sign;
      current_pos++;
    } else if (comp_str[current_pos] == 'z' || comp_str[current_pos] == 'Z') {
      op.rotation(row, 2) = sign;
      current_pos++;
    }
    // Check for translation part
    else if (isdigit(comp_str[current_pos])) {
      size_t next_pos;
      double num = std::stod(comp_str.substr(current_pos), &next_pos);
      current_pos += next_pos;
      if (current_pos < comp_str.length() && comp_str[current_pos] == '/') {
        current_pos++; // Skip '/'
        double den = std::stod(comp_str.substr(current_pos), &next_pos);
        current_pos += next_pos;
        op.translation[row] += sign * (num / den);
      } else {
        op.translation[row] += sign * num;
      }
    } else {
      // Should not happen with valid CIF
      current_pos++;
    }
  }
}

// Parses a full symmetry operation string like 'x, y, z+1/2'
SymmetryOp parseSymmetryString(const std::string &op_str) {
  SymmetryOp op;
  // Set rotation to zero initially
  op.rotation = linalg::Matrix3<double>({0, 0, 0}, {0, 0, 0}, {0, 0, 0});
  std::stringstream ss(op_str);
  std::string component;
  int row = 0;

  while (std::getline(ss, component, ',') && row < 3) {
    parseSymmetryComponent(component, row, op);
    row++;
  }
  return op;
}

// A more robust tokenizer that handles quoted strings with spaces
std::vector<std::string> tokenizeCifLine(const std::string &line) {
  std::vector<std::string> tokens;
  std::string current_token;
  char quote_char = '\0';
  bool in_token = false;

  for (char c : line) {
    if (quote_char != '\0') { // Inside a quoted string
      if (c == quote_char) {
        quote_char = '\0'; // End of quoted string
      } else {
        current_token += c;
      }
    } else { // Not inside a quoted string
      if (c == '\'' || c == '"') {
        if (!in_token) {
          quote_char = c;
          in_token = true;
        }
      } else if (isspace(c)) {
        if (in_token) {
          tokens.push_back(current_token);
          current_token.clear();
          in_token = false;
        }
      } else {
        if (!in_token) {
          in_token = true;
        }
        current_token += c;
      }
    }
  }
  if (in_token) {
    tokens.push_back(current_token);
  }
  return tokens;
}

} // namespace

namespace CifReader {

Cell read(const std::string &file_name) {
  std::ifstream file(file_name);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open CIF file: " + file_name);
  }

  Cell tempCell;
  std::string line;

  // --- Data storage during parsing ---
  std::map<std::string, std::string> cif_data;
  struct AsymmetricAtom {
    std::string symbol;
    linalg::Vector3<double> frac_pos;
  };
  std::vector<AsymmetricAtom> asymmetric_atoms;
  std::vector<SymmetryOp> symmetry_ops;

  // --- Parser state machine ---
  enum class ParseState { GLOBAL, LOOP_HEADER, LOOP_DATA };
  ParseState state = ParseState::GLOBAL;
  std::vector<std::string> loop_headers;

  // Read the CIF file line by line
  // We use a state machine approach to handle loops and global data blocks.
  while (std::getline(file, line)) {
    // Trim leading whitespace
    line.erase(0, line.find_first_not_of(" \t\n\r"));
    if (line.empty() || line[0] == '#')
      continue;

    // A new global tag or a new loop definition ends a previous loop's data
    // section
    if (state == ParseState::LOOP_DATA &&
        (line[0] == '_' || line.rfind("loop_", 0) == 0)) {
      state = ParseState::GLOBAL;
      loop_headers.clear();
    }

    if (line.rfind("loop_", 0) == 0) {
      state = ParseState::LOOP_HEADER;
      loop_headers.clear();
      continue;
    }

    if (state == ParseState::GLOBAL) {
      if (line[0] == '_') {
        std::stringstream ss(line);
        std::string key, value;
        ss >> key;
        std::getline(ss, value); // The rest of the line is the value
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
      // Map headers to their column index for efficient lookup
      std::map<std::string, int> header_map;
      for (size_t i = 0; i < loop_headers.size(); ++i) {
        header_map[loop_headers[i]] = i;
      }

      // Check which loop we're in by looking for key headers
      bool is_atom_loop = header_map.count("_atom_site_fract_x");
      bool is_symm_loop = header_map.count("_symmetry_equiv_pos_as_xyz") ||
                          header_map.count("_space_group_symop_operation_xyz");

      auto tokens = tokenizeCifLine(line);
      if (tokens.empty())
        continue;

      if (is_atom_loop) {
        if (tokens.size() < 3)
          continue; // Malformed line
        try {
          std::string element =
              cleanCifValue(tokens.at(header_map.at("_atom_site_type_symbol")));
          element.erase(
              std::remove_if(element.begin(), element.end(), ::isdigit),
              element.end());

          linalg::Vector3<double> pos = {
              std::stod(cleanCifValue(
                  tokens.at(header_map.at("_atom_site_fract_x")))),
              std::stod(cleanCifValue(
                  tokens.at(header_map.at("_atom_site_fract_y")))),
              std::stod(cleanCifValue(
                  tokens.at(header_map.at("_atom_site_fract_z"))))};
          asymmetric_atoms.push_back({element, pos});
        } catch (const std::out_of_range &oor) {
          throw std::runtime_error("CIF Error: Missing required atom site data "
                                   "(e.g., _atom_site_fract_x).");
        }
      } else if (is_symm_loop) {
        // The symmetry operation might be the only token on the line
        std::string op_key = header_map.count("_symmetry_equiv_pos_as_xyz")
                                 ? "_symmetry_equiv_pos_as_xyz"
                                 : "_space_group_symop_operation_xyz";
        symmetry_ops.push_back(
            parseSymmetryString(tokens.at(header_map.at(op_key))));
      }
    }
  }

  // --- Post-processing and Cell Construction ---

  // 1. Set Lattice Parameters
  try {
    std::array<double, 6> params = {
        std::stod(cif_data.at("_cell_length_a")),
        std::stod(cif_data.at("_cell_length_b")),
        std::stod(cif_data.at("_cell_length_c")),
        std::stod(cif_data.at("_cell_angle_alpha")),
        std::stod(cif_data.at("_cell_angle_beta")),
        std::stod(cif_data.at("_cell_angle_gamma"))};
    tempCell.setLatticeParameters(params);
  } catch (const std::out_of_range &oor) {
    throw std::runtime_error(
        "CIF Error: Missing one or more required cell parameters.");
  }

  // 2. Generate all atoms by applying symmetry
  if (symmetry_ops.empty()) { // If no symmetry specified, 'x,y,z' is implicit
    symmetry_ops.push_back(parseSymmetryString("x,y,z"));
  }

  std::vector<AsymmetricAtom> final_atoms;
  const double tolerance = 1e-4;
  for (const auto &atom : asymmetric_atoms) {
    for (const auto &op : symmetry_ops) {
      linalg::Vector3<double> frac_pos = op.apply(atom.frac_pos);

      // Normalize fractional coordinates to be within [0-epsilon, 1-epsilon)
      frac_pos.x() = std::fmod(frac_pos.x(), 1.0);
      if (frac_pos.x() < 0)
        frac_pos.x() += 1.0;
      frac_pos.y() = std::fmod(frac_pos.y(), 1.0);
      if (frac_pos.y() < 0)
        frac_pos.y() += 1.0;
      frac_pos.z() = std::fmod(frac_pos.z(), 1.0);
      if (frac_pos.z() < 0)
        frac_pos.z() += 1.0;

      // Avoid adding duplicate atoms
      bool exists = false;
      for (const auto &final_atom : final_atoms) {
        if (final_atom.symbol != atom.symbol)
          continue;

        linalg::Vector3<double> diff = frac_pos - final_atom.frac_pos;
        // Account for periodic boundary wrapping
        diff.x() = std::fmod(diff.x(), 1.0);
        if (std::abs(diff.x()) > 0.5)
          diff.x() -= std::copysign(1.0, diff.x());
        diff.y() = std::fmod(diff.y(), 1.0);
        if (std::abs(diff.y()) > 0.5)
          diff.y() -= std::copysign(1.0, diff.y());
        diff.z() = std::fmod(diff.z(), 1.0);
        if (std::abs(diff.z()) > 0.5)
          diff.z() -= std::copysign(1.0, diff.z());

        if (linalg::norm(diff) < tolerance) {
          exists = true;
          break;
        }
      }
      if (!exists) {
        final_atoms.push_back({atom.symbol, frac_pos});
      }
    }
  }

  // 3. Convert to Cartesian coordinates and add to the cell
  const auto &lattice = tempCell.latticeVectors();
  for (const auto &atom : final_atoms) {
    tempCell.addAtom(atom.symbol, lattice * atom.frac_pos);
  }

  return tempCell;
}

} // namespace CifReader
