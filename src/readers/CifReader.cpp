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

#include <array>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace correlation::readers {

// Automatic registration
const bool REGISTERED = ReaderFactory::registerTypeSafe<CifReader>("CifReader");

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
  [[nodiscard]] correlation::math::Vector3<real_t>
  apply(const correlation::math::Vector3<real_t> &pos) const {
    return rotation * pos + translation;
  }
};

// Helper to trim whitespace and remove trailing parentheses with uncertainties
[[nodiscard]] std::string_view cleanCifValue(std::string_view str_view) {
  // Trim leading whitespace
  size_t const start = str_view.find_first_not_of(" \t\n\r");
  if (start == std::string_view::npos) {
    return {};
  }
  str_view.remove_prefix(start);
  // Trim trailing whitespace
  size_t const end = str_view.find_last_not_of(" \t\n\r");
  str_view = str_view.substr(0, end + 1);

  // Remove uncertainty in parentheses, e.g., "1.234(5)" -> "1.234"
  size_t const p_pos = str_view.find('(');
  if (p_pos != std::string_view::npos) {
    str_view = str_view.substr(0, p_pos);
  }
  // Remove surrounding quotes
  if (str_view.size() >= 2 && ((str_view.front() == '\'' && str_view.back() == '\'') ||
                               (str_view.front() == '"' && str_view.back() == '"'))) {
    return str_view.substr(1, str_view.size() - 2);
  }
  return str_view;
}

template <typename T> [[nodiscard]] bool parseCifNumber(std::string_view str_view, T &out) {
  std::string_view const cleaned = cleanCifValue(str_view);
  if (cleaned.empty()) {
    return false;
  }
  std::array<char, 64> buf{};
  char const *c_str = nullptr;
  std::string fallback;
  if (cleaned.size() < buf.size()) {
    std::memcpy(buf.data(), cleaned.data(), cleaned.size());
    buf.at(cleaned.size()) = '\0';
    c_str = buf.data();
  } else {
    fallback = std::string(cleaned);
    c_str = fallback.c_str();
  }
  char *end = nullptr;
  errno = 0;
  double const val = std::strtod(c_str, &end);
  if (end != c_str + cleaned.size() || errno == ERANGE) {
    return false;
  }
  out = static_cast<T>(val);
  return true;
}

void parseRotationAxis(char axis_char, int row, real_t sign, SymmetryOp &sym_op) {
  if (axis_char == 'x' || axis_char == 'X') {
    sym_op.rotation(row, 0) = sign;
  } else if (axis_char == 'y' || axis_char == 'Y') {
    sym_op.rotation(row, 1) = sign;
  } else if (axis_char == 'z' || axis_char == 'Z') {
    sym_op.rotation(row, 2) = sign;
  }
}

void parseTranslationPart(std::string_view comp_str, size_t &current_pos, int row, real_t sign,
                          SymmetryOp &sym_op) {
  std::string_view const remaining = comp_str.substr(current_pos);
  std::array<char, 64> buf{};
  char const *c_str = nullptr;
  std::string fallback;
  if (remaining.size() < buf.size()) {
    std::memcpy(buf.data(), remaining.data(), remaining.size());
    buf.at(remaining.size()) = '\0';
    c_str = buf.data();
  } else {
    fallback = std::string(remaining);
    c_str = fallback.c_str();
  }

  char *end1 = nullptr;
  errno = 0;
  double const num = std::strtod(c_str, &end1);
  if (end1 == c_str || errno == ERANGE) {
    current_pos++;
    return;
  }
  auto const consumed1 = static_cast<size_t>(end1 - c_str);
  current_pos += consumed1;

  if (current_pos < comp_str.length() && comp_str[current_pos] == '/') {
    current_pos++; // Skip '/'
    std::string_view const den_remaining = comp_str.substr(current_pos);
    std::array<char, 64> den_buf{};
    char const *den_c_str = nullptr;
    std::string den_fallback;
    if (den_remaining.size() < den_buf.size()) {
      std::memcpy(den_buf.data(), den_remaining.data(), den_remaining.size());
      den_buf.at(den_remaining.size()) = '\0';
      den_c_str = den_buf.data();
    } else {
      den_fallback = std::string(den_remaining);
      den_c_str = den_fallback.c_str();
    }

    char *end2 = nullptr;
    errno = 0;
    double const den = std::strtod(den_c_str, &end2);
    if (end2 != den_c_str && errno != ERANGE && den != 0.0) {
      auto const consumed2 = static_cast<size_t>(end2 - den_c_str);
      current_pos += consumed2;
      sym_op.translation[row] += sign * static_cast<real_t>(num / den);
    }
  } else {
    sym_op.translation[row] += sign * static_cast<real_t>(num);
  }
}

// Parses a single component of a symmetry string like "-y+1/2"
void parseSymmetryComponent(std::string comp_str, int row, SymmetryOp &sym_op) {
  std::erase_if(comp_str, ::isspace);

  real_t sign = 1.0;
  size_t current_pos = 0;

  while (current_pos < comp_str.length()) {
    char const current_ch = comp_str[current_pos];
    if (current_ch == '+') {
      sign = 1.0;
      current_pos++;
    } else if (current_ch == '-') {
      sign = -1.0;
      current_pos++;
    } else if (isdigit(static_cast<unsigned char>(current_ch)) != 0) {
      parseTranslationPart(comp_str, current_pos, row, sign, sym_op);
    } else {
      parseRotationAxis(current_ch, row, sign, sym_op);
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

enum class TokenizerState : std::uint8_t { OutsideToken, InsideUnquoted, InsideQuoted };

void handleOutsideToken(char chr, char &quote_char, TokenizerState &state,
                        std::string &current_token) {
  if (chr == '\'' || chr == '"') {
    quote_char = chr;
    state = TokenizerState::InsideQuoted;
  } else if (::isspace(static_cast<unsigned char>(chr)) == 0) {
    current_token += chr;
    state = TokenizerState::InsideUnquoted;
  }
}

void handleInsideUnquoted(char chr, TokenizerState &state, std::string &current_token,
                          std::vector<std::string> &tokens) {
  if (::isspace(static_cast<unsigned char>(chr)) != 0) {
    tokens.push_back(current_token);
    current_token.clear();
    state = TokenizerState::OutsideToken;
  } else if (chr != '\'' && chr != '"') {
    current_token += chr;
  }
}

void handleInsideQuoted(char chr, char &quote_char, TokenizerState &state,
                        std::string &current_token) {
  if (chr == quote_char) {
    quote_char = '\0';
    state = TokenizerState::InsideUnquoted;
  } else {
    current_token += chr;
  }
}

// A more robust tokenizer that handles quoted strings with spaces
std::vector<std::string> tokenizeCifLine(const std::string &line) {
  std::vector<std::string> tokens;
  std::string current_token;
  TokenizerState state = TokenizerState::OutsideToken;
  char quote_char = '\0';

  for (char const chr : line) {
    switch (state) {
    case TokenizerState::OutsideToken:
      handleOutsideToken(chr, quote_char, state, current_token);
      break;
    case TokenizerState::InsideUnquoted:
      handleInsideUnquoted(chr, state, current_token, tokens);
      break;
    case TokenizerState::InsideQuoted:
      handleInsideQuoted(chr, quote_char, state, current_token);
      break;
    }
  }
  if (state != TokenizerState::OutsideToken) {
    tokens.push_back(current_token);
  }
  return tokens;
}

struct AsymmetricAtom {
  std::string symbol;
  correlation::math::Vector3<real_t> frac_pos;
};

enum class ParseState : std::uint8_t { Global, LoopHeader, LoopData };

void processLoopDataLine(const std::string &line, const std::vector<std::string> &loop_headers,
                         std::vector<AsymmetricAtom> &asymmetric_atoms,
                         std::vector<SymmetryOp> &symmetry_ops) {
  // Map headers to their column index for efficient lookup
  std::map<std::string, size_t> header_map;
  for (size_t i = 0; i < loop_headers.size(); ++i) {
    header_map[loop_headers[i]] = i;
  }

  // Check which loop we're in by looking for key headers
  bool const is_atom_loop = header_map.contains("_atom_site_fract_x");
  bool const is_symm_loop = header_map.contains("_symmetry_equiv_pos_as_xyz") ||
                            header_map.contains("_space_group_symop_operation_xyz");

  auto tokens = tokenizeCifLine(line);
  if (tokens.empty()) {
    return;
  }

  if (is_atom_loop) {
    if (tokens.size() < 3) {
      return; // Malformed line
    }
    try {
      std::string element(cleanCifValue(tokens.at(header_map.at("_atom_site_type_symbol"))));
      std::erase_if(element, ::isdigit);

      real_t x = 0;
      real_t y = 0;
      real_t z = 0;
      if (!parseCifNumber(tokens.at(header_map.at("_atom_site_fract_x")), x) ||
          !parseCifNumber(tokens.at(header_map.at("_atom_site_fract_y")), y) ||
          !parseCifNumber(tokens.at(header_map.at("_atom_site_fract_z")), z)) {
        throw std::runtime_error("CIF Error: Invalid atom site coordinates.");
      }

      correlation::math::Vector3<real_t> const pos = {x, y, z};
      asymmetric_atoms.push_back({
          .symbol = element,
          .frac_pos = pos,
      });
    } catch (const std::out_of_range &oor) {
      throw std::runtime_error("CIF Error: Missing required atom site data "
                               "(e.g., _atom_site_fract_x).");
    }
  } else if (is_symm_loop) {
    // The symmetry operation might be the only token on the line
    std::string const op_key =
        (static_cast<unsigned int>(header_map.contains("_symmetry_equiv_pos_as_xyz")) != 0U)
            ? "_symmetry_equiv_pos_as_xyz"
            : "_space_group_symop_operation_xyz";
    symmetry_ops.push_back(parseSymmetryString(tokens.at(header_map.at(op_key))));
  }
}

void processCifLine(const std::string &line, ParseState &state,
                    std::vector<std::string> &loop_headers,
                    std::map<std::string, std::string> &cif_data,
                    std::vector<AsymmetricAtom> &asymmetric_atoms,
                    std::vector<SymmetryOp> &symmetry_ops) {
  // A new global tag or a new loop definition ends a previous loop's data section
  if (state == ParseState::LoopData && (line[0] == '_' || line.starts_with("loop_"))) {
    state = ParseState::Global;
    loop_headers.clear();
  }

  if (line.starts_with("loop_")) {
    state = ParseState::LoopHeader;
    loop_headers.clear();
    return;
  }

  if (state == ParseState::Global) {
    if (line[0] == '_') {
      std::stringstream str_stream(line);
      std::string key;
      std::string value;
      str_stream >> key;
      std::getline(str_stream, value); // The rest of the line is the value
      cif_data[key] = std::string(cleanCifValue(value));
    }
  } else if (state == ParseState::LoopHeader) {
    if (line[0] == '_') {
      loop_headers.push_back(line);
    } else {
      state = ParseState::LoopData;
      // Fall through to process this line as the first data line
    }
  }

  if (state == ParseState::LoopData) {
    processLoopDataLine(line, loop_headers, asymmetric_atoms, symmetry_ops);
  }
}

void parseCifFile(std::ifstream &file, std::map<std::string, std::string> &cif_data,
                  std::vector<AsymmetricAtom> &asymmetric_atoms,
                  std::vector<SymmetryOp> &symmetry_ops) {
  std::string line;
  ParseState state = ParseState::Global;
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

void setupLatticeParameters(correlation::core::Cell &cell,
                            const std::map<std::string, std::string> &cif_data) {
  try {
    std::array<real_t, 6> params{};
    const std::array<const char *, 6> keys = {"_cell_length_a",   "_cell_length_b",
                                              "_cell_length_c",   "_cell_angle_alpha",
                                              "_cell_angle_beta", "_cell_angle_gamma"};

    for (size_t i = 0; i < 6; ++i) {
      if (!parseCifNumber(cif_data.at(keys.at(i)), params.at(i))) {
        throw std::runtime_error("Invalid cell parameter value for " + std::string(keys.at(i)));
      }
    }
    cell.setLatticeParameters(params);
  } catch (const std::exception &e) {
    throw std::runtime_error("CIF Error: Missing or invalid cell parameters: " +
                             std::string(e.what()));
  }
}

bool isDuplicateAtom(const correlation::math::Vector3<real_t> &pos1,
                     const correlation::math::Vector3<real_t> &pos2, real_t tolerance) {
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

std::vector<AsymmetricAtom>
generateSymmetryAtoms(const std::vector<AsymmetricAtom> &asymmetric_atoms,
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
      auto normalize_coord = [](real_t val) -> real_t {
        auto res = static_cast<real_t>(std::fmod(val, 1.0));
        if (res < 0) {
          res += 1.0;
        }
        return res;
      };
      frac_pos.x() = normalize_coord(frac_pos.x());
      frac_pos.y() = normalize_coord(frac_pos.y());
      frac_pos.z() = normalize_coord(frac_pos.z());

      // Avoid adding duplicate atoms
      bool exists = false;
      for (const auto &final_atom : final_atoms) {
        if (final_atom.symbol == atom.symbol &&
            isDuplicateAtom(frac_pos, final_atom.frac_pos, tolerance)) {
          exists = true;
          break;
        }
      }
      if (!exists) {
        final_atoms.push_back({
            .symbol = atom.symbol,
            .frac_pos = frac_pos,
        });
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

  correlation::core::Cell temp_cell;
  std::map<std::string, std::string> cif_data;
  std::vector<AsymmetricAtom> asymmetric_atoms;
  std::vector<SymmetryOp> symmetry_ops;

  // 1. Parse CIF file
  parseCifFile(file, cif_data, asymmetric_atoms, symmetry_ops);

  // 2. Set Lattice Parameters
  setupLatticeParameters(temp_cell, cif_data);

  // 3. Generate all atoms by applying symmetry and filtering duplicates
  const std::vector<AsymmetricAtom> final_atoms =
      generateSymmetryAtoms(asymmetric_atoms, symmetry_ops);

  // 4. Convert to Cartesian coordinates and add to the cell
  const auto &lattice = temp_cell.latticeVectors();
  for (const auto &atom : final_atoms) {
    temp_cell.addAtom(atom.symbol, lattice * atom.frac_pos);
  }

  return temp_cell;
}

} // namespace correlation::readers
