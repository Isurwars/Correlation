// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "../include/ReadFiles.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

//---------------------------------------------------------------------------//
//---------------------------- Helper Funcitons -----------------------------//
//---------------------------------------------------------------------------//
// Helper to send all UPPERCASE to lowercase
void toLower(std::string &s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
}

// Helper to trim whitespace and remove trailing parentheses with uncertainties
static std::string cleanToken(std::string s) {
  // Remove uncertainty part like "(1)" or ".123(45)" -> ".123"
  size_t paren_pos = s.find('(');
  if (paren_pos != std::string::npos) {
    s = s.substr(0, paren_pos);
  }
  // Trim leading/trailing whitespace
  s.erase(0, s.find_first_not_of(" \t\n\r"));
  s.erase(s.find_last_not_of(" \t\n\r") + 1);
  return s;
}

// Helper for robust substring replacement
static void replaceAll(std::string &str, const std::string &from,
                       const std::string &to) {
  size_t start_pos = 0;
  while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length(); // Move past the replaced instance
  }
}

// Helper to apply a symmetry operation string like '-x, y+1/2, z'
static Vector3D applySymmetry(const std::string &op, const Vector3D &v) {
  Vector3D result = {0.0, 0.0, 0.0};
  std::stringstream ss(op);
  std::string component;
  int i = 0;

  while (std::getline(ss, component, ',')) {
    double val = 0.0;
    double sign = 1.0;
    std::string term;
    // Replace x,y,z with their values to make it a simple expression
    replaceAll(component, "x", std::to_string(v[0]));
    replaceAll(component, "y", std::to_string(v[1]));
    replaceAll(component, "z", std::to_string(v[2]));

    // Pad operators with spaces to ensure correct tokenization
    replaceAll(component, "+", " + ");
    replaceAll(component, "-", " - ");

    // Simple expression evaluation
    std::stringstream term_ss(component);
    while (term_ss >> term) {
      if (term == "+") {
        sign = 1.0;
      } else if (term == "-") {
        sign = -1.0;
      } else {
        // Handle fractions like 1/2
        size_t slash_pos = term.find('/');
        if (slash_pos != std::string::npos) {
          double num = std::stod(term.substr(0, slash_pos));
          double den = std::stod(term.substr(slash_pos + 1));
          val += sign * (num / den);
        } else {
          val += sign * std::stod(term);
        }
        sign = 1.0; // Reset sign after processing a number
      }
    }
    result[i++] = std::fmod(val, 1.0);
  }
  return result;
}
//---------------------------------------------------------------------------//
//----------------------------- Read Bonds File -----------------------------//
//---------------------------------------------------------------------------//
std::vector<std::vector<double>> readBond(const std::string &file_name,
                                          const Cell &cell) {
  /*
   * This function reads the in_bond_file to populate the _bond_length_ Tensor.
   * The file should be in the format:
   * element_B element_A distance(in Angstroms)
   *
   * For example:
   *
   * Si Si 2.29
   * Mg Mg 2.85
   * C  C  1.55
   * C  Si 1.86
   * Si Mg 2.57
   * C  Mg 2.07
   *
   * Any missing pair of elements will use the bond_parameter as a default.
   */

  // Open the file and check for errors.
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to open bond file: " + file_name);
  }

  std::vector<std::string> elements = cell.elements();
  std::vector<std::vector<double>> bonds = cell.bond_length();

  std::unordered_map<std::string, int> element_map;
  for (int i = 0; i < elements.size(); ++i) {
    element_map[elements[i]] = i;
  }

  std::string line;
  std::string elemA;
  std::string elemB;
  double distance;
  // Read the file line by line.
  while (std::getline(myfile, line)) {
    // 2. Use stringstream for faster, simpler parsing.
    std::stringstream ss(line);
    if (ss >> elemA >> elemB >> distance) {
      // 3. Use the map to find indices quickly.
      auto itA = element_map.find(elemA);
      auto itB = element_map.find(elemB);

      // Ensure both elements were found in our map.
      if (itA != element_map.end() && itB != element_map.end()) {
        int i = itA->second;
        int j = itB->second;
        bonds[i][j] = distance;
        bonds[j][i] = distance;
      }
    }
  }
  return bonds;
} // readBond

Cell readCar(const std::string &file_name) {
  /*
   * This function reads a CAR file and returns a list of atoms objects
   * with the element, number and position inside. It also returns the
   * cell parameters (a, b, c, alpha, beta and gamma) for the periodic
   * repetition of the cell.
   */
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" +
                             std::strerror(errno) + ").");
  }

  Cell tempCell;
  std::string line;

  while (std::getline(myfile, line)) {

    // Ignore empty lines or comment lines
    if (line.empty() || line[0] == '!') {
      continue;
    }

    std::stringstream line_stream(line);
    std::string first_token;
    line_stream >> first_token;

    if (first_token == "PBC") {
      std::array<double, 6> lat;
      // The token "PBC" is consumed, so we read the 6 numbers that follow.
      if (line_stream >> lat[0] >> lat[1] >> lat[2] >> lat[3] >> lat[4] >>
          lat[5]) {
        tempCell.setLatticeParameters(lat);
        tempCell.calculateLatticeVectors();
      }
      continue;
    }

    if (first_token == "PBC=OFF") {
      std::array<double, 6> lat = {100.0, 100.0, 100.0, 90.0, 90.0, 90.0};
      tempCell.setLatticeParameters(lat);
      tempCell.calculateLatticeVectors();
      continue;
    }

    // Reset the stream to parse the full line as an atom entry.
    line_stream.clear();
    line_stream.seekg(0);

    // Declare variables for all 8 columns we need to read.
    std::string u1, u5, u6, u7, element;
    double x, y, z;

    // Read exactly 8 columns to get the element.
    if (line_stream >> u1 >> x >> y >> z >> u5 >> u6 >> u7 >> element) {
      Atom tempAtom(element, {x, y, z});
      tempCell.addAtom(tempAtom);
    }
  }

  return tempCell;
} // readCar

Cell readCell(const std::string &file_name) {
  /*
   * This function reads a CELL file and returns a list of atoms objects
   * with the element, number and position inside. It also returns the
   * cell parameters (a, b, c, alpha, beta and gamma) for the periodic
   * repetition of the cell.
   */
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" +
                             std::strerror(errno) + ").");
  }

  Cell tempCell;
  bool in_block = false;
  bool frac_flag = false;
  std::string current_block_type;
  std::string line;
  int lattice_row_count = 0;
  std::array<double, 6> lat;

  while (std::getline(myfile, line)) {
    /* Read line by line */
    std::stringstream line_stream(line);
    std::string token;
    line_stream >> token;
    toLower(token); // For case-insensitive matching

    if (token == "%block") {
      in_block = true;
      line_stream >> current_block_type;
      toLower(current_block_type);
      lattice_row_count = 0; // Reset for a new lattice block
      continue;
    }

    if (token == "%endblock") {
      in_block = false;
      current_block_type.clear();
      continue;
    }

    if (in_block) {
      if (current_block_type == "lattice_cart") {
        double v[3][3];
        line_stream.clear();
        line_stream.seekg(0); // Reread the full line
        if (line_stream >> v[lattice_row_count][0] >> v[lattice_row_count][1] >>
            v[lattice_row_count][2]) {
          lattice_row_count++;
          if (lattice_row_count == 3) {
            tempCell =
                Cell({v[0][0], v[0][1], v[0][2]}, {v[1][0], v[1][1], v[1][2]},
                     {v[2][0], v[2][1], v[2][2]});
          }
        }
      }

      if (current_block_type == "lattice_abc") {

        line_stream.clear();
        line_stream.seekg(0);
        if (line_stream >> lat[0 + 3 * lattice_row_count] >>
            lat[1 + 3 * lattice_row_count] >> lat[2 + 3 * lattice_row_count]) {
          lattice_row_count++;
          if (lattice_row_count == 2) {
            tempCell.setLatticeParameters(lat);
            tempCell.calculateLatticeVectors();
          }
        }
      }

      if (current_block_type == "positions_abs") {
        std::string element;
        double x, y, z;
        line_stream.clear();
        line_stream.seekg(0);
        if (line_stream >> element >> x >> y >> z) {
          Atom tempAtom(element, {x, y, z});
          tempCell.addAtom(tempAtom);
        }
      }

      if (current_block_type == "positions_frac") {
        std::string element;
        double x, y, z;

        frac_flag = true;
        line_stream.clear();
        line_stream.seekg(0);
        if (line_stream >> element >> x >> y >> z) {
          Atom tempAtom(element, {x, y, z});
          tempCell.addAtom(tempAtom);
        }
      }
    }
  }

  if (frac_flag)
    tempCell.correctFracPositions();
  return tempCell;
} // readCell

Cell readCif(const std::string &file_name) {
  std::ifstream file(file_name);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open CIF file: " + file_name);
  }

  Cell tempCell;
  std::string line;

  std::map<std::string, std::string> cif_data;
  std::vector<Atom> asymmetric_atoms;
  std::vector<std::string> symmetry_ops;

  enum class ParseState { GLOBAL, LOOP_HEADER, LOOP_DATA };
  ParseState state = ParseState::GLOBAL;
  std::vector<std::string> loop_headers;

  while (std::getline(file, line)) {
    // Trim whitespace
    line.erase(0, line.find_first_not_of(" \t\n\r"));
    line.erase(line.find_last_not_of(" \t\n\r") + 1);

    if (line.empty() || line[0] == '#')
      continue;

    std::string lower_line = line;
    std::transform(lower_line.begin(), lower_line.end(), lower_line.begin(),
                   ::tolower);

    if (lower_line.rfind("loop_", 0) == 0) {
      state = ParseState::LOOP_HEADER;
      loop_headers.clear();
      continue;
    }

    if (state == ParseState::GLOBAL) {
      if (line[0] == '_') {
        std::stringstream ss(line);
        std::string key, value;
        ss >> key;
        // The rest of the line is the value
        std::getline(ss, value);
        cif_data[key] = cleanToken(value);
      }
    } else if (state == ParseState::LOOP_HEADER) {
      if (line[0] == '_') {
        loop_headers.push_back(line);
      } else {
        // First data line after headers
        state = ParseState::LOOP_DATA;
        // Fall through to process this line as data
      }
    }

    if (state == ParseState::LOOP_DATA) {
      // Find which loop we are in
      bool is_atom_loop = std::find(loop_headers.begin(), loop_headers.end(),
                                    "_atom_site_fract_x") != loop_headers.end();
      bool is_symm_loop =
          std::find(loop_headers.begin(), loop_headers.end(),
                    "_symmetry_equiv_pos_as_xyz") != loop_headers.end();

      std::stringstream ss(line);
      if (is_atom_loop) {
        // Map headers to their column index
        auto get_idx = [&](const std::string &key) {
          auto it = std::find(loop_headers.begin(), loop_headers.end(), key);
          return (it == loop_headers.end())
                     ? -1
                     : std::distance(loop_headers.begin(), it);
        };

        int sym_idx = get_idx("_atom_site_type_symbol");
        int x_idx = get_idx("_atom_site_fract_x");
        int y_idx = get_idx("_atom_site_fract_y");
        int z_idx = get_idx("_atom_site_fract_z");

        if (x_idx == -1 || y_idx == -1 || z_idx == -1 || sym_idx == -1) {
          throw std::runtime_error("Incomplete atom site loop in CIF file.");
        }

        std::vector<std::string> tokens;
        std::string token;
        while (ss >> token)
          tokens.push_back(token);

        if (tokens.size() >= loop_headers.size()) {
          std::string element = cleanToken(tokens[sym_idx]);
          // Remove numbers from element labels like "C1", "Si2"
          element.erase(
              std::remove_if(element.begin(), element.end(), ::isdigit),
              element.end());

          Vector3D pos = {std::stod(cleanToken(tokens[x_idx])),
                          std::stod(cleanToken(tokens[y_idx])),
                          std::stod(cleanToken(tokens[z_idx]))};
          asymmetric_atoms.emplace_back(element, pos);
        } else {
          // End of loop data
          state = ParseState::GLOBAL;
        }

      } else if (is_symm_loop) {
        // The entire line (quotes removed) is the symmetry op
        line.erase(std::remove(line.begin(), line.end(), '\''), line.end());
        line.erase(std::remove(line.begin(), line.end(), '"'), line.end());
        symmetry_ops.push_back(cleanToken(line));
      } else {
        // Some other loop we don't care about, wait for it to end
        if (line[0] == '_')
          state = ParseState::LOOP_HEADER;
        if (line.rfind("loop_", 0) == 0)
          state = ParseState::LOOP_HEADER;
      }
    }
  }

  // --- Post-processing ---
  // 1. Set Lattice Parameters
  std::array<double, 6> params = {std::stod(cif_data.at("_cell_length_a")),
                                  std::stod(cif_data.at("_cell_length_b")),
                                  std::stod(cif_data.at("_cell_length_c")),
                                  std::stod(cif_data.at("_cell_angle_alpha")),
                                  std::stod(cif_data.at("_cell_angle_beta")),
                                  std::stod(cif_data.at("_cell_angle_gamma"))};
  tempCell.setLatticeParameters(params);

  // 2. Generate all atoms by applying symmetry
  std::vector<Atom> final_atoms;
  if (symmetry_ops.empty()) { // If no symmetry specified, 'x,y,z' is implicit
    symmetry_ops.push_back("x, y, z");
  }

  for (const auto &atom : asymmetric_atoms) {
    for (const auto &op : symmetry_ops) {
      Vector3D frac_pos = applySymmetry(op, atom.position());

      // Normalize fractional coordinates to be within [0, 1)
      frac_pos[0] = std::fmod(frac_pos[0], 1.0);
      if (frac_pos[0] < 0)
        frac_pos[0] += 1.0;
      frac_pos[1] = std::fmod(frac_pos[1], 1.0);
      if (frac_pos[1] < 0)
        frac_pos[1] += 1.0;
      frac_pos[2] = std::fmod(frac_pos[2], 1.0);
      if (frac_pos[2] < 0)
        frac_pos[2] += 1.0;

      // Avoid adding duplicate atoms with a small tolerance
      bool exists = false;
      for (const auto &final_atom : final_atoms) {
        Vector3D diff = {std::abs(final_atom.position()[0] - frac_pos[0]),
                         std::abs(final_atom.position()[1] - frac_pos[1]),
                         std::abs(final_atom.position()[2] - frac_pos[2])};
        // Check periodic boundary difference
        if (diff[0] > 0.5)
          diff[0] = 1.0 - diff[0];
        if (diff[1] > 0.5)
          diff[1] = 1.0 - diff[1];
        if (diff[2] > 0.5)
          diff[2] = 1.0 - diff[2];

        if (final_atom.element() == atom.element() && norm(diff) < 1e-4) {
          exists = true;
          break;
        }
      }
      if (!exists) {
        final_atoms.emplace_back(atom.element(), frac_pos);
      }
    }
  }

  // 3. Set atoms and convert to Cartesian coordinates
  tempCell.setAtoms(final_atoms);
  tempCell.correctFracPositions();

  return tempCell;
}

Cell readLammpsDump(const std::string &file_name) {
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" +
                             std::strerror(errno) + ").");
  }

  Cell tempCell;
  std::string line;
  int num_atoms = 0;

  // Read the file line by line
  while (std::getline(myfile, line)) {
    // Find the beginning of a new timestep block
    if (line.find("ITEM: TIMESTEP") != std::string::npos) {
      // A new timestep is starting, so we clear the old cell data
      // This ensures we only keep the data for the *last* timestep
      tempCell = Cell();

      // The next section will be the header for this timestep
      // Read NUMBER OF ATOMS
      std::getline(myfile, line); // Timestep number (we ignore it)
      std::getline(myfile, line); // "ITEM: NUMBER OF ATOMS"
      std::getline(myfile, line); // The actual number
      num_atoms = std::stoi(line);

      // Read BOX BOUNDS
      std::getline(myfile, line); // "ITEM: BOX BOUNDS..."
      double xlo, xhi, ylo, yhi, zlo, zhi;
      std::getline(myfile, line);
      std::stringstream(line) >> xlo >> xhi;
      std::getline(myfile, line);
      std::stringstream(line) >> ylo >> yhi;
      std::getline(myfile, line);
      std::stringstream(line) >> zlo >> zhi;

      // Create an orthorhombic cell from the box bounds
      tempCell = Cell({xhi - xlo, 0.0, 0.0}, {0.0, yhi - ylo, 0.0},
                      {0.0, 0.0, zhi - zlo});

      // Read ATOMS header
      std::getline(myfile, line); // "ITEM: ATOMS id type x y z"

      // Now, read the specified number of atom lines
      for (int i = 0; i < num_atoms; ++i) {
        std::getline(myfile, line);
        std::stringstream atom_stream(line);
        int id, type;
        double x, y, z;

        // We assume the format "id type x y z"
        if (atom_stream >> id >> type >> x >> y >> z) {
          // Use the numeric type as the element string
          std::string element = std::to_string(type);
          Atom tempAtom(element, {x, y, z});
          tempCell.addAtom(tempAtom);
        }
      }
    }
  }

  return tempCell;
} // readLammpsDump

Cell readOnetepDat(std::string file_name) {
  /*
   * This is a Stub for reading a ONETEP file and return and empty cell.
   */
  std::ifstream myfile(file_name);
  // std::string   line;
  Cell tempCell;
  // std::smatch match;

  return tempCell;
} // readOnetepDat
