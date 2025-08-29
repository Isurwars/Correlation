// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "../include/ReadFiles.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

//---------------------------------------------------------------------------//
//----------------------------- Helper Funciton -----------------------------//
//---------------------------------------------------------------------------//
void toLower(std::string &s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
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

Cell readCif(std::string file_name) {
  /*
   * This is a Stub for reading a CIF file and return and empty cell. Since
   * CIF files are mostly used for crystalline materials we don't only need
   * to read the file, but we also need to calculate the correct positions of
   * the images in the cell and to verify that no atom is multiply accounted
   * for, due to the different symmetries.
   *
   * This will be implemented in V1.2.
   */
  std::ifstream myfile(file_name);
  // std::string   line;
  Cell tempCell;
  // std::smatch match;

  return tempCell;
} // readCif

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
