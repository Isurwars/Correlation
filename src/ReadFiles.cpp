// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "../include/ReadFiles.hpp"

#include <array>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <iostream>
#include <ostream>
#include <regex>
#include <string>
#include <utility>
#include <vector>

#include "../include/Templates.hpp"

//---------------------------------------------------------------------------//
//----------------------------- Read Bonds File -----------------------------//
//---------------------------------------------------------------------------//
std::vector<std::vector<double>> readBond(std::string file_name, Cell cell) {
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
  std::ifstream myfile(file_name);
  std::smatch match;
  std::vector<std::string> elements = cell.elements();
  std::vector<std::vector<double>> bonds = cell.bond_length();
  /*
   * Every line should have two elements and a bond length separeted by spaces:
   *
   * element_A element_B _bond_length_
   */

  std::regex regex_bond("^([A-Z][a-z]?)"
                        "(\\s+)"
                        "([A-Z][a-z]?)"
                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)");

  if (myfile.is_open()) {
    /* Check if the file is open */
    std::string line;
    while (std::getline(myfile, line)) {
      /* Read line by line */
      if (std::regex_search(line, match, regex_bond)) {
        /* Bond found */
        int i, j;
        std::pair<bool, int> IdA =
            findInVector(elements, std::string(match.str(1).data()));
        if (IdA.first)
          i = IdA.second;
        std::pair<bool, int> IdB =
            findInVector(elements, std::string(match.str(1).data()));
        if (IdB.first)
          j = IdB.second;
        double dist = std::stof(match.str(4).data());
        if (IdA.first && IdB.first) {
          bonds[i][j] = dist;
          bonds[j][i] = dist;
        }
      }
    }
  }
  return bonds;
} // readBond

Cell readCar(std::string file_name) {
  /*
   * This function reads a CAR file and returns a list of atoms objects
   * with the element, number and position inside. It also returns the
   * cell parameters (a, b, c, alpha, beta and gamma) for the periodic
   * repetition of the cell.
   */
  std::ifstream myfile(file_name);
  Cell tempCell;
  std::smatch match;

  /*
   * CAR file is CASE SENSITIVE, and has a strict policy for position
   * and spacing, this simplifies the reading of the file by external codes.
   */

  /*
   * The parameters are in a line that start with PBC, followed by the
   * six lattice parameters: A, B, C, alpha, beta, gamma.
   */
  std::regex regex_parameters("(^!?PBC)"
                              "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                              "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                              "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                              "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                              "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                              "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)");

  /*
   * The atoms start with the element symbol, the first letter is UPPERCASE,
   * a second optional letter is lowercase, followed by a identifier.
   * Then the absolute coordinates (x, y, z) in Angstroms.
   */
  std::regex regex_atom("^([A-Z][a-z]?)"
                        "([A-Z0-9]*)"
                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                        "(\\s+[A-Za-z]+)"
                        "(\\s\\d+)"
                        "(\\s+[A-Za-z]+\\s+)"
                        "([A-Z][a-z]?)");

  if (myfile.is_open()) {
    std::array<double, 6> lat;
    std::string line;
    while (std::getline(myfile, line)) {
      if (std::regex_search(line, match, regex_parameters)) {
        lat[0] = std::stod(match.str(2).data());
        lat[1] = std::stod(match.str(4).data());
        lat[2] = std::stod(match.str(6).data());
        lat[3] = std::stod(match.str(8).data());
        lat[4] = std::stod(match.str(10).data());
        lat[5] = std::stod(match.str(12).data());
        tempCell.setLatticeParameters(lat);
        tempCell.calculateLatticeVectors();
      }
      if (std::regex_search(line, match, regex_atom)) {
        Atom tempAtom(match.str(12).data(), {std::stod(match.str(3).data()),
                                             std::stod(match.str(5).data()),
                                             std::stod(match.str(7).data())});
        tempCell.addAtom(tempAtom);
      }
    }
  } else {
    std::cout << "Unable to read file: " << file_name << " ("
              << std::strerror(errno) << ")." << std::endl;
    exit(1);
  }
  return tempCell;
} // readCar

Cell readCell(std::string file_name) {
  /*
   * This function reads a CELL file and returns a list of atoms objects
   * with the element, number and position inside. It also returns the
   * cell parameters (a, b, c, alpha, beta and gamma) for the periodic
   * repetition of the cell.
   */
  std::ifstream myfile(file_name);
  std::smatch match;

  Cell tempCell;
  bool frac_flag = false;

  /*
   * CELL files are CasE InSenSitiVE. The file is separeted in "BLOCKS".
   * Ordering the BLOCKS is not mandatory, and most of the BLOCKS
   * are optional, while some are mutally exclusive.
   * Inside the BLOCK the data is stored in a strict order, both in
   * position and separation, so we can use regex to read the data inside
   * the BLOCKS, and deal with the disorder between BLOCKS in a different
   * approach.
   */

  /*
   * BLOCK starts with %BLOCK "something", CasE InSenSitiVE
   */
  std::regex regex_block("^(%block)\\s+"
                         "([^\n\r]+)",
                         std::regex::icase);

  /*
   * BLOCK ends with %ENDBLOCK "something", CasE InSenSitiVE
   */
  std::regex regex_endblock("^(%endblock)\\s+"
                            "([^\n\r]+)",
                            std::regex::icase);

  /*
   * There are two lattice block, mutally exclusive:
   * - LATTICE_CART
   *   The three lattice vectors are provided in a 3x3 matrix:
   *   v_ai  v_aj  v_ak
   *   v_bi  v_bj  v_bk
   *   v_ci  v_cj  v_ck
   *
   * - LATTICE_ABC
   *   The 6 lattice parameters are provided in a 3x2 matrix:
   *   A     B     c
   *   alpha beta  gamma
   */

  std::regex regex_lattice_cart("(lattice_cart)", std::regex::icase);
  std::regex regex_lattice_abc("(lattice_abc)", std::regex::icase);
  std::regex regex_lattice("(\\s*[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                           "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                           "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)");

  /*
   * Atoms have two options: absolute coordinates and relative coordinates,
   * both of them have the following structure:
   * Element Symbol, X, Y, Z
   */
  std::regex regex_positions_abs("(positions_abs)", std::regex::icase);
  std::regex regex_positions_frac("(positions_frac)", std::regex::icase);
  std::regex regex_atom("([A-Z][a-z]?)"
                        "([a-z]?[0-9]?)"
                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)");

  if (myfile.is_open()) {
    /* Check if the file is open */
    std::string line;
    std::array<double, 6> lat;
    while (std::getline(myfile, line)) {
      /* Read line by line */
      if (std::regex_search(line, match, regex_block)) {
        /* BLOCK found */
        if (std::regex_search(line, match, regex_lattice_cart)) {
          /* LATTICE_CART case.*/
          int i = 0;
          double aux_v[3][3];
          while (!std::regex_search(line, match, regex_endblock)) {
            std::getline(myfile, line);
            if (std::regex_search(line, match, regex_lattice)) {
              aux_v[i][0] = std::stof(match.str(1).data());
              aux_v[i][1] = std::stof(match.str(3).data());
              aux_v[i][2] = std::stof(match.str(5).data());
              i++;
            }
          }
          tempCell = Cell(
              {aux_v[0][0], aux_v[0][1], aux_v[0][2]},
              {aux_v[1][0], aux_v[1][1], aux_v[1][2]},
              {aux_v[2][0], aux_v[2][1], aux_v[2][2]});
        }
        if (std::regex_search(line, match, regex_lattice_abc)) {
          /* LATTICE_ABC case.*/
          int i = 0;
          while (!std::regex_search(line, match, regex_endblock)) {
            std::getline(myfile, line);
            if (std::regex_search(line, match, regex_lattice)) {
              lat[0 + 3 * i] = std::stof(match.str(1).data());
              lat[1 + 3 * i] = std::stof(match.str(3).data());
              lat[2 + 3 * i] = std::stof(match.str(5).data());
              tempCell.setLatticeParameters(lat);
              i++;
            }
          }
          tempCell.calculateLatticeVectors();
        }
        if (std::regex_search(line, match, regex_positions_frac)) {
          /* POSITIONS_FRAC case.*/
          frac_flag = true;
          while (!std::regex_search(line, match, regex_endblock)) {
            std::getline(myfile, line);
            if (std::regex_search(line, match, regex_atom)) {
              Atom tempAtom(match.str(1).data(),
                            {std::stod(match.str(3).data()),
                             std::stod(match.str(5).data()),
                             std::stod(match.str(7).data())});
              tempCell.addAtom(tempAtom);
              if (!(findInVector(tempCell.elements(), tempAtom.element())
                        .first)) {
                tempCell.addElement(tempAtom.element());
              }
            }
          }
        }
        if (std::regex_search(line, match, regex_positions_abs)) {
          /* POSITIONS_ABS case.*/
          frac_flag = false;
          while (!std::regex_search(line, match, regex_endblock)) {
            std::getline(myfile, line);
            if (std::regex_search(line, match, regex_atom)) {
              Atom tempAtom(match.str(1).data(),
                            {std::stod(match.str(3).data()),
                             std::stod(match.str(5).data()),
                             std::stod(match.str(7).data())});
              tempCell.addAtom(tempAtom);
              if (!(findInVector(tempCell.elements(), tempAtom.element())
                        .first)) {
                tempCell.addElement(tempAtom.element());
              }
            }
          }
        }
      }
    }
  } else {
    std::cout << "Unable to read file: " << file_name << " ("
              << std::strerror(errno) << ")." << std::endl;
    exit(1);
  }
  if (frac_flag)
    tempCell.correctFracPositions();
  return tempCell;
} // readCell

Cell readOnetepDat(std::string file_name) {
  /*
   * This function reads a CELL file and returns a list of atoms objects
   * with the element, number and position inside. It also returns the
   * cell parameters (a, b, c, alpha, beta and gamma) for the periodic
   * repetition of the cell.
   */
  std::ifstream myfile(file_name);
  std::smatch match;
  Cell tempCell;
  bool frac_flag = false;

  /*
   * CELL files are CasE InSenSitiVE. The file is separeted in "BLOCKS".
   * Ordering the BLOCKS is not mandatory, and most of the BLOCKS
   * are optional, while some are mutally exclusive.
   * Inside the BLOCK the data is stored in a strict order, both in
   * position and separation, so we can use regex to read the data inside
   * the BLOCKS, and deal with the disorder between BLOCKS in a different
   * approach.
   */

  /*
   * BLOCK starts with %BLOCK "something", CasE InSenSitiVE
   */
  std::regex regex_block("^(%block)\\s+"
                         "([^\n\r]+)",
                         std::regex::icase);

  /*
   * BLOCK ends with %ENDBLOCK "something", CasE InSenSitiVE
   */
  std::regex regex_endblock("^(%endblock)\\s+"
                            "([^\n\r]+)",
                            std::regex::icase);

  /*
   * There are two lattice block, mutally exclusive:
   * - LATTICE_CART
   *   The three lattice vectors are provided in a 3x3 matrix:
   *   v_ai  v_aj  v_ak
   *   v_bi  v_bj  v_bk
   *   v_ci  v_cj  v_ck
   *
   * - LATTICE_ABC
   *   The 6 lattice parameters are provided in a 3x2 matrix:
   *   A     B     c
   *   alpha beta  gamma
   */

  std::regex regex_lattice_cart("(lattice_cart)", std::regex::icase);
  std::regex regex_lattice_abc("(lattice_abc)", std::regex::icase);
  std::regex regex_lattice("(\\s*[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                           "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                           "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)");

  /*
   * Onetep default distance units is a0, we need to know if we should
   * change from a.u. to Å.
   */
  std::regex regex_ang("(ang)", std::regex::icase);

  /*
   * Atoms have two options: absolute coordinates and relative coordinates,
   * both of them have the following structure:
   * Element Symbol, X, Y, Z
   */
  std::regex regex_positions_abs("(positions_abs)", std::regex::icase);
  std::regex regex_positions_frac("(positions_frac)", std::regex::icase);
  std::regex regex_atom("([A-Z][a-z]?)"
                        "([a-z]?[0-9]?)"
                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)");

  if (myfile.is_open()) {
    std::string line;
    std::array<double, 6> lat;
    double p_f = 0.52918;
    /* Check if the file is open */
    while (std::getline(myfile, line)) {
      /* Read line by line */
      if (std::regex_search(line, match, regex_block)) {
        /* BLOCK found */
        if (std::regex_search(line, match, regex_lattice_cart)) {
          /* LATTICE_CART case.*/
          int i = 0;
          double aux_v[3][3];
          while (!std::regex_search(line, match, regex_endblock)) {
            std::getline(myfile, line);
            if (std::regex_search(line, match, regex_ang)) {
              p_f = 1.0;
            }
            if (std::regex_search(line, match, regex_lattice)) {
              aux_v[i][0] = std::stof(match.str(1).data()) * p_f;
              aux_v[i][1] = std::stof(match.str(3).data()) * p_f;
              aux_v[i][2] = std::stof(match.str(5).data()) * p_f;
              i++;
            }
          }
          tempCell = Cell({aux_v[0][0], aux_v[0][1], aux_v[0][2]},
                          {aux_v[1][0], aux_v[1][1], aux_v[1][2]},
                          {aux_v[2][0], aux_v[2][1], aux_v[2][2]});
        }
        if (std::regex_search(line, match, regex_lattice_abc)) {
          /* LATTICE_ABC case.*/
          int i = 0;
          while (!std::regex_search(line, match, regex_endblock)) {
            std::getline(myfile, line);
            if (std::regex_search(line, match, regex_ang)) {
              p_f = 1.0;
            }
            if (std::regex_search(line, match, regex_lattice)) {
              lat[0 + 3 * i] = std::stof(match.str(1).data()) * p_f;
              lat[1 + 3 * i] = std::stof(match.str(3).data()) * p_f;
              lat[2 + 3 * i] = std::stof(match.str(5).data()) * p_f;
              tempCell.setLatticeParameters(lat);
              i++;
            }
          }
          tempCell.calculateLatticeVectors();
        }
        if (std::regex_search(line, match, regex_positions_frac)) {
          /* POSITIONS_FRAC case.*/
          frac_flag = true;
          while (!std::regex_search(line, match, regex_endblock)) {
            std::getline(myfile, line);
            if (std::regex_search(line, match, regex_ang)) {
              p_f = 1.0;
            }
            if (std::regex_search(line, match, regex_atom)) {
              Atom tempAtom(match.str(1).data(),
                            {std::stod(match.str(3).data()) * p_f,
                             std::stod(match.str(5).data()) * p_f,
                             std::stod(match.str(7).data()) * p_f});
              tempCell.addAtom(tempAtom);
              if (!(findInVector(tempCell.elements(), tempAtom.element())
                        .first)) {
                tempCell.addElement(tempAtom.element());
              }
            }
          }
        }
        if (std::regex_search(line, match, regex_positions_abs)) {
          /* POSITIONS_ABS case.*/
          frac_flag = false;
          while (!std::regex_search(line, match, regex_endblock)) {
            std::getline(myfile, line);
            if (std::regex_search(line, match, regex_ang)) {
              p_f = 1.0;
            }
            if (std::regex_search(line, match, regex_atom)) {
              Atom tempAtom(match.str(1).data(),
                            {std::stod(match.str(3).data()) * p_f,
                             std::stod(match.str(5).data()) * p_f,
                             std::stod(match.str(7).data()) * p_f});
              tempCell.addAtom(tempAtom);
              if (!(findInVector(tempCell.elements(), tempAtom.element())
                        .first)) {
                tempCell.addElement(tempAtom.element());
              }
            }
          }
        }
      }
    }
  } else {
    std::cout << "Unable to read file: " << file_name << " ("
              << std::strerror(errno) << ")." << std::endl;
    exit(1);
  }
  if (frac_flag)
    tempCell.correctFracPositions();
  return tempCell;
} // readOnetepDat

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

Cell readLammpsDump(std::string file_name) {
  /*
   * This is a Stub for reading a Dump file from LAMMPS. Since DUMP files
   * have custom format we are supposing the standard.
   * ID TYPE X Y Z
   *
   *
   */
  std::ifstream myfile(file_name);
  // std::string   line;
  Cell tempCell;
  // std::smatch match;

  return tempCell;
} // readCif
