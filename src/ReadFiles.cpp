/* ---------------------------------------------------------------------
 * Correlation: An Analysis Tool for Liquids and for Amorphous Solids
 * Copyright (c) 2013-2021 Isaías Rodríguez <isurwars@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the MIT License version as published in:
 * https://github.com/Isurwars/Correlation/blob/main/LICENSE
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * ----------------------------------------------------------------------
 */
#include "ReadFiles.h"

#include <array>
#include <fstream>
#include <iostream>
#include <list>
#include <regex>
#include <string.h>
#include <utility>
#include <vector>

#include "Atom.h"

/*
 * Generic function to find if an element of any type exists in vector,
 * if true, then return the index.
 */
template<typename T>
std::pair<bool, int> findInVector(const std::vector<T>& vecOfElements,
  const T                                             & element) {
  std::pair<bool, int> result;
  // Find given element in vector
  auto it = std::find(vecOfElements.begin(), vecOfElements.end(), element);

  if (it != vecOfElements.end()) {
    result.second = std::distance(vecOfElements.begin(), it);
    result.first  = true;
  } else {
    result.first  = false;
    result.second = -1;
  }
  return result;
}  // findInVector
Cell read_CAR(std::string file_name) {
  /*
   * This function reads a CAR file and returns a list of atoms objects
   * with the element, number and position inside. It also returns the
   * cell parameters (a, b, c, alpha, beta and gamma) for the periodic
   * repetition of the cell.
   */
  std::ifstream myfile(file_name);
  std::string   line;
  Cell        tempCell;
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
   * a second optional letter is lowercase, followed by a identifing integer.
   * Then the absolute coordinates (x, y, z) in Angstroms.
   */
  std::regex regex_atom("^([A-Z][a-z]?)"
                        "(\\d+)"
                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                        "(\\s+[A-Z]+)"
                        "(\\s\\d+)"
                        "(\\s+[a-z]+\\s+)"
                        "([A-Z][a-z]?)");

  if (myfile.is_open()) {
    while (std::getline(myfile, line)) {
      if (std::regex_search(line, match, regex_parameters)) {
        tempCell.lattice_parameters[0] = std::stod(match.str(2).data());
        tempCell.lattice_parameters[1] = std::stod(match.str(4).data());
        tempCell.lattice_parameters[2] = std::stod(match.str(6).data());
        tempCell.lattice_parameters[3] = std::stod(match.str(8).data());
        tempCell.lattice_parameters[4] = std::stod(match.str(10).data());
        tempCell.lattice_parameters[5] = std::stod(match.str(12).data());
        tempCell.SetLatticeVectors();
      }
      if (std::regex_search(line, match, regex_atom)) {
        Atom tempAtom(match.str(12).data(),
          { std::stod(match.str(3).data()),
            std::stod(match.str(5).data()),
            std::stod(match.str(7).data()) });
        tempCell.atoms.push_back(tempAtom);
        if (!(findInVector(tempCell.elements, tempAtom.element).first)) {
          tempCell.elements.push_back(tempAtom.element);
        }
      }
    }
  } else {
    std::cout << "Unable to read file: "
              << file_name
              << " ("
              << strerror(errno)
              << ")."
              << std::endl;
    exit(1);
  }
  return tempCell;
}  // read_CAR
Cell read_CELL(std::string file_name) {
  /*
   * This function reads a CELL file and returns a list of atoms objects
   * with the element, number and position inside. It also returns the
   * cell parameters (a, b, c, alpha, beta and gamma) for the periodic
   * repetition of the cell.
   */
  std::ifstream myfile(file_name);
  std::string   line;
  std::smatch   match;
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
                         "([^\n\r]+)", std::regex::icase);

  /*
   * BLOCK ends with %ENDBLOCK "something", CasE InSenSitiVE
   */
  std::regex regex_endblock("^(%endblock)\\s+"
                            "([^\n\r]+)", std::regex::icase);


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
  std::regex regex_atom("([A-Z][a-z]?)" "([a-z]?[0-9]?)"
                                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)");

  if (myfile.is_open()) {
    /* Check if the file is open */
    while (std::getline(myfile, line)) {
      /* Read line by line */
      if (std::regex_search(line, match, regex_block)) {
        /* BLOCK found */
        if (std::regex_search(line, match, regex_lattice_cart)) {
          /* LATTICE_CART case.*/
          int    index = 0;
          double aux_v[3][3];
          while (!std::regex_search(line, match, regex_endblock)) {
            std::getline(myfile, line);
            if (std::regex_search(line, match, regex_lattice)) {
              aux_v[index][0] = std::stof(match.str(1).data());
              aux_v[index][1] = std::stof(match.str(3).data());
              aux_v[index][2] = std::stof(match.str(5).data());
              index++;
            }
          }
          tempCell.SetFromVectors({ aux_v[0][0], aux_v[0][1], aux_v[0][2] },
            { aux_v[1][0], aux_v[1][1], aux_v[1][2] },
            { aux_v[2][0], aux_v[2][1], aux_v[2][2] });
        }
        if (std::regex_search(line, match, regex_lattice_abc)) {
          /* LATTICE_ABC case.*/
          int index = 0;
          while (!std::regex_search(line, match, regex_endblock)) {
            std::getline(myfile, line);
            if (std::regex_search(line, match, regex_lattice)) {
              tempCell.lattice_parameters[0 + 3 * index] = std::stof(match.str(1).data());
              tempCell.lattice_parameters[1 + 3 * index] = std::stof(match.str(3).data());
              tempCell.lattice_parameters[2 + 3 * index] = std::stof(match.str(5).data());
              index++;
            }
          }
          tempCell.SetLatticeVectors();
        }
        if (std::regex_search(line, match, regex_positions_frac)) {
          /* POSITIONS_FRAC case.*/
          frac_flag = true;
          while (!std::regex_search(line, match, regex_endblock)) {
            std::getline(myfile, line);
            if (std::regex_search(line, match, regex_atom)) {
              Atom tempAtom(match.str(1).data(),
                { std::stod(match.str(3).data()),
                  std::stod(match.str(5).data()),
                  std::stod(match.str(7).data()) });
              tempCell.atoms.push_back(tempAtom);
              if (!(findInVector(tempCell.elements, tempAtom.element).first)) {
                tempCell.elements.push_back(tempAtom.element);
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
                { std::stod(match.str(3).data()),
                  std::stod(match.str(5).data()),
                  std::stod(match.str(7).data()) });
              tempCell.atoms.push_back(tempAtom);
              if (!(findInVector(tempCell.elements, tempAtom.element).first)) {
                tempCell.elements.push_back(tempAtom.element);
              }
            }
          }
        }
      }
    }
  } else {
    std::cout << "Unable to read file: "
              << file_name
              << " ("
              << strerror(errno)
              << ")."
              << std::endl;
    exit(1);
  }
  if (frac_flag) tempCell.CorrectFracPositions();
  return tempCell;
}  // read_CELL
Cell read_ONETEP_DAT(std::string file_name) {
  /*
   * This function reads a CELL file and returns a list of atoms objects
   * with the element, number and position inside. It also returns the
   * cell parameters (a, b, c, alpha, beta and gamma) for the periodic
   * repetition of the cell.
   */
  std::ifstream myfile(file_name);
  std::string   line;
  std::smatch   match;
  Cell   tempCell;
  bool   frac_flag = false;
  double p_f       = 0.52918;

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
                         "([^\n\r]+)", std::regex::icase);

  /*
   * BLOCK ends with %ENDBLOCK "something", CasE InSenSitiVE
   */
  std::regex regex_endblock("^(%endblock)\\s+"
                            "([^\n\r]+)", std::regex::icase);


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
  std::regex regex_atom("([A-Z][a-z]?)" "([a-z]?[0-9]?)"
                                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
                                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)");

  if (myfile.is_open()) {
    /* Check if the file is open */
    while (std::getline(myfile, line)) {
      /* Read line by line */
      if (std::regex_search(line, match, regex_block)) {
        /* BLOCK found */
        if (std::regex_search(line, match, regex_lattice_cart)) {
          /* LATTICE_CART case.*/
          int    index = 0;
          double aux_v[3][3];
          while (!std::regex_search(line, match, regex_endblock)) {
            std::getline(myfile, line);
            if (std::regex_search(line, match, regex_ang)) {
              p_f = 1.0;
            }
            if (std::regex_search(line, match, regex_lattice)) {
              aux_v[index][0] = std::stof(match.str(1).data()) * p_f;
              aux_v[index][1] = std::stof(match.str(3).data()) * p_f;
              aux_v[index][2] = std::stof(match.str(5).data()) * p_f;
              index++;
            }
          }
          tempCell.SetFromVectors({ aux_v[0][0], aux_v[0][1], aux_v[0][2] },
            { aux_v[1][0], aux_v[1][1], aux_v[1][2] },
            { aux_v[2][0], aux_v[2][1], aux_v[2][2] });
        }
        if (std::regex_search(line, match, regex_lattice_abc)) {
          /* LATTICE_ABC case.*/
          int index = 0;
          while (!std::regex_search(line, match, regex_endblock)) {
            std::getline(myfile, line);
            if (std::regex_search(line, match, regex_ang)) {
              p_f = 1.0;
            }
            if (std::regex_search(line, match, regex_lattice)) {
              tempCell.lattice_parameters[0 + 3 * index] = std::stof(match.str(1).data()) * p_f;
              tempCell.lattice_parameters[1 + 3 * index] = std::stof(match.str(3).data()) * p_f;
              tempCell.lattice_parameters[2 + 3 * index] = std::stof(match.str(5).data()) * p_f;
              index++;
            }
          }
          tempCell.SetLatticeVectors();
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
                { std::stod(match.str(3).data()) * p_f,
                  std::stod(match.str(5).data()) * p_f,
                  std::stod(match.str(7).data()) * p_f });
              tempCell.atoms.push_back(tempAtom);
              if (!(findInVector(tempCell.elements, tempAtom.element).first)) {
                tempCell.elements.push_back(tempAtom.element);
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
                { std::stod(match.str(3).data()) * p_f,
                  std::stod(match.str(5).data()) * p_f,
                  std::stod(match.str(7).data()) * p_f });
              tempCell.atoms.push_back(tempAtom);
              if (!(findInVector(tempCell.elements, tempAtom.element).first)) {
                tempCell.elements.push_back(tempAtom.element);
              }
            }
          }
        }
      }
    }
  } else {
    std::cout << "Unable to read file: "
              << file_name
              << " ("
              << strerror(errno)
              << ")."
              << std::endl;
    exit(1);
  }
  if (frac_flag) tempCell.CorrectFracPositions();
  return tempCell;
}  // read_ONETEP_DAT
Cell read_CIF(std::string file_name) {
  /*
   * This is a Stub for reading a CIF file and return and empty cell.
   * Since CIF files are mostly used for crystalline materials we don't
   * only need to read the file, but we also need to calculate the correct
   * positions of the images in the cell and to verify that no atom is
   * multiply accounted for, due to the different symmetries.
   * This will be implemented AFTER the bond angle is finished.
   */
  std::ifstream myfile(file_name);
  std::string   line;
  Cell        tempCell;
  std::smatch match;


  return tempCell;
}  // read_CIF
