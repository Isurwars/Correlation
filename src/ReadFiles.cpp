#include <iostream>
#include <fstream>
#include <list>
#include <regex>

#include "Atom.h"
#include "Cell.h"

/*
 * Generic function to find if an element of any type exists in vector,
 * if true, then return the index.
 */
template <typename T>
std::pair<bool, int> findInVector(const std::vector<T> & vecOfElements, const T  & element)
{
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
}// findInVector

Cell read_CAR(std::string file_name)
{
    /*
     * This function reads a CAR file, it returns a list of atoms objects
     * with the element, number and position inside, it also returns the
     * cell parameters (a, b, c, alpha, beta and gamma) for the periodic
     * repetition of the cell.
     */
    std::ifstream myfile(file_name);
    std::string line;
    Cell tempCell;
    std::smatch match;

    /*
     * CAR file is CASE SENSITIVE, and has a strict policy for position
     * and spacing this simplify the reading of the file by external codes.
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
     * The atoms start with the element symbol, First letter in UPPERCASE,
     * second optional letter in lowercase followed by a identifing integer.
     * Then the absolute coordinates (x, y, z) in Angstroms.
     */
    std::regex regex_atom("^([A-Z][a-z]?)"
      "(\\d+)"
      "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
      "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
      "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
      "(\\s+[A-Z]+)"
      "(\\s\\d+)"
      "(\\s+[a-z]+)"
      "(\\s+[A-Z][a-z]?)");

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
    }
    return tempCell;
} // read_CAR

Cell read_CELL(std::string file_name)
{
    /*
     * This function reads a CELL file, it returns a list of atoms objects
     * with the element, number and position inside, it also returns the
     * cell parameters (a, b, c, alpha, beta and gamma) for the periodic
     * repetition of the cell.
     */
    std::ifstream myfile(file_name);
    std::string line;
    std::smatch match;
    Cell tempCell;
    bool frac_flag = false;

    /*
     * CELL files are CasE InSenSitiVE. The file is separeted in "BLOCKS".
     * The order of the BLOCKS are not mandatory, and most of the BLOCKS
     * are optional, while some are mutal exclusive.
     * Inside the block the data is stored in an strict order, both in
     * position and separation, so we can use regex to read the data inside
     * the BLOCKS, and deal with the disorder between BLOCKS in a different
     * approach.
     */

    /*
     * BLOCK start with %BLOCK "something", CasE InSenSitiVE
     */
    std::regex regex_block("^(%block)\\s+"
      "([^\n\r]+)", std::regex::icase);

    /*
     * BLOCK end with %ENDBLOCK "something", CasE InSenSitiVE
     */
    std::regex regex_endblock("^(%endblock)\\s+"
      "([^\n\r]+)", std::regex::icase);


    /*
     * There are two lattice block, mutal exclusive:
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
    std::regex regex_lattice("(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
      "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
      "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)");

    /*
     * Atoms has two options: absolute coordinates and relative coordinates,
     * both  of them has the following structure:
     * Element Symbol, x, Y, Z
     */
    std::regex regex_positions_abs("(positions_abs)", std::regex::icase);
    std::regex regex_positions_frac("(positions_frac)", std::regex::icase);
    std::regex regex_atom("([A-Z][a-z]?)"
      "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
      "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)"
      "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)");

    if (myfile.is_open()) {
        /* Check if tje file is open */
        while (std::getline(myfile, line)) {
            /* Begin reading line by line */
            if (std::regex_search(line, match, regex_block)) {
                /* BLOCK found */
                if (std::regex_search(line, match, regex_lattice_cart)) {
                    /* LATTICE_CART case.*/
                    int index = 0;
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
                              { std::stod(match.str(2).data()),
                                std::stod(match.str(4).data()),
                                std::stod(match.str(6).data()) });
                            tempCell.atoms.push_back(tempAtom);
                            if (!(findInVector(tempCell.elements, tempAtom.element).first)) {
                                tempCell.elements.push_back(tempAtom.element);
                            }
                        }
                    }
                }
                if (std::regex_search(line, match, regex_positions_abs)) {
                    /* POSITIONS_ABS case.*/
                    frac_flag = true;
                    while (!std::regex_search(line, match, regex_endblock)) {
                        std::getline(myfile, line);
                        if (std::regex_search(line, match, regex_atom)) {
                            Atom tempAtom(match.str(1).data(),
                              { std::stod(match.str(2).data()),
                                std::stod(match.str(4).data()),
                                std::stod(match.str(6).data()) });
                            tempCell.atoms.push_back(tempAtom);
                            if (!(findInVector(tempCell.elements, tempAtom.element).first)) {
                                tempCell.elements.push_back(tempAtom.element);
                            }
                        }
                    }
                }
            }
        }
    }
    if (frac_flag) tempCell.CorrectFracPositions();
    return tempCell;
} // read_CELL

Cell read_CIF(std::string file_name)
{
    /*
     * This is a Stub for reading a CIF file and return and empthy cell.
     * Since most of the CIF files are used for crystalline materials
     * we not only need to read the file, but also calculate the correct
     * positions of the images in the cell and that no atom is duplicated
     * by the different symmetries.
     * This will be implemented AFTER the bond angle is finished.
     */
    std::ifstream myfile(file_name);
    std::string line;
    Cell tempCell;
    std::smatch match;


    return tempCell;
} // read_CIF
