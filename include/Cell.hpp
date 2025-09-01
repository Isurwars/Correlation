#ifndef INCLUDE_CELL_HPP_
#define INCLUDE_CELL_HPP_
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include <array>
#include <string>
#include <vector>

#include "../include/Atom.hpp"
#include "../include/LinearAlgebra.hpp"

class Cell {
  // This object contains the lattice parameters,
  // and the atoms that compose the material

private:
  // Cell Vectors private modification, public read with getter
  linalg::Vector3<double> v_a_;
  linalg::Vector3<double> v_b_;
  linalg::Vector3<double> v_c_;
  // Lattice Parameters
  std::array<double, 6> lattice_parameters_;
  // Volume of the cell
  double volume_;
  // List of atoms in the Cell
  std::vector<Atom> atoms_;
  // List of elements in Cell
  std::vector<std::string> elements_;
  // List with the number of atoms of each kind of elements in Cell
  std::vector<int> element_numbers_;
  // List of weights for partials funcions
  std::vector<double> w_ij_;
  // Matrix of Bond-lengths
  std::vector<std::vector<double>> bond_length_;
  // 3D Tensor of Distances
  std::vector<std::vector<std::vector<double>>> distances_;
  // 4D Tensor of Angles
  std::vector<std::vector<std::vector<std::vector<double>>>> angles_;

public:
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//

  explicit Cell(const std::array<double, 6> &);
  explicit Cell(const linalg::Vector3<double> &,
                const linalg::Vector3<double> &,
                const linalg::Vector3<double> &);
  explicit Cell();

  //-------------------------------------------------------------------------//
  //------------------------------- Accessors -------------------------------//
  //-------------------------------------------------------------------------//

  // Lattice Parameters
  const std::array<double, 6> &lattice_parameters() const {
    return lattice_parameters_;
  }
  void setLatticeParameters(std::array<double, 6>);

  // Lattice Vectors
  const linalg::Vector3<double> &v_a() const { return v_a_; };
  const linalg::Vector3<double> &v_b() const { return v_b_; };
  const linalg::Vector3<double> &v_c() const { return v_c_; };
  void calculateLatticeVectors();

  // Volume
  const double &volume() const { return volume_; }
  void calculateVolume();

  // Atoms
  const std::vector<Atom> &atoms() const { return atoms_; }
  void addAtom(Atom);
  void setAtoms(std::vector<Atom>);
  int getElementID(std::string);

  // Elements
  const std::vector<std::string> &elements() const { return elements_; }
  void setElements(const std::vector<std::string> &e) { elements_ = e; }
  void addElement(const std::string &);

  // Element Numbers
  const std::vector<int> &element_numbers() const { return element_numbers_; }
  void setElementsNumbers(const std::vector<int> &e) { element_numbers_ = e; }
  void calculateElementNumbers();

  // Weigth factos
  const std::vector<double> &w_ij() const { return w_ij_; }
  void setWeightFactors(const std::vector<double> &w) { w_ij_ = w; }

  // Distances
  const std::vector<std::vector<std::vector<double>>> &distances() const {
    return distances_;
  }

  // Angles
  const std::vector<std::vector<std::vector<std::vector<double>>>> &
  angles() const {
    return angles_;
  }

  // Bonds
  const std::vector<std::vector<double>> &bond_length() const {
    return bond_length_;
  }
  void setBondLength(const std::vector<std::vector<double>> &bonds) {
    bond_length_ = bonds;
  }

  //-------------------------------------------------------------------------//
  //------------------------------- Methods ---------------------------------//
  //-------------------------------------------------------------------------//

  // element_id
  const int element_id(const std::string &);
  // Populate element_id in all  atoms
  void populateElementID();
  // Populate the Bond length Matrix
  void populateBondLength(double);
  // Correct atom positions
  void correctPositions();
  void correctFracPositions();
  // Calculate Distances (max distance between atoms)
  void distancePopulation(double, bool);
  // Coordination Numbers Calculation
  void coordinationNumber();
  // Bond-Angle Calulation (degrees=true, radian=false)
  void planeAnglePopulation(bool = true);
};
#endif // INCLUDE_CELL_HPP_
