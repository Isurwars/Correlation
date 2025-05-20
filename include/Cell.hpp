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
  Vector3D v_a_;
  Vector3D v_b_;
  Vector3D v_c_;
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
  // Matrix of g(r) Histograms
  std::vector<std::vector<double>> g_;
  // Matrix of J(r) Histograms
  std::vector<std::vector<double>> J_;
  // Matrix of G(r) Histograms
  std::vector<std::vector<double>> G_;
  // Matrix of f(theta) Histograms
  std::vector<std::vector<double>> F_;
  // Matrix of S(Q) Histograms
  std::vector<std::vector<double>> S_;
  // Matrix of calculateXRD Histograms
  std::vector<std::vector<double>> X_;
  // Matrix of Coordination Number Histograms
  std::vector<std::vector<int>> Z_;
  // Matrix of J(r) Smoothed Histograms
  std::vector<std::vector<double>> J_smoothed_;
  // Matrix of g(r) Smoothed Histograms
  std::vector<std::vector<double>> g_smoothed_;
  // Matrix of G(r) Smoothed Histograms
  std::vector<std::vector<double>> G_smoothed_;
  // Matrix of f(theta) Smoothed Histograms
  std::vector<std::vector<double>> F_smoothed_;
  // Matrix of S(Q) Smoothed Histograms
  std::vector<std::vector<double>> S_smoothed_;
  // Matrix of calculateXRD Smoothed Histograms
  std::vector<std::vector<double>> X_smoothed_;

public:
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//

  explicit Cell(const std::array<double, 6> &);
  explicit Cell(const Vector3D &, const Vector3D &, const Vector3D &);
  explicit Cell();

  //-------------------------------------------------------------------------//
  //------------------------------- Accessors -------------------------------//
  //-------------------------------------------------------------------------//

  // Lattice Parameters
  std::array<double, 6> lattice_parameters() { return lattice_parameters_; }
  void setLatticeParameters(std::array<double, 6>);

  // Lattice Vectors
  const Vector3D &v_a() const { return v_a_; };
  const Vector3D &v_b() const { return v_b_; };
  const Vector3D &v_c() const { return v_c_; };
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

  // Distribution Functions
  const std::vector<std::vector<double>> &g() const { return g_; }
  const std::vector<std::vector<double>> &J() const { return J_; }
  const std::vector<std::vector<double>> &G() const { return G_; }
  const std::vector<std::vector<double>> &F() const { return F_; }
  const std::vector<std::vector<double>> &S() const { return S_; }
  const std::vector<std::vector<double>> &X() const { return X_; }
  const std::vector<std::vector<int>> &Z() const { return Z_; }
  const std::vector<std::vector<double>> &J_smoothed() const {
    return J_smoothed_;
  }
  const std::vector<std::vector<double>> &g_smoothed() const {
    return g_smoothed_;
  }
  const std::vector<std::vector<double>> &G_smoothed() const {
    return G_smoothed_;
  }
  const std::vector<std::vector<double>> &F_smoothed() const {
    return F_smoothed_;
  }
  const std::vector<std::vector<double>> &S_smoothed() const {
    return S_smoothed_;
  }
  const std::vector<std::vector<double>> &X_smoothed() const {
    return X_smoothed_;
  }

  //-------------------------------------------------------------------------//
  //------------------------------- Methods ---------------------------------//
  //-------------------------------------------------------------------------//
  // Populate element_id in all  atoms
  void populateElementID();
  // In Cell corrected positions
  void correctPositions();
  void correctFracPositions();
  // Populate the Bond length Matrix
  void populateBondLength(double);
  // Calculate Distances (max distance between atoms)
  void distancePopulation(double, bool);
  // Coordination Numbers Calculation
  void coordinationNumber();
  // Bond-Angle Calulation (degrees=true, radian=false)
  void planeAnglePopulation(bool = true);
  // RDF Histograms  (max distance between atoms, bin width, normalize)
  void calculateRDF(double = 20.0, double = 0.05, bool = false);
  // PAD Histograms (max angle, bin width)
  void calculatePAD(double = 180.0, double = 1.0);
  // Structure Factor Calculation
  void calculateSQ(double = 25.0, double = 0.0, bool = false);
  // calculateXRD Calculation
  void calculateXRD(double = 1.5406, double = 5.0, double = 90.0, double = 1.0);
  // RDF Smoothing (sigma)
  void Smoothing(double, int = 1);
  // VoronoiIndex()
  void voronoiIndex();
};
#endif // INCLUDE_CELL_HPP_
