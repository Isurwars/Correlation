#ifndef SRC_CELL_H_
#define SRC_CELL_H_
/* ----------------------------------------------------------------------------
 * Correlation: An Analysis Tool for Liquids and for Amorphous Solids
 * Copyright (c) 2013-2025 Isaías Rodríguez <isurwars@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the MIT License version as published in:
 * https://github.com/Isurwars/Correlation/blob/main/LICENSE
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 * ----------------------------------------------------------------------------
 */
#include <array>
#include <list>
#include <string>
#include <vector>

#include "Atom.hpp"
#include "Templates.hpp" // Include template definitions

class Cell {
  // This object contains the lattice parameters,
  // and the atoms that compose the material

private:
  // Cell Vectors private modification, public read with getter
  std::array<double, 3> _v_a_;
  std::array<double, 3> _v_b_;
  std::array<double, 3> _v_c_;
  // Lattice Parameters
  std::array<double, 6> _lattice_parameters_;
  // List of atoms in the Cell
  std::vector<Atom> _atoms_;
  // List of elements in Cell
  std::vector<std::string> _elements_;
  // List with the number of atoms of each kind of elements in Cell
  std::vector<int> _element_numbers_;
  // List of weights for partials funcions
  std::vector<double> _w_ij_;
  // Matrix of Bond-lengths
  std::vector<std::vector<double>> _bond_length_;
  // 3D Tensor of Distances
  std::vector<std::vector<std::vector<double>>> _distances_;
  // 4D Tensor of Angles
  std::vector<std::vector<std::vector<std::vector<double>>>> _angles_;
  // Matrix of g(r) Histograms
  std::vector<std::vector<double>> _g_;
  // Matrix of J(r) Histograms
  std::vector<std::vector<double>> _J_;
  // Matrix of G(r) Histograms
  std::vector<std::vector<double>> _G_;
  // Matrix of f(theta) Histograms
  std::vector<std::vector<double>> _F_;
  // Matrix of S(Q) Histograms
  std::vector<std::vector<double>> _S_;
  // Matrix of XRD Histograms
  std::vector<std::vector<double>> _X_;
  // Matrix of Coordination Number Histograms
  std::vector<std::vector<int>> _Z_;
  // Matrix of J(r) Smoothed Histograms
  std::vector<std::vector<double>> _J_smoothed_;
  // Matrix of g(r) Smoothed Histograms
  std::vector<std::vector<double>> _g_smoothed_;
  // Matrix of G(r) Smoothed Histograms
  std::vector<std::vector<double>> _G_smoothed_;
  // Matrix of f(theta) Smoothed Histograms
  std::vector<std::vector<double>> _F_smoothed_;
  // Matrix of S(Q) Smoothed Histograms
  std::vector<std::vector<double>> _S_smoothed_;
  // Matrix of XRD Smoothed Histograms
  std::vector<std::vector<double>> _X_smoothed_;
  // Volume of the cell
  double volume;

public:
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//

  explicit Cell(std::array<double, 6>);
  Cell();

  //-------------------------------------------------------------------------//
  //--------------------------- Setters & Getters ---------------------------//
  //-------------------------------------------------------------------------//
  // Lattice Parameters
  std::array<double, 6> lattice_parameters() {
    return this->_lattice_parameters_;
  }
  void setFromVectors(std::vector<double>, std::vector<double>,
		      std::vector<double>);
  // Lattice Vectors
  const std::array<double, 3> &v_a() { return _v_a_; };
  const std::array<double, 3> &v_b() { return _v_b_; };
  const std::array<double, 3> &v_c() { return _v_c_; };
  void setLatticeVectors();
  void setLatticeParameters(std::array<double, 6>);
  // Volume
  double getVolume() { return this->volume; }
  // Atoms
  std::vector<Atom> atoms() { return this->_atoms_; }
  void addAtom(Atom);
  void setAtoms(std::vector<Atom>);
  int getElementID(std::string);
  // Elements
  std::vector<std::string> elements() { return this->_elements_; }
  std::vector<int> element_numbers() { return this->_element_numbers_; }
  void setElements(const std::vector<std::string> &ele) {
    this->_elements_ = ele;
  }
  void setElementsNumbers(const std::vector<int> &ele) {
    this->_element_numbers_ = ele;
  }
  void addElement(std::string ele) { this->_elements_.push_back(ele); }
  void addElementNumber(int ele) { this->_element_numbers_.push_back(ele); }
  // Distances
  int distancesSize() { return this->_distances_.size(); }
  // Bonds
  std::vector<std::vector<double>> bond_length() { return this->_bond_length_; }
  void setBondLength(const std::vector<std::vector<double>> &bonds) {
    this->_bond_length_ = bonds;
  }
  // Distribution Functions
  std::vector<std::vector<double>> g() { return this->_g_; }
  std::vector<std::vector<double>> J() { return this->_J_; }
  std::vector<std::vector<double>> G() { return this->_G_; }
  std::vector<std::vector<double>> F() { return this->_F_; }
  std::vector<std::vector<double>> S() { return this->_S_; }
  std::vector<std::vector<double>> X() { return this->_X_; }
  std::vector<std::vector<int>> Z() { return this->_Z_; }
  std::vector<std::vector<double>> J_smoothed() { return this->_J_smoothed_; }
  std::vector<std::vector<double>> g_smoothed() { return this->_g_smoothed_; }
  std::vector<std::vector<double>> G_smoothed() { return this->_G_smoothed_; }
  std::vector<std::vector<double>> F_smoothed() { return this->_F_smoothed_; }
  std::vector<std::vector<double>> S_smoothed() { return this->_S_smoothed_; }
  std::vector<std::vector<double>> X_smoothed() { return this->_X_smoothed_; }

  //-------------------------------------------------------------------------//
  //------------------------------- Methods ---------------------------------//
  //-------------------------------------------------------------------------//
  // Populate element_id in all  atoms
  void populateElementID();
  // Populate element_numbers
  void populateElementNumbers();
  // In Cell corrected positions
  void correctPositions();
  void correctFracPositions();
  // Populate the Bond length Matrix
  void populateBondLength(double);
  // Update Progress Bar
  void updateProgressBar(double);
  // Calculate Distances (max distance between atoms)
  void distancePopulation(double, bool);
  // Coordination Numbers Calculation
  void coordinationNumber();
  // Bond-Angle Calulation (degrees=true, radian=false)
  void planeAnglePopulation(bool = true);
  // RDF Histograms  (max distance between atoms, bin width, normalize)
  void radialDistributionFunctions(double = 20.0, double = 0.05, bool = false);
  // PAD Histograms (max angle, bin width)
  void planeAngleDistribution(double = 180.0, double = 1.0);
  // Structure Factor Calculation
  void SQ(double = 25.0, double = 0.0, bool = false);
  // void SQExact(double = 0.05, double = 0.05, bool = false);
  // XRD Calculation
  void XRD(double = 1.5406, double = 5.0, double = 90.0, double = 1.0);
  // RDF Smoothing (sigma)
  void Smoothing(double, int = 1);
  // VoronoiIndex()
  void voronoiIndex();
};
#endif // SRC_CELL_H_
