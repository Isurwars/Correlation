#ifndef SRC_CELL_H_
#define SRC_CELL_H_
/* ----------------------------------------------------------------------------
 * Correlation: An Analysis Tool for Liquids and for Amorphous Solids
 * Copyright (c) 2013-2024 Isaías Rodríguez <isurwars@gmail.com>
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
  std::array<double, 3> v_a_;
  std::array<double, 3> v_b_;
  std::array<double, 3> v_c_;
  // Lattice Parameters
  std::array<double, 6> _lattice_parameters_;
  // List of atoms in the Cell
  std::vector<Atom> _atoms_;
  // List of elements in Cell
  std::vector<std::string> _elements_;
  // List with the number of atoms of each kind of elements in Cell
  std::vector<int> _element_numbers_;
  // List of weights for partials funcions
  std::vector<double> w_ij;
  // Matrix of Bond-lengths
  std::vector<std::vector<double>> bond_length;
  // 3D Tensor of Distances
  std::vector<std::vector<std::vector<double>>> distances;
  // 3D Tensor of Coordination Numbers
  std::vector<std::vector<std::vector<int>>> coordination;
  // 4D Tensor of Angles
  std::vector<std::vector<std::vector<std::vector<double>>>> angles;
  // Matrix of g(r) Histograms
  std::vector<std::vector<double>> _g_;
  // Matrix of J(r) Histograms
  std::vector<std::vector<double>> _J_;
  // Matrix of G(r) Histograms
  std::vector<std::vector<double>> _G_;
  // Matrix of S(Q) Histograms
  std::vector<std::vector<double>> _S_;
  // Matrix of XRD Histograms
  std::vector<std::vector<double>> _X_;
  // Matrix of f(theta) Histograms
  std::vector<std::vector<double>> _f_theta_;
  // Matrix of g(r) Histograms
  std::vector<std::vector<double>> _J_smoothed_;
  // Matrix of g(r) Histograms
  std::vector<std::vector<double>> _g_smoothed_;
  // Matrix of g(r) Histograms
  std::vector<std::vector<double>> _G_smoothed_;
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
  void SetLatticeParameters(std::array<double, 6> lat) {
    this->_lattice_parameters_ = lat;
  }
  void SetFromVectors(std::vector<double>, std::vector<double>,
                      std::vector<double>);
  // Lattice Vectors
  const std::array<double, 3> &v_a() { return v_a_; };
  const std::array<double, 3> &v_b() { return v_b_; };
  const std::array<double, 3> &v_c() { return v_c_; };
  void SetLatticeVectors();
  // Atoms
  std::vector<Atom> atoms() { return this->_atoms_; }
  void SetAtoms(const std::vector<Atom> &ats) { this->_atoms_ = ats; }
  void AddAtom(Atom at) { this->_atoms_.push_back(at); }
  // Elements
  std::vector<std::string> elements() { return this->_elements_; }
  std::vector<int> element_numbers() { return this->_element_numbers_; }
  void SetElements(const std::vector<std::string> &ele) {
    this->_elements_ = ele;
  }
  void SetElementsNumbers(const std::vector<int> &ele) {
    this->_element_numbers_ = ele;
  }
  void AddElement(std::string ele) { this->_elements_.push_back(ele); }
  void AddElementNumber(int ele) { this->_element_numbers_.push_back(ele); }
  // Distribution Functions
  std::vector<std::vector<double>> g() { return this->_g_; }
  std::vector<std::vector<double>> J() { return this->_J_; }
  std::vector<std::vector<double>> G() { return this->_G_; }
  std::vector<std::vector<double>> S() { return this->_S_; }
  std::vector<std::vector<double>> X() { return this->_X_; }
  std::vector<std::vector<double>> f_theta() { return this->_f_theta_; }
  std::vector<std::vector<double>> J_smoothed() { return this->_J_smoothed_; }
  std::vector<std::vector<double>> g_smoothed() { return this->_g_smoothed_; }
  std::vector<std::vector<double>> G_smoothed() { return this->_G_smoothed_; }

  //-------------------------------------------------------------------------//
  //------------------------------- Methods ---------------------------------//
  //-------------------------------------------------------------------------//

  // In Cell corrected positions
  void CorrectPositions();
  void CorrectFracPositions();
  // Populate the Bond length Matrix
  void PopulateBondLength(double);
  //  Read Bond File
  void ReadBOND(std::string);
  // Update Progress Bar
  void UpdateProgressBar(double);
  // Calculate Distances (max distance between atoms)
  void DistancePopulation(double, bool);
  // Coordination Numbers Calculation
  void CoordinationNumber();
  // Bond-Angle Calulation ()
  void PAD(bool = true);
  // RDF Histograms  (max distance between atoms, bin width)
  void RDFHistogram(std::string, double = 20.0, double = 0.05, bool = false);
  // RDF Smoothing (sigma)
  void RDFSmoothing(std::string, double, int = 1);
  // Nc Histograms ()
  void CoordinationNumberHistogram(std::string);
  // VoronoiIndex()
  void VoronoiIndex();
  // Structure Factor Calculation
  void SQ(std::string, double = 0.1571, double = 0.05, bool = false);
  // XRD Calculation
  void XRD(std::string, double = 1.5406, double = 5.0, double = 90.0,
           double = 1.0);
  // PAD Histograms ()
  void PADHistogram(std::string, double = 180.0, double = 1.0);
};

#endif // SRC_CELL_H_
