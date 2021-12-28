#ifndef SRC_CELL_H_
#define SRC_CELL_H_
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
#include <array>
#include <list>
#include <string>
#include <vector>

#include "Atom.h"

class Cell {
  // This object contains the lattice parameters,
  // and the atoms that compose the material

 private:
  // Cell Vectors private modification, public read with getter
  std::array < double, 3 > v_a_;
  std::array < double, 3 > v_b_;
  std::array < double, 3 > v_c_;

 public:
  // Lattice Parameters
  std::array < double, 6 > lattice_parameters;
  // List of atoms in the Cell
  std::list < Atom > atoms;
  // List of elements in Cell
  std::vector < std::string > elements;
  // List with the number of atoms of each kind of elements in Cell
  std::vector < int > element_numbers;
  // List of weights for partials funcions
  std::vector < double > w_ij;
  // Matrix of Bond-lengths
  std::vector < std::vector < double >> bond_length;
  // 3D Tensor of Distances
  std::vector < std::vector < std::vector < double >> > distances;
  // 3D Tensor of Coordination Numbers
  std::vector < std::vector < std::vector < int >> > coordination;
  // 4D Tensor of Angles
  std::vector < std::vector < std::vector < std::vector < double >> >> angles;
  // Matrix of g(r) Histograms
  std::vector < std::vector < double >> g;
  // Matrix of J(r) Histograms
  std::vector < std::vector < double >> J;
  // Matrix of G(r) Histograms
  std::vector < std::vector < double >> G;
  // Matrix of S(Q) Histograms
  std::vector < std::vector < double >> S;
  // Matrix of XRD Histograms
  std::vector < std::vector < double >> X;
  // Matrix of f(theta) Histograms
  std::vector < std::vector < double >> f_theta;
  // Matrix of g(r) Histograms
  std::vector < std::vector < double >> g_smoothed;
  // Volume of the cell
  double volume;

  // Standar Constructors
  explicit Cell(std::array < double, 6 >);
  Cell();

  // Lattice Vectors
  void SetFromVectors(std::vector < double >,
    std::vector < double >,
    std::vector < double >);
  void SetLatticeVectors();
  // In Cell corrected positions
  void CorrectPositions();
  void CorrectFracPositions();
  // Populate the Bond length Matrix
  void PopulateBondLength(double);
  //  Read Bond File
  void ReadBOND(std::string);
  // Update Progress Bar
  void UpdateProgressBar(double);
  // Calculate Distances MultiThreading
  void DistancePopulationMultiThreading(double);
  // Calculate Distances (max distance between atoms)
  void DistancePopulation(double, bool);
  // Coordination Numbers Calculation
  void CoordinationNumber();
  // Bond-Angle Calulation ()
  void PAD(bool = true);
  // RDF Histograms  (max distance between atoms, bin width)
  void RDFHistogram(std::string,
    double = 20.0,
    double = 0.05,
    bool   = false);
  // RDF Smoothing (sigma)
  void RDFSmoothing(std::string, double, int = 1);
  // Nc Histograms ()
  void CoordinationNumberHistogram(std::string);
  // VoronoiIndex()
  void VoronoiIndex(std::string);
  // Structure Factor Calculation
  void SQ(std::string,
    double = 0.1571,
    double = 0.05,
    double = 20.0,
    bool   = false);
  // XRD Calculation
  void XRD(std::string,
    double = 1.5406,
    double = 5.0,
    double = 90.0,
    double = 1.0);
  // PAD Histograms ()
  void PADHistogram(std::string, double = 180.0, double = 1.0);
  // Read Only Lattive Vectors
  std::array < double, 3 > v_a() {
    return v_a_;
  };
  std::array < double, 3 > v_b() {
    return v_b_;
  };
  std::array < double, 3 > v_c() {
    return v_c_;
  };
};

#endif  // SRC_CELL_H_
