#ifndef INCLUDE_DISTRIBUTION_FUNCTIONS_HPP_
#define INCLUDE_DISTRIBUTION_FUNCTIONS_HPP_
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "../include/Cell.hpp"

class DistributionFunctions {
private:
  // Cell
  Cell cell_;
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
  std::vector<std::vector<double>> Z_;
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

  explicit DistributionFunctions(const Cell &);

  //-------------------------------------------------------------------------//
  //------------------------------- Accessors -------------------------------//
  //-------------------------------------------------------------------------//

  // Cell
  const Cell &cell() const { return cell_; }
  void setCell(const Cell &c) { cell_ = c; }

  // Distribution Functions
  const std::vector<std::vector<double>> &g() const { return g_; }
  const std::vector<std::vector<double>> &J() const { return J_; }
  const std::vector<std::vector<double>> &G() const { return G_; }
  const std::vector<std::vector<double>> &F() const { return F_; }
  const std::vector<std::vector<double>> &S() const { return S_; }
  const std::vector<std::vector<double>> &X() const { return X_; }
  const std::vector<std::vector<double>> &Z() const { return Z_; }
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

  // Coordination Number
  void coordinationNumber();
  // RDF Histograms  (max distance between atoms, bin width, normalize)
  void calculateRDF(double = 20.0, double = 0.05, bool = false);
  // PAD Histograms (max angle, bin width)
  void calculatePAD(double = 180.0, double = 1.0);
  // Structure Factor Calculation
  void calculateSQ(double = 25.0, double = 0.0, double = 8.0, bool = false);
  // Structure Factor Calculation
  void calculateSQDirac(double = 25.0, double = 0.0, bool = false);
  // calculateXRD Calculation
  void calculateXRD(double = 1.5406, double = 5.0, double = 90.0, double = 1.0);
  // RDF Smoothing (sigma)
  void Smoothing(double, int = 1);
  // VoronoiIndex()
  void voronoiIndex();
};
#endif // INCLUDE_DISTRIBUTION_FUNCTIONS_HPP_
