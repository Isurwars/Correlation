// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "../include/StructureFactor.hpp"

#include <cmath>

#include "../include/Constants.hpp"

//---------------------------------------------------------------------------//
//---------------------------------- Sinc -----------------------------------//
//---------------------------------------------------------------------------//
double sinc(double x) {
  if (x == 0.0) {
    return 1.0; // sinc(0) = 1
  }
  return sin(x) / x;
}

//---------------------------------------------------------------------------//
//--------------------------------- Gaussian --------------------------------//
//---------------------------------------------------------------------------//
double gaussian(double x, double a, double b) {
  double aux = x / (4 * constants::pi);
  return a * exp(-1.0 * b * aux * aux);
}

//---------------------------------------------------------------------------//
//--------------------------- Alpha Weight Factor ---------------------------//
//---------------------------------------------------------------------------//
double alphaWeightFactor(double x, double a, double b) {
  double aux = x / (4 * constants::pi);
  return a * exp(-1.0 * b * aux * aux);
}

//---------------------------------------------------------------------------//
//--------------------------- Atomic Form Factor  ---------------------------//
//---------------------------------------------------------------------------//
double atomicFormFactor(double q, std::string &element) {
  std::vector<double> param = atomicFormFactorParameters(element);
  double sum = 0.0;
  for (int i = 0; i < 4; ++i) {
    sum += gaussian(q, param[i * 2], param[i * 2 + 1]);
  }
  return sum + param[8];
}

//---------------------------------------------------------------------------//
//-------------------- Average Atomic Scattering Factor ---------------------//
//---------------------------------------------------------------------------//
double avrgAtomicScattering(const std::vector<double> &q_values,
          std::string &element) {
  const size_t num_points = q_values.size();
  if (num_points == 0)
    return 0.0; // Handle empty input

  const std::vector<double> parameters = atomicFormFactorParameters(element);

  double total_sum = 0.0;

  for (double q : q_values) {
    double sum = parameters[8]; // c parameter

    // Calculate Gaussian terms (a1-4 * exp(-b1-4 * q²))
    for (size_t i = 0; i < 4; ++i) {
      const double a = parameters[i * 2];
      const double b = parameters[i * 2 + 1];
      sum += gaussian(q, a, b);
    }

    total_sum += sum;
  }

  return total_sum / num_points;
}

//---------------------------------------------------------------------------//
//------------------------ Average Scattering Factor ------------------------//
//---------------------------------------------------------------------------//
double avrgScatteringFactor(std::vector<double> q_,
          std::vector<std::string> elements) {
  int n_ = elements.size();
  double aux = 0.0;
  for (int i = 0; i < n_; ++i) {
    // aux += avrgAtomicScattering(q_, elements[i]);
  }
  return aux / n_;
}
