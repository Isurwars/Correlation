#ifndef INCLUDE_STRUCTURE_FACTOR_HPP_
#define INCLUDE_STRUCTURE_FACTOR_HPP_
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <string>
#include <vector>

//---------------------------------------------------------------------------//
//------------------------------- Methods -----------------------------------//
//---------------------------------------------------------------------------//

// Sinc Function, sin(x)/x
double sinc(double);
// Gaussian Function a*exp (-b * (x / 4Pi)**2)
double gaussian(double, double, double);
// Atomic Form Factor
double atomicFormFactor(double, std::string_view);
// Average Scattering Factor
double avrgAtomicScattering(const std::vector<double>, std::string_view);
double avrgScatteringFactor(std::vector<double>, std::vector<std::string>);
#endif // INCLUDE_STRUCTURE_FACTOR_HPP_
