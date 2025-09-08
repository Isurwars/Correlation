#ifndef INCLUDE_SMOOTHING_HPP_
#define INCLUDE_SMOOTHING_HPP_
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <vector>

//---------------------------------------------------------------------------//
//-------------------------------- Methods ----------------------------------//
//---------------------------------------------------------------------------//

std::vector<double> GaussianKernel(const std::vector<double> &r_vals, double r0,
                                   double sigma);

std::vector<double> BumpKernel(const std::vector<double> &r_vals, double r0,
                               double radius);

std::vector<double> TriweightKernel(const std::vector<double> &r_vals,
                                    double r0, double radius);

std::vector<double> KernelSmoothing(const std::vector<double> &r,
                                    const std::vector<double> &y, double sigma,
                                    int _kernel_);

#endif // INCLUDE_SMOOTHING_HPP
