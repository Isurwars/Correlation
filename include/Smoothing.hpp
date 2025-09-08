#ifndef INCLUDE_SMOOTHING_HPP_
#define INCLUDE_SMOOTHING_HPP_
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <vector>

// Enum class for type-safe selection of kernel types.
// This prevents passing invalid integer values.
enum class KernelType { Gaussian, Bump, Triweight };

/**
 * @brief Applies kernel smoothing to a data series.
 *
 * This function convolves a given kernel with the input data `y` to produce a
 * smoothed version.
 *
 * @param r The independent variable values (x-axis).
 * @param y The dependent variable values (y-axis) to be smoothed.
 * @param sigma The bandwidth (standard deviation for Gaussian) of the kernel.
 * @param type The type of kernel to use for smoothing.
 * @return A vector containing the smoothed data.
 */
std::vector<double> KernelSmoothing(const std::vector<double> &r,
                                    const std::vector<double> &y, double sigma,
                                    KernelType type = KernelType::Gaussian);

#endif // INCLUDE_SMOOTHING_HPP
