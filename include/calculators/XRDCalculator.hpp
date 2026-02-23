// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "DistributionFunctions.hpp"
#include <map>
#include <string>

class XRDCalculator {
public:
  static Histogram
  calculate(const Histogram &g_r_hist, const Cell &cell,
            const std::map<std::string, double> &ashcroft_weights,
            double lambda, double theta_min, double theta_max,
            double bin_width);
};
