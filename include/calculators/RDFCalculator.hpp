// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "DistributionFunctions.hpp"
#include <map>
#include <string>

class RDFCalculator {
public:
  // Returns a map containing J(r), g(r), and G(r) histograms.
  static std::map<std::string, Histogram>
  calculate(const Cell &cell, const StructureAnalyzer *neighbors,
            const std::map<std::string, double> &ashcroft_weights, double r_max,
            double r_bin_width);
};
