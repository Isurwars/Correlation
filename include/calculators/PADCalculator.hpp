// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "DistributionFunctions.hpp"

class PADCalculator {
public:
  static Histogram calculate(const Cell &cell,
                             const StructureAnalyzer *neighbors,
                             double bin_width);
};
