// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "../Cell.hpp"
#include "../NeighborGraph.hpp"
#include <vector>

namespace calculators {

using AngleTensor = std::vector<std::vector<std::vector<std::vector<double>>>>;

class AngleCalculator {
public:
  static void compute(const Cell &cell, const NeighborGraph &graph,
                      AngleTensor &out_angles);
};

} // namespace calculators
