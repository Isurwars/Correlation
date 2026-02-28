// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "../Cell.hpp"
#include "../NeighborGraph.hpp"
#include <vector>

namespace calculators {

// [Element A] [Element B] [Element C] [Element D] -> List of Angles
using DihedralTensor =
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>;

class DihedralCalculator {
public:
  static void compute(const Cell &cell, const NeighborGraph &graph,
                      DihedralTensor &out_dihedrals);
};

} // namespace calculators
