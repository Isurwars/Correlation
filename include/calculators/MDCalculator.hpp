// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "DistributionFunctions.hpp"

class MDCalculator {
public:
  static Histogram calculate(const NeighborGraph &graph, size_t max_ring_size);
};
