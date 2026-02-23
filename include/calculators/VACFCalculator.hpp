// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "DistributionFunctions.hpp"
#include "Trajectory.hpp"
#include <map>
#include <string>

class VACFCalculator {
public:
  static std::map<std::string, Histogram> calculate(const Trajectory &traj,
                                                    int max_correlation_frames,
                                                    size_t start_frame,
                                                    size_t end_frame);
};
