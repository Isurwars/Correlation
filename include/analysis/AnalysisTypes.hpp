#pragma once

#include "math/Precision.hpp"
#include <cstddef>
#include <vector>

namespace correlation::analysis {

struct MaxFrames {
  int value;
};

struct StartFrame {
  size_t value;
};

struct EndFrame {
  size_t value;
};

/**
 * @struct BondCutoffRange
 * @brief Represents the squared minimum and maximum bond distance bounds for an element pair.
 */
struct BondCutoffRange {
  real_t min_sq{static_cast<real_t>(0.36)}; ///< Squared minimum bond distance (default 0.60 Å -> 0.36 Å²).
  real_t max_sq{static_cast<real_t>(0.0)};  ///< Squared maximum bond distance (Å²).
};

/** @brief Type definition for a 2D matrix of element-pair bond cutoff ranges. */
using BondCutoffMatrix = std::vector<std::vector<BondCutoffRange>>;

} // namespace correlation::analysis
