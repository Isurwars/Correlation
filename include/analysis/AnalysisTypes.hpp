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
  real_t min_sq{static_cast<real_t>(0.0)}; ///< Squared minimum bond distance (Å²).
  real_t max_sq{static_cast<real_t>(0.0)}; ///< Squared maximum bond distance (Å²).

  [[nodiscard]] constexpr bool isValid() const noexcept { return min_sq >= 0.0 && max_sq >= 0.0 && min_sq <= max_sq; }
};

/** @brief Type definition for a 2D matrix of element-pair bond cutoff ranges. */
using BondCutoffMatrix = std::vector<std::vector<BondCutoffRange>>;

} // namespace correlation::analysis
