#pragma once

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <math/Precision.hpp>

using correlation::real_t;

namespace correlation::testing {

inline constexpr real_t kTestTolerance =
    correlation::is_single_precision ? static_cast<real_t>(1e-4) : static_cast<real_t>(1e-12);

inline auto IsRealEq(real_t expected) {
  if constexpr (correlation::is_single_precision) {
    return ::testing::FloatNear(expected, kTestTolerance);
  } else {
    return ::testing::DoubleNear(expected, kTestTolerance);
  }
}

} // namespace correlation::testing