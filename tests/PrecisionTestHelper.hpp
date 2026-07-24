#pragma once

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <math/Precision.hpp>

using correlation::real_t;

namespace correlation::testing {

#ifdef CORRELATION_USE_SINGLE_PRECISION
constexpr float kTestTolerance = 1e-4F;
inline auto IsRealEq(real_t expected) { return ::testing::FloatNear(expected, kTestTolerance); }
#else
constexpr double kTestTolerance = 1e-12;
inline auto IsRealEq(real_t expected) { return ::testing::DoubleNear(expected, kTestTolerance); }
#endif

} // namespace correlation::testing