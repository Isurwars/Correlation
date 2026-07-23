#pragma once
#include <gtest/gtest.h>
#include <math/Precision.hpp>

using correlation::real_t;

namespace correlation::testing {

#ifdef CORRELATION_USE_SINGLE_PRECISION
constexpr double kTestTolerance = 1e-4;
#else
constexpr double kTestTolerance = 1e-12;
#endif

} // namespace correlation::testing

#ifdef CORRELATION_USE_SINGLE_PRECISION
#undef EXPECT_DOUBLE_EQ
#define EXPECT_DOUBLE_EQ(val1, val2) EXPECT_NEAR(val1, val2, 1e-4) // NOLINT(cppcoreguidelines-macro-usage)
#endif