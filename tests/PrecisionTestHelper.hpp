#pragma once
#include <math/Precision.hpp>
#include <gtest/gtest.h>

using correlation::real_t;

#ifdef CORRELATION_USE_SINGLE_PRECISION
#undef EXPECT_DOUBLE_EQ
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define EXPECT_DOUBLE_EQ(val1, val2) EXPECT_NEAR(val1, val2, 1e-4)
#endif
