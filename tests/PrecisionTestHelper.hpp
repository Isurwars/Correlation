#pragma once
#include <gtest/gtest.h>

#ifdef CORRELATION_USE_SINGLE_PRECISION
#undef EXPECT_DOUBLE_EQ
#define EXPECT_DOUBLE_EQ(val1, val2) EXPECT_NEAR(val1, val2, 1e-4)
#endif
