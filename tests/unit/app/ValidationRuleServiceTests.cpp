// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "app/ValidationRuleService.hpp"
#include <gtest/gtest.h>

namespace {

using correlation::app::ValidationRuleService;

TEST(ValidationRuleServiceTests, ParsePositiveFloat) {
  auto res1 = ValidationRuleService::parsePositiveFloat("1.25");
  ASSERT_TRUE(res1.has_value());
  EXPECT_FLOAT_EQ(*res1, 1.25F);

  auto res2 = ValidationRuleService::parsePositiveFloat("0.0");
  EXPECT_FALSE(res2.has_value());

  auto res3 = ValidationRuleService::parsePositiveFloat("-2.5");
  EXPECT_FALSE(res3.has_value());

  auto res4 = ValidationRuleService::parsePositiveFloat("abc");
  EXPECT_FALSE(res4.has_value());

  auto res5 = ValidationRuleService::parsePositiveFloat("1.5abc");
  EXPECT_FALSE(res5.has_value());

  auto res6 = ValidationRuleService::parsePositiveFloat("");
  EXPECT_FALSE(res6.has_value());
}

TEST(ValidationRuleServiceTests, ParseNonNegativeFloat) {
  auto res1 = ValidationRuleService::parseNonNegativeFloat("0.0");
  ASSERT_TRUE(res1.has_value());
  EXPECT_FLOAT_EQ(*res1, 0.0F);

  auto res2 = ValidationRuleService::parseNonNegativeFloat("3.14");
  ASSERT_TRUE(res2.has_value());
  EXPECT_FLOAT_EQ(*res2, 3.14F);

  auto res3 = ValidationRuleService::parseNonNegativeFloat("-0.01");
  EXPECT_FALSE(res3.has_value());
}

TEST(ValidationRuleServiceTests, ParsePositiveInt) {
  auto res1 = ValidationRuleService::parsePositiveInt("42");
  ASSERT_TRUE(res1.has_value());
  EXPECT_EQ(*res1, 42);

  auto res2 = ValidationRuleService::parsePositiveInt("0");
  EXPECT_FALSE(res2.has_value());

  auto res3 = ValidationRuleService::parsePositiveInt("-5");
  EXPECT_FALSE(res3.has_value());

  auto res4 = ValidationRuleService::parsePositiveInt("1.5");
  EXPECT_FALSE(res4.has_value());
}

TEST(ValidationRuleServiceTests, ParseMinFrame) {
  EXPECT_EQ(ValidationRuleService::parseMinFrame("start", 10).value(), 0);
  EXPECT_EQ(ValidationRuleService::parseMinFrame("Start", 10).value(), 0);
  EXPECT_EQ(ValidationRuleService::parseMinFrame("", 10).value(), 0);
  EXPECT_EQ(ValidationRuleService::parseMinFrame("end", 10).value(), 9);
  EXPECT_EQ(ValidationRuleService::parseMinFrame("1", 10).value(), 0);
  EXPECT_EQ(ValidationRuleService::parseMinFrame("5", 10).value(), 4);

  // Out of bounds (> total_frames)
  auto err = ValidationRuleService::parseMinFrame("15", 10);
  EXPECT_FALSE(err.has_value());
  EXPECT_EQ(err.error(), "Must be ≤ total frames (10)");

  // Invalid text
  auto err2 = ValidationRuleService::parseMinFrame("invalid", 10);
  EXPECT_FALSE(err2.has_value());
}

TEST(ValidationRuleServiceTests, ParseMaxFrame) {
  EXPECT_EQ(ValidationRuleService::parseMaxFrame("end", 10).value(), 9);
  EXPECT_EQ(ValidationRuleService::parseMaxFrame("End", 10).value(), 9);
  EXPECT_EQ(ValidationRuleService::parseMaxFrame("", 10).value(), 9);
  EXPECT_EQ(ValidationRuleService::parseMaxFrame("end", 0).value(), -1);
  EXPECT_EQ(ValidationRuleService::parseMaxFrame("start", 10).value(), 0);
  EXPECT_EQ(ValidationRuleService::parseMaxFrame("8", 10).value(), 7);

  // Out of bounds
  auto err = ValidationRuleService::parseMaxFrame("12", 10);
  EXPECT_FALSE(err.has_value());
  EXPECT_EQ(err.error(), "Must be ≤ total frames (10)");
}

TEST(ValidationRuleServiceTests, ParseFrameStride) {
  EXPECT_EQ(ValidationRuleService::parseFrameStride("").value(), 1);
  EXPECT_EQ(ValidationRuleService::parseFrameStride("1").value(), 1);
  EXPECT_EQ(ValidationRuleService::parseFrameStride("5").value(), 5);

  EXPECT_FALSE(ValidationRuleService::parseFrameStride("0").has_value());
  EXPECT_FALSE(ValidationRuleService::parseFrameStride("-2").has_value());
  EXPECT_FALSE(ValidationRuleService::parseFrameStride("xyz").has_value());
}

TEST(ValidationRuleServiceTests, ValidateBinWithinMax) {
  EXPECT_TRUE(ValidationRuleService::validateBinWithinMax(0.02F, 20.0F, "r_max").has_value());
  auto err = ValidationRuleService::validateBinWithinMax(25.0F, 20.0F, "r_max");
  EXPECT_FALSE(err.has_value());
  EXPECT_EQ(err.error(), "Must be ≤ r_max");
}

TEST(ValidationRuleServiceTests, ValidateAngleDegrees) {
  EXPECT_TRUE(ValidationRuleService::validateAngleDegrees(45.0F).has_value());
  EXPECT_TRUE(ValidationRuleService::validateAngleDegrees(180.0F).has_value());
  auto err = ValidationRuleService::validateAngleDegrees(180.1F);
  EXPECT_FALSE(err.has_value());
  EXPECT_EQ(err.error(), "Must be ≤ 180°");
}

TEST(ValidationRuleServiceTests, ValidateXrdTheta) {
  EXPECT_TRUE(ValidationRuleService::validateXrdTheta(10.0F, 80.0F).has_value());

  auto err1 = ValidationRuleService::validateXrdTheta(180.0F, 190.0F);
  EXPECT_FALSE(err1.has_value());
  EXPECT_EQ(err1.error(), "Must be < 180°");

  auto err2 = ValidationRuleService::validateXrdTheta(20.0F, 190.0F);
  EXPECT_FALSE(err2.has_value());
  EXPECT_EQ(err2.error(), "Must be ≤ 180°");

  auto err3 = ValidationRuleService::validateXrdTheta(50.0F, 40.0F);
  EXPECT_FALSE(err3.has_value());
  EXPECT_EQ(err3.error(), "Must be > Min 2θ");
}

} // namespace
