/**
 * @file ColormapServiceTests.cpp
 * @brief Unit tests for scientific colormap engine and precomputed LUTs.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "plotters/ColormapService.hpp"

#include <gtest/gtest.h>
#include <limits>

namespace correlation::plotters::testing {

namespace {

[[nodiscard]] double computeLuminance(const RgbaColor &color) noexcept {
  return 0.2126 * static_cast<double>(color.r) + 0.7152 * static_cast<double>(color.g) +
         0.0722 * static_cast<double>(color.b);
}

struct SampleTestCase {
  double value;
  size_t expected_index;
};

void verifyBoundaryClamping(ColormapType type) {
  const auto lut = ColormapService::getLut(type);
  constexpr std::array<SampleTestCase, 7> K_CASES = {
      SampleTestCase{.value = -1.0, .expected_index = 0},
      SampleTestCase{.value = -0.0001, .expected_index = 0},
      SampleTestCase{.value = 0.0, .expected_index = 0},
      SampleTestCase{.value = 0.5, .expected_index = 128},
      SampleTestCase{.value = 1.0, .expected_index = 255},
      SampleTestCase{.value = 1.0001, .expected_index = 255},
      SampleTestCase{.value = 100.0, .expected_index = 255},
  };
  for (const auto &test_case : K_CASES) {
    EXPECT_EQ(ColormapService::sample(test_case.value, type), lut[test_case.expected_index]);
  }
  const double nan_val = std::numeric_limits<double>::quiet_NaN();
  EXPECT_EQ(ColormapService::sample(nan_val, type), lut[0]);
}

void verifyLutSpan(ColormapType type) {
  const auto lut = ColormapService::getLut(type);
  EXPECT_EQ(lut.size(), 256U);
  for (const auto &color : lut) {
    EXPECT_EQ(color.a, 255U);
  }
}

} // namespace

TEST(ColormapServiceTests, BoundaryChecksClampCorrectly) {
  verifyBoundaryClamping(ColormapType::Viridis);
  verifyBoundaryClamping(ColormapType::Plasma);
  verifyBoundaryClamping(ColormapType::Coolwarm);
  verifyBoundaryClamping(ColormapType::Turbo);
}

TEST(ColormapServiceTests, LutSpansHaveValidLengthAndAlpha) {
  verifyLutSpan(ColormapType::Viridis);
  verifyLutSpan(ColormapType::Plasma);
  verifyLutSpan(ColormapType::Coolwarm);
  verifyLutSpan(ColormapType::Turbo);
}

TEST(ColormapServiceTests, ViridisReferenceColorsMatchExpected) {
  const auto color_min = ColormapService::sample(0.0, ColormapType::Viridis);
  EXPECT_EQ(color_min.r, 68U);
  EXPECT_EQ(color_min.g, 1U);
  EXPECT_EQ(color_min.b, 84U);

  const auto color_mid = ColormapService::sample(0.5, ColormapType::Viridis);
  EXPECT_EQ(color_mid.r, 33U);
  EXPECT_EQ(color_mid.g, 145U);
  EXPECT_EQ(color_mid.b, 140U);

  const auto color_max = ColormapService::sample(1.0, ColormapType::Viridis);
  EXPECT_EQ(color_max.r, 253U);
  EXPECT_EQ(color_max.g, 231U);
  EXPECT_EQ(color_max.b, 37U);
}

TEST(ColormapServiceTests, PlasmaReferenceColorsMatchExpected) {
  const auto color_min = ColormapService::sample(0.0, ColormapType::Plasma);
  EXPECT_EQ(color_min.r, 13U);
  EXPECT_EQ(color_min.g, 8U);
  EXPECT_EQ(color_min.b, 135U);

  const auto color_mid = ColormapService::sample(0.5, ColormapType::Plasma);
  EXPECT_EQ(color_mid.r, 204U);
  EXPECT_EQ(color_mid.g, 71U);
  EXPECT_EQ(color_mid.b, 120U);

  const auto color_max = ColormapService::sample(1.0, ColormapType::Plasma);
  EXPECT_EQ(color_max.r, 240U);
  EXPECT_EQ(color_max.g, 249U);
  EXPECT_EQ(color_max.b, 33U);
}

TEST(ColormapServiceTests, CoolwarmReferenceColorsMatchExpected) {
  const auto color_min = ColormapService::sample(0.0, ColormapType::Coolwarm);
  EXPECT_EQ(color_min.r, 59U);
  EXPECT_EQ(color_min.g, 76U);
  EXPECT_EQ(color_min.b, 192U);

  const auto color_mid = ColormapService::sample(0.5, ColormapType::Coolwarm);
  EXPECT_EQ(color_mid.r, 221U);
  EXPECT_EQ(color_mid.g, 220U);
  EXPECT_EQ(color_mid.b, 220U);

  const auto color_max = ColormapService::sample(1.0, ColormapType::Coolwarm);
  EXPECT_EQ(color_max.r, 180U);
  EXPECT_EQ(color_max.g, 4U);
  EXPECT_EQ(color_max.b, 38U);
}

TEST(ColormapServiceTests, TurboReferenceColorsMatchExpected) {
  const auto color_min = ColormapService::sample(0.0, ColormapType::Turbo);
  EXPECT_EQ(color_min.r, 48U);
  EXPECT_EQ(color_min.g, 18U);
  EXPECT_EQ(color_min.b, 59U);

  const auto color_mid = ColormapService::sample(0.5, ColormapType::Turbo);
  EXPECT_EQ(color_mid.r, 164U);
  EXPECT_EQ(color_mid.g, 252U);
  EXPECT_EQ(color_mid.b, 60U);

  const auto color_max = ColormapService::sample(1.0, ColormapType::Turbo);
  EXPECT_EQ(color_max.r, 122U);
  EXPECT_EQ(color_max.g, 4U);
  EXPECT_EQ(color_max.b, 3U);
}

TEST(ColormapServiceTests, MonotonicLuminanceInViridisAndPlasma) {
  for (const auto type : {ColormapType::Viridis, ColormapType::Plasma}) {
    const auto lut = ColormapService::getLut(type);
    const double lum_start = computeLuminance(lut.front());
    const double lum_end = computeLuminance(lut.back());

    // Both Viridis and Plasma transition from dark to bright
    EXPECT_GT(lum_end, lum_start + 150.0);

    // Verify monotonicity across entries (tolerance for 8-bit quantization)
    for (size_t idx = 1; idx < lut.size(); ++idx) {
      const double prev_lum = computeLuminance(lut[idx - 1]);
      const double curr_lum = computeLuminance(lut[idx]);
      EXPECT_GE(curr_lum, prev_lum - 0.75)
          << "Luminance dip detected at index " << idx << " for " << ColormapService::name(type);
    }
  }
}

TEST(ColormapServiceTests, CoolwarmDivergingBehavior) {
  const auto color_min = ColormapService::sample(0.0, ColormapType::Coolwarm);
  EXPECT_GT(color_min.b, color_min.r);
  EXPECT_GT(color_min.b, color_min.g);

  const auto color_mid = ColormapService::sample(0.5, ColormapType::Coolwarm);
  EXPECT_NEAR(static_cast<double>(color_mid.r), static_cast<double>(color_mid.g), 2.0);
  EXPECT_NEAR(static_cast<double>(color_mid.r), static_cast<double>(color_mid.b), 2.0);

  const auto color_max = ColormapService::sample(1.0, ColormapType::Coolwarm);
  EXPECT_GT(color_max.r, color_max.b);
  EXPECT_GT(color_max.r, color_max.g);
}

TEST(ColormapServiceTests, NameAndFromStringRoundtrip) {
  EXPECT_EQ(ColormapService::name(ColormapType::Viridis), "Viridis");
  EXPECT_EQ(ColormapService::name(ColormapType::Plasma), "Plasma");
  EXPECT_EQ(ColormapService::name(ColormapType::Coolwarm), "Coolwarm");
  EXPECT_EQ(ColormapService::name(ColormapType::Turbo), "Turbo");

  EXPECT_EQ(ColormapService::fromString("Viridis"), ColormapType::Viridis);
  EXPECT_EQ(ColormapService::fromString("viridis"), ColormapType::Viridis);
  EXPECT_EQ(ColormapService::fromString("VIRIDIS"), ColormapType::Viridis);

  EXPECT_EQ(ColormapService::fromString("Plasma"), ColormapType::Plasma);
  EXPECT_EQ(ColormapService::fromString("plasma"), ColormapType::Plasma);

  EXPECT_EQ(ColormapService::fromString("Coolwarm"), ColormapType::Coolwarm);
  EXPECT_EQ(ColormapService::fromString("COOLWARM"), ColormapType::Coolwarm);

  EXPECT_EQ(ColormapService::fromString("Turbo"), ColormapType::Turbo);
  EXPECT_EQ(ColormapService::fromString("turbo"), ColormapType::Turbo);

  EXPECT_EQ(ColormapService::fromString(""), std::nullopt);
  EXPECT_EQ(ColormapService::fromString("unknown"), std::nullopt);
  EXPECT_EQ(ColormapService::fromString("jet"), std::nullopt);
}

} // namespace correlation::plotters::testing
