// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "app/SettingsManager.hpp"
#include <filesystem>
#include <gtest/gtest.h>

namespace {

class SettingsManagerTests : public ::testing::Test {
protected:
  void SetUp() override {
    test_path_ = correlation::app::SettingsManager::settingsFilePath();
    if (std::filesystem::exists(test_path_)) {
      original_settings_ = correlation::app::SettingsManager::load();
    }
  }

  void TearDown() override {
    // Restore original settings if existed
    if (original_settings_.window_width > 0) {
      correlation::app::SettingsManager::save(original_settings_);
    }
  }

private:
  std::filesystem::path test_path_;
  correlation::app::AppSettings original_settings_;
};

TEST_F(SettingsManagerTests, DefaultSettingsFallback) {
  correlation::app::AppSettings settings;
  EXPECT_EQ(settings.window_width, 1180U);
  EXPECT_EQ(settings.window_height, 795U);
  EXPECT_FLOAT_EQ(settings.left_col_width, 260.0F);
  EXPECT_FLOAT_EQ(settings.middle_col_width, 260.0F);
}

TEST_F(SettingsManagerTests, RoundTripSerialization) {
  correlation::app::AppSettings custom;
  custom.window_width = 1440U;
  custom.window_height = 900U;
  custom.left_col_width = 315.5F;
  custom.middle_col_width = 280.0F;

  ASSERT_TRUE(correlation::app::SettingsManager::save(custom));

  correlation::app::AppSettings loaded = correlation::app::SettingsManager::load();
  EXPECT_EQ(loaded.window_width, 1440U);
  EXPECT_EQ(loaded.window_height, 900U);
  EXPECT_FLOAT_EQ(loaded.left_col_width, 315.5F);
  EXPECT_FLOAT_EQ(loaded.middle_col_width, 280.0F);
}

} // namespace
