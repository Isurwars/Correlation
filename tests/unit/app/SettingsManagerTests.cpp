// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "app/core/SettingsManager.hpp"
#include <chrono>
#include <filesystem>
#include <gtest/gtest.h>

namespace {

class SettingsManagerTests : public ::testing::Test {
protected:
  void SetUp() override {
    const auto *test_info = ::testing::UnitTest::GetInstance()->current_test_info();
    const std::string suffix =
        (test_info != nullptr ? std::string(test_info->name()) : std::string("test")) + "_" +
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
    temp_dir_ = std::filesystem::temp_directory_path() / ("corr_test_settings_" + suffix);
    std::filesystem::create_directories(temp_dir_);
    test_path_ = temp_dir_ / "settings.json";
#ifdef _WIN32
    _putenv_s("CORRELATION_SETTINGS_PATH", test_path_.string().c_str());
#else
    setenv("CORRELATION_SETTINGS_PATH", test_path_.c_str(), 1);
#endif
  }

  void TearDown() override {
#ifdef _WIN32
    _putenv_s("CORRELATION_SETTINGS_PATH", "");
#else
    unsetenv("CORRELATION_SETTINGS_PATH");
#endif
    std::error_code error_code;
    std::filesystem::remove_all(temp_dir_, error_code);
  }

private:
  std::filesystem::path temp_dir_;
  std::filesystem::path test_path_;
};

TEST_F(SettingsManagerTests, DefaultSettingsFallback) {
  const correlation::app::AppSettings settings;
  EXPECT_EQ(settings.window_width, 1180U);
  EXPECT_EQ(settings.window_height, 795U);
  EXPECT_FLOAT_EQ(settings.left_col_width, 220.0F);
  EXPECT_FLOAT_EQ(settings.middle_col_width, 220.0F);
  EXPECT_FLOAT_EQ(settings.right_col_width, 320.0F);
  EXPECT_FLOAT_EQ(settings.dual_split_ratio, 0.5F);
  EXPECT_TRUE(settings.dark_theme);
}

TEST_F(SettingsManagerTests, RoundTripSerialization) {
  correlation::app::AppSettings custom;
  custom.window_width = 1440U;
  custom.window_height = 900U;
  custom.left_col_width = 315.5F;
  custom.middle_col_width = 280.0F;
  custom.right_col_width = 390.0F;
  custom.dual_split_ratio = 0.65F;
  custom.dark_theme = false;

  ASSERT_TRUE(correlation::app::SettingsManager::save(custom));

  const correlation::app::AppSettings loaded = correlation::app::SettingsManager::load();
  EXPECT_EQ(loaded.window_width, 1440U);
  EXPECT_EQ(loaded.window_height, 900U);
  EXPECT_FLOAT_EQ(loaded.left_col_width, 315.5F);
  EXPECT_FLOAT_EQ(loaded.middle_col_width, 280.0F);
  EXPECT_FLOAT_EQ(loaded.right_col_width, 390.0F);
  EXPECT_FLOAT_EQ(loaded.dual_split_ratio, 0.65F);
  EXPECT_FALSE(loaded.dark_theme);
}

TEST_F(SettingsManagerTests, ClampingBounds) {
  correlation::app::AppSettings out_of_bounds;
  out_of_bounds.window_width = 500U;
  out_of_bounds.window_height = 400U;
  out_of_bounds.left_col_width = 120.0F;
  out_of_bounds.middle_col_width = 500.0F;
  out_of_bounds.right_col_width = 800.0F;
  out_of_bounds.dual_split_ratio = 0.99F;

  ASSERT_TRUE(correlation::app::SettingsManager::save(out_of_bounds));

  const correlation::app::AppSettings loaded = correlation::app::SettingsManager::load();
  EXPECT_EQ(loaded.window_width, 800U);
  EXPECT_EQ(loaded.window_height, 600U);
  EXPECT_FLOAT_EQ(loaded.left_col_width, 180.0F);
  EXPECT_FLOAT_EQ(loaded.middle_col_width, 450.0F);
  EXPECT_FLOAT_EQ(loaded.right_col_width, 520.0F);
  EXPECT_FLOAT_EQ(loaded.dual_split_ratio, 0.85F);
}

} // namespace
