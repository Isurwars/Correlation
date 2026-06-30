// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "app/PresetManager.hpp"
#include <filesystem>
#include <gtest/gtest.h>

namespace {

class PresetManagerTests : public ::testing::Test {
protected:
  void SetUp() override {
    // Ensure presets folder is clean or exists
    std::filesystem::create_directories(correlation::app::PresetManager::presetsDirectory());
  }

  void TearDown() override {
    // Clean up test presets
    correlation::app::PresetManager::remove("TestUnitPreset");
    correlation::app::PresetManager::remove("TestUnitPreset2");
    correlation::app::PresetManager::remove("TestUnitPreset_Special");
    correlation::app::PresetManager::remove("TestUnitPreset_Sanitize!@#");
    correlation::app::PresetManager::remove("TestUnitPreset_A");
    correlation::app::PresetManager::remove("TestUnitPreset_B");
    correlation::app::PresetManager::remove("TestUnitPreset_C");
  }
};

TEST_F(PresetManagerTests, RoundTripSerialization) {
  correlation::app::Preset preset;
  preset.name = "TestUnitPreset";
  preset.description = "Unit testing preset description with quotes \"escaped\"";
  preset.options.r_max = 18.5;
  preset.options.r_bin_width = 0.05;
  preset.options.smoothing = true;
  preset.options.active_calculators = {{"RDF", true}, {"S_K", false}};

  // Serialize
  std::string const json = correlation::app::PresetManager::toJson(preset);

  // Deserialize
  correlation::app::Preset loaded = correlation::app::PresetManager::fromJson(json);

  // Verify
  EXPECT_EQ(loaded.name, preset.name);
  EXPECT_EQ(loaded.description, preset.description);
  EXPECT_DOUBLE_EQ(loaded.options.r_max, preset.options.r_max);
  EXPECT_DOUBLE_EQ(loaded.options.r_bin_width, preset.options.r_bin_width);
  EXPECT_EQ(loaded.options.smoothing, preset.options.smoothing);

  auto it1 = loaded.options.active_calculators.find("RDF");
  ASSERT_NE(it1, loaded.options.active_calculators.end());
  EXPECT_TRUE(it1->second);

  auto it2 = loaded.options.active_calculators.find("S_K");
  ASSERT_NE(it2, loaded.options.active_calculators.end());
  EXPECT_FALSE(it2->second);
}

TEST_F(PresetManagerTests, SerializationWithSpecialCharacters) {
  correlation::app::Preset preset;
  preset.name = "TestUnitPreset_Special";
  preset.description = "Line 1\nLine 2\tTabbed\\Backslash\rCarriageReturn\"Quote\"";
  preset.options.r_max = 15.0;

  std::string const json = correlation::app::PresetManager::toJson(preset);
  correlation::app::Preset loaded = correlation::app::PresetManager::fromJson(json);

  EXPECT_EQ(loaded.name, preset.name);
  EXPECT_EQ(loaded.description, preset.description);
}

TEST_F(PresetManagerTests, MissingKeysFallback) {
  std::string const json = R"({
    "name": "TestUnitPreset"
  })";

  correlation::app::Preset loaded = correlation::app::PresetManager::fromJson(json);

  EXPECT_EQ(loaded.name, "TestUnitPreset");
  EXPECT_EQ(loaded.description, ""); // default string fallback
  EXPECT_DOUBLE_EQ(loaded.options.r_max, correlation::app::AppDefaults::R_MAX);
  EXPECT_DOUBLE_EQ(loaded.options.r_bin_width, correlation::app::AppDefaults::R_BIN_WIDTH);
  EXPECT_EQ(loaded.options.smoothing, true); // fallback is true
  EXPECT_TRUE(loaded.options.active_calculators.empty());
}

TEST_F(PresetManagerTests, MalformedJsonHandling) {
  // Empty or garbage JSON should not crash and should return default values
  correlation::app::Preset loaded = correlation::app::PresetManager::fromJson("");
  EXPECT_EQ(loaded.name, "");
  EXPECT_DOUBLE_EQ(loaded.options.r_max, correlation::app::AppDefaults::R_MAX);

  loaded = correlation::app::PresetManager::fromJson("{ invalid json");
  EXPECT_EQ(loaded.name, "");
  EXPECT_DOUBLE_EQ(loaded.options.r_max, correlation::app::AppDefaults::R_MAX);
}

TEST_F(PresetManagerTests, SaveAndLoadAll) {
  correlation::app::Preset preset;
  preset.name = "TestUnitPreset";
  preset.description = "Desc";
  preset.options.r_max = 12.0;

  // Save
  correlation::app::PresetManager::save(preset);

  // Load all
  auto all = correlation::app::PresetManager::loadAll();
  bool found = false;
  for (const auto &preset : all) {
    if (preset.name == "TestUnitPreset") {
      found = true;
      EXPECT_DOUBLE_EQ(preset.options.r_max, 12.0);
    }
  }
  EXPECT_TRUE(found);

  // Remove
  correlation::app::PresetManager::remove("TestUnitPreset");

  // Load again
  all = correlation::app::PresetManager::loadAll();
  found = false;
  for (const auto &preset : all) {
    if (preset.name == "TestUnitPreset") {
      found = true;
    }
  }
  EXPECT_FALSE(found);
}

TEST_F(PresetManagerTests, PresetNameSanitization) {
  correlation::app::Preset preset;
  preset.name = "TestUnitPreset_Sanitize!@#";
  preset.description = "Sanitized name test";
  preset.options.r_max = 14.0;

  // Save
  correlation::app::PresetManager::save(preset);

  // Check that the file was created with sanitized name
  std::filesystem::path const dir = correlation::app::PresetManager::presetsDirectory();
  std::filesystem::path const expected_file = dir / "TestUnitPreset_Sanitize___.json";
  EXPECT_TRUE(std::filesystem::exists(expected_file));

  // Load all to see if it reads it back with original name
  auto all = correlation::app::PresetManager::loadAll();
  bool found = false;
  for (const auto &preset : all) {
    if (preset.name == "TestUnitPreset_Sanitize!@#") {
      found = true;
      EXPECT_DOUBLE_EQ(preset.options.r_max, 14.0);
    }
  }
  EXPECT_TRUE(found);

  // Remove using original name
  correlation::app::PresetManager::remove("TestUnitPreset_Sanitize!@#");
  EXPECT_FALSE(std::filesystem::exists(expected_file));
}

TEST_F(PresetManagerTests, LoadAllSorting) {
  correlation::app::Preset preset_c;
  preset_c.name = "TestUnitPreset_C";
  preset_c.options.r_max = 1.0;

  correlation::app::Preset preset_a;
  preset_a.name = "TestUnitPreset_A";
  preset_a.options.r_max = 2.0;

  correlation::app::Preset preset_b;
  preset_b.name = "TestUnitPreset_B";
  preset_b.options.r_max = 3.0;

  // Save in non-alphabetical order
  correlation::app::PresetManager::save(preset_c);
  correlation::app::PresetManager::save(preset_a);
  correlation::app::PresetManager::save(preset_b);

  // Load all
  auto all = correlation::app::PresetManager::loadAll();

  // Find indices of our test presets in the result
  std::vector<std::string> names;
  for (const auto &preset : all) {
    if (preset.name.starts_with("TestUnitPreset_")) {
      names.push_back(preset.name);
    }
  }

  // They must be sorted: A, B, C
  ASSERT_GE(names.size(), 3);
  auto it_a = std::ranges::find(names, "TestUnitPreset_A");
  auto it_b = std::ranges::find(names, "TestUnitPreset_B");
  auto it_c = std::ranges::find(names, "TestUnitPreset_C");

  ASSERT_NE(it_a, names.end());
  ASSERT_NE(it_b, names.end());
  ASSERT_NE(it_c, names.end());

  EXPECT_TRUE(it_a < it_b);
  EXPECT_TRUE(it_b < it_c);

  // Clean up
  correlation::app::PresetManager::remove("TestUnitPreset_A");
  correlation::app::PresetManager::remove("TestUnitPreset_B");
  correlation::app::PresetManager::remove("TestUnitPreset_C");
}

} // namespace