// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "app/PresetManager.hpp"
#include <gtest/gtest.h>
#include <filesystem>

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
    for (const auto &p : all) {
        if (p.name == "TestUnitPreset") {
            found = true;
            EXPECT_DOUBLE_EQ(p.options.r_max, 12.0);
        }
    }
    EXPECT_TRUE(found);

    // Remove
    correlation::app::PresetManager::remove("TestUnitPreset");

    // Load again
    all = correlation::app::PresetManager::loadAll();
    found = false;
    for (const auto &p : all) {
        if (p.name == "TestUnitPreset") {
            found = true;
        }
    }
    EXPECT_FALSE(found);
}
