/**
 * @file DomainPresetRegistryTests.cpp
 * @brief Unit tests for DomainPresetRegistry material templates and option application.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/services/DomainPresetRegistry.hpp"

#include <gtest/gtest.h>

using namespace correlation::app;

TEST(DomainPresetRegistryTests, ReturnsAllFourDomainPresets) {
  const auto &presets = DomainPresetRegistry::getPresets();
  ASSERT_EQ(presets.size(), 4);

  EXPECT_EQ(presets[0].id, DomainPresetId::SilicateBorateGlasses);
  EXPECT_EQ(presets[0].name, "Silicate & Borate Oxide Glasses");
  EXPECT_EQ(presets[0].badge_text, "Oxide Glass");

  EXPECT_EQ(presets[1].id, DomainPresetId::AqueousMolecularLiquids);
  EXPECT_EQ(presets[1].name, "Aqueous & Molecular Liquids");
  EXPECT_EQ(presets[1].badge_text, "Aqueous Liquid");

  EXPECT_EQ(presets[2].id, DomainPresetId::MoltenSaltsIonicLiquids);
  EXPECT_EQ(presets[2].name, "Molten Salts & Ionic Liquids");
  EXPECT_EQ(presets[2].badge_text, "Molten Salt");

  EXPECT_EQ(presets[3].id, DomainPresetId::MetallicGlassesHEAs);
  EXPECT_EQ(presets[3].name, "Metallic Glasses & High-Entropy Alloys");
  EXPECT_EQ(presets[3].badge_text, "Metallic Glass / HEA");
}

TEST(DomainPresetRegistryTests, FindsPresetById) {
  const auto *preset = DomainPresetRegistry::getPreset(DomainPresetId::SilicateBorateGlasses);
  ASSERT_NE(preset, nullptr);
  EXPECT_EQ(preset->name, "Silicate & Borate Oxide Glasses");

  const auto *invalid = DomainPresetRegistry::getPreset(DomainPresetId::Invalid);
  EXPECT_EQ(invalid, nullptr);
}

TEST(DomainPresetRegistryTests, AppliesSilicatePresetToOptions) {
  ProgramOptions opts;
  const bool applied =
      DomainPresetRegistry::applyPreset(DomainPresetId::SilicateBorateGlasses, opts);
  EXPECT_TRUE(applied);

  EXPECT_NEAR(opts.r_max, 15.0, 1e-4);
  EXPECT_NEAR(opts.r_bin_width, 0.02, 1e-4);
  EXPECT_NEAR(opts.q_max, 20.0, 1e-4);
  EXPECT_NEAR(opts.q_bin_width, 0.02, 1e-4);
  EXPECT_NEAR(opts.angle_bin_width, 0.25, 1e-4);
  EXPECT_EQ(opts.max_ring_size, 12);
  EXPECT_EQ(opts.material_type, 0);

  EXPECT_TRUE(opts.active_calculators.at("RDF"));
  EXPECT_TRUE(opts.active_calculators.at("PAD"));
  EXPECT_TRUE(opts.active_calculators.at("RD"));
  EXPECT_TRUE(opts.active_calculators.at("CN"));
  EXPECT_TRUE(opts.active_calculators.at("S(Q)"));
}

TEST(DomainPresetRegistryTests, AppliesAqueousPresetToOptions) {
  ProgramOptions opts;
  const bool applied =
      DomainPresetRegistry::applyPreset(DomainPresetId::AqueousMolecularLiquids, opts);
  EXPECT_TRUE(applied);

  EXPECT_NEAR(opts.r_max, 10.0, 1e-4);
  EXPECT_NEAR(opts.r_bin_width, 0.05, 1e-4);
  EXPECT_EQ(opts.material_type, 1);
  EXPECT_TRUE(opts.active_calculators.at("RDF"));
  EXPECT_TRUE(opts.active_calculators.at("Hydrogen Bond"));
  EXPECT_TRUE(opts.active_calculators.at("MSD"));
  EXPECT_TRUE(opts.active_calculators.at("VACF"));
}

TEST(DomainPresetRegistryTests, AppliesMoltenSaltsPresetToOptions) {
  ProgramOptions opts;
  const bool applied =
      DomainPresetRegistry::applyPreset(DomainPresetId::MoltenSaltsIonicLiquids, opts);
  EXPECT_TRUE(applied);

  EXPECT_NEAR(opts.r_max, 12.0, 1e-4);
  EXPECT_NEAR(opts.r_bin_width, 0.04, 1e-4);
  EXPECT_EQ(opts.material_type, 1);
  EXPECT_TRUE(opts.active_calculators.at("RDF"));
  EXPECT_TRUE(opts.active_calculators.at("CN"));
  EXPECT_TRUE(opts.active_calculators.at("MSD"));
  EXPECT_TRUE(opts.active_calculators.at("S(Q)"));
}

TEST(DomainPresetRegistryTests, AppliesMetallicGlassesPresetToOptions) {
  ProgramOptions opts;
  const bool applied = DomainPresetRegistry::applyPreset(DomainPresetId::MetallicGlassesHEAs, opts);
  EXPECT_TRUE(applied);

  EXPECT_NEAR(opts.r_max, 16.0, 1e-4);
  EXPECT_NEAR(opts.r_bin_width, 0.02, 1e-4);
  EXPECT_EQ(opts.material_type, 0);
  EXPECT_TRUE(opts.active_calculators.at("Voronoi"));
  EXPECT_TRUE(opts.active_calculators.at("CNA"));
  EXPECT_TRUE(opts.active_calculators.at("Local Entropy"));
  EXPECT_TRUE(opts.active_calculators.at("RDF"));
}

TEST(DomainPresetRegistryTests, ReturnsFalseOnInvalidPresetId) {
  ProgramOptions opts;
  const bool applied = DomainPresetRegistry::applyPreset(DomainPresetId::Invalid, opts);
  EXPECT_FALSE(applied);
}
