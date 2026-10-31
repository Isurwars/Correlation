/**
 * @file Structure3DControllerPhaseSearchTests.cpp
 * @brief Unit tests for Phase Search and Steinhardt Halo Visualization in Structure3DController.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "AppWindow.h"
#include "app/core/AppOptions.hpp"
#include "app/viewmodel/Structure3DController.hpp"
#include "core/Cell.hpp"

#include "../../CrystalTestHelper.hpp"

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <optional>
#include <stdexcept>

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#endif

namespace correlation::app {
namespace {

class Structure3DControllerPhaseSearchTests : public ::testing::Test {
public:
  [[nodiscard]] AppWindow &window() {
    if (!window_.has_value()) {
      throw std::runtime_error("Window is not initialized");
    }
    return **window_;
  }

protected:
  void SetUp() override {
#ifndef _WIN32
    setenv("SLINT_BACKEND", "software", 1);
#else
    _putenv_s("SLINT_BACKEND", "software");
#endif
    window_.emplace(AppWindow::create());
  }

private:
  std::optional<slint::ComponentHandle<AppWindow>> window_;
};

TEST_F(Structure3DControllerPhaseSearchTests, DefaultPhaseSearchSettings) {
  const ProgramOptions options;
  const Structure3DController controller(window(), options);

  EXPECT_FALSE(controller.isPhaseSearchEnabled());
  EXPECT_EQ(controller.phaseDescriptor(), PhaseDescriptor::Q6Bar);
  EXPECT_EQ(controller.crystalPreset(), CrystalPreset::Custom);
  EXPECT_FLOAT_EQ(controller.phaseMinThreshold(), 0.46F);
  EXPECT_FLOAT_EQ(controller.phaseMaxThreshold(), 0.58F);
  EXPECT_EQ(controller.phaseColormap(), plotters::ColormapType::Viridis);
  EXPECT_NEAR(controller.phaseHaloIntensity(), 0.85F, 1e-4F);
  EXPECT_NEAR(controller.phaseHaloRadiusRatio(), 1.35F, 1e-4F);
  EXPECT_FALSE(controller.isPhaseDimNonMatching());
}

TEST_F(Structure3DControllerPhaseSearchTests, CrystalPresetApplication) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  // Apply FCC preset
  controller.setCrystalPreset(CrystalPreset::FCC);
  EXPECT_EQ(controller.crystalPreset(), CrystalPreset::FCC);
  EXPECT_EQ(controller.phaseDescriptor(), PhaseDescriptor::Q6Bar);
  EXPECT_FLOAT_EQ(controller.phaseMinThreshold(), 0.46F);
  EXPECT_FLOAT_EQ(controller.phaseMaxThreshold(), 0.58F);

  // Apply Liquid preset
  controller.setCrystalPreset(CrystalPreset::Liquid);
  EXPECT_EQ(controller.crystalPreset(), CrystalPreset::Liquid);
  EXPECT_EQ(controller.phaseDescriptor(), PhaseDescriptor::Q6Bar);
  EXPECT_FLOAT_EQ(controller.phaseMinThreshold(), 0.10F);
  EXPECT_FLOAT_EQ(controller.phaseMaxThreshold(), 0.28F);

  // Apply AmorphousGlass preset (1D)
  controller.setCrystalPreset(CrystalPreset::AmorphousGlass);
  EXPECT_EQ(controller.crystalPreset(), CrystalPreset::AmorphousGlass);
  EXPECT_EQ(controller.phaseDescriptor(), PhaseDescriptor::LocalEntropy);
  EXPECT_FLOAT_EQ(controller.phaseMinThreshold(), -8.0F);
  EXPECT_FLOAT_EQ(controller.phaseMaxThreshold(), -4.5F);

  // Apply AmorphousNetwork preset (1D)
  controller.setCrystalPreset(CrystalPreset::AmorphousNetwork);
  EXPECT_EQ(controller.crystalPreset(), CrystalPreset::AmorphousNetwork);
  EXPECT_EQ(controller.phaseDescriptor(), PhaseDescriptor::QTet);
  EXPECT_FLOAT_EQ(controller.phaseMinThreshold(), 0.75F);
  EXPECT_FLOAT_EQ(controller.phaseMaxThreshold(), 1.00F);

  // Apply Icosahedral preset
  controller.setCrystalPreset(CrystalPreset::Icosahedral);
  EXPECT_EQ(controller.crystalPreset(), CrystalPreset::Icosahedral);
  EXPECT_EQ(controller.phaseDescriptor(), PhaseDescriptor::W6Hat);
  EXPECT_FLOAT_EQ(controller.phaseMinThreshold(), -0.17F);
  EXPECT_FLOAT_EQ(controller.phaseMaxThreshold(), -0.14F);

  // Manual threshold adjustment resets preset mode to Custom
  controller.setPhaseThresholds(0.30F, 0.70F);
  EXPECT_EQ(controller.crystalPreset(), CrystalPreset::Custom);
  EXPECT_FLOAT_EQ(controller.phaseMinThreshold(), 0.30F);
  EXPECT_FLOAT_EQ(controller.phaseMaxThreshold(), 0.70F);

  // 2D mode preset application
  controller.setPhase2DMode(true);
  controller.setCrystalPreset(CrystalPreset::Liquid);
  EXPECT_EQ(controller.crystalPreset(), CrystalPreset::Liquid);
  EXPECT_EQ(controller.phaseDescriptor(), PhaseDescriptor::Q6Bar);
  EXPECT_FLOAT_EQ(controller.phaseMinThreshold(), 0.05F);
  EXPECT_FLOAT_EQ(controller.phaseMaxThreshold(), 0.28F);
  EXPECT_EQ(controller.phaseDescriptorY(), PhaseDescriptor::LocalEntropy);
  EXPECT_FLOAT_EQ(controller.phaseMinThresholdY(), -3.5F);
  EXPECT_FLOAT_EQ(controller.phaseMaxThresholdY(), -1.0F);

  controller.setCrystalPreset(CrystalPreset::AmorphousGlass);
  EXPECT_EQ(controller.crystalPreset(), CrystalPreset::AmorphousGlass);
  EXPECT_EQ(controller.phaseDescriptor(), PhaseDescriptor::Q6Bar);
  EXPECT_FLOAT_EQ(controller.phaseMinThreshold(), 0.05F);
  EXPECT_FLOAT_EQ(controller.phaseMaxThreshold(), 0.32F);
  EXPECT_EQ(controller.phaseDescriptorY(), PhaseDescriptor::LocalEntropy);
  EXPECT_FLOAT_EQ(controller.phaseMinThresholdY(), -8.0F);
  EXPECT_FLOAT_EQ(controller.phaseMaxThresholdY(), -4.5F);

  controller.setPhaseThresholds2D(0.10F, 0.40F, -5.0F, -2.0F);
  EXPECT_EQ(controller.crystalPreset(), CrystalPreset::Custom);
  EXPECT_FLOAT_EQ(controller.phaseMinThreshold(), 0.10F);
  EXPECT_FLOAT_EQ(controller.phaseMaxThreshold(), 0.40F);
  EXPECT_FLOAT_EQ(controller.phaseMinThresholdY(), -5.0F);
  EXPECT_FLOAT_EQ(controller.phaseMaxThresholdY(), -2.0F);
}

TEST_F(Structure3DControllerPhaseSearchTests, ColormapAndIntensityConfiguration) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  controller.setPhaseColormap(plotters::ColormapType::Plasma);
  EXPECT_EQ(controller.phaseColormap(), plotters::ColormapType::Plasma);

  controller.setPhaseColormap(plotters::ColormapType::Turbo);
  EXPECT_EQ(controller.phaseColormap(), plotters::ColormapType::Turbo);

  controller.setPhaseHaloIntensity(0.5F);
  EXPECT_FLOAT_EQ(controller.phaseHaloIntensity(), 0.5F);

  controller.setPhaseHaloRadiusRatio(3.8F);
  EXPECT_FLOAT_EQ(controller.phaseHaloRadiusRatio(), 3.8F);

  controller.setAtomScale(0.75F);
  EXPECT_FLOAT_EQ(controller.atomScale(), 0.75F);
  EXPECT_FLOAT_EQ(window().get_structure_3d_atom_scale(), 0.75F);
}

TEST_F(Structure3DControllerPhaseSearchTests, LargeCellAdaptiveHaloRadiusBoostAndAtomScale) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  // 8x8x8 FCC cell of Cu (a = 3.61 * 8 = 28.88 A, diagonal ~50.0 A > 40.0 A threshold)
  const auto large_cell = testing::crystals::createFCCCell(3.61, "Cu", 8, 8, 8);
  controller.setCell(large_cell);
  controller.setAtomScale(0.60F);
  ASSERT_NE(controller.scene(), nullptr);
  ASSERT_FALSE(controller.scene()->atoms().empty());

  const auto &first_atom = controller.scene()->atoms().front();
  EXPECT_NEAR(first_atom.radius, first_atom.covalent_radius * 0.60F, 1e-4F);

  controller.setCrystalPreset(CrystalPreset::FCC);
  controller.setPhaseHaloRadiusRatio(2.5F);
  controller.setPhaseSearchEnabled(true);
  static_cast<void>(controller.getOrComputeSteinhardt(0));
  controller.setPhaseThresholds(0.46F, 0.58F);

  // Large cell diagonal > 40 A receives an adaptive halo radius boost > 1.0x
  EXPECT_GT(controller.scene()->atoms().front().halo_radius_ratio, 2.5F);
}

void assertFCCQ6BarValues(const std::vector<real_t> &q6_bar) {
  ASSERT_FALSE(q6_bar.empty());
  for (const auto val : q6_bar) {
    EXPECT_NEAR(val, 0.575, 0.05);
  }
}

TEST_F(Structure3DControllerPhaseSearchTests, SteinhardtEvaluationAndLRUCacheOnFCCCell) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  // Setup FCC cell
  const auto cell = testing::crystals::createFCCCell(3.61, "Cu", 2, 2, 2);
  controller.setCell(cell);
  ASSERT_NE(controller.scene(), nullptr);
  EXPECT_GT(controller.scene()->atoms().size(), 0UL);

  // Initial cache state
  EXPECT_FALSE(controller.lruCache().contains(0));

  // Compute / get Steinhardt order parameters
  const auto params = controller.getOrComputeSteinhardt(0);
  ASSERT_NE(params, nullptr);
  EXPECT_EQ(params->Q6_bar.size(), cell.atoms().size());
  EXPECT_TRUE(controller.lruCache().contains(0));

  // For ideal FCC, Q6_bar is ~0.575, matching FCC preset [0.46, 0.58]
  assertFCCQ6BarValues(params->Q6_bar);

  // Second call must be a cache hit returning identical pointer
  const auto cached_params = controller.getOrComputeSteinhardt(0);
  EXPECT_EQ(params.get(), cached_params.get());
}

void assertAtomsHaveActiveHalos(const std::vector<render3d::AtomInstance> &atoms) {
  ASSERT_FALSE(atoms.empty());
  for (const auto &atom : atoms) {
    EXPECT_GT(atom.halo_color.a, 0);
    EXPECT_FLOAT_EQ(atom.dim_factor, 1.0F);
  }
}

void assertAtomsHaveNoHalos(const std::vector<render3d::AtomInstance> &atoms, float expected_dim) {
  ASSERT_FALSE(atoms.empty());
  for (const auto &atom : atoms) {
    EXPECT_EQ(atom.halo_color.a, 0);
    EXPECT_FLOAT_EQ(atom.dim_factor, expected_dim);
  }
}

TEST_F(Structure3DControllerPhaseSearchTests, PhaseHalosAndDimmingOnScene) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  const auto cell = testing::crystals::createFCCCell(3.61, "Cu", 1, 1, 1);
  controller.setCell(cell);
  controller.setCrystalPreset(CrystalPreset::FCC);

  // Enable phase search
  controller.setPhaseSearchEnabled(true);
  EXPECT_TRUE(controller.isPhaseSearchEnabled());

  // Force synchronous computation to apply visuals to scene
  static_cast<void>(controller.getOrComputeSteinhardt(0));
  controller.setPhaseThresholds(0.46F, 0.58F);
  assertAtomsHaveActiveHalos(controller.scene()->atoms());

  // If we set thresholds outside the FCC range (e.g. [0.10, 0.20]), halos are removed
  controller.setPhaseThresholds(0.10F, 0.20F);
  assertAtomsHaveNoHalos(controller.scene()->atoms(), 1.0F);

  // If dim_non_matching is enabled, non-matching atoms should have dim_factor = 0.3F
  controller.setPhaseDimNonMatching(true);
  EXPECT_TRUE(controller.isPhaseDimNonMatching());
  assertAtomsHaveNoHalos(controller.scene()->atoms(), 0.3F);

  // Disabling phase search restores dim_factor = 1.0F and halo_color.a = 0
  controller.setPhaseSearchEnabled(false);
  assertAtomsHaveNoHalos(controller.scene()->atoms(), 1.0F);
}

TEST_F(Structure3DControllerPhaseSearchTests, TelemetryHistogramGeneration) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  const auto cell = testing::crystals::createFCCCell(3.61, "Cu", 2, 2, 2);
  controller.setCell(cell);
  controller.setCrystalPreset(CrystalPreset::FCC);
  controller.setPhaseSearchEnabled(true);

  // Compute Steinhardt parameters and refresh histogram
  static_cast<void>(controller.getOrComputeSteinhardt(0));
  controller.setPhaseThresholds(0.46F, 0.58F);

  const auto hist = controller.phaseHistogram();
  EXPECT_EQ(hist.counts.size(), 50UL);
  EXPECT_EQ(hist.bin_centers.size(), 50UL);
  EXPECT_EQ(hist.total_atom_count, cell.atoms().size());
  EXPECT_GT(hist.matching_atom_count, 0UL);
  EXPECT_LE(hist.matching_atom_count, hist.total_atom_count);
}

TEST_F(Structure3DControllerPhaseSearchTests, SlintWindowBindingSynchronization) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  controller.setPhaseSearchEnabled(true);
  controller.setCrystalPreset(CrystalPreset::FCC);
  controller.setPhaseHaloIntensity(0.92F);
  controller.setPhaseHaloRadiusRatio(1.65F);
  controller.setPhaseDimNonMatching(true);

  EXPECT_TRUE(window().get_structure_3d_phase_enabled());
  EXPECT_EQ(window().get_structure_3d_phase_preset(), static_cast<int>(CrystalPreset::FCC));
  EXPECT_EQ(window().get_structure_3d_phase_descriptor(), static_cast<int>(PhaseDescriptor::Q6Bar));
  EXPECT_FLOAT_EQ(window().get_structure_3d_phase_min_threshold(), 0.46F);
  EXPECT_FLOAT_EQ(window().get_structure_3d_phase_max_threshold(), 0.58F);
  EXPECT_NEAR(window().get_structure_3d_phase_halo_intensity(), 0.92F, 1e-4F);
  EXPECT_NEAR(window().get_structure_3d_phase_halo_radius_ratio(), 1.65F, 1e-4F);
  EXPECT_TRUE(window().get_structure_3d_phase_dim_non_matching());
  EXPECT_EQ(window().get_structure_3d_phase_descriptor_name(), "q̄₆");
}

TEST_F(Structure3DControllerPhaseSearchTests, Phase2DModeAndThresholdConfiguration) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  EXPECT_FALSE(controller.isPhase2DMode());
  EXPECT_EQ(controller.phaseDescriptorY(), PhaseDescriptor::Q4Bar);

  controller.setPhase2DMode(true);
  EXPECT_TRUE(controller.isPhase2DMode());
  EXPECT_TRUE(window().get_structure_3d_phase_is_2d_mode());

  controller.setPhaseDescriptorY(PhaseDescriptor::QTet);
  EXPECT_EQ(controller.phaseDescriptorY(), PhaseDescriptor::QTet);
  EXPECT_EQ(window().get_structure_3d_phase_descriptor_y(),
            static_cast<int>(PhaseDescriptor::QTet));
  EXPECT_EQ(window().get_structure_3d_phase_descriptor_name_y(), "q_tet (Tetrahedral)");

  controller.setPhaseThresholds2D(0.40F, 0.60F, 0.20F, 0.80F);
  EXPECT_FLOAT_EQ(controller.phaseMinThreshold(), 0.40F);
  EXPECT_FLOAT_EQ(controller.phaseMaxThreshold(), 0.60F);
  EXPECT_FLOAT_EQ(controller.phaseMinThresholdY(), 0.20F);
  EXPECT_FLOAT_EQ(controller.phaseMaxThresholdY(), 0.80F);
  EXPECT_FLOAT_EQ(window().get_structure_3d_phase_min_threshold_y(), 0.20F);
  EXPECT_FLOAT_EQ(window().get_structure_3d_phase_max_threshold_y(), 0.80F);
}

void verifyCrystalMarkersQ6Q4(const std::vector<CrystalMarkerData> &markers) {
  EXPECT_FALSE(markers.empty());
  float fcc_x_ratio = -1.0F;
  for (const auto &marker : markers) {
    if (marker.label == "FCC") {
      fcc_x_ratio = marker.x_ratio;
    }
  }
  EXPECT_GT(fcc_x_ratio, 0.4F);
  EXPECT_LT(fcc_x_ratio, 0.7F);
}

void verifyHeatmapPngExport(Structure3DController &controller) {
  const auto test_png = std::filesystem::temp_directory_path() / "test_heatmap_export.png";
  std::error_code error_code;
  std::filesystem::remove(test_png, error_code);

  const bool exported = controller.exportHeatmapPng(test_png.string());
  EXPECT_TRUE(exported);

  std::ifstream check_file(test_png, std::ios::binary | std::ios::ate);
  ASSERT_TRUE(check_file.is_open());
  EXPECT_GT(check_file.tellg(), 0);
  check_file.close();
  EXPECT_TRUE(std::filesystem::remove(test_png, error_code));
}

TEST_F(Structure3DControllerPhaseSearchTests, PhaseHeatmap2DAndPngExport) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  const auto cell = testing::crystals::createFCCCell(3.61, "Cu", 2, 2, 2);
  controller.setCell(cell);
  controller.setPhaseSearchEnabled(true);
  controller.setPhase2DMode(true);
  controller.setPhaseDescriptor(PhaseDescriptor::Q6Bar);
  controller.setPhaseDescriptorY(PhaseDescriptor::Q4Bar);

  // Compute descriptors and 2D heatmap
  static_cast<void>(controller.getOrComputeSteinhardt(0));
  controller.setPhaseThresholds2D(0.45F, 0.60F, 0.10F, 0.25F);

  const auto heatmap = controller.phaseHeatmap();
  EXPECT_EQ(heatmap.counts.size(), 64UL * 64UL);
  EXPECT_EQ(heatmap.total_atom_count, cell.atoms().size());
  EXPECT_GT(heatmap.matching_atom_count, 0UL);
  EXPECT_LE(heatmap.matching_atom_count, heatmap.total_atom_count);

  const auto markers = getCrystalReferenceMarkers(PhaseDescriptor::Q6Bar, PhaseDescriptor::Q4Bar);
  verifyCrystalMarkersQ6Q4(markers);
  verifyHeatmapPngExport(controller);
}

TEST_F(Structure3DControllerPhaseSearchTests, SetPhaseDescriptorYResetsPresetToCustom) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  controller.setCrystalPreset(CrystalPreset::FCC);
  EXPECT_EQ(controller.crystalPreset(), CrystalPreset::FCC);

  // Changing Y descriptor must reset preset to Custom
  controller.setPhaseDescriptorY(PhaseDescriptor::W6Hat);
  EXPECT_EQ(controller.crystalPreset(), CrystalPreset::Custom);
  EXPECT_EQ(controller.phaseDescriptorY(), PhaseDescriptor::W6Hat);
}

TEST_F(Structure3DControllerPhaseSearchTests, CrystalPreset2DConfigurationConsistency) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  // Configure FCC preset while in 1D mode
  controller.setPhase2DMode(false);
  controller.setCrystalPreset(CrystalPreset::FCC);

  // Switch to 2D mode, verify Y descriptor and bounds were calibrated by the preset
  controller.setPhase2DMode(true);
  const auto settings = controller.phaseSearchSettings();
  EXPECT_EQ(settings.descriptor, PhaseDescriptor::Q6Bar);
  EXPECT_EQ(settings.descriptor_y, PhaseDescriptor::Q4Bar);
  EXPECT_FLOAT_EQ(settings.min_threshold_y, 0.15F);
  EXPECT_FLOAT_EQ(settings.max_threshold_y, 0.24F);
}

TEST_F(Structure3DControllerPhaseSearchTests, OffThreadSetCellAndTrajectorySafety) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  const auto cell = testing::crystals::createFCCCell(3.61, "Cu", 1, 1, 1);

  // Call setCell from a background thread
  {
    const std::jthread bg_thread([&controller, &cell] { controller.setCell(cell); });
  }

  // Pump Slint event loop to process queued property mutations
  slint::invoke_from_event_loop([] { slint::quit_event_loop(); });
  slint::run_event_loop();

  EXPECT_EQ(window().get_structure_3d_total_frames(), 1);
  EXPECT_EQ(window().get_structure_3d_current_frame(), 0);
}

} // namespace
} // namespace correlation::app
