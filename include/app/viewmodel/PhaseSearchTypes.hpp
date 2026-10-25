/**
 * @file PhaseSearchTypes.hpp
 * @brief Data types, descriptor enums, crystal presets, and telemetry models for Phase Search.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "plotters/ColormapService.hpp"

#include <cstdint>
#include <string_view>
#include <vector>

namespace correlation::app {

/**
 * @enum PhaseDescriptor
 * @brief Bond-orientational order parameter descriptors available for phase filtering.
 */
enum class PhaseDescriptor : std::uint8_t {
  Q6Bar,        ///< Lechner-Dellago averaged q̄₆ (primary crystal symmetry descriptor).
  Q4Bar,        ///< Lechner-Dellago averaged q̄₄.
  Q6,           ///< Standard single-atom q₆.
  Q4,           ///< Standard single-atom q₄.
  W6Hat,        ///< Normalized third-order invariant ŵ₆ (distinguishes icosahedral symmetry).
  W4Hat,        ///< Normalized third-order invariant ŵ₄.
  LocalEntropy, ///< Local Entropic Fingerprint (s_i / k_B).
  TwfNconn,     ///< Ten Wolde-Frenkel coherent bond count (N_conn).
  QTet,         ///< Chau-Hardwick tetrahedral order parameter (q_tet).
};

/**
 * @enum CrystalPreset
 * @brief Canonical crystal structure and disordered state presets.
 */
enum class CrystalPreset : std::uint8_t {
  Custom,           ///< User-defined manual threshold boundaries.
  Liquid,           ///< Disordered / equilibrium liquid state.
  AmorphousGlass,   ///< Dense random-packed metallic / structural glass.
  AmorphousNetwork, ///< Tetrahedral network glass (e.g., a-Si, a-Ge, LDA water).
  FCC,              ///< Face-Centered Cubic lattice.
  HCP,              ///< Hexagonal Close-Packed lattice.
  BCC,              ///< Body-Centered Cubic lattice.
  Icosahedral,      ///< Icosahedral quasi-crystalline / cluster symmetry.
  Diamond,          ///< Diamond / tetrahedral network lattice.
};

/**
 * @struct PhasePresetConfig
 * @brief Nominal descriptor and threshold window for a canonical crystal structure or disordered
 * preset.
 */
struct PhasePresetConfig {
  PhaseDescriptor descriptor{PhaseDescriptor::Q6Bar}; ///< 1D order descriptor used by preset.
  float min_threshold{0.46F};                         ///< 1D lower scalar bound.
  float max_threshold{0.58F};                         ///< 1D upper scalar bound.
  bool is_2d{false};                                  ///< True if preset defines 2D joint bounds.
  PhaseDescriptor descriptor_2d_x{PhaseDescriptor::Q6Bar}; ///< 2D X-axis order descriptor.
  float min_threshold_2d_x{0.0F};                          ///< 2D X-axis lower scalar bound.
  float max_threshold_2d_x{1.0F};                          ///< 2D X-axis upper scalar bound.
  PhaseDescriptor descriptor_2d_y{PhaseDescriptor::Q4Bar}; ///< 2D Y-axis order descriptor.
  float min_threshold_2d_y{0.0F};                          ///< 2D Y-axis lower scalar bound.
  float max_threshold_2d_y{0.3F};                          ///< 2D Y-axis upper scalar bound.
};

/**
 * @struct PhaseSearchSettings
 * @brief Comprehensive configuration bundle for real-time 3D phase transition visualization.
 */
struct PhaseSearchSettings {
  bool enabled{false};                                ///< Master toggle for phase halos.
  PhaseDescriptor descriptor{PhaseDescriptor::Q6Bar}; ///< Active order descriptor.
  CrystalPreset preset{CrystalPreset::Custom};        ///< Active preset mode.
  float min_threshold{0.46F};                         ///< Active scalar filter lower bound.
  float max_threshold{0.58F};                         ///< Active scalar filter upper bound.
  plotters::ColormapType colormap{
      ///< Colormap for scalar gradient.
      plotters::ColormapType::Viridis,
  };
  float halo_intensity{0.85F};    ///< Alpha glow strength [0.1, 1.0].
  float halo_radius_ratio{1.35F}; ///< Outer glow radius multiplier.
  bool dim_non_matching{false};   ///< Dim core atoms outside threshold.
  bool hide_non_matching{false};  ///< Hide/cull core atoms outside threshold.
  float dim_factor{0.3F};         ///< Luminance factor for non-matches.
  bool is_2d_mode{false};         ///< True if 2D joint distribution filtering is active.
  PhaseDescriptor descriptor_y{PhaseDescriptor::Q4Bar}; ///< Y-axis order descriptor for 2D mode.
  float min_threshold_y{0.0F};                          ///< Y-axis lower scalar bound.
  float max_threshold_y{0.3F};                          ///< Y-axis upper scalar bound.
};

/**
 * @struct PhaseHistogram
 * @brief 50-bin distribution histogram of the active descriptor across all atoms in the frame.
 */
struct PhaseHistogram {
  std::vector<float> bin_centers;  ///< Center scalar value for each of the 50 bins.
  std::vector<uint32_t> counts;    ///< Atom frequency count per bin.
  float min_val{0.0F};             ///< Nominal domain minimum.
  float max_val{1.0F};             ///< Nominal domain maximum.
  float bin_width{0.02F};          ///< Bin delta spacing.
  uint32_t matching_atom_count{0}; ///< Number of atoms within [min_threshold, max_threshold].
  uint32_t total_atom_count{0};    ///< Total atom count in configuration.
};

/**
 * @struct PhaseHeatmap2D
 * @brief 64x64 joint distribution histogram across two order parameters (X and Y).
 */
struct PhaseHeatmap2D {
  static constexpr std::size_t K_GRID_SIZE = 64;
  std::vector<uint32_t> counts; ///< 64x64 grid cells (row 0 is max_y, row 63 is min_y).
  float min_x{0.0F};            ///< Domain minimum on X-axis.
  float max_x{1.0F};            ///< Domain maximum on X-axis.
  float min_y{0.0F};            ///< Domain minimum on Y-axis.
  float max_y{1.0F};            ///< Domain maximum on Y-axis.
  float bin_width_x{0.015625F}; ///< Bin width along X.
  float bin_width_y{0.015625F}; ///< Bin width along Y.
  uint32_t max_count{0};        ///< Peak bin count across the 64x64 grid.
  uint32_t matching_atom_count{
      0,
  }; ///< Atoms falling inside 2D bounding box [Xmin, Xmax] x [Ymin, Ymax].
  uint32_t total_atom_count{0}; ///< Total atom count.
};

/**
 * @struct CrystalMarkerData
 * @brief Canonical crystal reference coordinate normalized to [0, 1] on 2D heatmap.
 */
struct CrystalMarkerData {
  std::string_view label; ///< Crystal phase abbreviation (FCC, HCP, BCC, Ico, Liq, Gl, a-Net).
  float x_ratio{0.0F};    ///< Normalized X position [0.0, 1.0].
  float y_ratio{0.0F};    ///< Normalized Y position [0.0, 1.0] (0 = top, 1 = bottom).
};

/**
 * @brief Converts PhaseDescriptor to human-readable string representation.
 * @param[in] descriptor Descriptor enum value.
 * @return String view of descriptor name.
 */
[[nodiscard]] constexpr std::string_view descriptorToString(PhaseDescriptor descriptor) noexcept {
  switch (descriptor) {
  case PhaseDescriptor::Q6Bar:
    return "q̄₆";
  case PhaseDescriptor::Q4Bar:
    return "q̄₄";
  case PhaseDescriptor::Q6:
    return "q₆";
  case PhaseDescriptor::Q4:
    return "q₄";
  case PhaseDescriptor::W6Hat:
    return "ŵ₆";
  case PhaseDescriptor::W4Hat:
    return "ŵ₄";
  case PhaseDescriptor::LocalEntropy:
    return "s₂ (LEF)";
  case PhaseDescriptor::TwfNconn:
    return "N_conn (TWF)";
  case PhaseDescriptor::QTet:
    return "q_tet (Tetrahedral)";
  }
  return "Unknown";
}

/**
 * @brief Converts CrystalPreset to human-readable string representation.
 * @param[in] preset Preset enum value.
 * @return String view of preset name.
 */
[[nodiscard]] constexpr std::string_view presetToString(CrystalPreset preset) noexcept {
  switch (preset) {
  case CrystalPreset::Custom:
    return "Custom";
  case CrystalPreset::Liquid:
    return "Liquid";
  case CrystalPreset::AmorphousGlass:
    return "Glass";
  case CrystalPreset::AmorphousNetwork:
    return "a-Net";
  case CrystalPreset::FCC:
    return "FCC";
  case CrystalPreset::HCP:
    return "HCP";
  case CrystalPreset::BCC:
    return "BCC";
  case CrystalPreset::Icosahedral:
    return "Icosahedral";
  case CrystalPreset::Diamond:
    return "Diamond";
  }
  return "Unknown";
}

/**
 * @brief Retrieves nominal descriptor and threshold bounds for canonical crystal presets.
 * @param[in] preset Preset enum value.
 * @return PhasePresetConfig populated with standard Steinhardt bounds.
 */
[[nodiscard]] constexpr PhasePresetConfig getPresetConfig(CrystalPreset preset) noexcept {
  switch (preset) {
  case CrystalPreset::Liquid:
    return {
        .descriptor = PhaseDescriptor::Q6Bar,
        .min_threshold = 0.10F,
        .max_threshold = 0.28F,
        .is_2d = true,
        .descriptor_2d_x = PhaseDescriptor::Q6Bar,
        .min_threshold_2d_x = 0.05F,
        .max_threshold_2d_x = 0.28F,
        .descriptor_2d_y = PhaseDescriptor::LocalEntropy,
        .min_threshold_2d_y = -3.5F,
        .max_threshold_2d_y = -1.0F,
    };
  case CrystalPreset::AmorphousGlass:
    return {
        .descriptor = PhaseDescriptor::LocalEntropy,
        .min_threshold = -8.0F,
        .max_threshold = -4.5F,
        .is_2d = true,
        .descriptor_2d_x = PhaseDescriptor::Q6Bar,
        .min_threshold_2d_x = 0.05F,
        .max_threshold_2d_x = 0.32F,
        .descriptor_2d_y = PhaseDescriptor::LocalEntropy,
        .min_threshold_2d_y = -8.0F,
        .max_threshold_2d_y = -4.5F,
    };
  case CrystalPreset::AmorphousNetwork:
    return {
        .descriptor = PhaseDescriptor::QTet,
        .min_threshold = 0.75F,
        .max_threshold = 1.00F,
        .is_2d = true,
        .descriptor_2d_x = PhaseDescriptor::QTet,
        .min_threshold_2d_x = 0.75F,
        .max_threshold_2d_x = 1.00F,
        .descriptor_2d_y = PhaseDescriptor::LocalEntropy,
        .min_threshold_2d_y = -7.5F,
        .max_threshold_2d_y = -4.0F,
    };
  case CrystalPreset::FCC:
    return {
        .descriptor = PhaseDescriptor::Q6Bar,
        .min_threshold = 0.46F,
        .max_threshold = 0.58F,
        .is_2d = true,
        .descriptor_2d_x = PhaseDescriptor::Q6Bar,
        .min_threshold_2d_x = 0.46F,
        .max_threshold_2d_x = 0.58F,
        .descriptor_2d_y = PhaseDescriptor::Q4Bar,
        .min_threshold_2d_y = 0.15F,
        .max_threshold_2d_y = 0.24F,
    };
  case CrystalPreset::HCP:
    return {
        .descriptor = PhaseDescriptor::Q6Bar,
        .min_threshold = 0.44F,
        .max_threshold = 0.52F,
        .is_2d = true,
        .descriptor_2d_x = PhaseDescriptor::Q6Bar,
        .min_threshold_2d_x = 0.44F,
        .max_threshold_2d_x = 0.52F,
        .descriptor_2d_y = PhaseDescriptor::Q4Bar,
        .min_threshold_2d_y = 0.06F,
        .max_threshold_2d_y = 0.14F,
    };
  case CrystalPreset::BCC:
    return {
        .descriptor = PhaseDescriptor::Q6Bar,
        .min_threshold = 0.35F,
        .max_threshold = 0.44F,
        .is_2d = true,
        .descriptor_2d_x = PhaseDescriptor::Q6Bar,
        .min_threshold_2d_x = 0.35F,
        .max_threshold_2d_x = 0.44F,
        .descriptor_2d_y = PhaseDescriptor::Q4Bar,
        .min_threshold_2d_y = 0.01F,
        .max_threshold_2d_y = 0.06F,
    };
  case CrystalPreset::Icosahedral:
    return {
        .descriptor = PhaseDescriptor::W6Hat,
        .min_threshold = -0.17F,
        .max_threshold = -0.14F,
        .is_2d = true,
        .descriptor_2d_x = PhaseDescriptor::Q6Bar,
        .min_threshold_2d_x = 0.15F,
        .max_threshold_2d_x = 0.25F,
        .descriptor_2d_y = PhaseDescriptor::W6Hat,
        .min_threshold_2d_y = -0.17F,
        .max_threshold_2d_y = -0.14F,
    };
  case CrystalPreset::Diamond:
    return {
        .descriptor = PhaseDescriptor::QTet,
        .min_threshold = 0.85F,
        .max_threshold = 1.00F,
        .is_2d = true,
        .descriptor_2d_x = PhaseDescriptor::QTet,
        .min_threshold_2d_x = 0.85F,
        .max_threshold_2d_x = 1.00F,
        .descriptor_2d_y = PhaseDescriptor::Q6Bar,
        .min_threshold_2d_y = 0.35F,
        .max_threshold_2d_y = 0.55F,
    };
  case CrystalPreset::Custom:
  default:
    return {.descriptor = PhaseDescriptor::Q6Bar, .min_threshold = 0.0F, .max_threshold = 1.0F};
  }
}

/**
 * @brief Returns domain bounds [min, max] for a given order descriptor.
 * @param[in] descriptor Target order parameter.
 * @return Pair of floats representing the nominal [min, max] range.
 */
[[nodiscard]] constexpr std::pair<float, float>
getDescriptorDomain(PhaseDescriptor descriptor) noexcept {
  switch (descriptor) {
  case PhaseDescriptor::W6Hat:
  case PhaseDescriptor::W4Hat:
    return {-0.2F, 0.2F};
  case PhaseDescriptor::LocalEntropy:
    return {-15.0F, 0.0F};
  case PhaseDescriptor::TwfNconn:
    return {0.0F, 16.0F};
  case PhaseDescriptor::QTet:
  case PhaseDescriptor::Q6Bar:
  case PhaseDescriptor::Q4Bar:
  case PhaseDescriptor::Q6:
  case PhaseDescriptor::Q4:
  default:
    return {0.0F, 1.0F};
  }
}

namespace detail {

[[nodiscard]] inline std::vector<CrystalMarkerData> markersQ4Q6() noexcept {
  return {
      CrystalMarkerData{.label = "FCC", .x_ratio = 0.191F, .y_ratio = 0.425F},
      CrystalMarkerData{.label = "HCP", .x_ratio = 0.097F, .y_ratio = 0.515F},
      CrystalMarkerData{.label = "BCC", .x_ratio = 0.036F, .y_ratio = 0.489F},
      CrystalMarkerData{.label = "Ico", .x_ratio = 0.000F, .y_ratio = 0.801F},
      CrystalMarkerData{.label = "Liq", .x_ratio = 0.025F, .y_ratio = 0.880F},
  };
}

[[nodiscard]] inline std::vector<CrystalMarkerData> markersQ6Q4() noexcept {
  return {
      CrystalMarkerData{.label = "FCC", .x_ratio = 0.575F, .y_ratio = 0.809F},
      CrystalMarkerData{.label = "HCP", .x_ratio = 0.485F, .y_ratio = 0.903F},
      CrystalMarkerData{.label = "BCC", .x_ratio = 0.511F, .y_ratio = 0.964F},
      CrystalMarkerData{.label = "Ico", .x_ratio = 0.199F, .y_ratio = 1.000F},
      CrystalMarkerData{.label = "Liq", .x_ratio = 0.120F, .y_ratio = 0.975F},
  };
}

[[nodiscard]] inline std::vector<CrystalMarkerData> markersQ6W6() noexcept {
  return {
      CrystalMarkerData{.label = "FCC", .x_ratio = 0.575F, .y_ratio = 0.500F},
      CrystalMarkerData{.label = "HCP", .x_ratio = 0.485F, .y_ratio = 0.530F},
      CrystalMarkerData{.label = "BCC", .x_ratio = 0.511F, .y_ratio = 0.467F},
      CrystalMarkerData{.label = "Ico", .x_ratio = 0.199F, .y_ratio = 0.923F},
      CrystalMarkerData{.label = "Liq", .x_ratio = 0.120F, .y_ratio = 0.500F},
  };
}

[[nodiscard]] inline std::vector<CrystalMarkerData> markersW6Q6() noexcept {
  return {
      CrystalMarkerData{.label = "FCC", .x_ratio = 0.500F, .y_ratio = 0.425F},
      CrystalMarkerData{.label = "HCP", .x_ratio = 0.470F, .y_ratio = 0.515F},
      CrystalMarkerData{.label = "BCC", .x_ratio = 0.533F, .y_ratio = 0.489F},
      CrystalMarkerData{.label = "Ico", .x_ratio = 0.077F, .y_ratio = 0.801F},
      CrystalMarkerData{.label = "Liq", .x_ratio = 0.500F, .y_ratio = 0.880F},
  };
}

[[nodiscard]] inline std::vector<CrystalMarkerData> markersQ4W4() noexcept {
  return {
      CrystalMarkerData{.label = "FCC", .x_ratio = 0.191F, .y_ratio = 0.898F},
      CrystalMarkerData{.label = "HCP", .x_ratio = 0.097F, .y_ratio = 0.165F},
      CrystalMarkerData{.label = "BCC", .x_ratio = 0.036F, .y_ratio = 0.332F},
      CrystalMarkerData{.label = "Ico", .x_ratio = 0.000F, .y_ratio = 0.500F},
  };
}

[[nodiscard]] inline std::vector<CrystalMarkerData> markersQ6LocalEntropy() noexcept {
  return {
      CrystalMarkerData{.label = "FCC", .x_ratio = 0.575F, .y_ratio = 0.700F},
      CrystalMarkerData{.label = "HCP", .x_ratio = 0.485F, .y_ratio = 0.700F},
      CrystalMarkerData{.label = "BCC", .x_ratio = 0.511F, .y_ratio = 0.633F},
      CrystalMarkerData{.label = "Ico", .x_ratio = 0.199F, .y_ratio = 0.600F},
      CrystalMarkerData{.label = "Gl", .x_ratio = 0.160F, .y_ratio = 0.400F},
      CrystalMarkerData{.label = "Liq", .x_ratio = 0.120F, .y_ratio = 0.133F},
  };
}

[[nodiscard]] inline std::vector<CrystalMarkerData> markersLocalEntropyQ6() noexcept {
  return {
      CrystalMarkerData{.label = "FCC", .x_ratio = 0.300F, .y_ratio = 0.425F},
      CrystalMarkerData{.label = "HCP", .x_ratio = 0.300F, .y_ratio = 0.515F},
      CrystalMarkerData{.label = "BCC", .x_ratio = 0.367F, .y_ratio = 0.489F},
      CrystalMarkerData{.label = "Ico", .x_ratio = 0.400F, .y_ratio = 0.801F},
      CrystalMarkerData{.label = "Gl", .x_ratio = 0.600F, .y_ratio = 0.840F},
      CrystalMarkerData{.label = "Liq", .x_ratio = 0.867F, .y_ratio = 0.880F},
  };
}

[[nodiscard]] inline std::vector<CrystalMarkerData> markersQTetLocalEntropy() noexcept {
  return {
      CrystalMarkerData{.label = "Dia", .x_ratio = 0.950F, .y_ratio = 0.733F},
      CrystalMarkerData{.label = "a-Net", .x_ratio = 0.850F, .y_ratio = 0.367F},
      CrystalMarkerData{.label = "Liq", .x_ratio = 0.320F, .y_ratio = 0.133F},
  };
}

[[nodiscard]] inline std::vector<CrystalMarkerData> markersQTetQ6() noexcept {
  return {
      CrystalMarkerData{.label = "Dia", .x_ratio = 0.950F, .y_ratio = 0.550F},
      CrystalMarkerData{.label = "a-Net", .x_ratio = 0.850F, .y_ratio = 0.850F},
      CrystalMarkerData{.label = "Liq", .x_ratio = 0.320F, .y_ratio = 0.880F},
  };
}

} // namespace detail

/**
 * @brief Computes normalized reference coordinates for canonical crystal phases in 2D space.
 * @param[in] desc_x X-axis order parameter descriptor.
 * @param[in] desc_y Y-axis order parameter descriptor.
 * @return Vector of CrystalMarkerData coordinates normalized to canvas [0, 1].
 */
[[nodiscard]] inline std::vector<CrystalMarkerData>
getCrystalReferenceMarkers(PhaseDescriptor desc_x, PhaseDescriptor desc_y) noexcept {
  const auto is_q6 = [](PhaseDescriptor descriptor) noexcept {
    return descriptor == PhaseDescriptor::Q6Bar || descriptor == PhaseDescriptor::Q6;
  };
  const auto is_q4 = [](PhaseDescriptor descriptor) noexcept {
    return descriptor == PhaseDescriptor::Q4Bar || descriptor == PhaseDescriptor::Q4;
  };

  if (is_q4(desc_x) && is_q6(desc_y)) {
    return detail::markersQ4Q6();
  }
  if (is_q6(desc_x) && is_q4(desc_y)) {
    return detail::markersQ6Q4();
  }
  if (is_q6(desc_x) && desc_y == PhaseDescriptor::W6Hat) {
    return detail::markersQ6W6();
  }
  if (desc_x == PhaseDescriptor::W6Hat && is_q6(desc_y)) {
    return detail::markersW6Q6();
  }
  if (is_q4(desc_x) && desc_y == PhaseDescriptor::W4Hat) {
    return detail::markersQ4W4();
  }
  if (is_q6(desc_x) && desc_y == PhaseDescriptor::LocalEntropy) {
    return detail::markersQ6LocalEntropy();
  }
  if (desc_x == PhaseDescriptor::LocalEntropy && is_q6(desc_y)) {
    return detail::markersLocalEntropyQ6();
  }
  if (desc_x == PhaseDescriptor::QTet && desc_y == PhaseDescriptor::LocalEntropy) {
    return detail::markersQTetLocalEntropy();
  }
  if (desc_x == PhaseDescriptor::QTet && is_q6(desc_y)) {
    return detail::markersQTetQ6();
  }

  return {};
}

} // namespace correlation::app
