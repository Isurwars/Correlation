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
  Displacement, ///< Cumulative PBC-unwrapped atomic displacement |Δr| relative to frame 0 (Å).
  D2Min,        ///< Falk-Langer non-affine squared displacement D²_min relative to frame 0 (Å²).
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
  float min_threshold{0.52F};                         ///< Active scalar filter lower bound.
  float max_threshold{0.62F};                         ///< Active scalar filter upper bound.
  plotters::ColormapType colormap{
      ///< Colormap for scalar gradient.
      plotters::ColormapType::Viridis,
  };
  float halo_intensity{0.85F};    ///< Alpha glow strength [0.1, 1.0].
  float halo_radius_ratio{2.70F}; ///< Outer glow radius multiplier.
  float atom_scale{0.35F};        ///< Uniform ball radius scaling factor [0.05, 2.0].
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
 * @brief 128x128 joint distribution histogram across two order parameters (X and Y).
 */
struct PhaseHeatmap2D {
  static constexpr std::size_t K_GRID_SIZE = 128;
  std::vector<uint32_t> counts;  ///< 128x128 grid cells (row 0 is max_y, row 127 is min_y).
  float min_x{0.0F};             ///< Domain minimum on X-axis.
  float max_x{1.0F};             ///< Domain maximum on X-axis.
  float min_y{0.0F};             ///< Domain minimum on Y-axis.
  float max_y{1.0F};             ///< Domain maximum on Y-axis.
  float bin_width_x{0.0078125F}; ///< Bin width along X.
  float bin_width_y{0.0078125F}; ///< Bin width along Y.
  uint32_t max_count{0};         ///< Peak bin count across the 128x128 grid.
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
  std::string_view
      label; ///< Crystal phase abbreviation (FCC, HCP, BCC, Ico, Liq, Gl, a-Net, LDA, HDA).
  float x_ratio{0.0F};      ///< Normalized X position [0.0, 1.0].
  float y_ratio{0.0F};      ///< Normalized Y position [0.0, 1.0] (0 = top, 1 = bottom).
  float width_ratio{0.0F};  ///< Normalized width extent [0.0, 1.0].
  float height_ratio{0.0F}; ///< Normalized height extent [0.0, 1.0].
  bool is_area{false};      ///< True if marker represents an extended phase region.
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
  case PhaseDescriptor::Displacement:
    return "|Δr| (Disp)";
  case PhaseDescriptor::D2Min:
    return "D²_min (Strain)";
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
        .min_threshold = 0.06F,
        .max_threshold = 0.22F,
        .is_2d = true,
        .descriptor_2d_x = PhaseDescriptor::Q6Bar,
        .min_threshold_2d_x = 0.06F,
        .max_threshold_2d_x = 0.22F,
        .descriptor_2d_y = PhaseDescriptor::LocalEntropy,
        .min_threshold_2d_y = -3.5F,
        .max_threshold_2d_y = -1.5F,
    };
  case CrystalPreset::AmorphousGlass:
    return {
        .descriptor = PhaseDescriptor::LocalEntropy,
        .min_threshold = -5.5F,
        .max_threshold = -3.5F,
        .is_2d = true,
        .descriptor_2d_x = PhaseDescriptor::Q6Bar,
        .min_threshold_2d_x = 0.10F,
        .max_threshold_2d_x = 0.26F,
        .descriptor_2d_y = PhaseDescriptor::LocalEntropy,
        .min_threshold_2d_y = -5.5F,
        .max_threshold_2d_y = -3.5F,
    };
  case CrystalPreset::AmorphousNetwork:
    return {
        .descriptor = PhaseDescriptor::QTet,
        .min_threshold = 0.75F,
        .max_threshold = 0.95F,
        .is_2d = true,
        .descriptor_2d_x = PhaseDescriptor::QTet,
        .min_threshold_2d_x = 0.75F,
        .max_threshold_2d_x = 0.95F,
        .descriptor_2d_y = PhaseDescriptor::LocalEntropy,
        .min_threshold_2d_y = -5.5F,
        .max_threshold_2d_y = -3.5F,
    };
  case CrystalPreset::FCC:
    return {
        .descriptor = PhaseDescriptor::Q6Bar,
        .min_threshold = 0.52F,
        .max_threshold = 0.62F,
        .is_2d = true,
        .descriptor_2d_x = PhaseDescriptor::Q6Bar,
        .min_threshold_2d_x = 0.52F,
        .max_threshold_2d_x = 0.62F,
        .descriptor_2d_y = PhaseDescriptor::Q4Bar,
        .min_threshold_2d_y = 0.14F,
        .max_threshold_2d_y = 0.24F,
    };
  case CrystalPreset::HCP:
    return {
        .descriptor = PhaseDescriptor::Q6Bar,
        .min_threshold = 0.44F,
        .max_threshold = 0.53F,
        .is_2d = true,
        .descriptor_2d_x = PhaseDescriptor::Q6Bar,
        .min_threshold_2d_x = 0.44F,
        .max_threshold_2d_x = 0.53F,
        .descriptor_2d_y = PhaseDescriptor::Q4Bar,
        .min_threshold_2d_y = 0.05F,
        .max_threshold_2d_y = 0.14F,
    };
  case CrystalPreset::BCC:
    return {
        .descriptor = PhaseDescriptor::Q6Bar,
        .min_threshold = 0.47F,
        .max_threshold = 0.55F,
        .is_2d = true,
        .descriptor_2d_x = PhaseDescriptor::Q6Bar,
        .min_threshold_2d_x = 0.47F,
        .max_threshold_2d_x = 0.55F,
        .descriptor_2d_y = PhaseDescriptor::Q4Bar,
        .min_threshold_2d_y = 0.01F,
        .max_threshold_2d_y = 0.07F,
    };
  case CrystalPreset::Icosahedral:
    return {
        .descriptor = PhaseDescriptor::W6Hat,
        .min_threshold = -0.18F,
        .max_threshold = -0.15F,
        .is_2d = true,
        .descriptor_2d_x = PhaseDescriptor::Q6Bar,
        .min_threshold_2d_x = 0.10F,
        .max_threshold_2d_x = 0.20F,
        .descriptor_2d_y = PhaseDescriptor::W6Hat,
        .min_threshold_2d_y = -0.18F,
        .max_threshold_2d_y = -0.15F,
    };
  case CrystalPreset::Diamond:
    return {
        .descriptor = PhaseDescriptor::QTet,
        .min_threshold = 0.88F,
        .max_threshold = 1.05F,
        .is_2d = true,
        .descriptor_2d_x = PhaseDescriptor::QTet,
        .min_threshold_2d_x = 0.88F,
        .max_threshold_2d_x = 1.05F,
        .descriptor_2d_y = PhaseDescriptor::Q6Bar,
        .min_threshold_2d_y = 0.55F,
        .max_threshold_2d_y = 0.70F,
    };
  case CrystalPreset::Custom:
  default:
    return {.descriptor = PhaseDescriptor::Q6Bar, .min_threshold = 0.0F, .max_threshold = 1.0F};
  }
}

/**
 * @brief Returns domain bounds [min, max] for a given order descriptor.
 * @details Boundaries are expanded to full thermodynamic and physical limits (~20% zoom out).
 * @param[in] descriptor Target order parameter.
 * @return Pair of floats representing the nominal [min, max] range.
 */
[[nodiscard]] constexpr std::pair<float, float>
getDescriptorDomain(PhaseDescriptor descriptor) noexcept {
  switch (descriptor) {
  case PhaseDescriptor::W6Hat:
    return {-0.26F, 0.12F};
  case PhaseDescriptor::W4Hat:
    return {-0.26F, 0.26F};
  case PhaseDescriptor::LocalEntropy:
    return {-12.0F, 0.6F};
  case PhaseDescriptor::TwfNconn:
    return {0.0F, 20.0F};
  case PhaseDescriptor::QTet:
    return {-0.25F, 1.15F};
  case PhaseDescriptor::Displacement:
    return {0.0F, 5.0F};
  case PhaseDescriptor::D2Min:
    return {0.0F, 10.0F};
  case PhaseDescriptor::Q6Bar:
  case PhaseDescriptor::Q6:
    return {0.0F, 0.90F};
  case PhaseDescriptor::Q4Bar:
  case PhaseDescriptor::Q4:
  default:
    return {-0.05F, 1.05F};
  }
}

namespace detail {

struct PhysicalMarker {
  std::string_view label;
  float x{0.0F};
  float y{0.0F};
  float width{0.0F};
  float height{0.0F};
  bool is_area{true};
};

[[nodiscard]] constexpr CrystalMarkerData
toNormalizedMarker(const PhysicalMarker &marker, std::pair<float, float> domain_x,
                   std::pair<float, float> domain_y) noexcept {
  const float span_x = std::max(domain_x.second - domain_x.first, 1e-6F);
  const float span_y = std::max(domain_y.second - domain_y.first, 1e-6F);
  const float x_ratio = std::clamp((marker.x - domain_x.first) / span_x, 0.0F, 1.0F);
  const float y_ratio = std::clamp((domain_y.second - marker.y) / span_y, 0.0F, 1.0F);
  const float width_ratio = std::clamp(marker.width / span_x, 0.0F, 1.0F);
  const float height_ratio = std::clamp(marker.height / span_y, 0.0F, 1.0F);
  return CrystalMarkerData{
      .label = marker.label,
      .x_ratio = x_ratio,
      .y_ratio = y_ratio,
      .width_ratio = width_ratio,
      .height_ratio = height_ratio,
      .is_area = marker.is_area,
  };
}

[[nodiscard]] constexpr PhysicalMarker swapAxes(const PhysicalMarker &marker) noexcept {
  return PhysicalMarker{
      .label = marker.label,
      .x = marker.y,
      .y = marker.x,
      .width = marker.height,
      .height = marker.width,
      .is_area = marker.is_area,
  };
}

[[nodiscard]] inline std::vector<CrystalMarkerData>
normalizeMarkers(std::span<const PhysicalMarker> physical_markers, PhaseDescriptor desc_x,
                 PhaseDescriptor desc_y) noexcept {
  const auto domain_x = getDescriptorDomain(desc_x);
  const auto domain_y = getDescriptorDomain(desc_y);
  std::vector<CrystalMarkerData> result;
  result.reserve(physical_markers.size());
  for (const auto &phy_marker : physical_markers) {
    result.push_back(toNormalizedMarker(phy_marker, domain_x, domain_y));
  }
  return result;
}

[[nodiscard]] inline std::vector<CrystalMarkerData>
normalizeSwappedMarkers(std::span<const PhysicalMarker> physical_markers, PhaseDescriptor desc_x,
                        PhaseDescriptor desc_y) noexcept {
  const auto domain_x = getDescriptorDomain(desc_x);
  const auto domain_y = getDescriptorDomain(desc_y);
  std::vector<CrystalMarkerData> result;
  result.reserve(physical_markers.size());
  for (const auto &phy_marker : physical_markers) {
    result.push_back(toNormalizedMarker(swapAxes(phy_marker), domain_x, domain_y));
  }
  return result;
}

inline constexpr std::array<PhysicalMarker, 8> K_MARKERS_Q6_Q4{{
    PhysicalMarker{.label = "FCC", .x = 0.575F, .y = 0.191F, .width = 0.10F, .height = 0.10F},
    PhysicalMarker{.label = "HCP", .x = 0.485F, .y = 0.097F, .width = 0.09F, .height = 0.08F},
    PhysicalMarker{.label = "BCC", .x = 0.511F, .y = 0.036F, .width = 0.08F, .height = 0.06F},
    PhysicalMarker{.label = "Dia", .x = 0.628F, .y = 0.509F, .width = 0.10F, .height = 0.10F},
    PhysicalMarker{.label = "SC", .x = 0.354F, .y = 0.764F, .width = 0.08F, .height = 0.08F},
    PhysicalMarker{.label = "HD", .x = 0.445F, .y = 0.108F, .width = 0.08F, .height = 0.08F},
    PhysicalMarker{.label = "Ico", .x = 0.199F, .y = 0.000F, .width = 0.08F, .height = 0.06F},
    PhysicalMarker{.label = "Liq", .x = 0.120F, .y = 0.025F, .width = 0.14F, .height = 0.08F},
}};

inline constexpr std::array<PhysicalMarker, 8> K_MARKERS_Q6_W6{{
    PhysicalMarker{.label = "FCC", .x = 0.575F, .y = -0.013F, .width = 0.10F, .height = 0.03F},
    PhysicalMarker{.label = "HCP", .x = 0.485F, .y = -0.012F, .width = 0.09F, .height = 0.03F},
    PhysicalMarker{.label = "BCC", .x = 0.511F, .y = 0.013F, .width = 0.08F, .height = 0.03F},
    PhysicalMarker{.label = "Dia", .x = 0.628F, .y = -0.013F, .width = 0.10F, .height = 0.03F},
    PhysicalMarker{.label = "SC", .x = 0.354F, .y = 0.013F, .width = 0.08F, .height = 0.03F},
    PhysicalMarker{.label = "HD", .x = 0.445F, .y = -0.012F, .width = 0.08F, .height = 0.03F},
    PhysicalMarker{.label = "Ico", .x = 0.150F, .y = -0.165F, .width = 0.08F, .height = 0.04F},
    PhysicalMarker{.label = "Liq", .x = 0.120F, .y = 0.000F, .width = 0.14F, .height = 0.04F},
}};

inline constexpr std::array<PhysicalMarker, 6> K_MARKERS_Q4_W4{{
    PhysicalMarker{.label = "FCC", .x = 0.191F, .y = -0.159F, .width = 0.10F, .height = 0.05F},
    PhysicalMarker{.label = "HCP", .x = 0.097F, .y = 0.134F, .width = 0.08F, .height = 0.04F},
    PhysicalMarker{.label = "BCC", .x = 0.036F, .y = 0.159F, .width = 0.08F, .height = 0.04F},
    PhysicalMarker{.label = "Dia", .x = 0.509F, .y = -0.159F, .width = 0.10F, .height = 0.05F},
    PhysicalMarker{.label = "SC", .x = 0.764F, .y = 0.159F, .width = 0.08F, .height = 0.04F},
    PhysicalMarker{.label = "Ico", .x = 0.000F, .y = 0.000F, .width = 0.08F, .height = 0.05F},
}};

inline constexpr std::array<PhysicalMarker, 6> K_MARKERS_Q6_LOCAL_ENTROPY{{
    PhysicalMarker{.label = "FCC", .x = 0.575F, .y = -8.00F, .width = 0.10F, .height = 2.0F},
    PhysicalMarker{.label = "HCP", .x = 0.485F, .y = -8.00F, .width = 0.09F, .height = 2.0F},
    PhysicalMarker{.label = "BCC", .x = 0.511F, .y = -7.00F, .width = 0.09F, .height = 2.0F},
    PhysicalMarker{.label = "Ico", .x = 0.200F, .y = -6.50F, .width = 0.09F, .height = 2.0F},
    PhysicalMarker{.label = "Gl", .x = 0.178F, .y = -4.50F, .width = 0.16F, .height = 2.0F},
    PhysicalMarker{.label = "Liq", .x = 0.140F, .y = -2.50F, .width = 0.16F, .height = 2.0F},
}};

// Calibrated against Errington & Debenedetti (Nature 2001), Martelli et al. (PRB 2017 / PNAS
// 2018), and Piaggi & Parrinello (JCP 2017) / Gartner et al. (PRL 2020) for (q_tet, s_2):
//   Dia/Ice: q_tet ~ 0.94–1.00, s_2 ~ -9.1 to -6.9 k_B
//   LDA:     q_tet ~ 0.80–0.95 (mean ~0.88–0.91), s_2 ~ -7.3 to -5.1 k_B
//   HDA:     q_tet ~ 0.62–0.78 (mean ~0.70–0.74), s_2 ~ -6.0 to -4.0 k_B
//   Liq:     q_tet ~ 0.45–0.71 (ambient/supercooled water), s_2 ~ -4.1 to -1.9 k_B
inline constexpr std::array<PhysicalMarker, 4> K_MARKERS_QTET_LOCAL_ENTROPY{{
    PhysicalMarker{.label = "Dia", .x = 0.965F, .y = -8.00F, .width = 0.14F, .height = 2.2F},
    PhysicalMarker{.label = "LDA", .x = 0.880F, .y = -6.20F, .width = 0.18F, .height = 2.2F},
    PhysicalMarker{.label = "HDA", .x = 0.700F, .y = -5.00F, .width = 0.18F, .height = 2.0F},
    PhysicalMarker{.label = "Liq", .x = 0.580F, .y = -3.00F, .width = 0.26F, .height = 2.2F},
}};

inline constexpr std::array<PhysicalMarker, 8> K_MARKERS_QTET_Q6{{
    PhysicalMarker{.label = "Dia", .x = 0.965F, .y = 0.628F, .width = 0.14F, .height = 0.10F},
    PhysicalMarker{.label = "HD", .x = 0.960F, .y = 0.445F, .width = 0.12F, .height = 0.08F},
    PhysicalMarker{.label = "a-Net", .x = 0.850F, .y = 0.150F, .width = 0.20F, .height = 0.12F},
    PhysicalMarker{.label = "FCC", .x = 0.000F, .y = 0.575F, .width = 0.10F, .height = 0.08F},
    PhysicalMarker{.label = "HCP", .x = 0.000F, .y = 0.485F, .width = 0.10F, .height = 0.08F},
    PhysicalMarker{.label = "BCC", .x = 0.000F, .y = 0.511F, .width = 0.10F, .height = 0.08F},
    PhysicalMarker{.label = "SC", .x = 0.000F, .y = 0.354F, .width = 0.10F, .height = 0.08F},
    PhysicalMarker{.label = "Liq", .x = 0.350F, .y = 0.120F, .width = 0.20F, .height = 0.12F},
}};

[[nodiscard]] inline std::vector<CrystalMarkerData>
getSteinhardtPairMarkers(PhaseDescriptor desc_x, PhaseDescriptor desc_y) noexcept {
  const auto is_q6 = [](PhaseDescriptor descriptor) noexcept {
    return descriptor == PhaseDescriptor::Q6Bar || descriptor == PhaseDescriptor::Q6;
  };
  const auto is_q4 = [](PhaseDescriptor descriptor) noexcept {
    return descriptor == PhaseDescriptor::Q4Bar || descriptor == PhaseDescriptor::Q4;
  };

  if (is_q4(desc_x) && is_q6(desc_y)) {
    return normalizeSwappedMarkers(K_MARKERS_Q6_Q4, desc_x, desc_y);
  }
  if (is_q6(desc_x) && is_q4(desc_y)) {
    return normalizeMarkers(K_MARKERS_Q6_Q4, desc_x, desc_y);
  }
  if (is_q6(desc_x) && desc_y == PhaseDescriptor::W6Hat) {
    return normalizeMarkers(K_MARKERS_Q6_W6, desc_x, desc_y);
  }
  if (desc_x == PhaseDescriptor::W6Hat && is_q6(desc_y)) {
    return normalizeSwappedMarkers(K_MARKERS_Q6_W6, desc_x, desc_y);
  }
  if (is_q4(desc_x) && desc_y == PhaseDescriptor::W4Hat) {
    return normalizeMarkers(K_MARKERS_Q4_W4, desc_x, desc_y);
  }
  if (desc_x == PhaseDescriptor::W4Hat && is_q4(desc_y)) {
    return normalizeSwappedMarkers(K_MARKERS_Q4_W4, desc_x, desc_y);
  }
  return {};
}

[[nodiscard]] inline std::vector<CrystalMarkerData>
getDisorderedPairMarkers(PhaseDescriptor desc_x, PhaseDescriptor desc_y) noexcept {
  const auto is_q6 = [](PhaseDescriptor descriptor) noexcept {
    return descriptor == PhaseDescriptor::Q6Bar || descriptor == PhaseDescriptor::Q6;
  };

  if (is_q6(desc_x) && desc_y == PhaseDescriptor::LocalEntropy) {
    return normalizeMarkers(K_MARKERS_Q6_LOCAL_ENTROPY, desc_x, desc_y);
  }
  if (desc_x == PhaseDescriptor::LocalEntropy && is_q6(desc_y)) {
    return normalizeSwappedMarkers(K_MARKERS_Q6_LOCAL_ENTROPY, desc_x, desc_y);
  }
  if (desc_x == PhaseDescriptor::QTet && desc_y == PhaseDescriptor::LocalEntropy) {
    return normalizeMarkers(K_MARKERS_QTET_LOCAL_ENTROPY, desc_x, desc_y);
  }
  if (desc_x == PhaseDescriptor::LocalEntropy && desc_y == PhaseDescriptor::QTet) {
    return normalizeSwappedMarkers(K_MARKERS_QTET_LOCAL_ENTROPY, desc_x, desc_y);
  }
  if (desc_x == PhaseDescriptor::QTet && is_q6(desc_y)) {
    return normalizeMarkers(K_MARKERS_QTET_Q6, desc_x, desc_y);
  }
  if (is_q6(desc_x) && desc_y == PhaseDescriptor::QTet) {
    return normalizeSwappedMarkers(K_MARKERS_QTET_Q6, desc_x, desc_y);
  }
  return {};
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
  auto steinhardt_markers = detail::getSteinhardtPairMarkers(desc_x, desc_y);
  if (!steinhardt_markers.empty()) {
    return steinhardt_markers;
  }
  return detail::getDisorderedPairMarkers(desc_x, desc_y);
}

} // namespace correlation::app
