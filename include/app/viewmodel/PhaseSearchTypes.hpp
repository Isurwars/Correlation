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
};

/**
 * @enum CrystalPreset
 * @brief Canonical crystal structure and disordered state presets.
 */
enum class CrystalPreset : std::uint8_t {
  Custom,      ///< User-defined manual threshold boundaries.
  Liquid,      ///< Disordered / liquid state.
  FCC,         ///< Face-Centered Cubic lattice.
  HCP,         ///< Hexagonal Close-Packed lattice.
  BCC,         ///< Body-Centered Cubic lattice.
  Icosahedral, ///< Icosahedral quasi-crystalline / cluster symmetry.
};

/**
 * @struct PhasePresetConfig
 * @brief Nominal descriptor and threshold window for a canonical crystal structure preset.
 */
struct PhasePresetConfig {
  PhaseDescriptor descriptor{PhaseDescriptor::Q6Bar}; ///< Order descriptor used by preset.
  float min_threshold{0.46F};                         ///< Lower scalar bound.
  float max_threshold{0.58F};                         ///< Upper scalar bound.
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
  case CrystalPreset::FCC:
    return "FCC";
  case CrystalPreset::HCP:
    return "HCP";
  case CrystalPreset::BCC:
    return "BCC";
  case CrystalPreset::Icosahedral:
    return "Icosahedral";
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
    return {.descriptor = PhaseDescriptor::Q6Bar, .min_threshold = 0.10F, .max_threshold = 0.28F};
  case CrystalPreset::FCC:
    return {.descriptor = PhaseDescriptor::Q6Bar, .min_threshold = 0.46F, .max_threshold = 0.58F};
  case CrystalPreset::HCP:
    return {.descriptor = PhaseDescriptor::Q6Bar, .min_threshold = 0.44F, .max_threshold = 0.52F};
  case CrystalPreset::BCC:
    return {.descriptor = PhaseDescriptor::Q6Bar, .min_threshold = 0.35F, .max_threshold = 0.44F};
  case CrystalPreset::Icosahedral:
    return {
        .descriptor = PhaseDescriptor::W6Hat,
        .min_threshold = -0.17F,
        .max_threshold = -0.14F,
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
  case PhaseDescriptor::Q6Bar:
  case PhaseDescriptor::Q4Bar:
  case PhaseDescriptor::Q6:
  case PhaseDescriptor::Q4:
  default:
    return {0.0F, 1.0F};
  }
}

} // namespace correlation::app
