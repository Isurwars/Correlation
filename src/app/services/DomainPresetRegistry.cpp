/**
 * @file DomainPresetRegistry.cpp
 * @brief Implementation of DomainPresetRegistry.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/services/DomainPresetRegistry.hpp"

namespace correlation::app {
namespace {

const std::vector<DomainPresetDefinition> &buildDomainPresets() {
  static const std::vector<DomainPresetDefinition> PRESETS = {
      DomainPresetDefinition{
          .id = DomainPresetId::SilicateBorateGlasses,
          .name = "Silicate & Borate Oxide Glasses",
          .description = "Covalent network glasses with ring statistics and tetrahedral bonding.",
          .badge_text = "Oxide Glass",
          .recommended_elements = "Si, B, O, Na, Ca, Al",
          .r_max = 15.0,
          .r_bin_width = 0.02,
          .q_max = 20.0,
          .q_bin_width = 0.02,
          .angle_bin_width = 0.25,
          .max_ring_size = 12,
          .smoothing = false,
          .smoothing_sigma = 0.02,
          .enabled_calculators = {"RDF", "PAD", "RD", "CN", "S(Q)"},
      },
      DomainPresetDefinition{
          .id = DomainPresetId::AqueousMolecularLiquids,
          .name = "Aqueous & Molecular Liquids",
          .description = "Hydrogen-bonded liquids, diffusion dynamics, and velocity correlation.",
          .badge_text = "Aqueous Liquid",
          .recommended_elements = "H, O, C, N",
          .r_max = 10.0,
          .r_bin_width = 0.05,
          .q_max = 15.0,
          .q_bin_width = 0.05,
          .angle_bin_width = 0.5,
          .max_ring_size = 8,
          .smoothing = true,
          .smoothing_sigma = 0.08,
          .enabled_calculators = {"RDF", "PAD", "MSD", "VACF", "Hydrogen Bond", "CN"},
      },
      DomainPresetDefinition{
          .id = DomainPresetId::MoltenSaltsIonicLiquids,
          .name = "Molten Salts & Ionic Liquids",
          .description = "Coulombic screening, charge alternation, and intermediate-range order.",
          .badge_text = "Molten Salt",
          .recommended_elements = "Li, Na, K, Cl, F, Al",
          .r_max = 12.0,
          .r_bin_width = 0.04,
          .q_max = 18.0,
          .q_bin_width = 0.03,
          .angle_bin_width = 0.5,
          .max_ring_size = 8,
          .smoothing = true,
          .smoothing_sigma = 0.05,
          .enabled_calculators = {"RDF", "S(Q)", "CN", "MSD"},
      },
      DomainPresetDefinition{
          .id = DomainPresetId::MetallicGlassesHEAs,
          .name = "Metallic Glasses & High-Entropy Alloys",
          .description = "Dense random packing, Voronoi polyhedra, and local entropy fingerprints.",
          .badge_text = "Metallic Glass / HEA",
          .recommended_elements = "Zr, Cu, Ni, Ti, Fe, Co, Cr",
          .r_max = 16.0,
          .r_bin_width = 0.02,
          .q_max = 22.0,
          .q_bin_width = 0.02,
          .angle_bin_width = 0.25,
          .max_ring_size = 10,
          .smoothing = false,
          .smoothing_sigma = 0.02,
          .enabled_calculators = {"RDF", "S(Q)", "CN", "CNA", "Voronoi", "Local Entropy"},
      },
  };
  return PRESETS;
}

} // namespace

auto DomainPresetRegistry::getPresets() -> const std::vector<DomainPresetDefinition> & {
  return buildDomainPresets();
}

auto DomainPresetRegistry::getPreset(DomainPresetId preset_id) -> const DomainPresetDefinition * {
  const auto &presets = buildDomainPresets();
  const auto iter = std::ranges::find_if(
      presets, [preset_id](const auto &preset) { return preset.id == preset_id; });
  if (iter != presets.end()) {
    return &(*iter);
  }
  return nullptr;
}

auto DomainPresetRegistry::applyPreset(DomainPresetId preset_id, ProgramOptions &options) -> bool {
  const auto *preset = getPreset(preset_id);
  if (preset == nullptr) {
    return false;
  }

  options.r_max = preset->r_max;
  options.r_bin_width = preset->r_bin_width;
  options.q_max = preset->q_max;
  options.q_bin_width = preset->q_bin_width;
  options.angle_bin_width = preset->angle_bin_width;
  options.dihedral_bin_width = preset->angle_bin_width;
  options.max_ring_size = preset->max_ring_size;
  options.smoothing = preset->smoothing;
  options.smoothing_sigma = preset->smoothing_sigma;

  // Set appropriate material type index:
  // 0: Amorphous, 1: Liquid, 2: Crystalline
  if (preset_id == DomainPresetId::AqueousMolecularLiquids ||
      preset_id == DomainPresetId::MoltenSaltsIonicLiquids) {
    options.material_type = 1; // Liquid
  } else {
    options.material_type = 0; // Amorphous
  }

  // Update calculator toggles
  options.active_calculators.clear();
  for (const auto &calc_name : preset->enabled_calculators) {
    options.active_calculators[calc_name] = true;
  }

  return true;
}

} // namespace correlation::app
