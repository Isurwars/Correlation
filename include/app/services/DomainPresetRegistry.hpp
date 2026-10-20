/**
 * @file DomainPresetRegistry.hpp
 * @brief Registry and templates for scientific material domain presets.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "app/core/AppOptions.hpp"
#include "math/Precision.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace correlation::app {

/**
 * @enum DomainPresetId
 * @brief Pre-packaged material domain template identifiers.
 */
enum class DomainPresetId : std::uint8_t {
  SilicateBorateGlasses = 0,   ///< Oxide glasses with covalent network and ring statistics.
  AqueousMolecularLiquids = 1, ///< Hydrogen-bonded liquids, diffusion, and velocity dynamics.
  MoltenSaltsIonicLiquids = 2, ///< Molten salts, charge alternation, and screening radii.
  MetallicGlassesHEAs = 3,     ///< Metallic glasses, Voronoi polyhedra, and local entropy.
  Invalid = 255,               ///< Sentinel for invalid preset identifiers.
};

/**
 * @struct DomainPresetDefinition
 * @brief Specification for a scientific domain preset.
 */
struct DomainPresetDefinition {
  DomainPresetId id = DomainPresetId::SilicateBorateGlasses; ///< Unique domain identifier.
  std::string name;                                          ///< Human-readable domain name.
  std::string description;                      ///< Scientific description of analysis workflow.
  std::string badge_text;                       ///< Compact category tag.
  std::string recommended_elements;             ///< Typical element archetypes.
  real_t r_max = 15.0;                          ///< Max radial distance in Å.
  real_t r_bin_width = 0.02;                    ///< RDF histogram bin width in Å.
  real_t q_max = 20.0;                          ///< Max scattering momentum transfer in Å^-1.
  real_t q_bin_width = 0.02;                    ///< S(Q) histogram bin width in Å^-1.
  real_t angle_bin_width = 0.25;                ///< Angular distribution bin width in degrees.
  size_t max_ring_size = 12;                    ///< Topological ring cutoff.
  bool smoothing = false;                       ///< Whether Gaussian smoothing is enabled.
  real_t smoothing_sigma = 0.02;                ///< Gaussian kernel sigma.
  std::vector<std::string> enabled_calculators; ///< Calculator IDs to enable for this domain.
};

/**
 * @class DomainPresetRegistry
 * @brief Registry of pre-configured material domain workflow templates.
 */
class DomainPresetRegistry {
public:
  DomainPresetRegistry() = delete;

  /**
   * @brief Retrieves all available domain preset templates.
   * @return Constant reference to vector of domain preset definitions.
   */
  [[nodiscard]] static auto getPresets() -> const std::vector<DomainPresetDefinition> &;

  /**
   * @brief Retrieves a specific preset by its identifier.
   * @param[in] preset_id Domain preset identifier.
   * @return Pointer to definition if found, nullptr otherwise.
   */
  [[nodiscard]] static auto getPreset(DomainPresetId preset_id) -> const DomainPresetDefinition *;

  /**
   * @brief Applies a domain preset configuration to the specified ProgramOptions.
   * @param[in] preset_id Domain preset identifier.
   * @param[in,out] options Target program options to update.
   * @return True if preset was applied, false if identifier was invalid.
   */
  static auto applyPreset(DomainPresetId preset_id, ProgramOptions &options) -> bool;
};

} // namespace correlation::app
