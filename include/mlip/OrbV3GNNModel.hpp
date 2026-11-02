/**
 * @file OrbV3GNNModel.hpp
 * @brief Native C++23 ORB-v3 Graph Neural Network + LDOSHead inference engine.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "mlip/MLIPInterface.hpp"

#include <array>
#include <cstdint>
#include <memory>
#include <string>

namespace correlation::mlip {

/**
 * @struct OrbV3ModelHeader
 * @brief 64-byte binary header for serialized Correlation ORB-v3 + LDOSHead weight files.
 */
struct OrbV3ModelHeader {
  std::array<char, 8> magic{};      /**< File magic identifier ("CORRORB3"). */
  std::uint32_t version{0};         /**< Binary format version (1). */
  std::uint32_t latent_dim{0};      /**< Node/edge latent embedding dimension (256). */
  std::uint32_t mlp_hidden_dim{0};  /**< Backbone MLP hidden width (1024). */
  std::uint32_t head_hidden_dim{0}; /**< Physical readout MLP hidden width (256). */
  std::uint32_t ldos_hidden_dim{0}; /**< LDOS readout MLP hidden width (1024). */
  std::uint32_t num_stacks{0};      /**< Number of AttentionInteractionNetwork stacks (5). */
  std::uint32_t num_rbf{0};         /**< Number of Bessel radial basis functions (8). */
  std::uint32_t l_max{0};           /**< Maximum spherical harmonics degree (3). */
  std::uint32_t ldos_bins{0};       /**< Number of continuous energy grid points (256). */
  std::uint32_t max_neighbors{0};   /**< Maximum nearest periodic neighbors per atom (120). */
  float r_max{0.0F};                /**< Radial neighbor cutoff radius in Angstroms (6.0). */
  float ldos_emin{0.0F};            /**< Lower bound of LDOS energy grid in eV (-30.0). */
  float ldos_emax{0.0F};            /**< Upper bound of LDOS energy grid in eV (15.0). */
  float ldos_sigma{0.0F};           /**< Gaussian broadening width in eV (0.2). */
};

static_assert(sizeof(OrbV3ModelHeader) == 64, "OrbV3ModelHeader must be exactly 64 bytes.");

/**
 * @class OrbV3GNNModel
 * @brief Zero-dependency native C++23 inference engine for Cellulose ORB-v3 + LDOSHead potentials.
 *
 * Evaluates total/per-atom potential energy, atomic forces, virial stress, and
 * 256-bin atom-projected Local Density of States (LDOS) directly from atomic
 * positions (`core::Cell`) and `PeriodicGraphBuilder::buildOrbGraph` tensors.
 */
class OrbV3GNNModel : public MLIPInterface {
public:
  /**
   * @brief Constructs an OrbV3GNNModel from a flat binary weight file (or auto-discovers default).
   * @param[in] model_path Path to `cellulose_ldos.bin` (if empty, searches default paths).
   */
  explicit OrbV3GNNModel(std::string model_path = "");

  ~OrbV3GNNModel() override;

  OrbV3GNNModel(OrbV3GNNModel &&) noexcept;
  OrbV3GNNModel &operator=(OrbV3GNNModel &&) noexcept;
  OrbV3GNNModel(const OrbV3GNNModel &) = delete;
  OrbV3GNNModel &operator=(const OrbV3GNNModel &) = delete;

  /**
   * @brief Searches standard workspace and user directories for `cellulose_ldos.bin`.
   * @param[in] preferred_path Optional explicit user path to check first.
   * @return Resolved existing file path, or empty string if not found.
   */
  [[nodiscard]] static std::string findDefaultModelPath(const std::string &preferred_path = "");

  /**
   * @brief Returns the model architecture and loaded weight file identifier.
   * @return Model descriptor string.
   */
  [[nodiscard]] std::string getModelName() const override;

  /**
   * @brief Evaluates energy, forces, stress, and atom-projected LDOS on a simulation Cell.
   * @param[in] cell Atomic configuration cell.
   * @return MLIPOutput containing energy, forces, stress, and site-resolved LDOS spectra.
   * @throws std::runtime_error if unmapped numeric species (Z=0) are present or model is unloaded.
   */
  [[nodiscard]] MLIPOutput evaluate(const correlation::core::Cell &cell) const override;

  /**
   * @brief Checks whether the binary weight file was successfully mapped and validated.
   * @return True if ready for inference.
   */
  [[nodiscard]] bool isLoaded() const noexcept;

  /**
   * @brief Returns the parsed 64-byte header of the loaded model.
   * @return Const reference to OrbV3ModelHeader.
   */
  [[nodiscard]] const OrbV3ModelHeader &header() const noexcept;

  /**
   * @brief Returns the resolved filesystem path of the loaded weight binary.
   * @return Model file path string.
   */
  [[nodiscard]] const std::string &modelPath() const noexcept { return model_path_; }

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
  std::string model_path_;
};

} // namespace correlation::mlip
