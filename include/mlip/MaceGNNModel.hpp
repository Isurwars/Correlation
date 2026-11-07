/**
 * @file MaceGNNModel.hpp
 * @brief LibTorch C++ inference engine for MACE (Higher-Order Equivariant Message Passing) models.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "mlip/MLIPInterface.hpp"

#include <memory>
#include <string>

namespace correlation::mlip {

/**
 * @class MaceGNNModel
 * @brief High-performance TorchScript MACE inference engine for atomic systems.
 *
 * Ingests serialized MACE TorchScript models (`.pt` or `.model` exported via
 * `mace.cli.create_lammps_model`) and evaluates potential energy, atomic forces,
 * and virial stress tensors over periodic boundary condition structures.
 */
class MaceGNNModel : public MLIPInterface {
public:
  /**
   * @brief Constructs a MACE TorchScript model from a serialized .pt or .model file.
   * @param[in] model_path Path to the compiled MACE TorchScript model artifact.
   * @param[in] device Execution target device (e.g. "cpu", "cuda", "cuda:0").
   * @param[in] cutoff Radial neighbor cutoff radius in Angstroms (default: 5.0).
   */
  explicit MaceGNNModel(std::string model_path, std::string device = "cpu",
                        real_t cutoff = static_cast<real_t>(5.0));

  ~MaceGNNModel() override;

  // Move-only semantics for model weights
  MaceGNNModel(MaceGNNModel &&) noexcept;
  MaceGNNModel &operator=(MaceGNNModel &&) noexcept;
  MaceGNNModel(const MaceGNNModel &) = delete;
  MaceGNNModel &operator=(const MaceGNNModel &) = delete;

  /**
   * @brief Returns the model architecture and source path identifier.
   * @return Model descriptor string.
   */
  [[nodiscard]] std::string getModelName() const override;

  /**
   * @brief Evaluates potential properties (energy, forces, stress) on an atomic Cell.
   * @param[in] cell Atomic configuration cell.
   * @return MLIPOutput containing atomic forces and virial stress tensor.
   */
  [[nodiscard]] MLIPOutput evaluate(const correlation::core::Cell &cell) const override;

  /**
   * @brief Gets the radial neighbor cutoff radius.
   * @return Cutoff distance in Angstroms.
   */
  [[nodiscard]] real_t cutoff() const noexcept { return cutoff_; }

  /**
   * @brief Sets the radial neighbor cutoff radius.
   * @param[in] cutoff New cutoff distance in Angstroms.
   */
  void setCutoff(real_t cutoff) noexcept { cutoff_ = cutoff; }

  /**
   * @brief Checks if LibTorch backend is compiled and model successfully loaded.
   * @return True if model is ready for forward evaluation.
   */
  [[nodiscard]] bool isLoaded() const noexcept;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
  std::string model_path_;
  std::string device_str_;
  real_t cutoff_{static_cast<real_t>(5.0)};
};

} // namespace correlation::mlip

