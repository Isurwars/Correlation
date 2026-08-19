/**
 * @file TorchGNNModel.hpp
 * @brief LibTorch C++ inference engine for Graph Neural Networks predicting forces and LDOS.
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
 * @class TorchGNNModel
 * @brief High-performance TorchScript GNN inference engine for atomic systems.
 */
class TorchGNNModel : public MLIPInterface {
public:
  /**
   * @brief Constructs a TorchScript GNN model from a serialized .pt file.
   * @param[in] model_path Path to the compiled TorchScript model artifact (.pt).
   * @param[in] device Execution target device (e.g. "cpu", "cuda", "cuda:0").
   * @param[in] cutoff Radial neighbor cutoff radius in Angstroms (default: 5.0).
   */
  explicit TorchGNNModel(std::string model_path, std::string device = "cpu", real_t cutoff = static_cast<real_t>(5.0));

  ~TorchGNNModel() override;

  // Move-only semantics for model weights
  TorchGNNModel(TorchGNNModel &&) noexcept;
  TorchGNNModel &operator=(TorchGNNModel &&) noexcept;
  TorchGNNModel(const TorchGNNModel &) = delete;
  TorchGNNModel &operator=(const TorchGNNModel &) = delete;

  /**
   * @brief Returns the model architecture and source path identifier.
   * @return Model descriptor string.
   */
  [[nodiscard]] std::string getModelName() const override;

  /**
   * @brief Evaluates potential properties (energy, forces, LDOS) on an atomic Cell.
   * @param[in] cell Atomic configuration cell.
   * @return MLIPOutput containing atomic forces and site-resolved LDOS bins.
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
