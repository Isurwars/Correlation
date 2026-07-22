/**
 * @file MLIPInterface.hpp
 * @brief Abstract interface for Machine Learning Interatomic Potentials (MLIP).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/Cell.hpp"
#include "math/LinearAlgebra.hpp"

#include <string>
#include <vector>

namespace correlation::mlip {

/**
 * @struct MLIPOutput
 * @brief Output container for machine learning interatomic potential calculations.
 */
struct MLIPOutput {
  real_t total_energy{static_cast<real_t>(0.0)};          /**< Total energy of the system. */
  std::vector<real_t> per_atom_energy;                    /**< Site-resolved per-atom energy. */
  std::vector<correlation::math::Vector3<real_t>> forces; /**< Atomic forces acting on each atom. */
  correlation::math::Matrix3<real_t> stress;              /**< Virial stress tensor. */
};

/**
 * @class MLIPInterface
 * @brief Abstract interface for machine learning interatomic potential evaluation engines (ORB-v3, MACE, CHGNet).
 */
class MLIPInterface {
public:
  virtual ~MLIPInterface() = default;

  /**
   * @brief Returns the model architecture name (e.g. "ORB-v3", "ORB-mol", "MACE").
   * @return Model name string.
   */
  [[nodiscard]] virtual std::string getModelName() const = 0;

  /**
   * @brief Evaluates potential energy, forces, and stress for a given unit cell.
   * @param[in] cell The atomic configuration cell.
   * @return Structural MLIPOutput containing energy, forces, and stress tensor.
   */
  [[nodiscard]] virtual MLIPOutput evaluate(const correlation::core::Cell &cell) const = 0;
};

} // namespace correlation::mlip
