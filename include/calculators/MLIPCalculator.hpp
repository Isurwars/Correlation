/**
 * @file MLIPCalculator.hpp
 * @brief Calculator for machine learning interatomic potential energy, forces, and stress.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "mlip/MLIPInterface.hpp"

#include <string>

namespace correlation::calculators {

/**
 * @class MLIPCalculator
 * @brief Computes atomic energy, forces, and stress using ORB-v3 / ORB-mol model family.
 */
class MLIPCalculator : public BaseCalculator {
public:
  MLIPCalculator() = default;
  explicit MLIPCalculator(const correlation::mlip::MLIPInterface *model) noexcept : model_(model) {}

  [[nodiscard]] std::string_view getName() const override { return "ML Interatomic Potential (ORB-v3)"; }
  [[nodiscard]] std::string_view getShortName() const override { return "MLIP"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Machine Learning"; }
  [[nodiscard]] std::string_view getDescription() const override {
    return "Evaluates atomic energy, forces, and stress using ORB-v3 / ORB-mol model family.";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  /**
   * @brief Checks if a valid MLIP evaluation model is attached.
   * @return True if model is non-null; false otherwise.
   */
  [[nodiscard]] bool isConfigured() const override { return model_ != nullptr; }

  /**
   * @brief Sets or updates the attached MLIP model engine.
   * @param model Pointer to the MLIPInterface engine.
   */
  void setModel(const correlation::mlip::MLIPInterface *model) noexcept { model_ = model; }

  /**
   * @brief Gets the attached MLIP model engine.
   * @return Pointer to the MLIPInterface engine, or nullptr if unset.
   */
  [[nodiscard]] const correlation::mlip::MLIPInterface *getModel() const noexcept { return model_; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Evaluates potential properties for a single unit cell using an MLIP engine.
   * @param[in] cell The atomic unit cell.
   * @param[in] model Optional custom MLIPInterface engine pointer.
   * @return Generated MLIPOutput struct.
   */
  static correlation::mlip::MLIPOutput calculate(const correlation::core::Cell &cell,
                                                 const correlation::mlip::MLIPInterface *model = nullptr);

private:
  const correlation::mlip::MLIPInterface *model_{nullptr};
};

} // namespace correlation::calculators
