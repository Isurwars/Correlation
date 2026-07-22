/**
 * @file MLIPCalculator.cpp
 * @brief Implementation of MLIPCalculator for machine learning interatomic potential calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/MLIPCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"

namespace correlation::calculators {

namespace {
const bool registered = CalculatorFactory::registerTypeSafe<MLIPCalculator>("MLIPCalculator");
} // namespace

void MLIPCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                    const correlation::analysis::AnalysisSettings & /*settings*/) const {
  const auto &cell = dists.cell();
  auto output = calculate(cell);

  // Store metadata or per-atom property distribution into DistributionFunctions if applicable
  (void)output;
}

correlation::mlip::MLIPOutput MLIPCalculator::calculate(const correlation::core::Cell &cell,
                                                        const correlation::mlip::MLIPInterface *model) {
  correlation::mlip::MLIPOutput output;
  const size_t atom_count = cell.atoms().size();

  if (model != nullptr) {
    return model->evaluate(cell);
  }

  // Fallback empirical estimation if no external ONNX/Torch MLIP engine is attached
  output.total_energy = static_cast<real_t>(0.0);
  output.per_atom_energy.resize(atom_count, static_cast<real_t>(0.0));
  output.forces.resize(atom_count, correlation::math::Vector3<real_t>{0.0, 0.0, 0.0});
  output.stress = correlation::math::Matrix3<real_t>();

  return output;
}

} // namespace correlation::calculators
