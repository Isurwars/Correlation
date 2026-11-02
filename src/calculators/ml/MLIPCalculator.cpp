/**
 * @file MLIPCalculator.cpp
 * @brief Implementation of MLIPCalculator for machine learning interatomic potential calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/ml/MLIPCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "mlip/OrbV3GNNModel.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

namespace correlation::calculators {

bool MLIPCalculator::isConfigured() const {
  return model_ != nullptr || !correlation::mlip::OrbV3GNNModel::findDefaultModelPath().empty();
}

namespace {
const bool REGISTERED = CalculatorFactory::registerTypeSafe<MLIPCalculator>("MLIPCalculator");

constexpr size_t K_FORCE_BINS = 100;
constexpr real_t K_FORCE_MIN = static_cast<real_t>(0.0);
constexpr real_t K_FORCE_MAX = static_cast<real_t>(10.0);

[[nodiscard]] correlation::analysis::Histogram
buildForceHistogram(const correlation::core::Cell &cell,
                    const correlation::mlip::MLIPOutput &output) {
  correlation::analysis::Histogram hist;
  const auto &atoms = cell.atoms();
  const size_t n_atoms = atoms.size();
  if (n_atoms == 0 || output.forces.size() < n_atoms) {
    return hist;
  }

  hist.bins.resize(K_FORCE_BINS);
  const real_t delta = (K_FORCE_MAX - K_FORCE_MIN) / static_cast<real_t>(K_FORCE_BINS);
  for (size_t i = 0; i < K_FORCE_BINS; ++i) {
    hist.bins[i] = K_FORCE_MIN + (static_cast<real_t>(i) + static_cast<real_t>(0.5)) * delta;
  }

  hist.title = "MLIP Atomic Force Magnitude Distribution";
  hist.x_label = "|F|";
  hist.y_label = "Probability Density";
  hist.x_unit = "eV/A";
  hist.y_unit = "A/eV";
  hist.description = "Distribution of ORB-v3 MLIP atomic force magnitudes.";
  hist.file_suffix = "_mlip_forces";
  hist.compute_count = 1;

  auto &total = hist.partials["Total"];
  total.assign(K_FORCE_BINS, static_cast<real_t>(0.0));
  const real_t weight = static_cast<real_t>(1.0) / (static_cast<real_t>(n_atoms) * delta);

  for (size_t i = 0; i < n_atoms; ++i) {
    const real_t f_mag = correlation::math::norm(output.forces[i]);
    const auto idx =
        static_cast<size_t>(std::clamp((f_mag - K_FORCE_MIN) / delta, static_cast<real_t>(0.0),
                                       static_cast<real_t>(K_FORCE_BINS - 1)));
    total[idx] += weight;

    const std::string &sym = atoms[i].element().symbol;
    auto &spec = hist.partials[sym];
    if (spec.empty()) {
      spec.assign(K_FORCE_BINS, static_cast<real_t>(0.0));
    }
    spec[idx] += weight;
  }

  return hist;
}

} // namespace

void MLIPCalculator::calculateFrame(
    correlation::analysis::DistributionFunctions &dists,
    const correlation::analysis::AnalysisSettings & /*settings*/) const {
  if (model_ == nullptr) {
    return;
  }
  const auto &cell = dists.cell();
  const auto output = calculate(cell, model_);

  // Store metadata or per-atom property distribution into DistributionFunctions if applicable
  (void)output;
  auto force_hist = buildForceHistogram(cell, output);
  if (!force_hist.bins.empty()) {
    dists.addHistogram("MLIP_Forces", std::move(force_hist));
  }
}

correlation::mlip::MLIPOutput
MLIPCalculator::calculate(const correlation::core::Cell &cell,
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
