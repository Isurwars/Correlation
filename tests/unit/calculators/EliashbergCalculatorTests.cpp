/**
 * @file EliashbergCalculatorTests.cpp
 * @brief Unit tests and physical sum-rule verification for EliashbergCalculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/CalculatorFactory.hpp"
#include "calculators/ml/EliashbergCalculator.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "mlip/MLIPInterface.hpp"

#include <atomic>
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

using correlation::real_t;
using correlation::analysis::CouplingEstimationMode;
using correlation::analysis::Histogram;

namespace correlation::calculators {

namespace {

constexpr double K_TWO_PI = 6.283185307179586;

class PositionDependentMockMLIP final : public correlation::mlip::MLIPInterface {
public:
  explicit PositionDependentMockMLIP(real_t slope = static_cast<real_t>(0.5)) : slope_(slope) {}

  [[nodiscard]] std::string getModelName() const override { return "PositionDependentMockMLIP"; }

  [[nodiscard]] correlation::mlip::MLIPOutput
  evaluate(const correlation::core::Cell &cell) const override {
    correlation::mlip::MLIPOutput out;
    const size_t n_atoms = cell.atoms().size();
    constexpr size_t K_BINS = 41;
    out.ldos_bins = K_BINS;
    out.ldos_e_min = static_cast<real_t>(-10.0);
    out.ldos_e_max = static_cast<real_t>(10.0);
    out.ldos.resize(n_atoms * K_BINS, static_cast<real_t>(0.0));

    for (size_t i = 0; i < n_atoms; ++i) {
      const real_t pos_x = cell.atoms()[i].position().x();
      const real_t n_ef =
          std::max(static_cast<real_t>(0.2), static_cast<real_t>(1.5) + slope_ * pos_x);
      for (size_t bin_idx = 0; bin_idx < K_BINS; ++bin_idx) {
        out.ldos[i * K_BINS + bin_idx] = n_ef;
      }
    }
    return out;
  }

private:
  real_t slope_{static_cast<real_t>(0.5)};
};

[[nodiscard]] Histogram makeSyntheticVdos(const std::string &species, real_t peak_thz) {
  Histogram hist;
  constexpr size_t K_PTS = 201;
  hist.bins.resize(K_PTS);
  std::vector<real_t> total(K_PTS, static_cast<real_t>(0.0));
  for (size_t i = 0; i < K_PTS; ++i) {
    const double freq_nu = static_cast<double>(i) * 0.1;
    hist.bins[i] = static_cast<real_t>(freq_nu);
    if (i > 0) {
      const double diff = freq_nu - static_cast<double>(peak_thz);
      total[i] = static_cast<real_t>(std::exp(-0.5 * diff * diff) / std::sqrt(K_TWO_PI));
    }
  }
  hist.partials["Total"] = total;
  hist.partials[species] = total;
  return hist;
}

} // namespace

TEST(EliashbergCalculatorTests, DiscoveryInCalculatorFactory) {
  const auto &factory = CalculatorFactory::instance();
  const auto *calc = factory.getCalculator("Eliashberg Spectral Function");
  ASSERT_NE(calc, nullptr);
  EXPECT_EQ(calc->getShortName(), "a2F");
  EXPECT_EQ(calc->getGroup(), "Machine Learning");
  EXPECT_FALSE(calc->isFrameCalculator());
  EXPECT_TRUE(calc->isTrajectoryCalculator());
}

TEST(EliashbergCalculatorTests, AnalyticalEinsteinModeSumRules) {
  constexpr double K_OMEGA0 = 15.0;
  constexpr double K_AREA = 7.5;
  constexpr double K_SIGMA = 0.25;
  constexpr size_t K_PTS = 2001;

  std::vector<real_t> omega(K_PTS);
  std::vector<real_t> a2f(K_PTS);
  const double norm = K_AREA / (K_SIGMA * std::sqrt(K_TWO_PI));
  for (size_t i = 0; i < K_PTS; ++i) {
    const double freq_w = static_cast<double>(i) * 0.02;
    omega[i] = static_cast<real_t>(freq_w);
    const double z_score = (freq_w - K_OMEGA0) / K_SIGMA;
    a2f[i] = static_cast<real_t>(norm * std::exp(-0.5 * z_score * z_score));
  }

  const auto diag = EliashbergCalculator::computeFromSpectrum(omega, a2f, static_cast<real_t>(0.10),
                                                              static_cast<real_t>(0.5));

  EXPECT_TRUE(diag.has_eliashberg);
  EXPECT_NEAR(diag.lambda, static_cast<real_t>(1.0), static_cast<real_t>(0.01));
  EXPECT_NEAR(diag.omega_log_mev, static_cast<real_t>(K_OMEGA0), static_cast<real_t>(0.15));
  EXPECT_NEAR(diag.omega_2_mev, static_cast<real_t>(K_OMEGA0), static_cast<real_t>(0.15));
  EXPECT_GT(diag.tc_allen_dynes_k, diag.tc_mcmillan_k);
}

TEST(EliashbergCalculatorTests, AllenDynesLiteratureCalibration) {
  constexpr double K_OMEGA0 = 4.826;
  constexpr double K_LAMBDA_TARGET = 1.55;
  constexpr double K_AREA = 0.5 * K_LAMBDA_TARGET * K_OMEGA0;
  constexpr double K_SIGMA = 0.15;
  constexpr size_t K_PTS = 1501;

  std::vector<real_t> omega(K_PTS);
  std::vector<real_t> a2f(K_PTS);
  const double norm = K_AREA / (K_SIGMA * std::sqrt(K_TWO_PI));
  for (size_t i = 0; i < K_PTS; ++i) {
    const double freq_w = static_cast<double>(i) * 0.01;
    omega[i] = static_cast<real_t>(freq_w);
    const double z_score = (freq_w - K_OMEGA0) / K_SIGMA;
    a2f[i] = static_cast<real_t>(norm * std::exp(-0.5 * z_score * z_score));
  }

  const auto diag = EliashbergCalculator::computeFromSpectrum(
      omega, a2f, static_cast<real_t>(0.105), static_cast<real_t>(0.5));
  EXPECT_NEAR(diag.lambda, static_cast<real_t>(1.55), static_cast<real_t>(0.02));
  EXPECT_NEAR(diag.tc_allen_dynes_k, static_cast<real_t>(7.2), static_cast<real_t>(0.8));

  const std::vector<real_t> weak_a2f(K_PTS, static_cast<real_t>(1e-4));
  const auto weak_diag = EliashbergCalculator::computeFromSpectrum(
      omega, weak_a2f, static_cast<real_t>(0.15), static_cast<real_t>(0.5));
  EXPECT_EQ(weak_diag.tc_mcmillan_k, static_cast<real_t>(0.0));
  EXPECT_EQ(weak_diag.tc_allen_dynes_k, static_cast<real_t>(0.0));
}

TEST(EliashbergCalculatorTests, TrajectoryAndFiniteDisplacementCouplingModes) {
  correlation::core::Trajectory traj;
  for (int frame_idx = 0; frame_idx < 6; ++frame_idx) {
    correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
    const real_t delta_x = static_cast<real_t>(frame_idx) * static_cast<real_t>(0.05);
    cell.addAtom("Nb", {delta_x, static_cast<real_t>(0.0), static_cast<real_t>(0.0)});
    cell.addAtom("Nb", {static_cast<real_t>(5.0) - delta_x, static_cast<real_t>(5.0),
                        static_cast<real_t>(5.0)});
    traj.addFrame(cell);
  }

  const PositionDependentMockMLIP model(static_cast<real_t>(2.0));
  const auto vdos = makeSyntheticVdos("Nb", static_cast<real_t>(5.0));

  EliashbergParams traj_params;
  traj_params.model = &model;
  traj_params.coupling_mode = CouplingEstimationMode::TrajectoryFluctuation;
  const auto res_traj = EliashbergCalculator::calculate(traj, vdos, traj_params);

  ASSERT_TRUE(res_traj.diagnostics.has_eliashberg);
  ASSERT_FALSE(res_traj.spectral_function.bins.empty());
  EXPECT_GT(res_traj.diagnostics.lambda, static_cast<real_t>(0.0));
  EXPECT_GT(res_traj.coupling_weights.eta_by_species.at("Nb"), static_cast<real_t>(0.0));

  EliashbergParams fd_params;
  fd_params.model = &model;
  fd_params.coupling_mode = CouplingEstimationMode::FiniteDisplacement;
  fd_params.displacement_step = static_cast<real_t>(0.02);
  const auto res_fd = EliashbergCalculator::calculate(traj, vdos, fd_params);

  ASSERT_TRUE(res_fd.diagnostics.has_eliashberg);
  EXPECT_GT(res_fd.diagnostics.lambda, static_cast<real_t>(0.0));
  EXPECT_NEAR(res_traj.coupling_weights.eta_by_species.at("Nb"),
              res_fd.coupling_weights.eta_by_species.at("Nb"), static_cast<real_t>(0.2));

  const std::atomic<bool> cancel{true};
  const auto res_cancel =
      EliashbergCalculator::calculate(traj, vdos, traj_params, nullptr, &cancel);
  EXPECT_TRUE(res_cancel.spectral_function.bins.empty());
}

} // namespace correlation::calculators
