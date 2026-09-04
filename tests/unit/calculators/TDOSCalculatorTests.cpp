/**
 * @file TDOSCalculatorTests.cpp
 * @brief Unit tests for TDOSCalculator (Total Density of States from MLIP).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/CalculatorFactory.hpp"
#include "calculators/TDOSCalculator.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "mlip/MLIPInterface.hpp"

#include <atomic>
#include <functional>
#include <gtest/gtest.h>
#include <initializer_list>

namespace correlation::calculators {

namespace {

struct MockTDOSModelConfig {
  size_t bins{20};
  real_t base_val{static_cast<real_t>(1.0)};
};

/**
 * @class MockTDOSModel
 * @brief Synthetic MLIP model producing deterministic per-atom LDOS.
 */
class MockTDOSModel : public correlation::mlip::MLIPInterface {
public:
  using Config = MockTDOSModelConfig;

  explicit MockTDOSModel(Config config = {}) : bins_(config.bins), base_val_(config.base_val) {}

  [[nodiscard]] std::string getModelName() const override { return "MockTDOSModel"; }

  [[nodiscard]] correlation::mlip::MLIPOutput evaluate(const correlation::core::Cell &cell) const override {
    correlation::mlip::MLIPOutput out;
    const size_t n_atoms = cell.atoms().size();
    out.ldos_bins = bins_;
    out.ldos.resize(n_atoms, std::vector<real_t>(bins_, static_cast<real_t>(0.0)));

    for (size_t i = 0; i < n_atoms; ++i) {
      const std::string &symbol = cell.atoms()[i].element().symbol;
      auto multiplier = static_cast<real_t>(1.0);
      if (symbol == "C") {
        multiplier = static_cast<real_t>(2.0);
      } else if (symbol == "O") {
        multiplier = static_cast<real_t>(3.0);
      } else if (symbol == "H") {
        multiplier = static_cast<real_t>(0.5);
      }

      for (size_t bin_idx = 0; bin_idx < bins_; ++bin_idx) {
        out.ldos[i][bin_idx] = base_val_ * multiplier + static_cast<real_t>(bin_idx) * static_cast<real_t>(0.1);
      }
    }
    return out;
  }

  void setBaseVal(real_t val) noexcept { base_val_ = val; }

private:
  size_t bins_{20};
  real_t base_val_{static_cast<real_t>(1.0)};
};

void verifyConservation(const std::vector<real_t> &total,
                        std::initializer_list<std::reference_wrapper<const std::vector<real_t>>> partials) {
  const size_t num_bins = total.size();
  for (size_t bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
    auto sum_species = static_cast<real_t>(0.0);
    for (const auto &partial : partials) {
      sum_species += partial.get()[bin_idx];
    }
    EXPECT_NEAR(total[bin_idx], sum_species, static_cast<real_t>(1e-5));
  }
}

} // namespace

TEST(TDOSCalculatorTests, DiscoveryInCalculatorFactory) {
  const auto &factory = CalculatorFactory::instance();
  const auto *calc = factory.getCalculator("Total Density of States (TDOS)");

  ASSERT_NE(calc, nullptr);
  EXPECT_EQ(calc->getShortName(), "TDOS");
  EXPECT_EQ(calc->getGroup(), "Machine Learning");
  EXPECT_TRUE(calc->isFrameCalculator());
  EXPECT_TRUE(calc->isTrajectoryCalculator());
}

TEST(TDOSCalculatorTests, NullModelGracefulFallback) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  const auto hist = TDOSCalculator::calculate(cell, TDOSParams{.model = nullptr});
  EXPECT_TRUE(hist.bins.empty());

  correlation::core::Trajectory traj;
  traj.addFrame(cell);
  const auto traj_hist = TDOSCalculator::calculateTrajectory(traj, TDOSParams{.model = nullptr});
  EXPECT_TRUE(traj_hist.bins.empty());
}

TEST(TDOSCalculatorTests, EnergyGridVerification) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("C", {0.0, 0.0, 0.0});

  const size_t num_bins = 50;
  const MockTDOSModel model({.bins = num_bins});
  const TDOSParams params{.e_min = static_cast<real_t>(-15.0), .e_max = static_cast<real_t>(5.0), .model = &model};

  const auto hist = TDOSCalculator::calculate(cell, params);

  ASSERT_EQ(hist.bins.size(), num_bins);
  EXPECT_EQ(hist.title, "Total Density of States (TDOS)");
  EXPECT_EQ(hist.x_unit, "eV");
  EXPECT_EQ(hist.y_unit, "states/eV/atom");

  const real_t delta = (static_cast<real_t>(5.0) - static_cast<real_t>(-15.0)) / static_cast<real_t>(num_bins);
  for (size_t bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
    const real_t expected_energy =
        static_cast<real_t>(-15.0) + (static_cast<real_t>(bin_idx) + static_cast<real_t>(0.5)) * delta;
    EXPECT_NEAR(hist.bins[bin_idx], expected_energy, static_cast<real_t>(1e-5));
  }
}

TEST(TDOSCalculatorTests, SpeciesPartialsConservation) {
  correlation::core::Cell cell({15.0, 15.0, 15.0, 90.0, 90.0, 90.0});
  cell.addAtom("C", {0.0, 0.0, 0.0});
  cell.addAtom("O", {1.0, 1.0, 1.0});
  cell.addAtom("H", {2.0, 2.0, 2.0});

  const size_t num_bins = 10;
  const MockTDOSModel model({.bins = num_bins});
  const TDOSParams params{.e_min = static_cast<real_t>(-15.0), .e_max = static_cast<real_t>(5.0), .model = &model};

  const auto hist = TDOSCalculator::calculate(cell, params);

  ASSERT_EQ(hist.bins.size(), num_bins);
  ASSERT_TRUE(hist.partials.contains("total"));
  ASSERT_TRUE(hist.partials.contains("C"));
  ASSERT_TRUE(hist.partials.contains("O"));
  ASSERT_TRUE(hist.partials.contains("H"));

  const auto &total = hist.partials.at("total");
  const auto &p_c = hist.partials.at("C");
  const auto &p_o = hist.partials.at("O");
  const auto &p_h = hist.partials.at("H");

  verifyConservation(total, {p_c, p_o, p_h});
}

TEST(TDOSCalculatorTests, TrajectoryAveraging) {
  correlation::core::Cell cell1({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell1.addAtom("C", {0.0, 0.0, 0.0});

  correlation::core::Cell cell2({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell2.addAtom("C", {0.5, 0.5, 0.5});

  correlation::core::Trajectory traj;
  traj.addFrame(cell1);
  traj.addFrame(cell2);

  const size_t num_bins = 10;
  MockTDOSModel model({.bins = num_bins});
  const TDOSParams params{.e_min = static_cast<real_t>(-15.0), .e_max = static_cast<real_t>(5.0), .model = &model};

  const auto hist1 = TDOSCalculator::calculate(cell1, params);

  const auto traj_hist = TDOSCalculator::calculateTrajectory(traj, params);

  ASSERT_EQ(traj_hist.bins.size(), num_bins);
  EXPECT_EQ(traj_hist.compute_count, 2);

  const auto &total_frame = hist1.partials.at("total");
  const auto &total_traj = traj_hist.partials.at("total");
  for (size_t bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
    EXPECT_NEAR(total_traj[bin_idx], total_frame[bin_idx], static_cast<real_t>(1e-5));
  }
}

TEST(TDOSCalculatorTests, CancellationResponsiveness) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("C", {0.0, 0.0, 0.0});

  correlation::core::Trajectory traj;
  traj.addFrame(cell);
  traj.addFrame(cell);

  const MockTDOSModel model({.bins = 10});
  const TDOSParams params{.e_min = static_cast<real_t>(-15.0), .e_max = static_cast<real_t>(5.0), .model = &model};

  const std::atomic<bool> cancel_flag{true};
  const auto hist = TDOSCalculator::calculateTrajectory(traj, params, &cancel_flag);

  EXPECT_TRUE(hist.bins.empty());
}

TEST(TDOSCalculatorTests, DistributionFunctionsIntegration) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("C", {0.0, 0.0, 0.0});

  correlation::analysis::DistributionFunctions dists(cell);
  const correlation::analysis::AnalysisSettings settings;

  const TDOSCalculator calc;
  // Without model attached, calculateFrame does not crash and leaves dists safe
  EXPECT_NO_THROW(calc.calculateFrame(dists, settings));
}

} // namespace correlation::calculators
