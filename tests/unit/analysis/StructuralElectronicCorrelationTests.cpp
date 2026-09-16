/**
 * @file StructuralElectronicCorrelationTests.cpp
 * @brief Unit tests for StructuralElectronicCorrelation pipelines.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "analysis/DistributionFunctions.hpp"
#include "analysis/StructuralElectronicCorrelation.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "mlip/MLIPInterface.hpp"

#include <atomic>
#include <gtest/gtest.h>

namespace correlation::analysis {

namespace {

class MockElectronicModel : public correlation::mlip::MLIPInterface {
public:
  explicit MockElectronicModel(size_t bins = 20) : bins_(bins) {}

  [[nodiscard]] std::string getModelName() const override { return "MockElectronicModel"; }

  [[nodiscard]] correlation::mlip::MLIPOutput
  evaluate(const correlation::core::Cell &cell) const override {
    correlation::mlip::MLIPOutput out;
    const size_t n_atoms = cell.atomCount();
    out.ldos_bins = bins_;
    out.ldos.resize(n_atoms, std::vector<real_t>(bins_, static_cast<real_t>(0.0)));

    for (size_t i = 0; i < n_atoms; ++i) {
      for (size_t bin_idx = 0; bin_idx < bins_; ++bin_idx) {
        out.ldos[i][bin_idx] = static_cast<real_t>(1.0 + static_cast<double>(i) * 0.2 +
                                                   static_cast<double>(bin_idx) * 0.05);
      }
    }
    return out;
  }

private:
  size_t bins_{20};
};

correlation::core::Cell makeFccCell(real_t lattice_a = static_cast<real_t>(3.615)) {
  correlation::core::Cell cell({lattice_a, lattice_a, lattice_a, 90.0, 90.0, 90.0});
  cell.addAtom("Cu", {0.0, 0.0, 0.0});
  cell.addAtom("Cu",
               {0.0, lattice_a * static_cast<real_t>(0.5), lattice_a * static_cast<real_t>(0.5)});
  cell.addAtom("Cu",
               {lattice_a * static_cast<real_t>(0.5), 0.0, lattice_a * static_cast<real_t>(0.5)});
  cell.addAtom("Cu",
               {lattice_a * static_cast<real_t>(0.5), lattice_a * static_cast<real_t>(0.5), 0.0});
  return cell;
}

void verifyTDOSConservation(const MotifProjectedTDOS &res, size_t num_bins) {
  ASSERT_EQ(res.energies.size(), num_bins);
  ASSERT_EQ(res.total_tdos.size(), num_bins);
  EXPECT_EQ(res.frame_count, 1);

  for (size_t bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
    real_t sum_motifs = 0.0;
    for (const auto &[motif_name, partial] : res.motif_tdos) {
      sum_motifs += partial[bin_idx];
    }
    EXPECT_NEAR(res.total_tdos[bin_idx], sum_motifs, static_cast<real_t>(1e-5));
  }
}

} // anonymous namespace

TEST(StructuralElectronicCorrelationTests, NullModelOrEmptyTrajectoryReturnsEmpty) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Cu", {0.0, 0.0, 0.0});
  correlation::analysis::DistributionFunctions dists(cell);
  const correlation::core::Trajectory empty_traj;
  const correlation::calculators::TDOSParams null_params{.model = nullptr};

  auto res1 = StructuralElectronicCorrelation::correlateCNA(dists, empty_traj, null_params);
  EXPECT_TRUE(res1.energies.empty());
  EXPECT_TRUE(res1.total_tdos.empty());

  auto res2 = StructuralElectronicCorrelation::correlateSteinhardt(dists, empty_traj, null_params);
  EXPECT_TRUE(res2.energies.empty());
  EXPECT_TRUE(res2.total_tdos.empty());
}

TEST(StructuralElectronicCorrelationTests, CNACorrelationPreservesTotalConservation) {
  const auto cell = makeFccCell();
  correlation::core::Trajectory traj;
  traj.addFrame(cell);

  const size_t num_bins = 15;
  const MockElectronicModel model(num_bins);
  const correlation::calculators::TDOSParams params{
      .e_min = static_cast<real_t>(-10.0), .e_max = static_cast<real_t>(5.0), .model = &model};

  correlation::analysis::DistributionFunctions dists(cell);
  const auto res = StructuralElectronicCorrelation::correlateCNA(dists, traj, params);

  verifyTDOSConservation(res, num_bins);
  const auto &hist = dists.getHistogram("MotifProjectedTDOS_CNA");
  EXPECT_EQ(hist.bins.size(), num_bins);
}

TEST(StructuralElectronicCorrelationTests, SteinhardtCorrelationPreservesTotalConservation) {
  const auto cell = makeFccCell();
  correlation::core::Trajectory traj;
  traj.addFrame(cell);

  const size_t num_bins = 15;
  const MockElectronicModel model(num_bins);
  const correlation::calculators::TDOSParams params{
      .e_min = static_cast<real_t>(-10.0), .e_max = static_cast<real_t>(5.0), .model = &model};

  correlation::analysis::DistributionFunctions dists(cell);
  const auto res = StructuralElectronicCorrelation::correlateSteinhardt(dists, traj, params);

  verifyTDOSConservation(res, num_bins);
  const auto &hist = dists.getHistogram("MotifProjectedTDOS_Steinhardt");
  EXPECT_EQ(hist.bins.size(), num_bins);
}

TEST(StructuralElectronicCorrelationTests, CancellationResponsiveness) {
  const auto cell = makeFccCell();
  correlation::core::Trajectory traj;
  traj.addFrame(cell);
  traj.addFrame(cell);

  const MockElectronicModel model(10);
  const correlation::calculators::TDOSParams params{.model = &model};
  correlation::analysis::DistributionFunctions dists(cell);

  const std::atomic<bool> cancel_flag{true};
  const auto res = StructuralElectronicCorrelation::correlateCNA(dists, traj, params, &cancel_flag);
  EXPECT_TRUE(res.energies.empty());
}

} // namespace correlation::analysis
