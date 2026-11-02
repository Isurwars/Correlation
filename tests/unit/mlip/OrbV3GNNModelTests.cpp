/**
 * @file OrbV3GNNModelTests.cpp
 * @brief Unit tests for native C++23 OrbV3GNNModel and Cellulose eDoS inference.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "analysis/CorrelationEngine.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "mlip/OrbV3GNNModel.hpp"
#include "mlip/PeriodicGraphBuilder.hpp"

#include <array>
#include <cmath>
#include <filesystem>
#include <gtest/gtest.h>
#include <span>
#include <stdexcept>
#include <string>

namespace correlation::mlip {

namespace {

[[nodiscard]] std::string findTestFixturePath() {
  const std::array<std::filesystem::path, 3> candidates = {
      "tests/data/mlip/mini_ldos.bin",
      "../tests/data/mlip/mini_ldos.bin",
      "../../tests/data/mlip/mini_ldos.bin",
  };
  for (const auto &cand : candidates) {
    if (std::filesystem::exists(cand)) {
      return std::filesystem::canonical(cand).string();
    }
  }
  return "";
}

[[nodiscard]] correlation::core::Cell buildCelluloseSnippetCell() {
  correlation::core::Cell cell({8.0, 8.0, 8.0, 90.0, 90.0, 90.0});
  cell.addAtom("C", {1.0, 1.0, 1.0});
  cell.addAtom("O", {2.4, 1.0, 1.0});
  cell.addAtom("H", {0.4, 1.8, 1.0});
  cell.addAtom("H", {2.9, 1.7, 1.0});
  return cell;
}

void verifyZeroNetForce(std::span<const correlation::math::Vector3<real_t>> forces) {
  correlation::math::Vector3<real_t> net_force{0.0, 0.0, 0.0};
  bool all_finite = true;
  for (const auto &force : forces) {
    all_finite = all_finite && std::isfinite(force.x()) && std::isfinite(force.y()) &&
                 std::isfinite(force.z());
    net_force += force;
  }
  EXPECT_TRUE(all_finite);
  EXPECT_LT(correlation::math::norm(net_force), static_cast<real_t>(1e-4));
}

void verifyNonNegativeLdos(std::span<const real_t> ldos) {
  bool valid = true;
  for (const real_t val : ldos) {
    valid = valid && std::isfinite(val) && (val >= static_cast<real_t>(0.0));
  }
  EXPECT_TRUE(valid);
}

[[nodiscard]] real_t maxAbsDiff(std::span<const real_t> actual, std::span<const real_t> expected) {
  auto max_err = static_cast<real_t>(0.0);
  for (size_t idx = 0; idx < expected.size(); ++idx) {
    max_err = std::max(max_err, std::abs(actual[idx] - expected[idx]));
  }
  return max_err;
}

void verifyParityForces(std::span<const correlation::math::Vector3<real_t>> forces) {
  ASSERT_EQ(forces.size(), 4U);
  const correlation::math::Vector3<real_t> ref_f0{
      static_cast<real_t>(2.230072), static_cast<real_t>(-4.376476), static_cast<real_t>(0.043121)};
  const correlation::math::Vector3<real_t> ref_f1{static_cast<real_t>(-4.991976),
                                                  static_cast<real_t>(-5.749089),
                                                  static_cast<real_t>(0.045077)};
  EXPECT_LT(correlation::math::norm(forces[0] - ref_f0), static_cast<real_t>(1e-3));
  EXPECT_LT(correlation::math::norm(forces[1] - ref_f1), static_cast<real_t>(1e-3));
}

void verifyParityLdos(const MLIPOutput &out) {
  ASSERT_EQ(out.ldos_bins, 256U);
  const auto ldos_0 = out.getLdos(0);
  const auto ldos_1 = out.getLdos(1);
  ASSERT_EQ(ldos_0.size(), 256U);
  ASSERT_EQ(ldos_1.size(), 256U);

  const std::array<real_t, 5> ref_0 = {
      static_cast<real_t>(0.001807), static_cast<real_t>(0.002277), static_cast<real_t>(0.001721),
      static_cast<real_t>(0.002137), static_cast<real_t>(0.002075),
  };
  const std::array<real_t, 5> ref_1 = {
      static_cast<real_t>(0.072283), static_cast<real_t>(0.056413), static_cast<real_t>(0.069160),
      static_cast<real_t>(0.055084), static_cast<real_t>(0.032377),
  };

  EXPECT_LT(maxAbsDiff(ldos_0.subspan(0, 5), ref_0), static_cast<real_t>(1e-4));
  EXPECT_LT(maxAbsDiff(ldos_1.subspan(128, 5), ref_1), static_cast<real_t>(1e-4));
}

} // namespace

TEST(OrbV3GNNModelTests, LoadMiniFixtureAndValidateHeader) {
  const std::string fixture = findTestFixturePath();
  ASSERT_FALSE(fixture.empty()) << "mini_ldos.bin test fixture not found.";

  OrbV3GNNModel const model(fixture);
  ASSERT_TRUE(model.isLoaded());

  const auto &hdr = model.header();
  EXPECT_EQ(hdr.version, 1U);
  EXPECT_EQ(hdr.latent_dim, 16U);
  EXPECT_EQ(hdr.mlp_hidden_dim, 32U);
  EXPECT_EQ(hdr.num_stacks, 2U);
  EXPECT_EQ(hdr.num_rbf, 8U);
  EXPECT_EQ(hdr.l_max, 3U);
  EXPECT_EQ(hdr.ldos_bins, 32U);
  EXPECT_EQ(hdr.max_neighbors, 120U);
  EXPECT_NEAR(hdr.r_max, 6.0F, 1e-5F);
  EXPECT_NEAR(hdr.ldos_emin, -30.0F, 1e-5F);
  EXPECT_NEAR(hdr.ldos_emax, 15.0F, 1e-5F);
}

TEST(OrbV3GNNModelTests, EvaluateCelluloseSnippetProducesFiniteEnergyForcesStressAndLDOS) {
  const std::string fixture = findTestFixturePath();
  ASSERT_FALSE(fixture.empty());

  OrbV3GNNModel const model(fixture);
  const auto cell = buildCelluloseSnippetCell();

  const auto out = model.evaluate(cell);
  EXPECT_TRUE(std::isfinite(out.total_energy));
  ASSERT_EQ(out.per_atom_energy.size(), 4U);
  ASSERT_EQ(out.forces.size(), 4U);
  EXPECT_EQ(out.ldos_bins, 32U);
  ASSERT_EQ(out.ldos.size(), 4U * 32U);
  EXPECT_NEAR(out.ldos_e_min, static_cast<real_t>(-30.0), 1e-5);
  EXPECT_NEAR(out.ldos_e_max, static_cast<real_t>(15.0), 1e-5);

  verifyZeroNetForce(out.forces);
  verifyNonNegativeLdos(out.ldos);
}

TEST(OrbV3GNNModelTests, RejectsUnmappedNumericSpeciesWithActionableError) {
  const std::string fixture = findTestFixturePath();
  ASSERT_FALSE(fixture.empty());

  OrbV3GNNModel const model(fixture);
  correlation::core::Cell unmapped_cell({8.0, 8.0, 8.0, 90.0, 90.0, 90.0});
  unmapped_cell.addAtom("1", {1.0, 1.0, 1.0});
  unmapped_cell.addAtom("2", {2.5, 1.0, 1.0});

  const auto report = PeriodicGraphBuilder::validateElementCompatibility(unmapped_cell);
  EXPECT_TRUE(report.has_unmapped_elements);
  ASSERT_EQ(report.unmapped_symbols.size(), 2U);

  EXPECT_THROW((void)model.evaluate(unmapped_cell), std::runtime_error);
}

TEST(OrbV3GNNModelTests, AllowsNonCelluloseValidElementsWithCalibrationReport) {
  const std::string fixture = findTestFixturePath();
  ASSERT_FALSE(fixture.empty());

  OrbV3GNNModel const model(fixture);
  correlation::core::Cell si_cell({8.0, 8.0, 8.0, 90.0, 90.0, 90.0});
  si_cell.addAtom("Si", {1.0, 1.0, 1.0});
  si_cell.addAtom("O", {2.6, 1.0, 1.0});

  const auto report = PeriodicGraphBuilder::validateElementCompatibility(si_cell);
  EXPECT_FALSE(report.has_unmapped_elements);
  EXPECT_TRUE(report.has_out_of_domain_elements);
  ASSERT_EQ(report.out_of_domain_symbols.size(), 1U);
  EXPECT_EQ(report.out_of_domain_symbols[0], "Si");

  const auto out = model.evaluate(si_cell);
  EXPECT_EQ(out.ldos.size(), 2U * 32U);
}

TEST(OrbV3GNNModelTests, EndToEndTDOSAndCorrelationEnginePipeline) {
  const std::string fixture = findTestFixturePath();
  ASSERT_FALSE(fixture.empty());

  correlation::core::Trajectory traj;
  traj.addFrame(buildCelluloseSnippetCell());

  correlation::analysis::CorrelationEngineConfig config;
  config.settings.mlip_model_path = fixture;
  config.settings.active_calculators["TDOS"] = true;
  config.settings.active_calculators["MLIP"] = true;
  config.settings.smoothing = false;

  auto res = correlation::analysis::CorrelationEngine::runAnalysis(traj, config);
  ASSERT_TRUE(res.has_value()) << res.error();
  const auto &dists = *res.value();

  ASSERT_TRUE(dists.hasHistogram("TDOS"));
  const auto &tdos = dists.getHistogram("TDOS");
  EXPECT_EQ(tdos.bins.size(), 32U);
  EXPECT_TRUE(tdos.partials.contains("total"));
  EXPECT_TRUE(tdos.partials.contains("C"));
  EXPECT_TRUE(tdos.partials.contains("O"));
  EXPECT_TRUE(tdos.partials.contains("H"));

  ASSERT_TRUE(dists.hasHistogram("MLIP_Forces"));
}

TEST(OrbV3GNNModelTests, NumericalParityAgainstPyTorchCelluloseModel) {
  const std::string full_model_path = OrbV3GNNModel::findDefaultModelPath();
  if (full_model_path.empty()) {
    GTEST_SKIP() << "models/cellulose_ldos.bin not present; skipping full checkpoint parity test.";
  }

  OrbV3GNNModel const model(full_model_path);
  ASSERT_TRUE(model.isLoaded());

  correlation::core::Cell cell({15.0, 15.0, 15.0, 90.0, 90.0, 90.0});
  cell.addAtom("C", {1.0, 1.0, 1.0});
  cell.addAtom("O", {2.4, 1.0, 1.0});
  cell.addAtom("H", {0.4, 1.8, 1.0});
  cell.addAtom("H", {2.9, 1.7, 1.0});

  const auto out = model.evaluate(cell);

  // Reference values computed via Python ASEFullHeadCalculator on cellulose_full_heads.pt
  EXPECT_NEAR(out.total_energy, static_cast<real_t>(-624.14917), static_cast<real_t>(1e-2));
  verifyParityForces(out.forces);
  verifyParityLdos(out);
}

} // namespace correlation::mlip
