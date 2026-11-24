/**
 * @file LDOSPostProcessorTests.cpp
 * @brief Unit tests for electronic LDOS/TDOS post-processing utilities.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/scattering/LDOSPostProcessor.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <vector>

using correlation::real_t;
using correlation::analysis::LDOSBroadeningKernel;

namespace correlation::calculators {

namespace {

[[nodiscard]] std::vector<real_t> makeUniformGrid(real_t e_min, real_t e_max, size_t count) {
  std::vector<real_t> grid(count);
  const real_t step = (e_max - e_min) / static_cast<real_t>(count - 1);
  for (size_t i = 0; i < count; ++i) {
    grid[i] = e_min + static_cast<real_t>(i) * step;
  }
  return grid;
}

[[nodiscard]] double trapezoidalArea(const std::vector<real_t> &grid,
                                     const std::vector<real_t> &values) {
  double area = 0.0;
  for (size_t i = 0; i + 1 < grid.size(); ++i) {
    area += 0.5 * (static_cast<double>(values[i]) + static_cast<double>(values[i + 1])) *
            (static_cast<double>(grid[i + 1]) - static_cast<double>(grid[i]));
  }
  return area;
}

} // namespace

TEST(LDOSPostProcessorTests, GaussianAndLorentzianBroadeningConserveArea) {
  const auto energies = makeUniformGrid(static_cast<real_t>(-10.0), static_cast<real_t>(10.0), 401);
  std::vector<real_t> dos(energies.size(), static_cast<real_t>(0.0));
  for (size_t i = 0; i < energies.size(); ++i) {
    const auto energy = static_cast<double>(energies[i]);
    dos[i] = static_cast<real_t>(std::exp(-2.0 * (energy + 1.5) * (energy + 1.5)) +
                                 0.6 * std::exp(-3.0 * (energy - 2.0) * (energy - 2.0)));
  }

  const double orig_area = trapezoidalArea(energies, dos);
  const auto gauss = LDOSPostProcessor::broadenSpectrum(
      energies, dos, LDOSBroadeningKernel::Gaussian, static_cast<real_t>(0.4));
  const auto lorentz = LDOSPostProcessor::broadenSpectrum(
      energies, dos, LDOSBroadeningKernel::Lorentzian, static_cast<real_t>(0.3));

  ASSERT_EQ(gauss.size(), dos.size());
  ASSERT_EQ(lorentz.size(), dos.size());
  EXPECT_NEAR(trapezoidalArea(energies, gauss), orig_area, 1e-3);
  EXPECT_NEAR(trapezoidalArea(energies, lorentz), orig_area, 1e-3);
}

TEST(LDOSPostProcessorTests, RebinSpectrumPreservesLinearFunction) {
  const auto src_e = makeUniformGrid(static_cast<real_t>(-5.0), static_cast<real_t>(5.0), 101);
  std::vector<real_t> src_dos(src_e.size());
  for (size_t i = 0; i < src_e.size(); ++i) {
    src_dos[i] = static_cast<real_t>(2.0) + static_cast<real_t>(0.2) * src_e[i];
  }

  const auto dst_e = makeUniformGrid(static_cast<real_t>(-4.0), static_cast<real_t>(4.0), 41);
  const auto rebinned = LDOSPostProcessor::rebinSpectrum(src_e, src_dos, dst_e);
  ASSERT_EQ(rebinned.size(), dst_e.size());
  for (size_t i = 0; i < dst_e.size(); ++i) {
    const real_t expected = static_cast<real_t>(2.0) + static_cast<real_t>(0.2) * dst_e[i];
    EXPECT_NEAR(rebinned[i], expected, static_cast<real_t>(1e-4));
  }
}

TEST(LDOSPostProcessorTests, FermiDosExactAndThermalSmearing) {
  const auto energies = makeUniformGrid(static_cast<real_t>(-5.0), static_cast<real_t>(5.0), 201);
  const std::vector<real_t> dos(energies.size(), static_cast<real_t>(1.75));

  const real_t exact_ef =
      LDOSPostProcessor::evaluateFermiDos(energies, dos, static_cast<real_t>(0.0));
  const real_t thermal_ef =
      LDOSPostProcessor::evaluateFermiDos(energies, dos, static_cast<real_t>(300.0));

  EXPECT_NEAR(exact_ef, static_cast<real_t>(1.75), static_cast<real_t>(1e-4));
  EXPECT_NEAR(thermal_ef, static_cast<real_t>(1.75), static_cast<real_t>(1e-3));
}

TEST(LDOSPostProcessorTests, BandGapDetectionMetallic) {
  const auto energies = makeUniformGrid(static_cast<real_t>(-6.0), static_cast<real_t>(6.0), 241);
  const std::vector<real_t> metal_dos(energies.size(), static_cast<real_t>(0.8));
  const auto metal_gap =
      LDOSPostProcessor::detectBandGap(energies, metal_dos, static_cast<real_t>(1e-3));
  EXPECT_TRUE(metal_gap.is_metallic);
  EXPECT_NEAR(metal_gap.band_gap, static_cast<real_t>(0.0), static_cast<real_t>(1e-6));
}

TEST(LDOSPostProcessorTests, BandGapDetectionSemiconductor) {
  const auto energies = makeUniformGrid(static_cast<real_t>(-6.0), static_cast<real_t>(6.0), 241);
  std::vector<real_t> semi_dos(energies.size(), static_cast<real_t>(0.0));
  for (size_t i = 0; i < energies.size(); ++i) {
    if (energies[i] <= static_cast<real_t>(-0.75) || energies[i] >= static_cast<real_t>(1.25)) {
      semi_dos[i] = static_cast<real_t>(1.5);
    }
  }
  const auto semi_gap =
      LDOSPostProcessor::detectBandGap(energies, semi_dos, static_cast<real_t>(1e-3));
  EXPECT_FALSE(semi_gap.is_metallic);
  EXPECT_NEAR(semi_gap.vbm_energy, static_cast<real_t>(-0.75), static_cast<real_t>(0.06));
  EXPECT_NEAR(semi_gap.cbm_energy, static_cast<real_t>(1.25), static_cast<real_t>(0.06));
  EXPECT_NEAR(semi_gap.band_gap, static_cast<real_t>(2.0), static_cast<real_t>(0.12));
}

TEST(LDOSPostProcessorTests, IntegratedDosMonotonicAndSumRule) {
  const auto energies = makeUniformGrid(static_cast<real_t>(0.0), static_cast<real_t>(4.0), 81);
  const std::vector<real_t> dos(energies.size(), static_cast<real_t>(2.5));

  const auto idos = LDOSPostProcessor::computeIntegratedDos(energies, dos);
  ASSERT_EQ(idos.size(), energies.size());
  EXPECT_NEAR(idos.front(), static_cast<real_t>(0.0), static_cast<real_t>(1e-6));
  EXPECT_NEAR(idos.back(), static_cast<real_t>(10.0), static_cast<real_t>(1e-4));
  for (size_t i = 1; i < idos.size(); ++i) {
    EXPECT_GE(idos[i], idos[i - 1]);
  }
}

} // namespace correlation::calculators
