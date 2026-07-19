// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/VoronoiCalculator.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::analysis {
namespace {

class VoronoiCalculatorTests : public ::testing::Test {
protected:
  // Helper to extract the peak bin index for a histogram
  static double getPeakValue(const Histogram &hist) {
    const auto &vals = hist.partials.at("Total");
    double max_val = -1.0;
    size_t max_bin = 0;
    for (size_t i = 0; i < vals.size(); ++i) {
      if (vals[i] > max_val) {
        max_val = vals[i];
        max_bin = i;
      }
    }
    return hist.bins[max_bin];
  }
};

TEST_F(VoronoiCalculatorTests, SimpleCubic) {
  // 1. Simple Cubic (SC) Structure
  // Box length 12.0 with SC lattice spacing 4.0 (27 atoms)
  correlation::core::Cell cell_sc({12.0, 12.0, 12.0, 90.0, 90.0, 90.0});
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        cell_sc.addAtom("Ar",
                        {static_cast<double>(i * 4.0), static_cast<double>(j * 4.0), static_cast<double>(k * 4.0)});
      }
    }
  }

  auto hists_sc = correlation::calculators::VoronoiCalculator::calculate(cell_sc);

  // Theoretical SC cell volume = 4.0^3 = 64.0 A^3
  double const vol_sc = getPeakValue(hists_sc.at("Voronoi Volume"));
  EXPECT_NEAR(vol_sc, 64.0, 0.5);

  // SC has coordination number of 6 (6 faces)
  double const cn_sc = getPeakValue(hists_sc.at("Voronoi Coordination Number"));
  EXPECT_NEAR(cn_sc, 6.0, 1e-4);

  // SC sphericity Psi = (pi^(1/3) * (6 * 64)^(2/3)) / 96 = pi^(1/3) * 6^(2/3) / 6 approx 0.806
  double const sph_sc = getPeakValue(hists_sc.at("Voronoi Sphericity"));
  EXPECT_NEAR(sph_sc, 0.806, 0.01);

  // SC polyhedral signature is (0, 6, 0, 0)
  EXPECT_TRUE(hists_sc.at("Voronoi Signatures").description.contains("(0, 6, 0, 0)"));
}

TEST_F(VoronoiCalculatorTests, BodyCenteredCubic) {
  // 2. Body-Centered Cubic (BCC) Structure
  // Box length 10.0 with BCC lattice spacing 5.0 (16 atoms)
  correlation::core::Cell cell_bcc({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        cell_bcc.addAtom("Ar",
                         {static_cast<double>(i * 5.0), static_cast<double>(j * 5.0), static_cast<double>(k * 5.0)});
        cell_bcc.addAtom("Ar", {i * 5.0 + 2.5, j * 5.0 + 2.5, k * 5.0 + 2.5});
      }
    }
  }

  auto hists_bcc = correlation::calculators::VoronoiCalculator::calculate(cell_bcc);

  // Theoretical BCC cell volume = 1000 / 16 = 62.5 A^3
  double const vol_bcc = getPeakValue(hists_bcc.at("Voronoi Volume"));
  EXPECT_NEAR(vol_bcc, 62.5, 0.5);

  // BCC has coordination number of 14 (14 faces: 8 nearest + 6 next-nearest)
  double const cn_bcc = getPeakValue(hists_bcc.at("Voronoi Coordination Number"));
  EXPECT_NEAR(cn_bcc, 14.0, 1e-4);

  // BCC sphericity Psi approx 0.905
  double const sph_bcc = getPeakValue(hists_bcc.at("Voronoi Sphericity"));
  EXPECT_NEAR(sph_bcc, 0.905, 0.01);

  // BCC signature is (0, 6, 0, 8)
  EXPECT_TRUE(hists_bcc.at("Voronoi Signatures").description.contains("(0, 6, 0, 8)"));
}

TEST_F(VoronoiCalculatorTests, FaceCenteredCubic) {
  // 3. Face-Centered Cubic (FCC) Structure
  // Box length 10.0 with FCC lattice spacing 5.0 (32 atoms)
  correlation::core::Cell cell_fcc({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        cell_fcc.addAtom("Ar", {i * 5.0, j * 5.0, k * 5.0});
        cell_fcc.addAtom("Ar", {i * 5.0 + 2.5, j * 5.0 + 2.5, k * 5.0});
        cell_fcc.addAtom("Ar", {i * 5.0 + 2.5, j * 5.0, k * 5.0 + 2.5});
        cell_fcc.addAtom("Ar", {i * 5.0, j * 5.0 + 2.5, k * 5.0 + 2.5});
      }
    }
  }

  auto hists_fcc = correlation::calculators::VoronoiCalculator::calculate(cell_fcc);

  // Theoretical FCC cell volume = 1000 / 32 = 31.25 A^3
  double const vol_fcc = getPeakValue(hists_fcc.at("Voronoi Volume"));
  EXPECT_NEAR(vol_fcc, 31.25, 0.5);

  // FCC has coordination number of 12 (12 faces)
  double const cn_fcc = getPeakValue(hists_fcc.at("Voronoi Coordination Number"));
  EXPECT_NEAR(cn_fcc, 12.0, 1e-4);

  // FCC sphericity Psi approx 0.905
  double const sph_fcc = getPeakValue(hists_fcc.at("Voronoi Sphericity"));
  EXPECT_NEAR(sph_fcc, 0.905, 0.01);

  // FCC signature is (0, 12, 0, 0)
  EXPECT_TRUE(hists_fcc.at("Voronoi Signatures").description.contains("(0, 12, 0, 0)"));
}

TEST_F(VoronoiCalculatorTests, HexagonalClosePacked) {
  // 4. Hexagonal Close-Packed (HCP) Structure
  // Supercell size 3x3x2 using lattice parameters (36 atoms)
  correlation::core::Cell cell_hcp({9.0, 9.0, 9.79795897, 90.0, 90.0, 120.0});
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 2; ++k) {
        double fx1 = (static_cast<double>(i) + 0.0) / 3.0;
        double fy1 = (static_cast<double>(j) + 0.0) / 3.0;
        double fz1 = (static_cast<double>(k) + 0.0) / 2.0;

        double fx2 = (static_cast<double>(i) + 1.0 / 3.0) / 3.0;
        double fy2 = (static_cast<double>(j) + 2.0 / 3.0) / 3.0;
        double fz2 = (static_cast<double>(k) + 0.5) / 2.0;

        auto pos1 = cell_hcp.latticeVectors() * correlation::math::Vector3<double>(fx1, fy1, fz1);
        auto pos2 = cell_hcp.latticeVectors() * correlation::math::Vector3<double>(fx2, fy2, fz2);

        cell_hcp.addAtom("Ar", pos1);
        cell_hcp.addAtom("Ar", pos2);
      }
    }
  }

  auto hists_hcp = correlation::calculators::VoronoiCalculator::calculate(cell_hcp);

  // Theoretical HCP cell volume = volume / 36 approx 19.0914 A^3
  double const vol_hcp = getPeakValue(hists_hcp.at("Voronoi Volume"));
  EXPECT_NEAR(vol_hcp, cell_hcp.volume() / 36.0, 0.5);

  // HCP has coordination number of 12 (12 faces)
  double const cn_hcp = getPeakValue(hists_hcp.at("Voronoi Coordination Number"));
  EXPECT_NEAR(cn_hcp, 12.0, 1e-4);

  // HCP sphericity Psi approx 0.905
  double const sph_hcp = getPeakValue(hists_hcp.at("Voronoi Sphericity"));
  EXPECT_NEAR(sph_hcp, 0.905, 0.01);

  // HCP signature is (0, 12, 0, 0) in double precision, or (0, 2, 4, 6) in single precision
  if (correlation::is_single_precision) {
    EXPECT_TRUE(hists_hcp.at("Voronoi Signatures").description.contains("(0, 2, 4, 6)"));
  } else {
    EXPECT_TRUE(hists_hcp.at("Voronoi Signatures").description.contains("(0, 12, 0, 0)"));
  }
}

} // namespace
} // namespace correlation::analysis
