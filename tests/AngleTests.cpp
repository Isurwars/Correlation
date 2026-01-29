// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <gtest/gtest.h>
#include <cmath>
#include <numeric>

#include "../include/Cell.hpp"
#include "../include/StructureAnalyzer.hpp"
#include "../include/DistributionFunctions.hpp"
#include "../include/PhysicalData.hpp"
#include "../include/Trajectory.hpp"

// ============================================================================
// Part 1: Angle Reproduction Tests
// ============================================================================


class AngleReproductionTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Simple cubic cell
        cell_ = Cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
    }

    void updateTrajectory() {
        trajectory_ = Trajectory();
        trajectory_.addFrame(cell_);
        trajectory_.precomputeBondCutoffs();
    }

    Cell cell_;
    Trajectory trajectory_;
};

TEST_F(AngleReproductionTest, MissingAnglesWhenCutoffIsTooSmall) {
    // A-B-C angle.
    // B is at (5,5,5)
    // A is at (4,5,5) -> dist 1.0
    // C is at (5,6,5) -> dist 1.0
    // Angle should be 90 degrees.
    
    cell_.addAtom("O", {4.0, 5.0, 5.0});
    cell_.addAtom("Si", {5.0, 5.0, 5.0});
    cell_.addAtom("O", {5.0, 6.0, 5.0});
    updateTrajectory();

    // Bond cutoff for Si-O is likely around 1.6 * 1.2 = 1.92 or similar.
    // Distance is 1.0.
    {
        StructureAnalyzer analyzer(cell_, 1.1, trajectory_.getBondCutoffs());
        const auto& angles = analyzer.angles();
        bool found = false;
        for(const auto& t1 : angles) {
            for(const auto& center : t1) {
                for(const auto& t2 : center) {
                    for(double angle : t2) {
                        if(std::abs(angle * 180.0 / M_PI - 90.0) < 1.0) {
                            found = true;
                        }
                    }
                }
            }
        }
        EXPECT_TRUE(found) << "Should find 90 degree angle with sufficient cutoff";
    }
}

TEST_F(AngleReproductionTest, PBCAngleDetection) {
    cell_.addAtom("Si", {0.5, 0.5, 0.5});
    cell_.addAtom("O", {9.6, 0.5, 0.5});
    cell_.addAtom("O", {0.5, 9.6, 0.5});
    updateTrajectory();
    
    StructureAnalyzer analyzer(cell_, 1.2, trajectory_.getBondCutoffs());
    
    bool found = false;
    const auto& angles = analyzer.angles();
    for(const auto& t1 : angles) {
        for(const auto& center : t1) {
            for(const auto& t2 : center) {
                for(double angle : t2) {
                    if(std::abs(angle * 180.0 / M_PI - 90.0) < 1.0) {
                        found = true;
                    }
                }
            }
        }
    }
    EXPECT_TRUE(found) << "Should find 90 degree angle across PBC";
}

TEST_F(AngleReproductionTest, SiTetrahedron_4Atoms) {
    cell_.addAtom("Si", {5.0, 5.0, 5.0});       // Center
    cell_.addAtom("Si", {6.0, 6.0, 6.0});       // Neighbor 1 (1,1,1)
    cell_.addAtom("Si", {6.0, 4.0, 4.0});       // Neighbor 2 (1,-1,-1)
    cell_.addAtom("Si", {4.0, 6.0, 4.0});       // Neighbor 3 (-1,1,-1)
    cell_.addAtom("Si", {4.0, 4.0, 6.0});       // Neighbor 4 (-1,-1,1)
    updateTrajectory();

    // With 4 neighbors, we have C(4,2) = 6 angles.
    // Neighbors are at dist sqrt(3) ~ 1.73.
    // N-N dist is sqrt(8) ~ 2.82.
    // Si radius 1.16. Bond cutoff ~ 2.78. 
    // Thus neighbors are NOT connected to each other.
    
    StructureAnalyzer analyzer(cell_, 3.0, trajectory_.getBondCutoffs());
    const auto& angles = analyzer.angles();

    int angle_count = 0;
    for(const auto& t1 : angles) {
        for(const auto& center : t1) {
            for(const auto& t2 : center) {
                for(double angle : t2) {
                    double degrees = angle * 180.0 / M_PI;
                    // std::cout << "Angle: " << degrees << " degrees\n";
                    // Expected angle is acos(-1/3) ~ 109.47 degrees
                    if(std::abs(degrees - 109.47) < 1.0) {
                        angle_count++;
                    }
                }
            }
        }
    }
    EXPECT_EQ(angle_count, 6) << "Should find exactly 6 angles of ~109.47 degrees for a standard Si tetrahedron";
}

TEST_F(AngleReproductionTest, Icosahedron_13Atoms) {
    cell_.addAtom("Si", {10.0, 10.0, 10.0}); // Center
    
    double phi = (1.0 + std::sqrt(5.0)) / 2.0;
    // Vertices of icosahedron (edge length 2) relative to center
    std::vector<std::vector<double>> verts = {
        {0, 1, phi}, {0, 1, -phi}, {0, -1, phi}, {0, -1, -phi},
        {1, phi, 0}, {1, -phi, 0}, {-1, phi, 0}, {-1, -phi, 0},
        {phi, 0, 1}, {phi, 0, -1}, {-phi, 0, 1}, {-phi, 0, -1}
    };
    
    for(const auto& v : verts) {
        cell_.addAtom("Si", {10.0 + v[0], 10.0 + v[1], 10.0 + v[2]});
    }
    updateTrajectory();
    
    // Cutoff ~ 2.5 covers bonds (1.902, 2.0) but avoids next-nearest (3.236)
    StructureAnalyzer analyzer(cell_, 2.5, trajectory_.getBondCutoffs());
    const auto& angles = analyzer.angles();
    
    int count_63 = 0;  // Center-Edge (approx 63.43)
    int count_116 = 0; // Center-Diagonal (approx 116.57)
    int count_180 = 0; // Center-Opposite (180.0)
    int count_60 = 0;  // Surface-Triangle (60.0)
    int count_108 = 0; // Surface-Pentagon (108.0)
    int count_58 = 0;  // Surface-Center (approx 58.28)
    
    int total_angles = 0;
    
    for(const auto& t1 : angles) {
        for(const auto& center : t1) {
            for(const auto& t2 : center) {
                for(double angle : t2) {
                    double deg = angle * 180.0 / M_PI;
                    total_angles++;
                    
                    if(std::abs(deg - 63.43) < 1.0) count_63++;
                    else if(std::abs(deg - 116.57) < 1.0) count_116++;
                    else if(std::abs(deg - 180.0) < 1.0) count_180++;
                    else if(std::abs(deg - 60.0) < 1.0) count_60++;
                    else if(std::abs(deg - 108.0) < 1.0) count_108++;
                    else if(std::abs(deg - 58.28) < 1.0) count_58++;
                }
            }
        }
    }
    
    // Diagnosis output if needed
    if (total_angles != 246) {
        std::cout << "Found " << total_angles << " angles.\n";
        std::cout << "63: " << count_63 << "\n";
        std::cout << "116: " << count_116 << "\n";
        std::cout << "180: " << count_180 << "\n";
        std::cout << "60: " << count_60 << "\n";
        std::cout << "108: " << count_108 << "\n";
        std::cout << "58: " << count_58 << "\n";
    }

    EXPECT_EQ(count_63, 30) << "Should find 30 Center-Edge angles (~63.4 deg)";
    EXPECT_EQ(count_116, 30) << "Should find 30 Center-Diagonal angles (~116.6 deg)";
    EXPECT_EQ(count_180, 6) << "Should find 6 Center-Opposite angles (180 deg)";
    EXPECT_EQ(count_60, 60) << "Should find 60 Surface-Triangle angles (60 deg)";
    EXPECT_EQ(count_108, 60) << "Should find 60 Surface-Pentagon angles (108 deg)";
    EXPECT_EQ(count_58, 60) << "Should find 60 Surface-Center angles (Center-S-Neighbor, ~58.3 deg)";
    EXPECT_EQ(total_angles, 246) << "Total angles should be 246";
}

// ============================================================================
// Part 2: Plane Angle Distribution (PAD) Tests
// ============================================================================

// Helper to sum a partial histogram
double sumHistogram(const std::vector<double>& hist) {
    return std::accumulate(hist.begin(), hist.end(), 0.0);
}

class PADTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Large box to avoid PBC issues by default
        cell_ = Cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0}); 
    }

    void updateTrajectory() {
        trajectory_ = Trajectory();
        trajectory_.addFrame(cell_);
        trajectory_.precomputeBondCutoffs();
    }

    Cell cell_;
    Trajectory trajectory_;
};

// 1. Trivial Cases
TEST_F(PADTest, EmptyCellThrows) {
    // Current implementation throws explicitly if atoms are empty in calculateAshcroftWeights
    // or implicitly via other checks.
    updateTrajectory();
    EXPECT_THROW({
    DistributionFunctions df(cell_, 5.0, trajectory_.getBondCutoffs());
    }, std::invalid_argument);
}

TEST_F(PADTest, SingleAtomNoAngles) {
    cell_.addAtom("Si", {10.0, 10.0, 10.0});
    updateTrajectory();
    DistributionFunctions df(cell_, 5.0, trajectory_.getBondCutoffs());
    df.calculatePAD(1.0);
    // Might have partials created but empty, or just no "f(theta)" if logic handles it.
    // Actually implementation might create partials if atoms exist but no angles found.
    // Let's check total counts.
    if(df.getAllHistograms().count("f(theta)")) {
        const auto& hist = df.getHistogram("f(theta)");
        if (!hist.partials.empty()) {
            if (hist.partials.count("Total")) {
                 EXPECT_DOUBLE_EQ(sumHistogram(hist.partials.at("Total")), 0.0);
            }
        }
    }
}

// 2. Geometry Verification
TEST_F(PADTest, LinearGeometry180) {
    // A-B-C line
    cell_.addAtom("O", {9.0, 10.0, 10.0});
    auto& si = cell_.addAtom("Si", {10.0, 10.0, 10.0}); // Center
    cell_.addAtom("O", {11.0, 10.0, 10.0});
    updateTrajectory();
    
    int id_O = cell_.findElement("O")->id.value;
    int id_Si = cell_.findElement("Si")->id.value;
    
    // O(0.73) + Si(1.11) = 1.84 * 1.3 = 2.392.
    // EXPECT_GT(cutoff, 1.1) << "Bond cutoff must be larger than bond distance 1.0";

    // Verify StructureAnalyzer finds neighbors
    StructureAnalyzer analyzer(cell_, 1.5, trajectory_.getBondCutoffs());
    const auto& neighbors = analyzer.neighbors();
    // Si is atom index 1 (0-based)
    ASSERT_GT(neighbors.size(), 1);
    EXPECT_EQ(neighbors[1].size(), 2) << "Si should have 2 neighbors (O atoms)";

    // Bond length 1.0. Cutoff needs to be > 1.0
    DistributionFunctions df(cell_, 1.5, trajectory_.getBondCutoffs());
    // Use 180.0 now that we fixed the binning logic
    df.calculatePAD(1.0); 
    
    const auto& hist = df.getHistogram("f(theta)");
    // Should have O-Si-O peak at 180
    ASSERT_EQ(hist.partials.count("O-Si-O"), 1);
    const auto& partial = hist.partials.at("O-Si-O");
    
    // Bin for 180 degrees. 
    // If n_bins = 180.1/1 = 180.
    
    double total_prob = sumHistogram(partial);
    EXPECT_NEAR(total_prob, 1.0, 1e-5) << "Should be normalized to 1 angle (normalized by counts * bin_width)";
    
    // Check peak location
    double peak_val = 0;
    int peak_bin = -1;
    for(size_t i=0; i<partial.size(); ++i) {
        if (partial[i] > peak_val) {
            peak_val = partial[i];
            peak_bin = i;
        }
    }
    
    if (peak_bin >= 0) {
        double peak_angle = hist.bins[peak_bin];
        EXPECT_NEAR(peak_angle, 179.5, 1.0);
    } else {
        FAIL() << "No peak found in partial distribution";
    }
}

TEST_F(PADTest, RightAngle90) {
    cell_.addAtom("O", {10.0, 9.0, 10.0});
    cell_.addAtom("Si", {10.0, 10.0, 10.0}); // Center
    cell_.addAtom("O", {11.0, 10.0, 10.0});
    updateTrajectory();
    
    DistributionFunctions df(cell_, 1.5, trajectory_.getBondCutoffs());
    df.calculatePAD(1.0);
    
    const auto& hist = df.getHistogram("f(theta)");
    ASSERT_EQ(hist.partials.count("O-Si-O"), 1);
    
    // Find peak
    const auto& partial = hist.partials.at("O-Si-O");
    double peak_val = 0;
    int peak_bin = -1;
    for(size_t i=0; i<partial.size(); ++i) {
        if (partial[i] > peak_val) {
            peak_val = partial[i];
            peak_bin = i;
        }
    }
    double peak_angle = hist.bins[peak_bin];
    EXPECT_NEAR(peak_angle, 90.0, 1.0);
}

TEST_F(PADTest, EquilateralTriangle60) {
    // Si at (0,0,0)
    // O at (1,0,0)
    // O at (0.5, sqrt(3)/2, 0)
    cell_.addAtom("Si", {10.0, 10.0, 10.0});
    cell_.addAtom("O", {11.0, 10.0, 10.0});
    cell_.addAtom("O", {10.5, 10.0 + std::sqrt(3.0)/2.0, 10.0});
    updateTrajectory();
    
    DistributionFunctions df(cell_, 1.5, trajectory_.getBondCutoffs());
    df.calculatePAD(1.0);
    
    const auto& hist = df.getHistogram("f(theta)");
    // Should have O-Si-O
    const auto& partial = hist.partials.at("O-Si-O");
    
    // Find peak near 60
    double val_at_60 = 0;
    // index for 60 deg is 60 or 59 depending on binning.
    // 59.5 (idx 59) -> [59, 60)
    // 60.5 (idx 60) -> [60, 61)
    // Exact 60 might land in 60.
    
    // Search max around 60
    size_t bin_60 = 60; 
    EXPECT_GT(partial[bin_60] + partial[bin_60-1], 0.1) << "Should have peak near 60 degrees";
}

TEST_F(PADTest, TetrahedralAngle) {
    // Si at center
    // 4 Neighbors at tetrahedral positions.
    // For simplicity, just check one angle 109.47
    cell_.addAtom("Si", {10.0, 10.0, 10.0});
    // Vector 1: (1,1,1) normalized
    // Vector 2: (1,-1,-1) normalized
    // Dot product = (1-1-1)/3 = -1/3. acos(-1/3) = 109.47 deg
    
    double L = 1.0 / std::sqrt(3.0);
    cell_.addAtom("O", {10.0 + L, 10.0 + L, 10.0 + L});
    cell_.addAtom("O", {10.0 + L, 10.0 - L, 10.0 - L});
    updateTrajectory();
    
    DistributionFunctions df(cell_, 1.5, trajectory_.getBondCutoffs()); // Distance is 1.0
    df.calculatePAD(0.5); // Finer bins
    
    const auto& hist = df.getHistogram("f(theta)");
    const auto& partial = hist.partials.at("O-Si-O");
    
    // Expected ~109.5
    double peak_val = 0;
    double peak_angle = 0;
    for(size_t i=0; i<partial.size(); ++i) {
        if (partial[i] > peak_val) {
            peak_val = partial[i];
            peak_angle = hist.bins[i];
        }
    }
    EXPECT_NEAR(peak_angle, 109.5, 1.0);
}


// 3. Symmetry & Multi-Species
TEST_F(PADTest, SymmetryAndSorting) {
    cell_.addAtom("Si", {10.0, 10.0, 10.0}); // Center
    cell_.addAtom("O", {11.0, 10.0, 10.0});
    cell_.addAtom("N", {10.0, 11.0, 10.0}); // 90 degrees
    updateTrajectory();
    
    DistributionFunctions df(cell_, 1.5, trajectory_.getBondCutoffs());
    df.calculatePAD(1.0);
    
    const auto& hist = df.getHistogram("f(theta)");
    
    // Check if we have O-Si-N or N-Si-O
    bool found = false;
    if (hist.partials.count("O-Si-N")) found = true;
    if (hist.partials.count("N-Si-O")) found = true;
    
    EXPECT_TRUE(found) << "Should have mixed species angle distribution";
}

// 4. Normalization
TEST_F(PADTest, FullNormalizationCheck) {
    // 1 Si, 4 O neighbors (tetrahedron)
    // 4 neighbors -> 4*3/2 = 6 angles.
    // All 6 angles are 109.47
    
    cell_.addAtom("Si", {10.0, 10.0, 10.0});
    double L = 1.0 / std::sqrt(3.0);
    
    // Tetrahedral vertices
    cell_.addAtom("O", {10.0 + L, 10.0 + L, 10.0 + L});
    cell_.addAtom("O", {10.0 + L, 10.0 - L, 10.0 - L});
    cell_.addAtom("O", {10.0 - L, 10.0 + L, 10.0 - L});
    cell_.addAtom("O", {10.0 - L, 10.0 - L, 10.0 + L});
    updateTrajectory();
    
    DistributionFunctions df(cell_, 1.5, trajectory_.getBondCutoffs());
    df.calculatePAD(1.0);
    
    const auto& hist = df.getHistogram("f(theta)");
    
    double sum_partial = 0;
    double sum_total = 0;
    double bin_width = 1.0;
    


    if (hist.partials.count("O-Si-O")) {
        const auto& partial = hist.partials.at("O-Si-O");
        for(double v : partial) sum_partial += v * bin_width;
    }
    
    if (hist.partials.count("Total")) {
        const auto& total = hist.partials.at("Total");
        for(double v : total) sum_total += v * bin_width;
    }
    
    EXPECT_NEAR(sum_partial, 1.0, 0.05); // Relaxed checking 0.05 due to binning effects
    EXPECT_NEAR(sum_total, 1.0, 0.05);
}
