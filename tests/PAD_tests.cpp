#include <gtest/gtest.h>
#include <cmath>
#include <numeric>
#include "../include/Cell.hpp"
#include "../include/StructureAnalyzer.hpp"
#include "../include/DistributionFunctions.hpp"

// Helper to sum a partial histogram
double sumHistogram(const std::vector<double>& hist) {
    return std::accumulate(hist.begin(), hist.end(), 0.0);
}

class PADTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Large box to avoid PBC issues by default
        cell_ = Cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0}); 
        cell_.setBondFactor(1.3); // Standard bond factor
    }
    Cell cell_;
};

// 1. Trivial Cases
TEST_F(PADTest, EmptyCellThrows) {
    // Current implementation throws explicitly if atoms are empty in calculateAshcroftWeights
    // or implicitly via other checks.
    EXPECT_THROW({
        DistributionFunctions df(cell_, 5.0, 1.3);
    }, std::invalid_argument);
}

TEST_F(PADTest, SingleAtomNoAngles) {
    cell_.addAtom("Si", {10.0, 10.0, 10.0});
    DistributionFunctions df(cell_, 5.0, 1.3);
    df.calculatePAD(180.0, 1.0);
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
    
    // Debug: Check bond cutoff is sufficient
    // We need element IDs.
    // Assuming 0=O, 1=Si (based on insertion order? No, map order).
    // Let's rely on cell internals or just trust if getBondCutoff works.
    int id_O = cell_.findElement("O")->id.value;
    int id_Si = cell_.findElement("Si")->id.value;
    double cutoff = cell_.getBondCutoff(id_O, id_Si);
    // O(0.73) + Si(1.11) = 1.84 * 1.3 = 2.392.
    EXPECT_GT(cutoff, 1.1) << "Bond cutoff must be larger than bond distance 1.0";

    // Verify StructureAnalyzer finds neighbors
    StructureAnalyzer analyzer(cell_, 1.5);
    const auto& neighbors = analyzer.neighbors();
    // Si is atom index 1 (0-based)
    ASSERT_GT(neighbors.size(), 1);
    EXPECT_EQ(neighbors[1].size(), 2) << "Si should have 2 neighbors (O atoms)";

    // Bond length 1.0. Cutoff needs to be > 1.0
    DistributionFunctions df(cell_, 1.5, 1.3);
    // Use 180.0 now that we fixed the binning logic
    df.calculatePAD(180.0, 1.0); 
    
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
    
    DistributionFunctions df(cell_, 1.5, 1.3);
    df.calculatePAD(180.0, 1.0);
    
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
    
    DistributionFunctions df(cell_, 1.5, 1.3);
    df.calculatePAD(180.0, 1.0);
    
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
    
    DistributionFunctions df(cell_, 1.5, 1.3); // Distance is 1.0
    df.calculatePAD(180.0, 0.5); // Finer bins
    
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
    // A-B-C should partial B-A-C (or A-B-C sorted) or similar? 
    // Implementation uses "Type1-Center-Type2" key.
    // If we have O-Si-N, is it stored as N-Si-O or O-Si-N?
    // It should be canonical.
    
    cell_.addAtom("Si", {10.0, 10.0, 10.0}); // Center
    cell_.addAtom("O", {11.0, 10.0, 10.0});
    // Add Nitrogen
    cell_.addAtom("N", {10.0, 11.0, 10.0}); // 90 degrees
    
    DistributionFunctions df(cell_, 1.5, 1.3);
    df.calculatePAD(180.0, 1.0);
    
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
    
    DistributionFunctions df(cell_, 1.5, 1.1);
    df.calculatePAD(180.0, 1.0);
    
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
