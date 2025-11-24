#include <gtest/gtest.h>
#include "../include/Cell.hpp"
#include "../include/StructureAnalyzer.hpp"
#include "../include/DistributionFunctions.hpp"
#include <cmath>

class AngleReproductionTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Simple cubic cell
        cell_ = Cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
    }
    Cell cell_;
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

    // Bond cutoff for Si-O is likely around 1.6 * 1.2 = 1.92 or similar.
    // Distance is 1.0.
    
    // If we set analyzer cutoff to 0.9, we should find nothing.
    {
        StructureAnalyzer analyzer(cell_, 0.9);
        const auto& angles = analyzer.angles();
        // angles[type][center][type]
        // We expect empty because cutoff < distance
        // But wait, StructureAnalyzer might throw if cutoff is too small? No.
        
        // Check if any angles are found
        bool found = false;
        for(const auto& t1 : angles) {
            for(const auto& center : t1) {
                for(const auto& t2 : center) {
                    if(!t2.empty()) found = true;
                }
            }
        }
        EXPECT_TRUE(found) << "Should find angles even if cutoff is smaller than bond distance (auto-correction)";

    }

    // If we set analyzer cutoff to 1.1, we should find the angle.
    {
        StructureAnalyzer analyzer(cell_, 1.1);
        const auto& angles = analyzer.angles();
        
        // Find Si (center)
        // We need to know element IDs.
        // Si is usually heavier, maybe index 1? 
        // Let's just iterate and find 90 degrees.
        
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
    // Atom A at 0.1
    // Atom B (center) at 9.9 (wrapped) -> effectively -0.1
    // Atom C at 9.9 (wrapped) but different axis?
    
    // Let's put Center at 0.5, 0.5, 0.5
    // Neighbor 1 at 9.6, 0.5, 0.5 (distance 0.9 across boundary)
    // Neighbor 2 at 0.5, 9.6, 0.5 (distance 0.9 across boundary)
    // Angle should be 90 degrees.
    
    cell_.addAtom("Si", {0.5, 0.5, 0.5});
    cell_.addAtom("O", {9.6, 0.5, 0.5});
    cell_.addAtom("O", {0.5, 9.6, 0.5});
    
    StructureAnalyzer analyzer(cell_, 1.5); // Cutoff > 0.9
    
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
