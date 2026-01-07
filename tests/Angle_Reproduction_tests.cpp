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
    
    StructureAnalyzer analyzer(cell_, 1.2);
    
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

    // With 4 neighbors, we have C(4,2) = 6 angles.
    // Neighbors are at dist sqrt(3) ~ 1.73.
    // N-N dist is sqrt(8) ~ 2.82.
    // Si radius 1.16. Bond cutoff ~ 2.78. 
    // Thus neighbors are NOT connected to each other.
    
    StructureAnalyzer analyzer(cell_, 3.0);
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
