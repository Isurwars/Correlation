#include <gtest/gtest.h>
#include <filesystem>
#include <vector>
#include <iostream>

#include "../include/DynamicsAnalyzer.hpp"
#include "../include/FileIO.hpp"
#include "../include/Trajectory.hpp"

// Assuming the test is run from the build directory or project root
// We need to locate the l-Bi.arc file
const std::string EXAMPLE_FILE = "examples/l-Bi/l-Bi.arc";

TEST(DynamicsAnalyzerTest, CalculatesVACFFromExampletraj) {
    // 1. Locate the file
    std::filesystem::path file_path = std::filesystem::absolute(EXAMPLE_FILE);
    
    // If running from build/, we might need to go up one level or adjust path logic
    // Usually standard CMake test layout puts us in build/tests or similar.
    // Let's try to find it relative to current path or project root.
    if (!std::filesystem::exists(file_path)) {
        // Try going up levels if not found
        file_path = std::filesystem::path("../" + EXAMPLE_FILE);
        if (!std::filesystem::exists(file_path)) {
            file_path = std::filesystem::path("../../" + EXAMPLE_FILE);
        }
    }
    
    ASSERT_TRUE(std::filesystem::exists(file_path)) << "Could not find example file: " << EXAMPLE_FILE;

    // 2. Read Trajectory
    Trajectory traj = FileIO::readTrajectory(file_path.string(), FileIO::FileType::Arc);
    ASSERT_GT(traj.getFrameCount(), 0) << "Trajectory should not be empty";
    
    // 3. Calculate Velocities
    traj.calculateVelocities();
    
    const auto& velocities = traj.getVelocities();
    ASSERT_EQ(velocities.size(), traj.getFrameCount());
    ASSERT_EQ(velocities[0].size(), traj.getFrames()[0].atomCount());
    
    // Check if we have some non-zero velocities (it's liquid Bi, particles move)
    double max_v_sq = 0.0;
    for (const auto& v : velocities[10]) { // Check some intermediate frame
        max_v_sq = std::max(max_v_sq, linalg::dot(v, v));
    }
    EXPECT_GT(max_v_sq, 0.0) << "Particles should be moving";

    // 4. Calculate VACF
    int max_lag = 50; // Calculate for 50 frames lag
    std::vector<double> vacf = DynamicsAnalyzer::calculateVACF(traj, max_lag);
    
    ASSERT_EQ(vacf.size(), max_lag + 1);
    
    // C(0) should be positive (autocorrelation at t=0 is <v^2>)
    EXPECT_GT(vacf[0], 0.0);
    
    // 5. Calculate Normalized VACF
    std::vector<double> norm_vacf = DynamicsAnalyzer::calculateNormalizedVACF(traj, max_lag);
    
    ASSERT_EQ(norm_vacf.size(), max_lag + 1);
    EXPECT_NEAR(norm_vacf[0], 1.0, 1e-5) << "Normalized VACF should start at 1.0";
    
    // Check decay (statistical nature means it won't strictly decrease, but should generally drop)
    // For a liquid, it eventually decorrelates.
    // Just print the first few values for sanity check in logs
    std::cout << "VACF[0-5]: ";
    for(int i=0; i<=5 && i < (int)norm_vacf.size(); ++i) {
        std::cout << norm_vacf[i] << " ";
    }
    std::cout << std::endl;
}
