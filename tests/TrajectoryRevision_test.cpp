#include <gtest/gtest.h>
#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#include <chrono>

#include "../include/Cell.hpp"
#include "../include/Trajectory.hpp"
#include "../include/FileIO.hpp"

// Test fixture for Trajectory tests
class TrajectoryRevisionTest : public ::testing::Test {
protected:
  void SetUp() override {
    // Create a temporary directory for test files
    std::filesystem::create_directory("test_data");
  }

  void TearDown() override {
    // Clean up test files
    std::filesystem::remove_all("test_data");
  }
};

TEST_F(TrajectoryRevisionTest, RemoveDuplicatedFrames) {
  // Create a few cells
  Cell c1;
  c1.addAtom("H", {0.0, 0.0, 0.0});
  c1.addAtom("O", {1.0, 0.0, 0.0});

  Cell c2; // Identical to c1
  c2.addAtom("H", {0.0, 0.0, 0.0});
  c2.addAtom("O", {1.0, 0.0, 0.0});

  Cell c3; // Different positions
  c3.addAtom("H", {0.1, 0.0, 0.0});
  c3.addAtom("O", {1.1, 0.0, 0.0});

  Cell c4; // Identical to c3
  c4.addAtom("H", {0.1, 0.0, 0.0});
  c4.addAtom("O", {1.1, 0.0, 0.0});

  Cell c5; // Different atoms (simulation: just partial difference or different position)
  // To satisfy Trajectory validation, it must have same atom count and symbols in order
  c5.addAtom("H", {0.2, 0.0, 0.0});
  c5.addAtom("O", {1.2, 0.0, 0.0});

  std::vector<Cell> frames = {c1, c2, c3, c4, c5};
  Trajectory traj(frames, 1.0);

  EXPECT_EQ(traj.getFrameCount(), 5);

  traj.removeDuplicatedFrames();

  EXPECT_EQ(traj.getFrameCount(), 3);
  // Expected frames: c1, c3, c5
  
  const auto& new_frames = traj.getFrames();
  EXPECT_EQ(new_frames[0].atoms()[0].position().x(), 0.0);
  EXPECT_EQ(new_frames[1].atoms()[0].position().x(), 0.1);
  EXPECT_EQ(new_frames[2].atoms()[0].position().x(), 0.2);
}

TEST_F(TrajectoryRevisionTest, ParseEnergyFromArc) {
  // Create a mock .arc file with energy
  std::string filename = "test_data/energy_test.arc";
  std::ofstream out(filename);
  out << "!BIOSYM archive 3\n";
  out << "PBC=ON\n";
  out << "                                      -12345.6789\n";
  out << "!DATE     Feb 13 17:21:41 2025\n";
  out << "PBC   10.0 10.0 10.0 90.0 90.0 90.0\n";
  out << "He     0.0 0.0 0.0 XXXX 1 xx He 0.0\n";
  out << "end\n";
  out << "end\n";
  out.close();

  auto traj = FileIO::readTrajectory(filename, FileIO::FileType::Arc);
  const auto& frames = traj.getFrames();
  ASSERT_EQ(frames.size(), 1);
  
  EXPECT_DOUBLE_EQ(frames[0].getEnergy(), -12345.6789);
}

TEST_F(TrajectoryRevisionTest, ParseMultipleFramesWithEnergy) {
    std::string filename = "test_data/multi_energy.arc";
    std::ofstream out(filename);
    
    // Frame 1
    out << "!BIOSYM archive 3\n";
    out << "PBC=ON\n";
    out << " -100.0\n"; 
    out << "!DATE\n";
    out << "PBC 10.0 10.0 10.0 90.0 90.0 90.0\n";
    out << "He 0.0 0.0 0.0 XXXX 1 xx He 0.0\n";
    out << "end\n";
    out << "end\n";

    // Frame 2
    out << "!BIOSYM archive 3\n";
    out << "PBC=ON\n";
    out << " -200.5\n"; 
    out << "!DATE\n";
    out << "PBC 10.0 10.0 10.0 90.0 90.0 90.0\n";
    out << "He 1.0 1.0 1.0 XXXX 1 xx He 0.0\n";
    out << "end\n";
    out << "end\n";
    out.close();

    auto traj2 = FileIO::readTrajectory(filename, FileIO::FileType::Arc);
    const auto& frames = traj2.getFrames();
    ASSERT_EQ(frames.size(), 2);
    EXPECT_DOUBLE_EQ(frames[0].getEnergy(), -100.0);
    EXPECT_DOUBLE_EQ(frames[1].getEnergy(), -200.5);
}

TEST_F(TrajectoryRevisionTest, RemoveDuplicatesFromRealFile) {
    // This file is known to have duplicates
    // Using absolute path to ensure checking the specific file requested
    std::string filename = "/home/isurwars/Projects/Correlation/examples/l-Bi/l-Bi.arc";
    
    // Check if file exists to avoid failing if example is missing (though user said it's there)
    if (!std::filesystem::exists(filename)) {
        GTEST_SKIP() << "Test file not found: " << filename;
    }

    auto traj = FileIO::readTrajectory(filename, FileIO::FileType::Arc);
    size_t initial_count = traj.getFrameCount();
    std::cout << "Initial frame count: " << initial_count << std::endl;

    auto time_start = std::chrono::high_resolution_clock::now();
    traj.removeDuplicatedFrames();
    auto time_end = std::chrono::high_resolution_clock::now();
    
    size_t final_count = traj.getFrameCount();
    std::cout << "Final frame count: " << final_count << std::endl;
    std::cout << "Removed: " << (initial_count - final_count) << " frames" << std::endl;
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count() << "ms" << std::endl;

    EXPECT_LT(final_count, initial_count);
    EXPECT_NEAR(initial_count - final_count, 100, 20); // Expect around 100 duplicates
}
