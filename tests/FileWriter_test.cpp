// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <algorithm> // For std::max_element
#include <cstdio>    // For std::remove
#include <fstream>
#include <gtest/gtest.h>
#include <iterator> // For std::distance

#include "../include/Cell.hpp"
#include "../include/DistributionFunctions.hpp"
#include "../include/FileIO.hpp"
#include "../include/FileWriter.hpp"
#include <highfive/highfive.hpp>

// Test fixture for FileWriter integration tests.
class FileWriterTest : public ::testing::Test {
protected:
  void SetUp() override {
    // Create a temporary CAR file for an 8-atom Silicon crystal.
    std::ofstream car_file("si_crystal.car");
    ASSERT_TRUE(car_file.is_open());
    car_file << "!BIOSYM archive 3\n";
    car_file << "PBC=ON\n";
    car_file << "PBC 5.431 5.431 5.431 90.0 90.0 90.0\n";
    car_file << "Si1   0.0000  0.0000  0.0000 XXXX 1    xx    Si   0.000\n";
    car_file << "Si2   0.0000  2.7155  2.7155 XXXX 1    xx    Si   0.000\n";
    car_file << "Si3   2.7155  0.0000  2.7155 XXXX 1    xx    Si   0.000\n";
    car_file << "Si4   2.7155  2.7155  0.0000 XXXX 1    xx    Si   0.000\n";
    car_file << "Si5   1.3578  1.3578  1.3578 XXXX 1    xx    Si   0.000\n";
    car_file << "Si6   1.3578  4.0733  4.0733 XXXX 1    xx    Si   0.000\n";
    car_file << "Si7   4.0733  1.3578  4.0733 XXXX 1    xx    Si   0.000\n";
    car_file << "Si8   4.0733  4.0733  1.3578 XXXX 1    xx    Si   0.000\n";
    car_file << "end\n";
    car_file << "end\n";
    car_file.close();
  }

  void TearDown() override {
    // Clean up all generated files.
    std::remove("si_crystal.car");
    std::remove("test_si_g.csv");
    std::remove("test_si_J.csv");
    std::remove("test_si_G.csv");
    std::remove("test_si_PAD.csv");
    std::remove("test_si_g_smoothed.csv");
    std::remove("test_si_J_smoothed.csv");
    std::remove("test_si_G_smoothed.csv");
    std::remove("test_si_PAD_smoothed.csv");
    std::remove("test_si.h5");
  }

  // Helper to check if a file exists and is not empty.
  bool fileExistsAndIsNotEmpty(const std::string &name) {
    if (std::ifstream f(name); f.good()) {
      return f.peek() != std::ifstream::traits_type::eof();
    }
    return false;
  }
};

TEST_F(FileWriterTest, CalculatesAndWritesSiliconDistributions) {
  // Arrange
  FileIO::FileType type = FileIO::determineFileType("si_crystal.car");
  Cell si_cell = FileIO::readStructure("si_crystal.car", type);
  DistributionFunctions df(si_cell, 20.0);

  // Act
  const double rdf_bin = 0.05;
  const double pad_bin = 1.0;
  df.calculateRDF(20.0, rdf_bin);
  df.calculatePAD(180.0, pad_bin);
  df.smoothAll(0.1);

  FileWriter writer(df);
  writer.writeAllCSVs("test_si", true);

  // Assert: Part 1 - Validate content of the calculated g(r) histogram.

  const auto &rdf_hist = df.getHistogram("g(r)");
  const auto &bins = rdf_hist.bins;
  const auto &si_si_rdf = rdf_hist.partials.at("Si-Si");

  // Helper to find the index of the maximum value in a given range.
  auto find_peak_idx = [&](size_t start, size_t end) {
    auto it =
        std::max_element(si_si_rdf.begin() + start, si_si_rdf.begin() + end);
    return std::distance(si_si_rdf.begin(), it);
  };

  // Find peaks in expected regions for crystalline silicon.
  // 1st neighbor shell: ~2.35 Å. Search from 2.0 to 3.0 Å.
  size_t first_peak_idx = find_peak_idx(2.0 / rdf_bin, 3.0 / rdf_bin);
  // 2nd neighbor shell: ~3.84 Å. Search from 3.5 to 4.2 Å.
  size_t second_peak_idx = find_peak_idx(3.5 / rdf_bin, 4.2 / rdf_bin);

  EXPECT_NEAR(bins[first_peak_idx], 2.35, rdf_bin * 2);
  EXPECT_NEAR(bins[second_peak_idx], 3.84, rdf_bin * 2);

  // Assert: Part 2 - Validate content of the calculated f(theta) histogram.

  const auto &pad_hist = df.getHistogram("f(theta)");
  const auto &pad_bins = pad_hist.bins;
  const auto &si_si_si_pad = pad_hist.partials.at("Si-Si-Si");

  auto max_it = std::max_element(si_si_si_pad.begin(), si_si_si_pad.end());
  size_t peak_index = std::distance(si_si_si_pad.begin(), max_it);

  // The dominant peak angle should be the tetrahedral angle in silicon.
  EXPECT_NEAR(pad_bins[peak_index], 109.5, pad_bin * 2.0);

  // Assert: Part 3 - Check that all expected files were created
  EXPECT_TRUE(fileExistsAndIsNotEmpty("test_si_g.csv"));
  EXPECT_TRUE(fileExistsAndIsNotEmpty("test_si_J.csv"));
  EXPECT_TRUE(fileExistsAndIsNotEmpty("test_si__G.csv"));
  EXPECT_TRUE(fileExistsAndIsNotEmpty("test_si_PAD.csv"));
}

TEST_F(FileWriterTest, WritesHDF5File) {
  // Arrange
  FileIO::FileType type = FileIO::determineFileType("si_crystal.car");
  Cell si_cell = FileIO::readStructure("si_crystal.car", type);
  DistributionFunctions df(si_cell, 5.0); // Use smaller r_max for faster test

  df.calculateRDF(5.0, 0.1);
  df.calculatePAD(180.0, 2.0);

  FileWriter writer(df);
  writer.writeHDF("test_si.h5");

  // Assert
  ASSERT_TRUE(fileExistsAndIsNotEmpty("test_si.h5"));

  // Verify content using HighFive
  HighFive::File file("test_si.h5", HighFive::File::ReadOnly);
  EXPECT_TRUE(file.exist("g(r)"));
  EXPECT_TRUE(file.exist("f(theta)"));

  HighFive::Group g_group = file.getGroup("g(r)");
  EXPECT_TRUE(g_group.exist("bins"));
  EXPECT_TRUE(g_group.exist("raw"));
  EXPECT_TRUE(g_group.getGroup("raw").exist("Si-Si"));

  // Check attribute
  EXPECT_TRUE(g_group.hasAttribute("bin_label"));
  std::string label;
  g_group.getAttribute("bin_label").read(label);
  EXPECT_EQ(label, "r (Å)");
}
