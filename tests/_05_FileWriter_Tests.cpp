// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <algorithm> // For std::max_element
#include <cstdio>    // For std::remove
#include <fstream>
#include <gtest/gtest.h>
#include <highfive/highfive.hpp>
#include <iterator> // For std::distance

#include "../include/Cell.hpp"
#include "../include/DistributionFunctions.hpp"
#include "../include/FileReader.hpp"
#include "../include/FileWriter.hpp"
#include "../include/Trajectory.hpp"

// Test fixture for FileWriter integration tests.
class _05_FileWriter_Tests : public ::testing::Test {
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
    std::remove("test_si__G.csv");
    std::remove("test_si_PAD.csv");
    std::remove("test_si_g_smoothed.csv");
    std::remove("test_si_J_smoothed.csv");
    std::remove("test_si_G_smoothed.csv");
    std::remove("test_si__G_smoothed.csv");
    std::remove("test_si_PAD_smoothed.csv");
    std::remove("test_si.h5");
    std::remove("test_vacf.h5");
    std::remove("test_vacf_new.h5");
    std::remove("test_vacf_vdos.h5");
    std::remove("test_vacf_vdos_VDOS.csv");
  }

  // Helper to check if a file exists and is not empty.
  bool fileExistsAndIsNotEmpty(const std::string &name) {
    if (std::ifstream f(name); f.good()) {
      return f.peek() != std::ifstream::traits_type::eof();
    }
    return false;
  }
};

TEST_F(_05_FileWriter_Tests, CalculatesAndWritesSiliconDistributions) {
  // Arrange
  FileReader::FileType type = FileReader::determineFileType("si_crystal.car");
  Cell si_cell = FileReader::readStructure("si_crystal.car", type);
  Trajectory trajectory;
  trajectory.addFrame(si_cell);
  trajectory.precomputeBondCutoffs();

  DistributionFunctions df(si_cell, 20.0, trajectory.getBondCutoffsSQ());

  // Act
  const double rdf_bin = 0.05;
  const double pad_bin = 1.0;
  df.calculateRDF(20.0, rdf_bin);
  df.calculatePAD(pad_bin);
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

TEST_F(_05_FileWriter_Tests, WritesHDF5File) {
  // Arrange
  FileReader::FileType type = FileReader::determineFileType("si_crystal.car");
  Cell si_cell = FileReader::readStructure("si_crystal.car", type);
  Trajectory trajectory;
  trajectory.addFrame(si_cell);
  trajectory.precomputeBondCutoffs();

  DistributionFunctions df(
      si_cell, 5.0,
      trajectory.getBondCutoffsSQ()); // Use smaller r_max for faster test

  df.calculateRDF(5.0, 0.1);
  df.calculatePAD(2.0);

  FileWriter writer(df);
  writer.writeHDF("test_si.h5");

  // Assert
  ASSERT_TRUE(fileExistsAndIsNotEmpty("test_si.h5"));

  // Verify content using HighFive
  HighFive::File file("test_si.h5", HighFive::File::ReadOnly);
  EXPECT_TRUE(file.exist("g_r"));
  EXPECT_TRUE(file.exist("f_theta"));

  HighFive::Group g_group = file.getGroup("g_r");
  // The old "data" dataset should not exist anymore
  EXPECT_FALSE(g_group.exist("data"));

  // Check group description
  EXPECT_TRUE(g_group.hasAttribute("description"));
  std::string description;
  g_group.getAttribute("description").read(description);
  EXPECT_EQ(description, "Radial Distribution Function");

  // Check for expected datasets (prefixed names)
  // For Si crystal, we expect "r (Å)" (bin) and "Si-Si" (raw)
  // Cleaned up names:
  // "r (Å)" -> "r__Å_" -> "00_r__Å_"
  // "Si-Si" -> "01_Si-Si"

  // Note: sanitize logic in FileWriter replaces '(', ')', '/', ' ' with '_'
  // "r (Å)" -> "r__Å_"
  std::string bin_ds_name = "00_r__Å_";
  std::string data_ds_name = "01_Si-Si";

  EXPECT_TRUE(g_group.exist(bin_ds_name));
  EXPECT_TRUE(g_group.exist(data_ds_name));

  // Verify bin dataset attributes
  HighFive::DataSet bin_ds = g_group.getDataSet(bin_ds_name);
  EXPECT_TRUE(bin_ds.hasAttribute("Units"));
  std::string bin_units;
  bin_ds.getAttribute("Units").read(bin_units);
  EXPECT_EQ(bin_units, "Å");

  EXPECT_TRUE(bin_ds.hasAttribute("Long Name"));
  std::string bin_label;
  bin_ds.getAttribute("Long Name").read(bin_label);
  EXPECT_EQ(bin_label, "r (Å)");

  EXPECT_TRUE(bin_ds.hasAttribute("Comments"));
  std::string bin_comment;
  bin_ds.getAttribute("Comments").read(bin_comment);
  EXPECT_EQ(bin_comment, "r (Å)");

  // Verify data dataset attributes
  HighFive::DataSet data_ds = g_group.getDataSet(data_ds_name);
  EXPECT_TRUE(data_ds.hasAttribute("Units"));
  std::string data_units;
  data_ds.getAttribute("Units").read(data_units);
  EXPECT_EQ(data_units, "Å^-1");

  EXPECT_TRUE(data_ds.hasAttribute("Long Name"));
  std::string data_label;
  data_ds.getAttribute("Long Name").read(data_label);
  EXPECT_EQ(data_label, "Si-Si");

  EXPECT_TRUE(data_ds.hasAttribute("Comments"));
  std::string data_comment;
  data_ds.getAttribute("Comments").read(data_comment);
  EXPECT_EQ(data_comment, "Si-Si");
}

TEST_F(_05_FileWriter_Tests, WritesVACFMetadata) {
  // Arrange
  FileReader::FileType type = FileReader::determineFileType("si_crystal.car");
  Cell frame1 = FileReader::readStructure("si_crystal.car", type);

  // Create frame2 with same lattice
  auto lv = frame1.latticeVectors();
  Cell frame2(lv[0], lv[1], lv[2]);

  // Modify frame2 slightly to create velocity
  // Use alternating signs to ensure zero net momentum (approx) for CoM removal
  // test
  int sign = 1;
  for (const auto &atom : frame1.atoms()) {
    auto pos = atom.position();
    pos.x() += 0.1 * sign;
    sign *= -1;
    frame2.addAtom(atom.element().symbol, pos);
  }

  Trajectory trajectory;
  trajectory.addFrame(frame1);
  trajectory.addFrame(frame2);
  trajectory.setTimeStep(1.0);        // 1 fs
  trajectory.precomputeBondCutoffs(); // Required for DistributionFunctions
  trajectory.calculateVelocities();

  DistributionFunctions df(frame1, 5.0, trajectory.getBondCutoffsSQ());
  df.calculateVACF(trajectory, 1);

  FileWriter writer(df);
  std::string filename = "test_vacf_new.h5";
  writer.writeHDF(filename);

  // Assert
  ASSERT_TRUE(fileExistsAndIsNotEmpty(filename));
  {
    HighFive::File file(filename, HighFive::File::ReadOnly);

    // Check VACF
    EXPECT_TRUE(file.exist("VACF"));
    HighFive::Group vacf_group = file.getGroup("VACF");

    // Description
    EXPECT_TRUE(vacf_group.hasAttribute("description"));
    std::string description;
    vacf_group.getAttribute("description").read(description);
    EXPECT_EQ(description, "Velocity Autocorrelation Function");

    // Check data dataset no longer exists
    EXPECT_FALSE(vacf_group.exist("data"));

    // Check new datasets
    // "Time" -> "00_Time"
    // "VACF" -> "Total" -> "01_Total"
    std::string time_ds_name = "00_Time";
    std::string vacf_ds_name = "01_Total";

    EXPECT_TRUE(vacf_group.exist(time_ds_name));
    EXPECT_TRUE(vacf_group.exist(vacf_ds_name));

    HighFive::DataSet time_ds = vacf_group.getDataSet(time_ds_name);
    EXPECT_TRUE(time_ds.hasAttribute("Units"));
    std::string bin_units;
    time_ds.getAttribute("Units").read(bin_units);
    EXPECT_EQ(bin_units, "fs");

    HighFive::DataSet vacf_ds = vacf_group.getDataSet(vacf_ds_name);
    EXPECT_TRUE(vacf_ds.hasAttribute("Units"));
    std::string data_units;
    vacf_ds.getAttribute("Units").read(data_units);
    EXPECT_EQ(data_units, "Å^2/fs^2");

    // Check Normalized VACF
    EXPECT_TRUE(file.exist("Normalized_VACF"));
    HighFive::Group norm_vacf_group = file.getGroup("Normalized_VACF");

    // Description
    EXPECT_TRUE(norm_vacf_group.hasAttribute("description"));
    std::string norm_desc;
    norm_vacf_group.getAttribute("description").read(norm_desc);
    EXPECT_EQ(norm_desc, "Normalized Velocity Autocorrelation Function");

    // Check new datasets
    // "Time" -> "00_Time"
    // "Normalized VACF" -> "Total" -> "01_Total"
    std::string norm_vacf_name = "01_Total";

    EXPECT_TRUE(norm_vacf_group.exist(norm_vacf_name));
    HighFive::DataSet norm_vacf_ds = norm_vacf_group.getDataSet(norm_vacf_name);

    // Data units
    EXPECT_TRUE(norm_vacf_ds.hasAttribute("Units"));
    std::string norm_data_units;
    norm_vacf_ds.getAttribute("Units").read(norm_data_units);
    EXPECT_EQ(norm_data_units, "normalized");
  } // Close file

  // Calculate and Check VDOS
  df.calculateVDOS();

  std::string vdos_filename = "test_vacf_vdos.h5";
  writer.writeHDF(vdos_filename);
  writer.writeAllCSVs("test_vacf_vdos", true);
  EXPECT_TRUE(fileExistsAndIsNotEmpty("test_vacf_vdos_VDOS.csv"));

  // Re-open to check VDOS
  HighFive::File file_vdos(vdos_filename, HighFive::File::ReadOnly);
  EXPECT_TRUE(file_vdos.exist("VDOS"));
  HighFive::Group vdos_group = file_vdos.getGroup("VDOS");

  EXPECT_TRUE(vdos_group.hasAttribute("description"));
  std::string vdos_desc;
  vdos_group.getAttribute("description").read(vdos_desc);
  EXPECT_EQ(vdos_desc, "Vibrational Density of States");

  std::string vdos_freq_name = "00_Frequency__THz_";
  std::string vdos_val_name = "03_Total";

  EXPECT_TRUE(vdos_group.exist(vdos_freq_name));
  EXPECT_TRUE(vdos_group.exist(vdos_val_name));

  // Check new units in HDF5
  std::string vdos_cm_name = "01_Frequency__cm-1_";
  std::string vdos_mev_name = "02_Frequency__meV_";

  EXPECT_TRUE(vdos_group.exist(vdos_cm_name));
  EXPECT_TRUE(vdos_group.exist(vdos_mev_name));

  HighFive::DataSet vdos_cm_ds = vdos_group.getDataSet(vdos_cm_name);
  EXPECT_TRUE(vdos_cm_ds.hasAttribute("Units"));
  std::string cm_units;
  vdos_cm_ds.getAttribute("Units").read(cm_units);
  EXPECT_EQ(cm_units, "cm^-1");

  HighFive::DataSet vdos_mev_ds = vdos_group.getDataSet(vdos_mev_name);
  EXPECT_TRUE(vdos_mev_ds.hasAttribute("Units"));
  std::string mev_units;
  vdos_mev_ds.getAttribute("Units").read(mev_units);
  EXPECT_EQ(mev_units, "meV");
}
