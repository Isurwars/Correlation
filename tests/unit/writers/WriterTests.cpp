// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DistributionFunctions.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "math/Precision.hpp"
#include "readers/FileReader.hpp"
#include "writers/FileWriter.hpp"

#include <algorithm> // For std::max_element
#include <cstdio>    // For std::remove
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#ifdef CORRELATION_USE_HDF5
#include <highfive/highfive.hpp>
#endif
#include <iterator> // For std::distance
#include <vector>

namespace {

// Helper to find the index of the maximum value in a given range.
size_t find_peak_idx(const std::vector<correlation::real_t> &vec, size_t start, size_t end) {
  auto const iterator = std::max_element(vec.begin() + static_cast<std::ptrdiff_t>(start),
                                         vec.begin() + static_cast<std::ptrdiff_t>(end));
  return static_cast<size_t>(std::distance(vec.begin(), iterator));
}

// Test fixture for FileWriter integration tests.
class FileWriterTests : public ::testing::Test {
protected:
  [[nodiscard]] std::string const &getDataDir() const { return data_dir_; }

private:
  std::string data_dir_;

protected:
  void SetUp() override {
    std::vector<std::string> const candidates = {
        "../../tests/data/",
        "../tests/data/",
        "tests/data/",
        "data/",
    };
    std::string base_dir = "../../tests/data/";
    for (const auto &dir : candidates) {
      if (std::filesystem::exists(dir + "car/si_crystal.car")) {
        base_dir = dir;
        break;
      }
    }
    data_dir_ = base_dir + "car/";
  }

  void TearDown() override {
    // Clean up all generated files.
    std::vector<std::string> const files_to_remove = {"test_si_g.csv",
                                                      "test_si_J.csv",
                                                      "test_si_G.csv",
                                                      "test_si_PAD.csv",
                                                      "test_si_DAD.csv",
                                                      "test_si_RD.csv",
                                                      "test_si_S.csv",
                                                      "test_si_XRD.csv",
                                                      "test_si_CN.csv",
                                                      "test_si_g_smoothed.csv",
                                                      "test_si_J_smoothed.csv",
                                                      "test_si_G_smoothed.csv",
                                                      "test_si_PAD_smoothed.csv",
                                                      "test_si.h5",
                                                      "test_vacf.h5",
                                                      "test_vacf_new.h5",
                                                      "test_vacf_new_VACF.csv",
                                                      "test_vacf_new_VACF_norm.csv",
                                                      "test_vacf_new_VDOS.csv",
                                                      "test_vacf_vdos.h5",
                                                      "test_vacf_vdos_VACF.csv",
                                                      "test_vacf_vdos_VACF_norm.csv",
                                                      "test_vacf_vdos_VDOS.csv",
                                                      "test_si_g.parquet",
                                                      "test_si_J.parquet",
                                                      "test_si_G.parquet",
                                                      "test_si_PAD.parquet",
                                                      "test_si_summary.txt",
                                                      "test_vacf_new_summary.txt",
                                                      "test_vacf_vdos_summary.txt"};

    for (const auto &file : files_to_remove) {
      std::error_code error_code;
      std::filesystem::remove(file, error_code);
    }
  }

  // Helper to check if a file exists and is not empty.
  static bool fileExistsAndIsNotEmpty(const std::string &name) {
    if (std::ifstream file(name); file.good()) {
      return file.peek() != std::ifstream::traits_type::eof();
    }
    return false;
  }
};
} // namespace

TEST_F(FileWriterTests, CalculatesAndWritesSiliconDistributions) {
  // Arrange
  std::string const path = getDataDir() + "si_crystal.car";
  correlation::readers::FileType const type = correlation::readers::determineFileType(path);
  correlation::core::Cell const si_cell = correlation::readers::readStructure(path, type);
  correlation::core::Trajectory trajectory;
  trajectory.addFrame(si_cell);
  trajectory.precomputeBondCutoffs();

  correlation::analysis::DistributionFunctions dists(si_cell, 20.0, trajectory.getBondCutoffsSQ());

  // Act
  const real_t rdf_bin = 0.05;
  const real_t pad_bin = 1.0;
  dists.calculateRDF({.r_max = 20.0, .r_bin_width = rdf_bin});
  dists.calculatePAD(pad_bin);
  dists.smoothAll(0.1);

  correlation::writers::FileWriter const writer(dists);
  writer.write("test_si", true, false, false, true);

  // Assert: Part 1 - Validate content of the calculated g_r histogram.

  const auto &rdf_hist = dists.getHistogram("g_r");
  const auto &bins = rdf_hist.bins;
  const auto &si_si_rdf = rdf_hist.partials.at("Si-Si");

  // Find peaks in expected regions for crystalline silicon.
  // 1st neighbor shell: ~2.35 Å. Search from 2.0 to 3.0 Å.
  size_t const first_peak_idx =
      find_peak_idx(si_si_rdf, static_cast<size_t>(2.0 / rdf_bin), static_cast<size_t>(3.0 / rdf_bin));
  // 2nd neighbor shell: ~3.84 Å. Search from 3.5 to 4.2 Å.
  size_t const second_peak_idx =
      find_peak_idx(si_si_rdf, static_cast<size_t>(3.5 / rdf_bin), static_cast<size_t>(4.2 / rdf_bin));

  EXPECT_NEAR(bins[first_peak_idx], 2.35, rdf_bin * 2);
  EXPECT_NEAR(bins[second_peak_idx], 3.84, rdf_bin * 2);

  // Assert: Part 2 - Validate content of the calculated BAD histogram.

  const auto &pad_hist = dists.getHistogram("BAD");
  const auto &pad_bins = pad_hist.bins;
  const auto &si_si_si_pad = pad_hist.partials.at("Si-Si-Si");

  auto max_it = std::max_element(si_si_si_pad.begin(), si_si_si_pad.end());
  size_t const peak_index = std::distance(si_si_si_pad.begin(), max_it);

  // The dominant peak angle should be the tetrahedral angle in silicon.
  EXPECT_NEAR(pad_bins[peak_index], 109.5, pad_bin * 2.0);

  // Assert: Part 3 - Check that all expected files were created
  EXPECT_TRUE(fileExistsAndIsNotEmpty("test_si_g.csv"));
  EXPECT_TRUE(fileExistsAndIsNotEmpty("test_si_J.csv"));
  EXPECT_TRUE(fileExistsAndIsNotEmpty("test_si_G.csv"));
  EXPECT_TRUE(fileExistsAndIsNotEmpty("test_si_PAD.csv"));
}

#ifdef CORRELATION_USE_HDF5
TEST_F(FileWriterTests, WritesHDF5File) {
  // Arrange
  std::string const path = getDataDir() + "si_crystal.car";
  correlation::readers::FileType type = correlation::readers::determineFileType(path);
  correlation::core::Cell si_cell = correlation::readers::readStructure(path, type);
  correlation::core::Trajectory trajectory;
  trajectory.addFrame(si_cell);
  trajectory.precomputeBondCutoffs();

  correlation::analysis::DistributionFunctions dists(
      si_cell, 5.0,
      trajectory.getBondCutoffsSQ()); // Use smaller r_max for faster test

  dists.calculateRDF(5.0, 0.1);
  dists.calculatePAD(2.0);

  correlation::writers::FileWriter writer(dists);
  writer.write("test_si", false, true, false, false);

  // Assert
  ASSERT_TRUE(fileExistsAndIsNotEmpty("test_si.h5"));

  // Verify content using HighFive
  HighFive::File file("test_si.h5", HighFive::File::ReadOnly);
  EXPECT_TRUE(file.exist("g_r"));
  EXPECT_TRUE(file.exist("BAD"));

  HighFive::Group g_group = file.getGroup("g_r");
  // The old "data" dataset should not exist anymore
  EXPECT_FALSE(g_group.exist("data"));

  // Check group description
  EXPECT_TRUE(g_group.hasAttribute("description"));
  std::string description;
  g_group.getAttribute("description").read(description);
  EXPECT_EQ(description, "Pair Distribution Function");

  // Note: sanitize logic in correlation::writers::FileWriter replaces '(', ')',
  // '/', ' ' with '_'
  std::string bin_ds_name = "00_r";
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
  EXPECT_EQ(bin_label, "r");

  EXPECT_TRUE(bin_ds.hasAttribute("Comments"));
  std::string bin_comment;
  bin_ds.getAttribute("Comments").read(bin_comment);
  EXPECT_EQ(bin_comment, "Pair Distribution Function");

  // Verify data dataset attributes
  HighFive::DataSet data_ds = g_group.getDataSet(data_ds_name);
  EXPECT_TRUE(data_ds.hasAttribute("Units"));
  std::string data_units;
  data_ds.getAttribute("Units").read(data_units);
  EXPECT_EQ(data_units, "Å⁻¹");

  EXPECT_TRUE(data_ds.hasAttribute("Long Name"));
  std::string data_label;
  data_ds.getAttribute("Long Name").read(data_label);
  EXPECT_EQ(data_label, "Si-Si");

  EXPECT_TRUE(data_ds.hasAttribute("Comments"));
  std::string data_comment;
  data_ds.getAttribute("Comments").read(data_comment);
  EXPECT_EQ(data_comment, "Si-Si");
}

TEST_F(FileWriterTests, WritesVACFMetadata) {
  // Arrange
  std::string const path = getDataDir() + "si_crystal.car";
  correlation::readers::FileType type = correlation::readers::determineFileType(path);
  correlation::core::Cell frame1 = correlation::readers::readStructure(path, type);

  // Create frame2 with same lattice
  auto lattice_vectors = frame1.latticeVectors();
  correlation::core::Cell frame2(lattice_vectors[0], lattice_vectors[1], lattice_vectors[2]);

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

  correlation::core::Trajectory trajectory;
  trajectory.addFrame(frame1);
  trajectory.addFrame(frame2);
  trajectory.setTimeStep(1.0);        // 1 fs
  trajectory.precomputeBondCutoffs(); // Required for DistributionFunctions
  trajectory.calculateVelocities();

  correlation::analysis::DistributionFunctions dists(frame1, 5.0, trajectory.getBondCutoffsSQ());
  dists.calculateVACF(trajectory, correlation::analysis::MaxFrames{1});

  correlation::writers::FileWriter writer(dists);
  std::string base_filename = "test_vacf_new";
  writer.write(base_filename, false, true, false, false);
  std::string filename = base_filename + ".h5";

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
    // "t (fs)" -> "00_t__fs_"
    // "VACF" -> "Total" -> "01_Total"
    std::string time_ds_name = "00_t";
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
    EXPECT_EQ(data_units, "Å² fs⁻²");

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
  dists.calculateVDOS();

  std::string vdos_base = "test_vacf_vdos";
  writer.write(vdos_base, true, true, false, true);
  std::string vdos_filename = vdos_base + ".h5";
  EXPECT_TRUE(fileExistsAndIsNotEmpty("test_vacf_vdos_VDOS.csv"));

  // Re-open to check VDOS
  HighFive::File file_vdos(vdos_filename, HighFive::File::ReadOnly);
  EXPECT_TRUE(file_vdos.exist("VDOS"));
  HighFive::Group vdos_group = file_vdos.getGroup("VDOS");

  EXPECT_TRUE(vdos_group.hasAttribute("description"));
  std::string vdos_desc;
  vdos_group.getAttribute("description").read(vdos_desc);
  EXPECT_EQ(vdos_desc, "Vibrational Density of States");

  std::string vdos_freq_name = "00_ν";
  std::string vdos_val_name = "03_Total";

  EXPECT_TRUE(vdos_group.exist(vdos_freq_name));
  EXPECT_TRUE(vdos_group.exist(vdos_val_name));

  // Check new units in HDF5
  std::string vdos_cm_name = "01_Frequency_cm_1";
  std::string vdos_mev_name = "02_Frequency_meV";

  EXPECT_TRUE(vdos_group.exist(vdos_cm_name));
  EXPECT_TRUE(vdos_group.exist(vdos_mev_name));

  HighFive::DataSet vdos_cm_ds = vdos_group.getDataSet(vdos_cm_name);
  EXPECT_TRUE(vdos_cm_ds.hasAttribute("Units"));
  std::string cm_units;
  vdos_cm_ds.getAttribute("Units").read(cm_units);
  EXPECT_EQ(cm_units, "arbitrary units");

  HighFive::DataSet vdos_mev_ds = vdos_group.getDataSet(vdos_mev_name);
  EXPECT_TRUE(vdos_mev_ds.hasAttribute("Units"));
  std::string mev_units;
  vdos_mev_ds.getAttribute("Units").read(mev_units);
  EXPECT_EQ(mev_units, "arbitrary units");
}
#endif

#ifdef CORRELATION_USE_ARROW
TEST_F(FileWriterTests, WritesParquetFiles) {
  // Arrange
  std::string const path = getDataDir() + "si_crystal.car";
  correlation::readers::FileType type = correlation::readers::determineFileType(path);
  correlation::core::Cell si_cell = correlation::readers::readStructure(path, type);
  correlation::core::Trajectory trajectory;
  trajectory.addFrame(si_cell);
  trajectory.precomputeBondCutoffs();

  correlation::analysis::DistributionFunctions dists(si_cell, 5.0, trajectory.getBondCutoffsSQ());

  dists.calculateRDF(5.0, 0.1);
  dists.calculatePAD(2.0);

  correlation::writers::FileWriter writer(dists);

  // Act
  // write(base_path, use_csv, use_hdf5, use_parquet, smoothing)
  writer.write("test_si", false, false, true, false);

  // Assert
  EXPECT_TRUE(fileExistsAndIsNotEmpty("test_si_g.parquet"));
  EXPECT_TRUE(fileExistsAndIsNotEmpty("test_si_J.parquet"));
  EXPECT_TRUE(fileExistsAndIsNotEmpty("test_si_G_reduced.parquet"));
  EXPECT_TRUE(fileExistsAndIsNotEmpty("test_si_PAD.parquet"));
  EXPECT_FALSE(fileExistsAndIsNotEmpty("test_si_S.parquet"));
}
#endif
