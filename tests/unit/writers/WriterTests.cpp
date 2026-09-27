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
#include <miniz.h>
#ifdef CORRELATION_USE_HDF5
#include <highfive/highfive.hpp>
#endif
#include <iterator> // For std::distance
#include <vector>

namespace {

// Helper to find the index of the maximum value in a given range.
size_t findPeakIdx(const std::vector<correlation::real_t> &vec, size_t start, size_t end) {
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

  void TearDown() override {}

  static void cleanPrefix(const std::string &prefix) {
    std::error_code error_code;
    for (const auto &entry : std::filesystem::directory_iterator(".", error_code)) {
      if (entry.is_regular_file()) {
        const auto filename = entry.path().filename().string();
        if (filename.rfind(prefix, 0) == 0) {
          std::filesystem::remove(entry.path(), error_code);
        }
      }
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
  const real_t dad_bin = 1.0;
  dists.calculateRDF({
      .r_max = 20.0,
      .r_bin_width = rdf_bin,
  });
  dists.calculatePAD(pad_bin);
  dists.calculateDAD(dad_bin);
  dists.smoothAll(0.1);

  const std::string prefix = "test_si_calc";
  cleanPrefix(prefix);

  correlation::writers::FileWriter const writer(dists);
  writer.write(prefix, true, false, false, true);

  // Assert: Part 1 - Validate content of the calculated g_r histogram.

  const auto &rdf_hist = dists.getHistogram("g_r");
  const auto &bins = rdf_hist.bins;
  const auto &si_si_rdf = rdf_hist.partials.at("Si-Si");

  // Find peaks in expected regions for crystalline silicon.
  // 1st neighbor shell: ~2.35 Å. Search from 2.0 to 3.0 Å.
  size_t const first_peak_idx = findPeakIdx(si_si_rdf, static_cast<size_t>(2.0 / rdf_bin),
                                            static_cast<size_t>(3.0 / rdf_bin));
  // 2nd neighbor shell: ~3.84 Å. Search from 3.5 to 4.2 Å.
  size_t const second_peak_idx = findPeakIdx(si_si_rdf, static_cast<size_t>(3.5 / rdf_bin),
                                             static_cast<size_t>(4.2 / rdf_bin));

  EXPECT_NEAR(bins[first_peak_idx], 2.35, rdf_bin * 2);
  EXPECT_NEAR(bins[second_peak_idx], 3.84, rdf_bin * 2);

  // Assert: Part 2 - Validate content of the calculated BAD histogram.

  const auto &pad_hist = dists.getHistogram("PAD");
  const auto &pad_bins = pad_hist.bins;
  const auto &si_si_si_pad = pad_hist.partials.at("Si-Si-Si");

  auto max_it = std::max_element(si_si_si_pad.begin(), si_si_si_pad.end());
  size_t const peak_index = std::distance(si_si_si_pad.begin(), max_it);

  // The dominant peak angle should be the tetrahedral angle in silicon.
  EXPECT_NEAR(pad_bins[peak_index], 109.5, pad_bin * 2.0);

  // Assert: Part 3 - Check that all expected files were created
  EXPECT_TRUE(fileExistsAndIsNotEmpty(prefix + "_g.csv"));
  EXPECT_TRUE(fileExistsAndIsNotEmpty(prefix + "_J.csv"));
  EXPECT_TRUE(fileExistsAndIsNotEmpty(prefix + "_G_reduced.csv"));
  EXPECT_TRUE(fileExistsAndIsNotEmpty(prefix + "_PAD.csv"));
  EXPECT_TRUE(fileExistsAndIsNotEmpty(prefix + "_PAD_raw.csv"));
  EXPECT_TRUE(fileExistsAndIsNotEmpty(prefix + "_DAD.csv"));
  EXPECT_TRUE(fileExistsAndIsNotEmpty(prefix + "_DAD_raw.csv"));

  cleanPrefix(prefix);
}

TEST_F(FileWriterTests, PADCsvContainsRawAndNormalizedAndSmoothedColumns) {
  // Arrange: single-element silicon crystal
  std::string const path = getDataDir() + "si_crystal.car";
  correlation::readers::FileType const type = correlation::readers::determineFileType(path);
  correlation::core::Cell const si_cell = correlation::readers::readStructure(path, type);
  correlation::core::Trajectory trajectory;
  trajectory.addFrame(si_cell);
  trajectory.precomputeBondCutoffs();

  correlation::analysis::DistributionFunctions dists(si_cell, 5.0, trajectory.getBondCutoffsSQ());
  dists.calculatePAD(2.0);
  dists.smoothAll(0.1);

  const std::string prefix = "test_si_pad_cols";
  cleanPrefix(prefix);

  correlation::writers::FileWriter const writer(dists);
  writer.write(prefix, true, false, false, true);

  // Act: read the PAD CSV header (line 1)
  std::ifstream pad_file(prefix + "_PAD.csv");
  ASSERT_TRUE(pad_file.good());
  std::string header_line;
  std::getline(pad_file, header_line);

  // Assert: header must contain raw companion columns, normalized, and smoothed
  EXPECT_NE(header_line.find("Si-Si-Si_raw"), std::string::npos)
      << "PAD CSV header missing raw companion column 'Si-Si-Si_raw'";
  EXPECT_NE(header_line.find("Total_raw"), std::string::npos)
      << "PAD CSV header missing raw companion column 'Total_raw'";
  EXPECT_NE(header_line.find("Si-Si-Si_smoothed"), std::string::npos)
      << "PAD CSV header missing smoothed column 'Si-Si-Si_smoothed'";
  // Verify the raw columns appear before the normalized columns
  auto raw_pos = header_line.find("Si-Si-Si_raw");
  auto norm_pos = header_line.find(",Si-Si-Si,");
  EXPECT_LT(raw_pos, norm_pos) << "Raw companion columns should appear before normalized columns";

  cleanPrefix(prefix);
}

TEST_F(FileWriterTests, PADRawCsvContainsRawCountsAndSmoothed) {
  // Arrange
  std::string const path = getDataDir() + "si_crystal.car";
  correlation::readers::FileType const type = correlation::readers::determineFileType(path);
  correlation::core::Cell const si_cell = correlation::readers::readStructure(path, type);
  correlation::core::Trajectory trajectory;
  trajectory.addFrame(si_cell);
  trajectory.precomputeBondCutoffs();

  correlation::analysis::DistributionFunctions dists(si_cell, 5.0, trajectory.getBondCutoffsSQ());
  dists.calculatePAD(2.0);
  dists.smoothAll(0.1);

  const std::string prefix = "test_si_pad_raw";
  cleanPrefix(prefix);

  correlation::writers::FileWriter const writer(dists);
  writer.write(prefix, true, false, false, true);

  // Assert: PAD_raw CSV exists and has smoothed columns
  std::ifstream raw_file(prefix + "_PAD_raw.csv");
  ASSERT_TRUE(raw_file.good());
  std::string header_line;
  std::getline(raw_file, header_line);

  // Raw file should have raw counts and smoothed raw counts
  EXPECT_NE(header_line.find("Si-Si-Si"), std::string::npos)
      << "PAD_raw CSV header missing raw partial 'Si-Si-Si'";
  EXPECT_NE(header_line.find("Si-Si-Si_smoothed"), std::string::npos)
      << "PAD_raw CSV header missing smoothed column 'Si-Si-Si_smoothed'";
  EXPECT_NE(header_line.find("Total"), std::string::npos)
      << "PAD_raw CSV header missing 'Total' partial";

  // Units row should show 'counts'
  std::string units_line;
  std::getline(raw_file, units_line);
  EXPECT_NE(units_line.find("counts"), std::string::npos)
      << "PAD_raw CSV units row should contain 'counts'";

  cleanPrefix(prefix);
}

TEST_F(FileWriterTests, DADCsvAndRawCsvAreCreated) {
  // Arrange
  std::string const path = getDataDir() + "si_crystal.car";
  correlation::readers::FileType const type = correlation::readers::determineFileType(path);
  correlation::core::Cell const si_cell = correlation::readers::readStructure(path, type);
  correlation::core::Trajectory trajectory;
  trajectory.addFrame(si_cell);
  trajectory.precomputeBondCutoffs();

  correlation::analysis::DistributionFunctions dists(si_cell, 5.0, trajectory.getBondCutoffsSQ());
  dists.calculateDAD(2.0);
  dists.smoothAll(0.1);

  const std::string prefix = "test_si_dad";
  cleanPrefix(prefix);

  correlation::writers::FileWriter const writer(dists);
  writer.write(prefix, true, false, false, true);

  // Both normalized and raw DAD CSVs must exist
  EXPECT_TRUE(fileExistsAndIsNotEmpty(prefix + "_DAD.csv"));
  EXPECT_TRUE(fileExistsAndIsNotEmpty(prefix + "_DAD_raw.csv"));

  // Verify DAD.csv header contains raw companion columns
  std::ifstream dad_file(prefix + "_DAD.csv");
  ASSERT_TRUE(dad_file.good());
  std::string header_line;
  std::getline(dad_file, header_line);
  EXPECT_NE(header_line.find("_raw"), std::string::npos)
      << "DAD CSV header should contain raw companion columns with '_raw' suffix";

  cleanPrefix(prefix);
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

  dists.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.1,
  });
  dists.calculatePAD(2.0);

  const std::string prefix = "test_si_h5";
  cleanPrefix(prefix);

  correlation::writers::FileWriter writer(dists);
  writer.write(prefix, false, true, false, false);

  // Assert
  ASSERT_TRUE(fileExistsAndIsNotEmpty(prefix + ".h5"));

  // Verify content using HighFive
  HighFive::File file(prefix + ".h5", HighFive::File::ReadOnly);
  EXPECT_TRUE(file.exist("g_r"));
  EXPECT_TRUE(file.exist("PAD"));

  HighFive::Group g_group = file.getGroup("g_r");
  // The old "data" dataset should not exist anymore
  EXPECT_FALSE(g_group.exist("data"));

  // Check group description
  EXPECT_TRUE(g_group.hasAttribute("description"));
  std::string description;
  g_group.getAttribute("description").read(description);
  EXPECT_EQ(description, "Radial Distribution Function");

  std::string bin_ds_name = "r";
  std::string data_ds_name = "Si-Si";

  EXPECT_TRUE(g_group.exist(bin_ds_name));
  EXPECT_TRUE(g_group.exist(data_ds_name));

  // Verify bin dataset attributes
  HighFive::DataSet bin_ds = g_group.getDataSet(bin_ds_name);
  EXPECT_TRUE(bin_ds.hasAttribute("units"));
  std::string bin_units;
  bin_ds.getAttribute("units").read(bin_units);
  EXPECT_EQ(bin_units, "Å");

  EXPECT_TRUE(bin_ds.hasAttribute("label"));
  std::string bin_label;
  bin_ds.getAttribute("label").read(bin_label);
  EXPECT_EQ(bin_label, "r");

  EXPECT_TRUE(bin_ds.hasAttribute("description"));
  std::string bin_comment;
  bin_ds.getAttribute("description").read(bin_comment);
  EXPECT_EQ(bin_comment, "Radial Distribution Function");

  // Verify data dataset attributes
  HighFive::DataSet data_ds = g_group.getDataSet(data_ds_name);
  EXPECT_TRUE(data_ds.hasAttribute("units"));
  std::string data_units;
  data_ds.getAttribute("units").read(data_units);
  EXPECT_EQ(data_units, "Å⁻¹");

  EXPECT_TRUE(data_ds.hasAttribute("label"));
  std::string data_label;
  data_ds.getAttribute("label").read(data_label);
  EXPECT_EQ(data_label, "Si-Si");
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
    pos.x() += static_cast<real_t>(0.1 * sign);
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
    std::string time_ds_name = "t";
    std::string vacf_ds_name = "Total";

    EXPECT_TRUE(vacf_group.exist(time_ds_name));
    EXPECT_TRUE(vacf_group.exist(vacf_ds_name));

    HighFive::DataSet time_ds = vacf_group.getDataSet(time_ds_name);
    EXPECT_TRUE(time_ds.hasAttribute("units"));
    std::string bin_units;
    time_ds.getAttribute("units").read(bin_units);
    EXPECT_EQ(bin_units, "fs");

    HighFive::DataSet vacf_ds = vacf_group.getDataSet(vacf_ds_name);
    EXPECT_TRUE(vacf_ds.hasAttribute("units"));
    std::string data_units;
    vacf_ds.getAttribute("units").read(data_units);
    EXPECT_EQ(data_units, "Å² fs⁻²");

    // Check Normalized VACF
    EXPECT_TRUE(file.exist("Normalized_VACF"));
    HighFive::Group norm_vacf_group = file.getGroup("Normalized_VACF");

    // Description
    EXPECT_TRUE(norm_vacf_group.hasAttribute("description"));
    std::string norm_desc;
    norm_vacf_group.getAttribute("description").read(norm_desc);
    EXPECT_EQ(norm_desc, "Normalized Velocity Autocorrelation Function");

    std::string norm_vacf_name = "Total";

    EXPECT_TRUE(norm_vacf_group.exist(norm_vacf_name));
    HighFive::DataSet norm_vacf_ds = norm_vacf_group.getDataSet(norm_vacf_name);

    // Data units
    EXPECT_TRUE(norm_vacf_ds.hasAttribute("units"));
    std::string norm_data_units;
    norm_vacf_ds.getAttribute("units").read(norm_data_units);
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

  std::string vdos_freq_name = "ν";
  std::string vdos_val_name = "Total";

  EXPECT_TRUE(vdos_group.exist(vdos_freq_name));
  EXPECT_TRUE(vdos_group.exist(vdos_val_name));
  HighFive::DataSet vdos_val_ds = vdos_group.getDataSet(vdos_val_name);
  EXPECT_TRUE(vdos_val_ds.hasAttribute("units"));
  std::string vdos_units;
  vdos_val_ds.getAttribute("units").read(vdos_units);
  EXPECT_EQ(vdos_units, "arbitrary units");
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

  dists.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.1,
  });
  dists.calculatePAD(2.0);

  const std::string prefix = "test_si_parquet";
  cleanPrefix(prefix);

  correlation::writers::FileWriter writer(dists);

  // Act
  // write(base_path, use_csv, use_hdf5, use_parquet, smoothing)
  writer.write(prefix, false, false, true, false);

  // Assert
  EXPECT_TRUE(fileExistsAndIsNotEmpty(prefix + "_g.parquet"));
  EXPECT_TRUE(fileExistsAndIsNotEmpty(prefix + "_J.parquet"));
  EXPECT_TRUE(fileExistsAndIsNotEmpty(prefix + "_G_reduced.parquet"));
  EXPECT_TRUE(fileExistsAndIsNotEmpty(prefix + "_PAD.parquet"));
  EXPECT_FALSE(fileExistsAndIsNotEmpty(prefix + "_S.parquet"));

  cleanPrefix(prefix);
}
#endif

TEST_F(FileWriterTests, WritesZipBundle) {
  // Arrange
  std::string const path = getDataDir() + "si_crystal.car";
  correlation::readers::FileType type = correlation::readers::determineFileType(path);
  correlation::core::Cell si_cell = correlation::readers::readStructure(path, type);
  correlation::core::Trajectory trajectory;
  trajectory.addFrame(si_cell);
  trajectory.precomputeBondCutoffs();

  correlation::analysis::DistributionFunctions dists(si_cell, 5.0, trajectory.getBondCutoffsSQ());

  dists.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.1,
  });
  dists.calculatePAD(2.0);

  const std::string zip_path = "test_si_bundle.zip";
  std::error_code error_code;
  std::filesystem::remove(zip_path, error_code);

  correlation::writers::FileWriter writer(dists);

  // Act: write bundle with CSV and SVG enabled
  auto res = writer.writeBundle(zip_path, /*use_csv=*/true, /*use_hdf5=*/false,
                                /*use_parquet=*/false, /*include_svg=*/true,
                                /*smoothing=*/false);

  // Assert
  ASSERT_TRUE(res.has_value()) << res.error();
  ASSERT_TRUE(std::filesystem::exists(zip_path));
  EXPECT_GT(std::filesystem::file_size(zip_path), 0U);

  // Verify contents using miniz
  mz_zip_archive zip;
  mz_zip_zero_struct(&zip);
  ASSERT_TRUE(mz_zip_reader_init_file(&zip, zip_path.c_str(), 0));

  EXPECT_GE(mz_zip_reader_locate_file(&zip, "summary.txt", nullptr, 0), 0);
  EXPECT_GE(mz_zip_reader_locate_file(&zip, "plots/g.svg", nullptr, 0), 0);
  EXPECT_GE(mz_zip_reader_locate_file(&zip, "plots/PAD.svg", nullptr, 0), 0);
  EXPECT_GE(mz_zip_reader_locate_file(&zip, "csv/g.csv", nullptr, 0), 0);
  EXPECT_GE(mz_zip_reader_locate_file(&zip, "csv/PAD.csv", nullptr, 0), 0);

  mz_zip_reader_end(&zip);
  std::filesystem::remove(zip_path, error_code);
}

TEST_F(FileWriterTests, WritesCategorizedFolder) {
  // Arrange
  std::string const path = getDataDir() + "si_crystal.car";
  correlation::readers::FileType type = correlation::readers::determineFileType(path);
  correlation::core::Cell si_cell = correlation::readers::readStructure(path, type);
  correlation::core::Trajectory trajectory;
  trajectory.addFrame(si_cell);
  trajectory.precomputeBondCutoffs();

  correlation::analysis::DistributionFunctions dists(si_cell, 5.0, trajectory.getBondCutoffsSQ());
  dists.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.1,
  });
  dists.calculatePAD(2.0);

  const std::string folder_path = "test_si_folder_export";
  std::error_code error_code;
  std::filesystem::remove_all(folder_path, error_code);

  correlation::writers::FileWriter writer(dists);

  // Act: write categorized folder with CSV and SVG enabled
  auto res = writer.writeFolder(folder_path, /*use_csv=*/true, /*use_hdf5=*/false,
                                /*use_parquet=*/false, /*include_svg=*/true,
                                /*smoothing=*/false);

  // Assert
  ASSERT_TRUE(res.has_value()) << res.error();
  ASSERT_TRUE(std::filesystem::is_directory(folder_path));

  EXPECT_TRUE(std::filesystem::exists(std::filesystem::path(folder_path) / "summary.txt"));
  EXPECT_TRUE(std::filesystem::exists(std::filesystem::path(folder_path) / "csv" / "g.csv"));
  EXPECT_TRUE(std::filesystem::exists(std::filesystem::path(folder_path) / "csv" / "PAD.csv"));
  EXPECT_TRUE(std::filesystem::exists(std::filesystem::path(folder_path) / "plots" / "g.svg"));
  EXPECT_TRUE(std::filesystem::exists(std::filesystem::path(folder_path) / "plots" / "PAD.svg"));

  std::filesystem::remove_all(folder_path, error_code);
}

TEST_F(FileWriterTests, WritesCategorizedFolderWithAlgorithmFiltering) {
  // Arrange
  std::string const path = getDataDir() + "si_crystal.car";
  correlation::readers::FileType type = correlation::readers::determineFileType(path);
  correlation::core::Cell si_cell = correlation::readers::readStructure(path, type);
  correlation::core::Trajectory trajectory;
  trajectory.addFrame(si_cell);
  trajectory.precomputeBondCutoffs();

  correlation::analysis::DistributionFunctions dists(si_cell, 5.0, trajectory.getBondCutoffsSQ());
  dists.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.1,
  });
  dists.calculatePAD(2.0);

  const std::string folder_path = "test_si_folder_filtered";
  std::error_code error_code;
  std::filesystem::remove_all(folder_path, error_code);

  correlation::writers::FileWriter writer(dists);

  // Filter only PAD (exclude RDF / g(r))
  std::vector<std::string> selected_algos = {"PAD"};

  // Act
  auto res = writer.writeFolder(folder_path, /*use_csv=*/true, /*use_hdf5=*/false,
                                /*use_parquet=*/false, /*include_svg=*/true,
                                /*smoothing=*/false, selected_algos);

  // Assert
  ASSERT_TRUE(res.has_value()) << res.error();
  ASSERT_TRUE(std::filesystem::is_directory(folder_path));

  EXPECT_TRUE(std::filesystem::exists(std::filesystem::path(folder_path) / "csv" / "PAD.csv"));
  EXPECT_TRUE(std::filesystem::exists(std::filesystem::path(folder_path) / "plots" / "PAD.svg"));
  EXPECT_FALSE(std::filesystem::exists(std::filesystem::path(folder_path) / "csv" / "g.csv"));
  EXPECT_FALSE(std::filesystem::exists(std::filesystem::path(folder_path) / "plots" / "g.svg"));

  std::filesystem::remove_all(folder_path, error_code);
}
