/**
 * @file ArrowWriterTests.cpp
 * @brief Unit tests for ArrowWriter Parquet export capability.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "writers/ArrowWriter.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include <gtest/gtest.h>
#include <filesystem>

namespace correlation::writers::testing {

TEST(ArrowWriterTests, MetadataAndExtensionInfo) {
  ArrowWriter writer;
  EXPECT_EQ(writer.getName(), "Parquet");
  auto exts = writer.getExtensions();
  ASSERT_EQ(exts.size(), 1U);
  EXPECT_EQ(exts[0], ".parquet");
}

TEST(ArrowWriterTests, WriteAllParquetCreatesOutputFile) {
  std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "correlation_arrow_tests";
  std::filesystem::create_directories(temp_dir);
  std::string base_output = (temp_dir / "test_sample").string();

  correlation::core::Cell cell(std::array<real_t, 6>{10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {2.0, 0.0, 0.0});

  correlation::core::Trajectory traj;
  traj.addFrame(cell);
  traj.precomputeBondCutoffs();

  correlation::analysis::DistributionFunctions dists(cell, 5.0, traj.getBondCutoffsSQ());
  dists.calculateRDF({.r_max = 5.0, .r_bin_width = 0.5});

  EXPECT_NO_THROW(ArrowWriter::writeAllParquet(base_output, dists, false));

  std::string expected_parquet = base_output + "_g.parquet";
  EXPECT_TRUE(std::filesystem::exists(expected_parquet));

  // Clean up
  std::filesystem::remove_all(temp_dir);
}

} // namespace correlation::writers::testing
