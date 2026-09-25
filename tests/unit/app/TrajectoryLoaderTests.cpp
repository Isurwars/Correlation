// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only

#include "app/services/TrajectoryLoader.hpp"

#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace {

std::string getTestDataDir() {
  std::vector<std::string> const candidates = {
      "../../tests/data/",
      "../tests/data/",
      "tests/data/",
      "data/",
  };
  for (const auto &dir : candidates) {
    if (std::filesystem::exists(dir + "xyz/clean.xyz")) {
      return dir;
    }
  }
  return "../../tests/data/";
}

TEST(TrajectoryLoaderTests, InitialStateIsEmpty) {
  correlation::app::TrajectoryLoader loader;

  EXPECT_EQ(loader.trajectory(), nullptr);
  EXPECT_EQ(loader.cell(), nullptr);
  EXPECT_EQ(loader.getFrameCount(), 0);
  EXPECT_EQ(loader.getTotalAtomCount(), 0);
  EXPECT_EQ(loader.getRemovedFrameCount(), 0);
  EXPECT_DOUBLE_EQ(loader.getTimeStep(), 1.0);
  EXPECT_TRUE(loader.getAtomCounts().empty());
}

TEST(TrajectoryLoaderTests, LoadsValidStructureFile) {
  correlation::app::TrajectoryLoader loader;
  std::string const filepath = getTestDataDir() + "xyz/clean.xyz";

  auto res = loader.loadFile(filepath);
  ASSERT_TRUE(res.has_value());
  EXPECT_NE(loader.trajectory(), nullptr);
  EXPECT_NE(loader.cell(), nullptr);
  EXPECT_GT(loader.getFrameCount(), 0);
  EXPECT_GT(loader.getTotalAtomCount(), 0);

  auto counts = loader.getAtomCounts();
  EXPECT_FALSE(counts.empty());
  EXPECT_GT(loader.getRecommendedTimeStep(), 0.0);
}

TEST(TrajectoryLoaderTests, RejectsNonExistentFile) {
  correlation::app::TrajectoryLoader loader;
  auto res = loader.loadFile("/non/existent/path/file.xyz");
  EXPECT_FALSE(res.has_value());
  EXPECT_EQ(loader.trajectory(), nullptr);
}

TEST(TrajectoryLoaderTests, ClearResetsState) {
  correlation::app::TrajectoryLoader loader;
  std::string const filepath = getTestDataDir() + "xyz/clean.xyz";

  auto res = loader.loadFile(filepath);
  ASSERT_TRUE(res.has_value());
  EXPECT_NE(loader.cell(), nullptr);

  loader.clear();
  EXPECT_EQ(loader.cell(), nullptr);
  EXPECT_EQ(loader.getFrameCount(), 0);
}

} // namespace
