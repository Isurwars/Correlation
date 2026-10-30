/**
 * @file NaturalSortTests.cpp
 * @brief Unit tests for NaturalSort and glob expansion utilities.
 */

#include "readers/NaturalSort.hpp"

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <string>
#include <vector>

using namespace correlation::readers;

TEST(NaturalSortTests, NaturalLessNumericComparison) {
  EXPECT_TRUE(naturalLess("dump.2", "dump.10"));
  EXPECT_FALSE(naturalLess("dump.10", "dump.2"));
  EXPECT_TRUE(naturalLess("dump.02", "dump.10"));
  EXPECT_TRUE(naturalLess("run1_part2.dump", "run1_part10.dump"));
  EXPECT_FALSE(naturalLess("file.dump", "file.dump"));
}

TEST(NaturalSortTests, NaturalSortOrdering) {
  std::vector<std::string> paths = {
      "dump.100.dump", "dump.2.dump", "dump.20.dump", "dump.1.dump", "dump.10.dump",
  };
  naturalSort(paths);

  const std::vector<std::string> expected = {
      "dump.1.dump", "dump.2.dump", "dump.10.dump", "dump.20.dump", "dump.100.dump",
  };
  EXPECT_EQ(paths, expected);
}

TEST(NaturalSortTests, WildcardMatchPatterns) {
  EXPECT_TRUE(detail::wildcardMatch("*.dump", "file.dump"));
  EXPECT_FALSE(detail::wildcardMatch("*.dump", "file.txt"));
  EXPECT_TRUE(detail::wildcardMatch("chunk_*.dump", "chunk_01.dump"));
  EXPECT_TRUE(detail::wildcardMatch("chunk_?.dump", "chunk_1.dump"));
  EXPECT_FALSE(detail::wildcardMatch("chunk_?.dump", "chunk_12.dump"));
  EXPECT_TRUE(detail::wildcardMatch("*", "anything"));
  EXPECT_TRUE(detail::wildcardMatch("exact.name", "exact.name"));
  EXPECT_FALSE(detail::wildcardMatch("exact.name", "other.name"));
}

TEST(NaturalSortTests, ExpandWildcardsNoWildcardPassthrough) {
  const auto res = expandWildcards("clean.dump");
  ASSERT_EQ(res.size(), 1);
  EXPECT_EQ(res[0], "clean.dump");
}

TEST(NaturalSortTests, ExpandWildcardsFilesystemMatching) {
  const std::filesystem::path temp_dir =
      std::filesystem::temp_directory_path() / "naturalsort_test_fs";
  std::filesystem::create_directories(temp_dir);

  const std::string f_1 = (temp_dir / "sample_10.txt").string();
  const std::string f_2 = (temp_dir / "sample_2.txt").string();
  {
    std::ofstream(f_1) << "10\n";
    std::ofstream(f_2) << "2\n";
  }

  auto matches = detail::expandWildcardsFilesystem((temp_dir / "sample_*.txt").string());
  ASSERT_EQ(matches.size(), 2);
  naturalSort(matches);

  EXPECT_EQ(matches[0], f_2);
  EXPECT_EQ(matches[1], f_1);

  std::filesystem::remove_all(temp_dir);
}

TEST(NaturalSortTests, ExpandWildcardsThrowsOnNoMatch) {
  EXPECT_THROW(expandWildcards("non_existent_dir_12345/*.xyz"), std::runtime_error);
  EXPECT_THROW(detail::expandWildcardsFilesystem("non_existent_dir_12345/*.xyz"),
               std::runtime_error);
}
