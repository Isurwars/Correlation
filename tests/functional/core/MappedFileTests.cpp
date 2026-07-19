// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/MappedFile.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <string>

namespace correlation::testing {

using namespace correlation::core;

namespace {
struct MappedFileFunctionalTests : public ::testing::Test {
  std::string test_dir_ = "mapped_file_functional_test_data";
  std::string file_a_path_ = test_dir_ + "/file_a.txt";
  std::string file_b_path_ = test_dir_ + "/file_b.txt";

  std::string content_a_ = "File A: Initial trajectory data with several frames and coordinates.";
  std::string content_b_ = "File B: Secondary configuration and parameters.";

  void SetUp() override {
    std::filesystem::create_directory(test_dir_);

    std::ofstream out_a(file_a_path_);
    out_a << content_a_;
    out_a.close();

    std::ofstream out_b(file_b_path_);
    out_b << content_b_;
    out_b.close();
  }

  void TearDown() override { std::filesystem::remove_all(test_dir_); }
};

TEST_F(MappedFileFunctionalTests, VerifyMoveAssignmentReleasesPreviousMapping) {
  // Construct two mapped files
  MappedFile mf_a(file_a_path_);
  MappedFile mf_b(file_b_path_);

  const char *ptr_b = mf_b.data();
  size_t size_b = mf_b.size();

  // Move-assign mf_b into mf_a. mf_a's original resources (file_a) must be released.
  mf_a = std::move(mf_b);

  // mf_a should now map file_b
  EXPECT_EQ(mf_a.data(), ptr_b);
  EXPECT_EQ(mf_a.size(), size_b);
  std::string read_content(mf_a.data(), mf_a.size());
  EXPECT_EQ(read_content, content_b_);

  // mf_b should be reset
  EXPECT_EQ(mf_b.data(), nullptr);
  EXPECT_EQ(mf_b.size(), 0);
}

TEST_F(MappedFileFunctionalTests, VerifyDataBufferIterationAndSearch) {
  // Test reading and parsing patterns in mapped memory using standard algorithms
  MappedFile mapped_file(file_a_path_);

  const char *begin = mapped_file.data();
  const char *end = mapped_file.end();

  // Find a specific word "trajectory" in the mapped data
  std::string search_target = "trajectory";
  const auto *iterator = std::search(begin, end, search_target.begin(), search_target.end());

  ASSERT_NE(iterator, end);

  // Verify characters match
  std::string found_word(iterator, iterator + search_target.size());
  EXPECT_EQ(found_word, search_target);

  // Count spaces in the mapped file
  long long space_count = std::count(begin, end, ' ');
  EXPECT_EQ(space_count, 9); // "File A: Initial trajectory data with several frames and coordinates." has 9 spaces
}

TEST_F(MappedFileFunctionalTests, VerifyPostCreationAppendsAreNotMapped) {
  // Open the file first
  MappedFile mapped_file(file_a_path_);
  size_t initial_size = mapped_file.size();

  // Append data to the file on disk
  std::ofstream out(file_a_path_, std::ios::app);
  out << " Extra appended text.";
  out.close();

  // The mapped file size and contents should remain unchanged (read-only view established at open)
  EXPECT_EQ(mapped_file.size(), initial_size);
  std::string current_content(mapped_file.data(), mapped_file.size());
  EXPECT_EQ(current_content, content_a_);
}

TEST_F(MappedFileFunctionalTests, VerifyEmptyFileDoesNotCrash) {
  std::string empty_file_path = test_dir_ + "/empty.txt";
  std::ofstream out(empty_file_path);
  out.close();

  EXPECT_NO_THROW({
    MappedFile mapped_file(empty_file_path);
    EXPECT_EQ(mapped_file.size(), 0);
    EXPECT_EQ(mapped_file.data(), nullptr);
  });
}

TEST_F(MappedFileFunctionalTests, VerifyThrowsOnNonExistentFile) {
  std::string non_existent_path = test_dir_ + "/does_not_exist.txt";
  EXPECT_THROW({
    MappedFile mapped_file(non_existent_path);
  }, std::runtime_error);
}
} // namespace
} // namespace correlation::testing
