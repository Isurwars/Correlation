// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/MappedFile.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <memory>
#include <string>

namespace correlation::testing {

using namespace correlation::core;

namespace {
class MappedFileFunctionalTests : public ::testing::Test {
protected:
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

  [[nodiscard]] const std::string &fileAPath() const { return file_a_path_; }
  [[nodiscard]] const std::string &fileBPath() const { return file_b_path_; }
  [[nodiscard]] const std::string &contentA() const { return content_a_; }
  [[nodiscard]] const std::string &contentB() const { return content_b_; }
  [[nodiscard]] const std::string &testDir() const { return test_dir_; }

private:
  std::string test_dir_ = "mapped_file_functional_test_data";
  std::string file_a_path_ = test_dir_ + "/file_a.txt";
  std::string file_b_path_ = test_dir_ + "/file_b.txt";

  std::string content_a_ = "File A: Initial trajectory data with several frames and coordinates.";
  std::string content_b_ = "File B: Secondary configuration and parameters.";
};

TEST_F(MappedFileFunctionalTests, VerifyMoveAssignmentReleasesPreviousMapping) {
  // Construct two mapped files
  MappedFile mf_a(fileAPath());
  auto p_mf_b = std::make_unique<MappedFile>(fileBPath());

  const char *ptr_b = p_mf_b->data();
  const size_t size_b = p_mf_b->size();

  // Move-assign mf_b into mf_a. mf_a's original resources (file_a) must be released.
  mf_a = std::move(*p_mf_b);

  // mf_a should now map file_b
  EXPECT_EQ(mf_a.data(), ptr_b);
  EXPECT_EQ(mf_a.size(), size_b);
  const std::string read_content(mf_a.data(), mf_a.size());
  EXPECT_EQ(read_content, contentB());

  // mf_b should be reset
  EXPECT_EQ(p_mf_b->data(), nullptr);
  EXPECT_EQ(p_mf_b->size(), 0);
}

TEST_F(MappedFileFunctionalTests, VerifyDataBufferIterationAndSearch) {
  // Test reading and parsing patterns in mapped memory using standard algorithms
  const MappedFile mapped_file(fileAPath());

  const char *begin = mapped_file.data();
  const char *end = mapped_file.end();

  // Find a specific word "trajectory" in the mapped data
  const std::string search_target = "trajectory";
  const auto *iterator = std::search(begin, end, search_target.begin(), search_target.end());

  ASSERT_NE(iterator, end);

  // Verify characters match
  const std::string found_word(iterator, iterator + search_target.size());
  EXPECT_EQ(found_word, search_target);

  // Count spaces in the mapped file
  const long long space_count = std::count(begin, end, ' ');
  EXPECT_EQ(
      space_count,
      9); // "File A: Initial trajectory data with several frames and coordinates." has 9 spaces
}

TEST_F(MappedFileFunctionalTests, VerifyPostCreationAppendsAreNotMapped) {
  // Open the file first
  const MappedFile mapped_file(fileAPath());
  const size_t initial_size = mapped_file.size();

  // Append data to the file on disk
  std::ofstream out(fileAPath(), std::ios::app);
  out << " Extra appended text.";
  out.close();

  // The mapped file size and contents should remain unchanged (read-only view established at open)
  EXPECT_EQ(mapped_file.size(), initial_size);
  const std::string current_content(mapped_file.data(), mapped_file.size());
  EXPECT_EQ(current_content, contentA());
}

TEST_F(MappedFileFunctionalTests, VerifyEmptyFileDoesNotCrash) {
  const std::string empty_file_path = testDir() + "/empty.txt";
  std::ofstream out(empty_file_path);
  out.close();

  EXPECT_NO_THROW({
    const MappedFile mapped_file(empty_file_path);
    EXPECT_EQ(mapped_file.size(), 0);
    EXPECT_EQ(mapped_file.data(), nullptr);
  });
}

TEST_F(MappedFileFunctionalTests, VerifyThrowsOnNonExistentFile) {
  const std::string non_existent_path = testDir() + "/does_not_exist.txt";
  EXPECT_THROW({ const MappedFile mapped_file(non_existent_path); }, std::runtime_error);
}
} // namespace
} // namespace correlation::testing
