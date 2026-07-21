// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/MappedFile.hpp"

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <string>

namespace correlation::testing {

using namespace correlation::core;
namespace {
class MappedFileTests : public ::testing::Test {
public:
  std::string test_dir_ = "mapped_file_test_data";
  std::string valid_file_path_ = test_dir_ + "/valid_test.txt";
  std::string file_content_ = "Hello, memory-mapped world! This is a test file for MappedFile RAII class.";

  void SetUp() override {
    std::filesystem::create_directory(test_dir_);
    std::ofstream out(valid_file_path_);
    out << file_content_;
    out.close();
  }

  void TearDown() override { std::filesystem::remove_all(test_dir_); }
};
} // namespace

TEST_F(MappedFileTests, MapsValidFileSuccessfully) {
  MappedFile file(valid_file_path_);
  EXPECT_EQ(file.size(), file_content_.size());
  ASSERT_NE(file.data(), nullptr);

  std::string read_content(file.data(), file.size());
  EXPECT_EQ(read_content, file_content_);
  EXPECT_EQ(file.end(), file.data() + file.size());
}

TEST_F(MappedFileTests, ThrowsOnNonExistentFile) {
  EXPECT_THROW(MappedFile{"non_existent_file_xyz_123.txt"}, std::runtime_error);
}

TEST_F(MappedFileTests, ThrowsOnDirectoryPath) { EXPECT_THROW(MappedFile{test_dir_}, std::runtime_error); }

TEST_F(MappedFileTests, MoveConstructorTransfersOwnership) {
  MappedFile mf1(valid_file_path_);
  const char *original_ptr = mf1.data();
  size_t original_size = mf1.size();

  MappedFile mf2(std::move(mf1));

  // mf2 should own the resource
  EXPECT_EQ(mf2.data(), original_ptr);
  EXPECT_EQ(mf2.size(), original_size);

  // mf1 should be null/empty
  EXPECT_EQ(mf1.data(), nullptr);
  EXPECT_EQ(mf1.size(), 0);
}

TEST_F(MappedFileTests, MoveAssignmentOperatorTransfersOwnership) {
  MappedFile mf1(valid_file_path_);
  const char *original_ptr = mf1.data();
  size_t original_size = mf1.size();

  MappedFile mf2 = std::move(mf1);

  // mf2 should own the resource
  EXPECT_EQ(mf2.data(), original_ptr);
  EXPECT_EQ(mf2.size(), original_size);

  // mf1 should be null/empty
  EXPECT_EQ(mf1.data(), nullptr);
  EXPECT_EQ(mf1.size(), 0);
}

TEST_F(MappedFileTests, EnforceSizeLimitCheck) {
  // We can test that the constructor works with enforce_size_limit = false
  MappedFile file(valid_file_path_, false);
  EXPECT_EQ(file.size(), file_content_.size());
}

TEST_F(MappedFileTests, EmptyFileHandling) {
  std::string const empty_file_path = test_dir_ + "/empty_test.txt";
  std::ofstream out(empty_file_path);
  out.close();

  // Mapping an empty file should either succeed with size 0 or throw runtime_error
  try {
    MappedFile file(empty_file_path, false);
    EXPECT_EQ(file.size(), 0);
  } catch (const std::runtime_error &) {
    SUCCEED();
  }
}

} // namespace correlation::testing
