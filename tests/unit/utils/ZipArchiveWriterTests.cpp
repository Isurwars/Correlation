/**
 * @file ZipArchiveWriterTests.cpp
 * @brief Unit tests for miniz RAII wrapper ZipArchiveWriter.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "utils/ZipArchiveWriter.hpp"

#include <gtest/gtest.h>
#include <miniz.h>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace {

class ZipArchiveWriterTest : public ::testing::Test {
protected:
  void SetUp() override {
    test_dir_ = std::filesystem::current_path() / "test_zip_out";
    std::filesystem::create_directories(test_dir_);
  }

  void TearDown() override {
    std::error_code ec;
    std::filesystem::remove_all(test_dir_, ec);
  }

  std::filesystem::path test_dir_;
};

TEST_F(ZipArchiveWriterTest, BasicCreationAndReadback) {
  const auto zip_path = test_dir_ / "test_basic.zip";

  {
    correlation::utils::ZipArchiveWriter writer;
    auto open_res = writer.open(zip_path);
    ASSERT_TRUE(open_res.has_value()) << open_res.error();
    EXPECT_TRUE(writer.isOpen());

    auto add_mem_res = writer.addFileFromString("hello.txt", "Hello, World!");
    EXPECT_TRUE(add_mem_res.has_value()) << add_mem_res.error();

    auto add_sub_res = writer.addFileFromString("nested/data.csv", "1,2,3\n4,5,6\n");
    EXPECT_TRUE(add_sub_res.has_value()) << add_sub_res.error();

    auto finalize_res = writer.finalize();
    EXPECT_TRUE(finalize_res.has_value()) << finalize_res.error();
    EXPECT_FALSE(writer.isOpen());
  }

  ASSERT_TRUE(std::filesystem::exists(zip_path));
  EXPECT_GT(std::filesystem::file_size(zip_path), 0u);

  // Verify using miniz reader API directly
  mz_zip_archive zip;
  mz_zip_zero_struct(&zip);
  mz_bool init_ok = mz_zip_reader_init_file(&zip, zip_path.string().c_str(), 0);
  ASSERT_TRUE(init_ok);

  mz_uint num_files = mz_zip_reader_get_num_files(&zip);
  EXPECT_EQ(num_files, 2u);

  int hello_idx = mz_zip_reader_locate_file(&zip, "hello.txt", nullptr, 0);
  EXPECT_GE(hello_idx, 0);

  size_t hello_size = 0;
  void *hello_buf = mz_zip_reader_extract_to_heap(&zip, hello_idx, &hello_size, 0);
  ASSERT_NE(hello_buf, nullptr);
  std::string hello_str(static_cast<const char *>(hello_buf), hello_size);
  mz_free(hello_buf);
  EXPECT_EQ(hello_str, "Hello, World!");

  int data_idx = mz_zip_reader_locate_file(&zip, "nested/data.csv", nullptr, 0);
  EXPECT_GE(data_idx, 0);

  size_t data_size = 0;
  void *data_buf = mz_zip_reader_extract_to_heap(&zip, data_idx, &data_size, 0);
  ASSERT_NE(data_buf, nullptr);
  std::string data_str(static_cast<const char *>(data_buf), data_size);
  mz_free(data_buf);
  EXPECT_EQ(data_str, "1,2,3\n4,5,6\n");

  mz_zip_reader_end(&zip);
}

TEST_F(ZipArchiveWriterTest, AddFromDisk) {
  const auto src_file = test_dir_ / "source.txt";
  {
    std::ofstream out(src_file);
    out << "disk file content\nsecond line";
  }

  const auto zip_path = test_dir_ / "test_disk.zip";
  {
    correlation::utils::ZipArchiveWriter writer;
    ASSERT_TRUE(writer.open(zip_path).has_value());
    auto add_res = writer.addFileFromDisk("from_disk/source.txt", src_file);
    EXPECT_TRUE(add_res.has_value()) << add_res.error();
    EXPECT_TRUE(writer.finalize().has_value());
  }

  // Verify
  mz_zip_archive zip;
  mz_zip_zero_struct(&zip);
  ASSERT_TRUE(mz_zip_reader_init_file(&zip, zip_path.string().c_str(), 0));

  int idx = mz_zip_reader_locate_file(&zip, "from_disk/source.txt", nullptr, 0);
  EXPECT_GE(idx, 0);

  size_t sz = 0;
  void *buf = mz_zip_reader_extract_to_heap(&zip, idx, &sz, 0);
  ASSERT_NE(buf, nullptr);
  std::string content(static_cast<const char *>(buf), sz);
  mz_free(buf);
  EXPECT_EQ(content, "disk file content\nsecond line");

  mz_zip_reader_end(&zip);
}

TEST_F(ZipArchiveWriterTest, RAIIAutoFinalizeDestruction) {
  const auto zip_path = test_dir_ / "test_raii.zip";
  {
    correlation::utils::ZipArchiveWriter writer;
    ASSERT_TRUE(writer.open(zip_path).has_value());
    ASSERT_TRUE(writer.addFileFromString("doc.txt", "abc").has_value());
    // Exiting scope without explicit finalize() should clean up cleanly without crashing
  }
}

} // namespace
