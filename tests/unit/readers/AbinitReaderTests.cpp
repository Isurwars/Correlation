/**
 * @file AbinitReaderTests.cpp
 * @brief Unit tests for AbinitReader structure and trajectory parsing.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/AbinitReader.hpp"
#include "readers/ReaderFactory.hpp"

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>

namespace correlation::readers {

namespace {

class AbinitReaderTests : public ::testing::Test {
protected:
  std::string test_file{"test_abinit_sample.abi"};

  void TearDown() override {
    if (std::filesystem::exists(test_file)) {
      std::filesystem::remove(test_file);
    }
  }
};

TEST_F(AbinitReaderTests, DiscoveryInReaderFactory) {
  auto &factory = ReaderFactory::instance();
  const auto *reader = factory.getReaderForExtension(".abi");

  ASSERT_NE(reader, nullptr);
  EXPECT_EQ(reader->getName(), "ABINIT Reader");
  EXPECT_TRUE(reader->isTrajectory());
}

TEST_F(AbinitReaderTests, ParsesAcellAndXangst) {
  std::ofstream out(test_file);
  out << "# ABINIT input sample\n";
  out << "acell 10.0 10.0 10.0\n";
  out << "rprim\n";
  out << "  1.0 0.0 0.0\n";
  out << "  0.0 1.0 0.0\n";
  out << "  0.0 0.0 1.0\n";
  out << "\n";
  out << "xangst\n";
  out << "  C  0.0 0.0 0.0\n";
  out << "  O  1.2 0.0 0.0\n";
  out << "\n";
  out.close();

  AbinitReader reader;
  auto cell = reader.readStructure(test_file);

  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_NEAR(cell.atoms()[1].position().x(), 1.2, 1e-4);
}

TEST_F(AbinitReaderTests, ThrowsOnNonExistentFile) {
  AbinitReader reader;
  EXPECT_THROW(reader.readStructure("non_existent_file.abi"), std::runtime_error);
}

} // namespace
} // namespace correlation::readers
