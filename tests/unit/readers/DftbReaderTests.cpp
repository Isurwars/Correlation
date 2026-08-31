/**
 * @file DftbReaderTests.cpp
 * @brief Unit tests for DftbReader GenFormat (.gen) structure parsing.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/DftbReader.hpp"
#include "readers/ReaderFactory.hpp"

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>

namespace correlation::readers {

namespace {

class DftbReaderTests : public ::testing::Test {
protected:
  std::string test_file{"test_dftb_sample.gen"};

  void TearDown() override {
    if (std::filesystem::exists(test_file)) {
      std::filesystem::remove(test_file);
    }
  }
};

TEST_F(DftbReaderTests, DiscoveryInReaderFactory) {
  auto &factory = ReaderFactory::instance();
  const auto *reader = factory.getReaderForExtension(".gen");

  ASSERT_NE(reader, nullptr);
  EXPECT_EQ(reader->getName(), "DFTB+ Reader");
  EXPECT_FALSE(reader->isTrajectory());
}

TEST_F(DftbReaderTests, ParsesClusterMode) {
  std::ofstream out(test_file);
  out << "3 C\n";
  out << "C O\n";
  out << "1 1 0.000 0.000 0.000\n";
  out << "2 2 1.160 0.000 0.000\n";
  out << "3 2 -1.160 0.000 0.000\n";
  out.close();

  DftbReader reader;
  auto cell = reader.readStructure(test_file);

  EXPECT_EQ(cell.atomCount(), 3);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "C");
  EXPECT_EQ(cell.atoms()[1].element().symbol, "O");
  EXPECT_EQ(cell.atoms()[2].element().symbol, "O");
  EXPECT_NEAR(cell.atoms()[1].position().x(), 1.16, 1e-4);
}

TEST_F(DftbReaderTests, ParsesPeriodicCartesianMode) {
  std::ofstream out(test_file);
  out << "2 S\n";
  out << "Si\n";
  out << "1 1 0.000 0.000 0.000\n";
  out << "2 1 2.500 2.500 2.500\n";
  out << "0.0 0.0 0.0\n";
  out << "5.0 0.0 0.0\n";
  out << "0.0 5.0 0.0\n";
  out << "0.0 0.0 5.0\n";
  out.close();

  DftbReader reader;
  auto cell = reader.readStructure(test_file);

  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "Si");
  EXPECT_NEAR(cell.volume(), 125.0, 1e-4);
}

TEST_F(DftbReaderTests, ThrowsOnNonExistentFile) {
  DftbReader reader;
  EXPECT_THROW(reader.readStructure("non_existent_file.gen"), std::runtime_error);
}

} // namespace
} // namespace correlation::readers
