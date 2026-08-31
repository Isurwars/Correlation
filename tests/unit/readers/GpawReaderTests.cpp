/**
 * @file GpawReaderTests.cpp
 * @brief Unit tests for GpawReader structure and trajectory parsing.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/GpawReader.hpp"
#include "readers/ReaderFactory.hpp"

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>

namespace correlation::readers {

namespace {

class GpawReaderTests : public ::testing::Test {
protected:
  std::string test_file{"test_gpaw_sample.txt"};

  void TearDown() override {
    if (std::filesystem::exists(test_file)) {
      std::filesystem::remove(test_file);
    }
  }
};

TEST_F(GpawReaderTests, DiscoveryInReaderFactory) {
  auto &factory = ReaderFactory::instance();
  const auto *reader = factory.getReaderForExtension(".gpaw");

  ASSERT_NE(reader, nullptr);
  EXPECT_EQ(reader->getName(), "GPAW Reader");
  EXPECT_TRUE(reader->isTrajectory());
}

TEST_F(GpawReaderTests, ParsesLatticeAndPositions) {
  std::ofstream out(test_file);
  out << "  ___ ___ ___ _ _ _  \n";
  out << " |   |   |   | | | | \n";
  out << " | | | | | | | | | | \n";
  out << " |__ |  _|___|_____| \n";
  out << " |___|_|             \n";
  out << " GPAW version 24.1.0\n";
  out << "\n";
  out << "Unit cell:\n";
  out << "    10.000000    0.000000    0.000000\n";
  out << "     0.000000   10.000000    0.000000\n";
  out << "     0.000000    0.000000   10.000000\n";
  out << "\n";
  out << "Positions:\n";
  out << "  0 Si   0.000000   0.000000   0.000000\n";
  out << "  1 Si   2.500000   2.500000   2.500000\n";
  out << "\n";
  out.close();

  GpawReader reader;
  auto cell = reader.readStructure(test_file);

  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_EQ(cell.elements().size(), 1);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "Si");
  EXPECT_NEAR(cell.volume(), 1000.0, 1e-4);
  EXPECT_NEAR(cell.atoms()[1].position().x(), 2.5, 1e-4);
}

TEST_F(GpawReaderTests, ThrowsOnNonExistentFile) {
  GpawReader reader;
  EXPECT_THROW(reader.readStructure("non_existent_file.gpaw"), std::runtime_error);
}

} // namespace
} // namespace correlation::readers
