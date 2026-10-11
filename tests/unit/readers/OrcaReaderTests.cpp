/**
 * @file OrcaReaderTests.cpp
 * @brief Unit tests for OrcaReader structure and optimization trajectory parsing.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/OrcaReader.hpp"
#include "readers/ReaderFactory.hpp"

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>

namespace correlation::readers {

namespace {

class OrcaReaderTests : public ::testing::Test {
protected:
  [[nodiscard]] const std::string &testFile() const noexcept { return test_file_; }

  void SetUp() override {
    const auto *test_info = ::testing::UnitTest::GetInstance()->current_test_info();
    const std::string name = (test_info != nullptr) ? test_info->name() : "orca";
    test_file_ = (std::filesystem::temp_directory_path() / ("test_orca_" + name + ".out")).string();
  }

  void TearDown() override {
    if (std::filesystem::exists(test_file_)) {
      std::filesystem::remove(test_file_);
    }
  }

private:
  std::string test_file_;
};

TEST_F(OrcaReaderTests, DiscoveryInReaderFactory) {
  auto &factory = ReaderFactory::instance();
  const auto *reader = factory.getReaderForExtension(".orca");

  ASSERT_NE(reader, nullptr);
  EXPECT_EQ(reader->getName(), "ORCA Reader");
  EXPECT_TRUE(reader->isTrajectory());
}

TEST_F(OrcaReaderTests, ParsesSingleStructureAngstrom) {
  std::ofstream out(testFile());
  out << "                                 *****************\n";
  out << "                                 * O   R   C   A *\n";
  out << "                                 *****************\n";
  out << "\n";
  out << "---------------------------------\n";
  out << "CARTESIAN COORDINATES (ANGSTROEM)\n";
  out << "---------------------------------\n";
  out << "  C      0.000000    0.000000    0.000000\n";
  out << "  O      1.160000    0.000000    0.000000\n";
  out << "  O     -1.160000    0.000000    0.000000\n";
  out << "\n";
  out.close();

  OrcaReader reader;
  auto cell = reader.readStructure(testFile());

  EXPECT_EQ(cell.atomCount(), 3);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "C");
  EXPECT_NEAR(cell.atoms()[1].position().x(), 1.16, 1e-4);
  EXPECT_NEAR(cell.atoms()[2].position().x(), -1.16, 1e-4);
}

TEST_F(OrcaReaderTests, ParsesOptimizationTrajectory) {
  std::ofstream out(testFile());
  out << "O   R   C   A\n";
  out << "CARTESIAN COORDINATES (ANGSTROEM)\n";
  out << "---------------------------------\n";
  out << "  H      0.000000    0.000000    0.000000\n";
  out << "  H      0.000000    0.000000    1.000000\n";
  out << "\n";
  out << "GEOMETRY OPTIMIZATION CYCLE   1\n";
  out << "CARTESIAN COORDINATES (ANGSTROEM)\n";
  out << "---------------------------------\n";
  out << "  H      0.000000    0.000000    0.000000\n";
  out << "  H      0.000000    0.000000    0.740000\n";
  out << "\n";
  out.close();

  OrcaReader reader;
  auto trajectory = reader.readTrajectory(testFile());

  EXPECT_EQ(trajectory.getFrameCount(), 2);
  EXPECT_NEAR(trajectory.getFrames()[1].atoms()[1].position().z(), 0.74, 1e-4);
}

TEST_F(OrcaReaderTests, ThrowsOnNonExistentFile) {
  OrcaReader reader;
  EXPECT_THROW(reader.readStructure("non_existent_file.out"), std::runtime_error);
}

} // namespace
} // namespace correlation::readers
