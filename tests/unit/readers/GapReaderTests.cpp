/**
 * @file GapReaderTests.cpp
 * @brief Unit tests for GAP (Gaussian Approximation Potentials / QUIP) trajectory reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */
#include "readers/GapReader.hpp"
#include "readers/ReaderFactory.hpp"

#include <filesystem>
#include <fstream>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {

struct TempFileInfo {
  std::string filename;
  std::string content;
};

class GapReaderTests : public ::testing::Test {
protected:
  void SetUp() override {
    temp_dir_ = std::filesystem::temp_directory_path() / "gap_reader_tests";
    std::filesystem::create_directories(temp_dir_);
  }

  void TearDown() override {
    std::error_code ec;
    std::filesystem::remove_all(temp_dir_, ec);
  }

  std::filesystem::path createTempFile(const TempFileInfo &info) {
    auto file_path = temp_dir_ / info.filename;
    std::ofstream ofs(file_path);
    ofs << info.content;
    return file_path;
  }

  std::filesystem::path temp_dir_;
};

TEST_F(GapReaderTests, ReadSingleFrameGap) {
  const std::string content =
      "2\n"
      "Lattice=\"5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0\" Properties=species:S:1:pos:R:3:gap_forces:R:3 "
      "gap_energy=-42.125\n"
      "C 0.0 0.0 0.0 0.1 -0.2 0.3\n"
      "C 1.5 1.5 1.5 -0.1 0.2 -0.3\n";

  auto file_path = createTempFile({.filename = "single_frame.gap", .content = content});

  correlation::readers::GapReader reader;
  correlation::core::Cell cell = reader.readStructure(file_path.string());

  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_THAT(cell.getEnergy(), correlation::testing::IsRealEq(-42.125));

  const auto lat = cell.lattice_parameters();
  EXPECT_THAT(lat[0], correlation::testing::IsRealEq(5.0));
  EXPECT_THAT(lat[1], correlation::testing::IsRealEq(5.0));
  EXPECT_THAT(lat[2], correlation::testing::IsRealEq(5.0));

  const auto &atoms = cell.atoms();
  EXPECT_EQ(atoms[0].element().symbol, "C");
  EXPECT_THAT(atoms[0].position().x(), correlation::testing::IsRealEq(0.0));
  EXPECT_THAT(atoms[0].velocity().x(), correlation::testing::IsRealEq(0.1));
  EXPECT_THAT(atoms[0].velocity().y(), correlation::testing::IsRealEq(-0.2));
  EXPECT_THAT(atoms[0].velocity().z(), correlation::testing::IsRealEq(0.3));

  EXPECT_EQ(atoms[1].element().symbol, "C");
  EXPECT_THAT(atoms[1].position().x(), correlation::testing::IsRealEq(1.5));
  EXPECT_THAT(atoms[1].velocity().x(), correlation::testing::IsRealEq(-0.1));
  EXPECT_THAT(atoms[1].velocity().y(), correlation::testing::IsRealEq(0.2));
  EXPECT_THAT(atoms[1].velocity().z(), correlation::testing::IsRealEq(-0.3));
}

TEST_F(GapReaderTests, ReadMultiFrameTrajectory) {
  const std::string content =
      "2\n"
      "Lattice=\"5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0\" Properties=species:S:1:pos:R:3:force:R:3 gap_energy=-10.0\n"
      "Si 0.0 0.0 0.0 0.01 0.02 0.03\n"
      "Si 2.0 2.0 2.0 -0.01 -0.02 -0.03\n"
      "2\n"
      "Lattice=\"5.1 0.0 0.0 0.0 5.1 0.0 0.0 0.0 5.1\" Properties=species:S:1:pos:R:3:force:R:3 gap_energy=-10.5\n"
      "Si 0.1 0.1 0.1 0.02 0.03 0.04\n"
      "Si 2.1 2.1 2.1 -0.02 -0.03 -0.04\n"
      "2\n"
      "Lattice=\"5.2 0.0 0.0 0.0 5.2 0.0 0.0 0.0 5.2\" Properties=species:S:1:pos:R:3:force:R:3 gap_energy=-11.0\n"
      "Si 0.2 0.2 0.2 0.03 0.04 0.05\n"
      "Si 2.2 2.2 2.2 -0.03 -0.04 -0.05\n";

  auto file_path = createTempFile({.filename = "trajectory.quip", .content = content});

  correlation::readers::GapReader reader;
  correlation::core::Trajectory traj = reader.readTrajectory(file_path.string());

  EXPECT_EQ(traj.getFrameCount(), 3);

  auto frame_0 = traj.getFrame(0);
  EXPECT_THAT(frame_0.getEnergy(), correlation::testing::IsRealEq(-10.0));
  EXPECT_THAT(frame_0.lattice_parameters()[0], correlation::testing::IsRealEq(5.0));
  EXPECT_THAT(frame_0.atoms()[0].velocity().x(), correlation::testing::IsRealEq(0.01));

  auto frame_1 = traj.getFrame(1);
  EXPECT_THAT(frame_1.getEnergy(), correlation::testing::IsRealEq(-10.5));
  EXPECT_THAT(frame_1.lattice_parameters()[0], correlation::testing::IsRealEq(5.1));
  EXPECT_THAT(frame_1.atoms()[0].velocity().x(), correlation::testing::IsRealEq(0.02));

  auto frame_2 = traj.getFrame(2);
  EXPECT_THAT(frame_2.getEnergy(), correlation::testing::IsRealEq(-11.0));
  EXPECT_THAT(frame_2.lattice_parameters()[0], correlation::testing::IsRealEq(5.2));
  EXPECT_THAT(frame_2.atoms()[0].velocity().x(), correlation::testing::IsRealEq(0.03));
}

TEST_F(GapReaderTests, ReadStructureReturnsLastFrame) {
  const std::string content =
      "1\n"
      "Lattice=\"4.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 4.0\" Properties=species:S:1:pos:R:3 gap_energy=-5.0\n"
      "C 0.0 0.0 0.0\n"
      "1\n"
      "Lattice=\"4.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 4.0\" Properties=species:S:1:pos:R:3 gap_energy=-6.0\n"
      "C 1.0 1.0 1.0\n";

  auto file_path = createTempFile({.filename = "two_frames.gap", .content = content});

  correlation::readers::GapReader reader;
  correlation::core::Cell cell = reader.readStructure(file_path.string());

  EXPECT_EQ(cell.atomCount(), 1);
  EXPECT_THAT(cell.getEnergy(), correlation::testing::IsRealEq(-6.0));
  EXPECT_THAT(cell.atoms()[0].position().x(), correlation::testing::IsRealEq(1.0));
}

TEST_F(GapReaderTests, MalformedHeaderThrows) {
  const std::string malformed_count =
      "not_a_number\n"
      "Properties=species:S:1:pos:R:3\n"
      "C 0.0 0.0 0.0\n";

  auto file1 = createTempFile({.filename = "malformed1.gap", .content = malformed_count});
  correlation::readers::GapReader reader;
  EXPECT_THROW(reader.readStructure(file1.string()), std::runtime_error);

  const std::string missing_comment = "1\n";
  auto file2 = createTempFile({.filename = "malformed2.gap", .content = missing_comment});
  EXPECT_THROW(reader.readStructure(file2.string()), std::runtime_error);
}

TEST_F(GapReaderTests, ReaderFactoryRegistration) {
  using correlation::readers::ReaderFactory;

  auto *reader_by_ext = ReaderFactory::instance().getReaderForExtension(".gap");
  ASSERT_NE(reader_by_ext, nullptr);
  EXPECT_EQ(reader_by_ext->getName(), "GAP Reader");

  auto *reader_by_quip = ReaderFactory::instance().getReaderForExtension(".quip");
  ASSERT_NE(reader_by_quip, nullptr);
  EXPECT_EQ(reader_by_quip->getName(), "GAP Reader");
}

} // namespace
