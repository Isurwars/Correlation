// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "readers/MaceReader.hpp"
#include "readers/ReaderFactory.hpp"

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>

namespace {

class MaceReaderTests : public ::testing::Test {
public:
  struct TempFileSpec {
    std::string filename;
    std::string content;
  };

  [[nodiscard]] std::filesystem::path createTempFile(const TempFileSpec &spec) const {
    auto file_path = test_dir_ / spec.filename;
    std::ofstream ofs(file_path);
    ofs << spec.content;
    return file_path;
  }

protected:
  void SetUp() override {
    test_dir_ = std::filesystem::temp_directory_path() / "correlation_mace_tests";
    std::filesystem::create_directories(test_dir_);
  }

  void TearDown() override {
    std::error_code err;
    std::filesystem::remove_all(test_dir_, err);
  }

private:
  std::filesystem::path test_dir_;
};

TEST_F(MaceReaderTests, ReadSingleFrameMace) {
  std::string const mace_content =
      "2\n"
      "Lattice=\"5.0 0.0 0.0 0.0 6.0 0.0 0.0 0.0 7.0\" Properties=species:S:1:pos:R:3:mace_forces:R:3 "
      "mace_energy=-25.432\n"
      "Si 0.1 0.2 0.3 0.01 -0.02 0.03\n"
      "O 1.5 2.5 3.5 -0.04 0.05 -0.06\n";

  auto file_path = createTempFile({.filename = "single_frame.extxyz", .content = mace_content});

  correlation::readers::MaceReader reader;
  auto cell = reader.readStructure(file_path.string());

  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_THAT(cell.getEnergy(), correlation::testing::IsRealEq(-25.432));

  const auto &params = cell.lattice_parameters();
  EXPECT_THAT(params[0], correlation::testing::IsRealEq(5.0));
  EXPECT_THAT(params[1], correlation::testing::IsRealEq(6.0));
  EXPECT_THAT(params[2], correlation::testing::IsRealEq(7.0));

  const auto &atoms = cell.atoms();
  EXPECT_EQ(atoms[0].element().symbol, "Si");
  EXPECT_THAT(atoms[0].position().x(), correlation::testing::IsRealEq(0.1));
  EXPECT_THAT(atoms[0].position().y(), correlation::testing::IsRealEq(0.2));
  EXPECT_THAT(atoms[0].position().z(), correlation::testing::IsRealEq(0.3));

  // Verify forces stored in velocity
  EXPECT_THAT(atoms[0].velocity().x(), correlation::testing::IsRealEq(0.01));
  EXPECT_THAT(atoms[0].velocity().y(), correlation::testing::IsRealEq(-0.02));
  EXPECT_THAT(atoms[0].velocity().z(), correlation::testing::IsRealEq(0.03));

  EXPECT_EQ(atoms[1].element().symbol, "O");
  EXPECT_THAT(atoms[1].position().x(), correlation::testing::IsRealEq(1.5));
  EXPECT_THAT(atoms[1].velocity().x(), correlation::testing::IsRealEq(-0.04));
}

TEST_F(MaceReaderTests, ReadMultiFrameTrajectory) {
  std::string const traj_content =
      "2\n"
      "Lattice=\"10.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 10.0\" Properties=species:S:1:pos:R:3:mace_forces:R:3 "
      "mace_energy=-100.5\n"
      "C 0.0 0.0 0.0 0.0 0.0 0.0\n"
      "C 1.5 0.0 0.0 0.1 0.0 0.0\n"
      "2\n"
      "Lattice=\"10.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 10.0\" Properties=species:S:1:pos:R:3:mace_forces:R:3 "
      "mace_energy=-102.3\n"
      "C 0.0 0.0 0.0 -0.05 0.0 0.0\n"
      "C 1.6 0.0 0.0 0.05 0.0 0.0\n";

  auto file_path = createTempFile({.filename = "trajectory.mace", .content = traj_content});

  correlation::readers::MaceReader reader;
  auto traj = reader.readTrajectory(file_path.string());

  EXPECT_EQ(traj.getFrameCount(), 2);

  auto frame_0 = traj.getFrame(0);
  EXPECT_EQ(frame_0.atomCount(), 2);
  EXPECT_THAT(frame_0.getEnergy(), correlation::testing::IsRealEq(-100.5));
  EXPECT_THAT(frame_0.atoms()[1].position().x(), correlation::testing::IsRealEq(1.5));

  auto frame_1 = traj.getFrame(1);
  EXPECT_EQ(frame_1.atomCount(), 2);
  EXPECT_THAT(frame_1.getEnergy(), correlation::testing::IsRealEq(-102.3));
  EXPECT_THAT(frame_1.atoms()[1].position().x(), correlation::testing::IsRealEq(1.6));
  EXPECT_THAT(frame_1.atoms()[1].velocity().x(), correlation::testing::IsRealEq(0.05));
}

TEST_F(MaceReaderTests, ReadStructureReturnsLastFrame) {
  std::string const traj_content =
      "1\n"
      "Properties=species:S:1:pos:R:3 energy=-10.0\n"
      "H 0.0 0.0 0.0\n"
      "1\n"
      "Properties=species:S:1:pos:R:3 energy=-15.0\n"
      "H 1.0 1.0 1.0\n";

  auto file_path = createTempFile({.filename = "last_frame.xyz", .content = traj_content});

  correlation::readers::MaceReader reader;
  auto cell = reader.readStructure(file_path.string());

  EXPECT_EQ(cell.atomCount(), 1);
  EXPECT_THAT(cell.getEnergy(), correlation::testing::IsRealEq(-15.0));
  EXPECT_THAT(cell.atoms()[0].position().x(), correlation::testing::IsRealEq(1.0));
}

TEST_F(MaceReaderTests, MalformedHeaderThrows) {
  std::string const malformed_content =
      "invalid_count\n"
      "Properties=species:S:1:pos:R:3\n"
      "H 0.0 0.0 0.0\n";

  auto file_path = createTempFile({.filename = "malformed.extxyz", .content = malformed_content});

  correlation::readers::MaceReader reader;
  EXPECT_THROW(reader.readStructure(file_path.string()), std::runtime_error);
}

TEST_F(MaceReaderTests, ReaderFactoryRegistration) {
  auto *reader = correlation::readers::ReaderFactory::instance().getReaderForExtension(".mace");
  ASSERT_NE(reader, nullptr);
  EXPECT_EQ(reader->getName(), "MACE Reader");
  EXPECT_TRUE(reader->isTrajectory());
}

} // namespace
