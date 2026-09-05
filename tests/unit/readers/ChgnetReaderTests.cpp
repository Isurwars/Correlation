// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "readers/ChgnetReader.hpp"
#include "readers/ReaderFactory.hpp"

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>

namespace {

class ChgnetReaderTests : public ::testing::Test {
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
    test_dir_ = std::filesystem::temp_directory_path() / "correlation_chgnet_tests";
    std::filesystem::create_directories(test_dir_);
  }

  void TearDown() override {
    std::error_code err;
    std::filesystem::remove_all(test_dir_, err);
  }

private:
  std::filesystem::path test_dir_;
};

TEST_F(ChgnetReaderTests, ReadSingleFrameChgnet) {
  std::string const chgnet_content =
      "2\n"
      "Lattice=\"4.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 4.0\" Properties=species:S:1:pos:R:3:chgnet_forces:R:3 "
      "chgnet_energy=-18.75\n"
      "Fe 0.0 0.0 0.0 0.02 -0.01 0.05\n"
      "O 1.2 1.2 1.2 -0.02 0.01 -0.05\n";

  auto file_path = createTempFile({.filename = "single_frame.chgnet", .content = chgnet_content});

  correlation::readers::ChgnetReader reader;
  auto cell = reader.readStructure(file_path.string());

  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_THAT(cell.getEnergy(), correlation::testing::IsRealEq(-18.75));

  const auto &params = cell.lattice_parameters();
  EXPECT_THAT(params[0], correlation::testing::IsRealEq(4.0));
  EXPECT_THAT(params[1], correlation::testing::IsRealEq(4.0));
  EXPECT_THAT(params[2], correlation::testing::IsRealEq(4.0));

  const auto &atoms = cell.atoms();
  EXPECT_EQ(atoms[0].element().symbol, "Fe");
  EXPECT_THAT(atoms[0].position().x(), correlation::testing::IsRealEq(0.0));
  EXPECT_THAT(atoms[0].position().y(), correlation::testing::IsRealEq(0.0));
  EXPECT_THAT(atoms[0].position().z(), correlation::testing::IsRealEq(0.0));

  // Verify forces stored in velocity
  EXPECT_THAT(atoms[0].velocity().x(), correlation::testing::IsRealEq(0.02));
  EXPECT_THAT(atoms[0].velocity().y(), correlation::testing::IsRealEq(-0.01));
  EXPECT_THAT(atoms[0].velocity().z(), correlation::testing::IsRealEq(0.05));

  EXPECT_EQ(atoms[1].element().symbol, "O");
  EXPECT_THAT(atoms[1].position().x(), correlation::testing::IsRealEq(1.2));
  EXPECT_THAT(atoms[1].velocity().x(), correlation::testing::IsRealEq(-0.02));
}

TEST_F(ChgnetReaderTests, ReadMultiFrameTrajectory) {
  std::string const traj_content =
      "2\n"
      "Lattice=\"8.0 0.0 0.0 0.0 8.0 0.0 0.0 0.0 8.0\" Properties=species:S:1:pos:R:3:chgnet_forces:R:3 "
      "chgnet_energy=-50.0\n"
      "Fe 0.0 0.0 0.0 0.0 0.0 0.0\n"
      "O 1.0 0.0 0.0 0.1 0.0 0.0\n"
      "2\n"
      "Lattice=\"8.0 0.0 0.0 0.0 8.0 0.0 0.0 0.0 8.0\" Properties=species:S:1:pos:R:3:chgnet_forces:R:3 "
      "chgnet_energy=-51.2\n"
      "Fe 0.0 0.0 0.0 -0.05 0.0 0.0\n"
      "O 1.1 0.0 0.0 0.05 0.0 0.0\n";

  auto file_path = createTempFile({.filename = "trajectory.chgnet", .content = traj_content});

  correlation::readers::ChgnetReader reader;
  auto traj = reader.readTrajectory(file_path.string());

  EXPECT_EQ(traj.getFrameCount(), 2);

  auto frame_0 = traj.getFrame(0);
  EXPECT_EQ(frame_0.atomCount(), 2);
  EXPECT_THAT(frame_0.getEnergy(), correlation::testing::IsRealEq(-50.0));
  EXPECT_THAT(frame_0.atoms()[1].position().x(), correlation::testing::IsRealEq(1.0));

  auto frame_1 = traj.getFrame(1);
  EXPECT_EQ(frame_1.atomCount(), 2);
  EXPECT_THAT(frame_1.getEnergy(), correlation::testing::IsRealEq(-51.2));
  EXPECT_THAT(frame_1.atoms()[1].position().x(), correlation::testing::IsRealEq(1.1));
  EXPECT_THAT(frame_1.atoms()[1].velocity().x(), correlation::testing::IsRealEq(0.05));
}

TEST_F(ChgnetReaderTests, ReadStructureReturnsLastFrame) {
  std::string const traj_content =
      "1\n"
      "Properties=species:S:1:pos:R:3 energy=-5.0\n"
      "Fe 0.0 0.0 0.0\n"
      "1\n"
      "Properties=species:S:1:pos:R:3 energy=-8.0\n"
      "Fe 0.5 0.5 0.5\n";

  auto file_path = createTempFile({.filename = "last_frame.chgnet", .content = traj_content});

  correlation::readers::ChgnetReader reader;
  auto cell = reader.readStructure(file_path.string());

  EXPECT_EQ(cell.atomCount(), 1);
  EXPECT_THAT(cell.getEnergy(), correlation::testing::IsRealEq(-8.0));
  EXPECT_THAT(cell.atoms()[0].position().x(), correlation::testing::IsRealEq(0.5));
}

TEST_F(ChgnetReaderTests, MalformedHeaderThrows) {
  std::string const malformed_content =
      "invalid_atoms\n"
      "Properties=species:S:1:pos:R:3\n"
      "Fe 0.0 0.0 0.0\n";

  auto file_path = createTempFile({.filename = "malformed.chgnet", .content = malformed_content});

  correlation::readers::ChgnetReader reader;
  EXPECT_THROW(reader.readStructure(file_path.string()), std::runtime_error);
}

TEST_F(ChgnetReaderTests, ReaderFactoryRegistration) {
  auto *reader = correlation::readers::ReaderFactory::instance().getReaderForExtension(".chgnet");
  ASSERT_NE(reader, nullptr);
  EXPECT_EQ(reader->getName(), "CHGNet Reader");
  EXPECT_TRUE(reader->isTrajectory());
}

} // namespace
