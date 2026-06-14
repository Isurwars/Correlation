#include "math/Constants.hpp"
#include "readers/OutmolReader.hpp"
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>

using namespace correlation::readers;

TEST(OutmolReaderTests, Properties) {
  OutmolReader const reader;
  EXPECT_EQ(reader.getName(), "Outmol");
  EXPECT_TRUE(reader.isTrajectory());
  auto exts = reader.getExtensions();
  EXPECT_EQ(exts.size(), 1);
  EXPECT_EQ(exts[0], "outmol");
}

TEST(OutmolReaderTests, ReadsTrajectoryFormat1) {
  // Format 1: $cell vectors + $coordinates
  std::ofstream out("temp_outmol_fmt1.outmol");
  out << "$cell vectors\n"
      << " 10.0 0.0 0.0\n"
      << " 0.0 10.0 0.0\n"
      << " 0.0 0.0 10.0\n"
      << "$coordinates\n"
      << " C 1.0 2.0 3.0\n"
      << " H 1.5 2.5 3.5\n"
      << "$end\n";
  out.close();

  OutmolReader reader;
  auto traj = reader.readTrajectory("temp_outmol_fmt1.outmol");
  auto cell = reader.readStructure("temp_outmol_fmt1.outmol");
  std::remove("temp_outmol_fmt1.outmol");

  EXPECT_EQ(traj.getFrameCount(), 1);
  const auto &f = traj.getFrame(0);
  EXPECT_EQ(f.atomCount(), 2);

  // Bohr to angstrom check
  EXPECT_NEAR(f.lattice_parameters()[0], 10.0 * correlation::math::bohr_to_angstrom, 1e-5);
  EXPECT_NEAR(f.atoms()[0].position().x(), 1.0 * correlation::math::bohr_to_angstrom, 1e-5);
  EXPECT_EQ(f.atoms()[0].element().symbol, "C");
  EXPECT_EQ(f.atoms()[1].element().symbol, "H");

  EXPECT_EQ(cell.atomCount(), 2);
}

TEST(OutmolReaderTests, ReadsTrajectoryFormat2) {
  // Format 2: $cell vectors + ATOMIC COORDINATES (au)
  std::ofstream out("temp_outmol_fmt2.outmol");
  out << "$cell vectors\n"
      << " 12.0 0.0 0.0\n"
      << " 0.0 12.0 0.0\n"
      << " 0.0 0.0 12.0\n"
      << "ATOMIC  COORDINATES (au)\n"
      << "  Header line\n"
      << "  df  Si  2.0 3.0 4.0\n"
      << "  df  O   2.5 3.5 4.5\n"
      << "  other line ends the block\n";
  out.close();

  OutmolReader reader;
  auto traj = reader.readTrajectory("temp_outmol_fmt2.outmol");
  std::remove("temp_outmol_fmt2.outmol");

  EXPECT_EQ(traj.getFrameCount(), 1);
  const auto &f = traj.getFrame(0);
  EXPECT_EQ(f.atomCount(), 2);

  EXPECT_NEAR(f.lattice_parameters()[0], 12.0 * correlation::math::bohr_to_angstrom, 1e-5);
  EXPECT_NEAR(f.atoms()[0].position().x(), 2.0 * correlation::math::bohr_to_angstrom, 1e-5);
  EXPECT_EQ(f.atoms()[0].element().symbol, "Si");
  EXPECT_EQ(f.atoms()[1].element().symbol, "O");
}

TEST(OutmolReaderTests, ThrowsOnInvalidFile) {
  OutmolReader reader;
  EXPECT_THROW(reader.readTrajectory("nonexistent_outmol.outmol"), std::runtime_error);
}
