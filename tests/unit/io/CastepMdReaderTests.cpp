#include "math/Constants.hpp"
#include "readers/CastepMdReader.hpp"
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>

using namespace correlation::readers;

TEST(CastepMdReaderTests, Properties) {
  CastepMdReader reader;
  EXPECT_EQ(reader.getName(), "CASTEP MD");
  EXPECT_TRUE(reader.isTrajectory());
  auto exts = reader.getExtensions();
  EXPECT_EQ(exts.size(), 1);
  EXPECT_EQ(exts[0], "md");
}

TEST(CastepMdReaderTests, ReadsTrajectory) {
  std::ofstream out("temp_md.md");
  // Write 2 frames
  out << " BEGIN header\n\n This is 8 atom cubic Si cell\n END header\n\n"
      // Frame 1
      << "               0.00000000E+000\n"
      << "              -3.18206146E+001     -3.18108683E+001      9.74270683E-003  <-- E\n"
      << "                9.27876841E-04                                            <-- T\n"
      << "                5.85402338E-06                                            <-- P\n"
      << "               10.00000000E+001      0.00000000E+000      0.00000000E+000  <-- h\n"
      << "               0.00000000E+000      11.00000000E+001      0.00000000E+000  <-- h\n"
      << "               0.00000000E+000      0.00000000E+000      12.00000000E+001  <-- h\n"
      << " Si     1      1.00000000E+000      2.00000000E+000      3.00000000E+000  <-- R\n"
      // Frame 2 (no new lattice <-- h, should reuse frame 1 lattice but update E and positions)
      << "               1.00000000E+000\n"
      << "              -4.00000000E+001     -3.99999999E+001      9.74270683E-003  <-- E\n"
      << "                9.27876841E-04                                            <-- T\n"
      << "                5.85402338E-06                                            <-- P\n"
      << " Si     1      1.10000000E+000      2.10000000E+000      3.10000000E+000  <-- R\n";
  out.close();

  CastepMdReader reader;
  auto traj = reader.readTrajectory("temp_md.md");
  auto struct_cell = reader.readStructure("temp_md.md");
  std::remove("temp_md.md");

  // Trajectory checks
  EXPECT_EQ(traj.getFrameCount(), 2);

  // Frame 1 check
  const auto &f1 = traj.getFrame(0);
  EXPECT_EQ(f1.atomCount(), 1);
  EXPECT_DOUBLE_EQ(f1.getEnergy(), -31.8206146);
  // Bohr to angstrom conversion checking:
  // lattice a = 100.0 bohr = 100.0 * 0.529177210903...
  EXPECT_NEAR(f1.lattice_parameters()[0], 100.0 * correlation::math::bohr_to_angstrom, 1e-5);
  EXPECT_NEAR(f1.atoms()[0].position().x(), 1.0 * correlation::math::bohr_to_angstrom, 1e-5);

  // Frame 2 check (reuses lattice but has different energy and positions)
  const auto &f2 = traj.getFrame(1);
  EXPECT_EQ(f2.atomCount(), 1);
  EXPECT_DOUBLE_EQ(f2.getEnergy(), -40.0);
  EXPECT_NEAR(f2.lattice_parameters()[0], 100.0 * correlation::math::bohr_to_angstrom, 1e-5);
  EXPECT_NEAR(f2.atoms()[0].position().x(), 1.1 * correlation::math::bohr_to_angstrom, 1e-5);

  // readStructure returns the first frame
  EXPECT_EQ(struct_cell.atomCount(), 1);
  EXPECT_DOUBLE_EQ(struct_cell.getEnergy(), -31.8206146);
}
