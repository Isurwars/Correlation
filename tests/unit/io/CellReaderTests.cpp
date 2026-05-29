#include "readers/CellReader.hpp"
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>

using namespace correlation::readers;

TEST(CellReaderTests, Properties) {
  CellReader reader;
  EXPECT_EQ(reader.getName(), "CASTEP CELL");
  EXPECT_FALSE(reader.isTrajectory());
  auto exts = reader.getExtensions();
  EXPECT_EQ(exts.size(), 1);
  EXPECT_EQ(exts[0], "cell");
}

TEST(CellReaderTests, ReadsStructureLatticeAbc) {
  std::ofstream out("temp_cell.cell");
  out << "%BLOCK LATTICE_ABC\n"
      << " 15.0 15.0 20.0\n"
      << " 90.0 90.0 120.0\n"
      << "%ENDBLOCK LATTICE_ABC\n\n"
      << "%BLOCK POSITIONS_ABS\n"
      << " C   1.1 2.2 3.3\n"
      << "%ENDBLOCK POSITIONS_ABS\n";
  out.close();

  CellReader reader;
  auto cell = reader.readStructure("temp_cell.cell");
  std::remove("temp_cell.cell");

  EXPECT_EQ(cell.atomCount(), 1);
  EXPECT_DOUBLE_EQ(cell.lattice_parameters()[0], 15.0);
  EXPECT_DOUBLE_EQ(cell.lattice_parameters()[5], 120.0);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "C");
}

TEST(CellReaderTests, ReadsStructureLatticeCartAndFrac) {
  std::ofstream out("temp_cell_cart.cell");
  out << "%block lattice_cart\n"
      << " 10.0 0.0 0.0\n"
      << " 0.0 10.0 0.0\n"
      << " 0.0 0.0 10.0\n"
      << "%endblock lattice_cart\n\n"
      << "%block positions_frac\n"
      << " Si  0.5 0.5 0.5\n"
      << " O   1.2 0.1 -0.3\n" // outside [0, 1), should wrap
      << "%endblock positions_frac\n";
  out.close();

  CellReader reader;
  auto cell = reader.readStructure("temp_cell_cart.cell");
  std::remove("temp_cell_cart.cell");

  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_DOUBLE_EQ(cell.lattice_parameters()[0], 10.0);
  EXPECT_DOUBLE_EQ(cell.lattice_parameters()[1], 10.0);
  EXPECT_DOUBLE_EQ(cell.lattice_parameters()[2], 10.0);

  EXPECT_EQ(cell.atoms()[0].element().symbol, "Si");
  EXPECT_EQ(cell.atoms()[1].element().symbol, "O");

  // Wrapped position check: Si should be (5.0, 5.0, 5.0)
  const auto &p1 = cell.atoms()[0].position();
  EXPECT_NEAR(p1.x(), 5.0, 1e-5);
  EXPECT_NEAR(p1.y(), 5.0, 1e-5);
  EXPECT_NEAR(p1.z(), 5.0, 1e-5);

  // Wrapped position check: O (1.2 0.1 -0.3) -> (0.2 0.1 0.7) -> (2.0, 1.0, 7.0)
  const auto &p2 = cell.atoms()[1].position();
  EXPECT_NEAR(p2.x(), 2.0, 1e-5);
  EXPECT_NEAR(p2.y(), 1.0, 1e-5);
  EXPECT_NEAR(p2.z(), 7.0, 1e-5);
}

TEST(CellReaderTests, ThrowsOnReadTrajectory) {
  CellReader reader;
  EXPECT_THROW(reader.readTrajectory("dummy.cell"), std::runtime_error);
}
