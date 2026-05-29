#include "readers/QEReader.hpp"
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>

using namespace correlation::readers;

class QEReaderTests : public ::testing::Test {
protected:
  void SetUp() override {
    std::ofstream out("test_qe.pwo");
    out << "CELL_PARAMETERS (angstrom)\n";
    out << "10.0 0.0 0.0\n";
    out << "0.0 10.0 0.0\n";
    out << "0.0 0.0 10.0\n";
    out << "ATOMIC_POSITIONS (angstrom)\n";
    out << "C 1.0 2.0 3.0\n";
    out << "O 4.0 5.0 6.0\n";
    out.close();
  }
  void TearDown() override { std::remove("test_qe.pwo"); }
};

TEST_F(QEReaderTests, ReadsSingleFrame) {
  QEReader reader;
  auto cell = reader.readStructure("test_qe.pwo");
  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "C");
  EXPECT_EQ(cell.lattice_parameters()[0], 10.0);
}

TEST_F(QEReaderTests, ReadsNonOrthogonalLattice) {
  std::ofstream out("test_qe_triclinic.pwo");
  out << "CELL_PARAMETERS (angstrom)\n";
  out << "10.0 0.0 0.0\n";
  out << "2.0 11.0 0.0\n";
  out << "1.0 3.0 12.0\n";
  out << "ATOMIC_POSITIONS (angstrom)\n";
  out << "C 1.0 2.0 3.0\n";
  out << "O 4.0 5.0 6.0\n";
  out.close();

  QEReader reader;
  auto cell = reader.readStructure("test_qe_triclinic.pwo");
  std::remove("test_qe_triclinic.pwo");

  EXPECT_EQ(cell.atomCount(), 2);

  // Check lattice vectors
  const auto &vectors = cell.latticeVectors();
  EXPECT_DOUBLE_EQ(vectors[0][0], 10.0);
  EXPECT_DOUBLE_EQ(vectors[0][1], 0.0);
  EXPECT_DOUBLE_EQ(vectors[0][2], 0.0);

  EXPECT_DOUBLE_EQ(vectors[1][0], 2.0);
  EXPECT_DOUBLE_EQ(vectors[1][1], 11.0);
  EXPECT_DOUBLE_EQ(vectors[1][2], 0.0);

  EXPECT_DOUBLE_EQ(vectors[2][0], 1.0);
  EXPECT_DOUBLE_EQ(vectors[2][1], 3.0);
  EXPECT_DOUBLE_EQ(vectors[2][2], 12.0);
}
