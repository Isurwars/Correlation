#include "readers/QEReader.hpp"

#include <gtest/gtest.h>

using namespace correlation::readers;

#include <filesystem>
#include <vector>

namespace {

std::string getTestDataDir() {
  std::vector<std::string> const candidates = {
      "../../tests/data/",
      "../tests/data/",
      "tests/data/",
      "data/",
  };
  for (const auto &dir : candidates) {
    if (std::filesystem::exists(dir + "qe/clean.pwo")) {
      return dir + "qe/";
    }
  }
  return "../../tests/data/qe/";
}

class QEReaderTests : public ::testing::Test {
public:
  std::string data_dir_;

protected:
  void SetUp() override { data_dir_ = getTestDataDir(); }
};
} // namespace

TEST_F(QEReaderTests, ReadsSingleFrame) {
  QEReader reader;
  auto cell = reader.readStructure(data_dir_ + "clean.pwo");
  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "C");
  EXPECT_EQ(cell.lattice_parameters()[0], 10.0);
}

TEST_F(QEReaderTests, ReadsNonOrthogonalLattice) {
  QEReader reader;
  auto cell = reader.readStructure(data_dir_ + "clean_triclinic.pwo");

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
