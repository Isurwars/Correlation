#include "readers/CP2KReader.hpp"

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
    if (std::filesystem::exists(dir + "cp2k/clean.restart")) {
      return dir + "cp2k/";
    }
  }
  return "../../tests/data/cp2k/";
}

class CP2KReaderTests : public ::testing::Test {
public:
  std::string data_dir_;

protected:
  void SetUp() override { data_dir_ = getTestDataDir(); }
};
} // namespace

TEST_F(CP2KReaderTests, ReadsSingleFrame) {
  CP2KReader reader;
  auto cell = reader.readStructure(data_dir_ + "clean.restart");
  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "H");
  EXPECT_EQ(cell.lattice_parameters()[0], 12.0);
}
