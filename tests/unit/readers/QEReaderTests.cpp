/**
 * @file QEReaderTests.cpp
 * @brief Unit tests for QEReader Quantum ESPRESSO output parsing.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "../../PrecisionTestHelper.hpp"
#include "readers/QEReader.hpp"

#include <filesystem>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace correlation::readers {

using correlation::real_t;
using correlation::testing::IsRealEq;

namespace {

[[nodiscard]] std::string getTestDataDir() {
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
protected:
  [[nodiscard]] const std::string &dataDir() const noexcept { return data_dir_; }

  void SetUp() override { data_dir_ = getTestDataDir(); }

private:
  std::string data_dir_;
};

} // namespace

TEST_F(QEReaderTests, ReadsSingleFrame) {
  QEReader reader;
  const auto cell = reader.readStructure(dataDir() + "clean.pwo");
  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "C");
  EXPECT_THAT(cell.latticeParameters()[0], IsRealEq(static_cast<real_t>(10.0)));
}

TEST_F(QEReaderTests, ReadsNonOrthogonalLattice) {
  QEReader reader;
  const auto cell = reader.readStructure(dataDir() + "clean_triclinic.pwo");

  EXPECT_EQ(cell.atomCount(), 2);

  // Check lattice vectors
  const auto &vectors = cell.latticeVectors();
  EXPECT_THAT(vectors[0][0], IsRealEq(static_cast<real_t>(10.0)));
  EXPECT_THAT(vectors[0][1], IsRealEq(static_cast<real_t>(0.0)));
  EXPECT_THAT(vectors[0][2], IsRealEq(static_cast<real_t>(0.0)));

  EXPECT_THAT(vectors[1][0], IsRealEq(static_cast<real_t>(2.0)));
  EXPECT_THAT(vectors[1][1], IsRealEq(static_cast<real_t>(11.0)));
  EXPECT_THAT(vectors[1][2], IsRealEq(static_cast<real_t>(0.0)));

  EXPECT_THAT(vectors[2][0], IsRealEq(static_cast<real_t>(1.0)));
  EXPECT_THAT(vectors[2][1], IsRealEq(static_cast<real_t>(3.0)));
  EXPECT_THAT(vectors[2][2], IsRealEq(static_cast<real_t>(12.0)));
}

} // namespace correlation::readers
