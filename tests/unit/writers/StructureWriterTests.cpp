/**
 * @file StructureWriterTests.cpp
 * @brief Unit and integration round-trip tests for atomic structure writers.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "core/Cell.hpp"
#include "readers/CarReader.hpp"
#include "readers/CellReader.hpp"
#include "readers/CifReader.hpp"
#include "readers/VaspReader.hpp"
#include "writers/structure/CarWriter.hpp"
#include "writers/structure/CellWriter.hpp"
#include "writers/structure/CifWriter.hpp"
#include "writers/structure/PoscarWriter.hpp"
#include "writers/structure/StructureWriterFactory.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>

namespace {

using namespace correlation::writers::structure;
using namespace correlation::readers;
using namespace correlation::core;

std::string findDataPath(const std::string &relative_path) {
  const std::vector<std::string> prefixes = {
      "../../tests/data/",
      "../tests/data/",
      "tests/data/",
  };
  for (const auto &prefix : prefixes) {
    std::string candidate = prefix + relative_path;
    if (std::filesystem::exists(candidate)) {
      return candidate;
    }
  }
  return "tests/data/" + relative_path;
}

class StructureWriterTests : public ::testing::Test {
protected:
  void SetUp() override {
    temp_dir_ = std::filesystem::temp_directory_path() / "correlation_writer_tests";
    std::filesystem::create_directories(temp_dir_);
  }

  void TearDown() override {
    std::error_code ec;
    std::filesystem::remove_all(temp_dir_, ec);
  }

  std::filesystem::path temp_dir_;
};

TEST_F(StructureWriterTests, FactoryRegistrationAndDiscovery) {
  auto &factory = StructureWriterFactory::instance();

  EXPECT_NE(factory.getWriterForExtension(".vasp"), nullptr);
  EXPECT_NE(factory.getWriterForExtension("poscar"), nullptr);
  EXPECT_NE(factory.getWriterForExtension(".cell"), nullptr);
  EXPECT_NE(factory.getWriterForExtension("cif"), nullptr);
  EXPECT_NE(factory.getWriterForExtension(".car"), nullptr);

  EXPECT_NE(factory.getWriterForFormat("poscar"), nullptr);
  EXPECT_NE(factory.getWriterForFormat("cell"), nullptr);
  EXPECT_NE(factory.getWriterForFormat("cif"), nullptr);
  EXPECT_NE(factory.getWriterForFormat("car"), nullptr);

  EXPECT_EQ(factory.getWriterForExtension("unknown_ext"), nullptr);
  EXPECT_EQ(factory.getWriterForFormat("non_existent_fmt"), nullptr);
}

TEST_F(StructureWriterTests, PoscarRoundTripFromCar) {
  std::string const source_car = findDataPath("car/si_crystal.car");
  ASSERT_TRUE(std::filesystem::exists(source_car));

  Cell initial_cell = CarReader::read(source_car);
  ASSERT_EQ(initial_cell.atomCount(), 8);

  std::filesystem::path const out_poscar = temp_dir_ / "silicon.vasp";
  PoscarWriter writer;
  writer.write(out_poscar.string(), initial_cell, {.fractional = true, .title = "Si Roundtrip"});
  ASSERT_TRUE(std::filesystem::exists(out_poscar));

  Cell parsed_cell = VaspReader::read(out_poscar.string());
  EXPECT_EQ(parsed_cell.atomCount(), initial_cell.atomCount());
  EXPECT_NEAR(parsed_cell.volume(), initial_cell.volume(), 1e-2);

  const auto &p_init = initial_cell.latticeParameters();
  const auto &p_parse = parsed_cell.latticeParameters();
  for (std::size_t i = 0; i < 6; ++i) {
    EXPECT_NEAR(p_init[i], p_parse[i], 1e-2);
  }
}

TEST_F(StructureWriterTests, CellRoundTripFromPoscar) {
  std::string const source_poscar = findDataPath("vasp/Si.poscar");
  ASSERT_TRUE(std::filesystem::exists(source_poscar));

  Cell initial_cell = VaspReader::read(source_poscar);
  ASSERT_EQ(initial_cell.atomCount(), 8);

  std::filesystem::path const out_cell = temp_dir_ / "silicon.cell";
  CellWriter writer;
  writer.write(out_cell.string(), initial_cell, {.fractional = true});
  ASSERT_TRUE(std::filesystem::exists(out_cell));

  Cell parsed_cell = CellReader::read(out_cell.string());
  EXPECT_EQ(parsed_cell.atomCount(), 8);
  EXPECT_NEAR(parsed_cell.volume(), initial_cell.volume(), 1e-2);
}

TEST_F(StructureWriterTests, CifRoundTrip) {
  std::string const source_poscar = findDataPath("vasp/Si.poscar");
  ASSERT_TRUE(std::filesystem::exists(source_poscar));

  Cell initial_cell = VaspReader::read(source_poscar);
  std::filesystem::path const out_cif = temp_dir_ / "silicon.cif";

  CifWriter writer;
  writer.write(out_cif.string(), initial_cell);
  ASSERT_TRUE(std::filesystem::exists(out_cif));

  Cell parsed_cell = CifReader::read(out_cif.string());
  EXPECT_EQ(parsed_cell.atomCount(), 8);
  EXPECT_NEAR(parsed_cell.volume(), initial_cell.volume(), 1e-2);
}

TEST_F(StructureWriterTests, CarRoundTrip) {
  std::string const source_car = findDataPath("car/si_crystal.car");
  ASSERT_TRUE(std::filesystem::exists(source_car));

  Cell initial_cell = CarReader::read(source_car);
  std::filesystem::path const out_car = temp_dir_ / "silicon_out.car";

  CarWriter writer;
  writer.write(out_car.string(), initial_cell);
  ASSERT_TRUE(std::filesystem::exists(out_car));

  Cell parsed_cell = CarReader::read(out_car.string());
  EXPECT_EQ(parsed_cell.atomCount(), 8);
  EXPECT_NEAR(parsed_cell.volume(), initial_cell.volume(), 1e-2);
}

TEST_F(StructureWriterTests, DegenerateCellThrowsException) {
  Cell empty_cell;
  std::filesystem::path const dummy_path = temp_dir_ / "fail.vasp";
  PoscarWriter writer;
  EXPECT_THROW(writer.write(dummy_path.string(), empty_cell), std::runtime_error);
}

} // namespace
