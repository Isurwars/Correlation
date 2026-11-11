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

#include <filesystem>
#include <gtest/gtest.h>
#include <span>

namespace {

using correlation::real_t;
using namespace correlation::writers::structure;
using namespace correlation::readers;
using namespace correlation::core;

[[nodiscard]] std::string findDataPath(const std::string &relative_path) {
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

void verifyLatticeParametersNear(const std::array<real_t, 6> &p_init,
                                 const std::array<real_t, 6> &p_parse, real_t tolerance) {
  std::span<const real_t, 6> const span_init(p_init);
  std::span<const real_t, 6> const span_parse(p_parse);
  for (size_t idx = 0; idx < span_init.size(); ++idx) {
    EXPECT_NEAR(span_init[idx], span_parse[idx], tolerance);
  }
}

class StructureWriterTests : public ::testing::Test {
protected:
  void SetUp() override {
    temp_dir_ = std::filesystem::temp_directory_path() / "correlation_writer_tests";
    std::filesystem::create_directories(temp_dir_);
  }

  void TearDown() override {
    std::error_code err_code;
    std::filesystem::remove_all(temp_dir_, err_code);
  }

  [[nodiscard]] const std::filesystem::path &tempDir() const noexcept { return temp_dir_; }

private:
  std::filesystem::path temp_dir_;
};

TEST_F(StructureWriterTests, FactoryRegistrationAndDiscovery) {
  const auto &factory = StructureWriterFactory::instance();

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

  const Cell initial_cell = CarReader::read(source_car);
  ASSERT_EQ(initial_cell.atomCount(), 8);

  std::filesystem::path const out_poscar = tempDir() / "silicon.vasp";
  const PoscarWriter writer;
  writer.write(out_poscar.string(), initial_cell, {.fractional = true, .title = "Si Roundtrip"});
  ASSERT_TRUE(std::filesystem::exists(out_poscar));

  const Cell parsed_cell = VaspReader::read(out_poscar.string());
  EXPECT_EQ(parsed_cell.atomCount(), initial_cell.atomCount());
  EXPECT_NEAR(parsed_cell.volume(), initial_cell.volume(), static_cast<real_t>(1e-2));

  verifyLatticeParametersNear(initial_cell.latticeParameters(), parsed_cell.latticeParameters(),
                              static_cast<real_t>(1e-2));
}

TEST_F(StructureWriterTests, CellRoundTripFromPoscar) {
  std::string const source_poscar = findDataPath("vasp/Si.poscar");
  ASSERT_TRUE(std::filesystem::exists(source_poscar));

  const Cell initial_cell = VaspReader::read(source_poscar);
  ASSERT_EQ(initial_cell.atomCount(), 8);

  std::filesystem::path const out_cell = tempDir() / "silicon.cell";
  const CellWriter writer;
  writer.write(out_cell.string(), initial_cell, {.fractional = true});
  ASSERT_TRUE(std::filesystem::exists(out_cell));

  const Cell parsed_cell = CellReader::read(out_cell.string());
  EXPECT_EQ(parsed_cell.atomCount(), 8);
  EXPECT_NEAR(parsed_cell.volume(), initial_cell.volume(), static_cast<real_t>(1e-2));
}

TEST_F(StructureWriterTests, CifRoundTrip) {
  std::string const source_poscar = findDataPath("vasp/Si.poscar");
  ASSERT_TRUE(std::filesystem::exists(source_poscar));

  const Cell initial_cell = VaspReader::read(source_poscar);
  std::filesystem::path const out_cif = tempDir() / "silicon.cif";

  const CifWriter writer;
  writer.write(out_cif.string(), initial_cell);
  ASSERT_TRUE(std::filesystem::exists(out_cif));

  const Cell parsed_cell = CifReader::read(out_cif.string());
  EXPECT_EQ(parsed_cell.atomCount(), 8);
  EXPECT_NEAR(parsed_cell.volume(), initial_cell.volume(), static_cast<real_t>(1e-2));
}

TEST_F(StructureWriterTests, CarRoundTrip) {
  std::string const source_car = findDataPath("car/si_crystal.car");
  ASSERT_TRUE(std::filesystem::exists(source_car));

  const Cell initial_cell = CarReader::read(source_car);
  std::filesystem::path const out_car = tempDir() / "silicon_out.car";

  const CarWriter writer;
  writer.write(out_car.string(), initial_cell);
  ASSERT_TRUE(std::filesystem::exists(out_car));

  const Cell parsed_cell = CarReader::read(out_car.string());
  EXPECT_EQ(parsed_cell.atomCount(), 8);
  EXPECT_NEAR(parsed_cell.volume(), initial_cell.volume(), static_cast<real_t>(1e-2));
}

TEST_F(StructureWriterTests, DegenerateCellThrowsException) {
  const Cell empty_cell;
  std::filesystem::path const dummy_path = tempDir() / "fail.vasp";
  const PoscarWriter writer;
  EXPECT_THROW(writer.write(dummy_path.string(), empty_cell), std::runtime_error);
}

} // namespace
