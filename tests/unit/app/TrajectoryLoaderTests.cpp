// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only

#include "app/services/TrajectoryLoader.hpp"
#include "readers/BaseReader.hpp"
#include "readers/contracts/IReaderRegistry.hpp"

#include <filesystem>
#include <gtest/gtest.h>
#include <string>
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
    if (std::filesystem::exists(dir + "xyz/clean.xyz")) {
      return dir;
    }
  }
  return "../../tests/data/";
}

TEST(TrajectoryLoaderTests, InitialStateIsEmpty) {
  const correlation::app::TrajectoryLoader loader;

  EXPECT_EQ(loader.trajectory(), nullptr);
  EXPECT_EQ(loader.cell(), nullptr);
  EXPECT_EQ(loader.getFrameCount(), 0);
  EXPECT_EQ(loader.getTotalAtomCount(), 0);
  EXPECT_EQ(loader.getRemovedFrameCount(), 0);
  EXPECT_DOUBLE_EQ(loader.getTimeStep(), 1.0);
  EXPECT_TRUE(loader.getAtomCounts().empty());
}

TEST(TrajectoryLoaderTests, LoadsValidStructureFile) {
  correlation::app::TrajectoryLoader loader;
  std::string const filepath = getTestDataDir() + "xyz/clean.xyz";

  const auto res = loader.loadFile(filepath);
  ASSERT_TRUE(res.has_value());
  EXPECT_NE(loader.trajectory(), nullptr);
  EXPECT_NE(loader.cell(), nullptr);
  EXPECT_GT(loader.getFrameCount(), 0);
  EXPECT_GT(loader.getTotalAtomCount(), 0);

  const auto counts = loader.getAtomCounts();
  EXPECT_FALSE(counts.empty());
  EXPECT_GT(loader.getRecommendedTimeStep(), 0.0);
}

TEST(TrajectoryLoaderTests, RejectsNonExistentFile) {
  correlation::app::TrajectoryLoader loader;
  const auto res = loader.loadFile("/non/existent/path/file.xyz");
  EXPECT_FALSE(res.has_value());
  EXPECT_EQ(loader.trajectory(), nullptr);
}

TEST(TrajectoryLoaderTests, ClearResetsState) {
  correlation::app::TrajectoryLoader loader;
  std::string const filepath = getTestDataDir() + "xyz/clean.xyz";

  const auto res = loader.loadFile(filepath);
  ASSERT_TRUE(res.has_value());
  EXPECT_NE(loader.cell(), nullptr);

  loader.clear();
  EXPECT_EQ(loader.cell(), nullptr);
  EXPECT_EQ(loader.getFrameCount(), 0);
}

class MockCustomReader : public correlation::readers::BaseReader {
public:
  [[nodiscard]] std::string getName() const override { return "MockCustomReader"; }
  [[nodiscard]] std::vector<std::string> getExtensions() const override { return {".mockext"}; }
  [[nodiscard]] bool isTrajectory() const override { return false; }
  correlation::core::Cell
  readStructure(const std::string & /*filename*/,
                correlation::readers::ProgressCallback /*progress_callback*/ = nullptr) override {
    correlation::core::Cell cell;
    cell.addAtom("Si", {0.0, 0.0, 0.0});
    return cell;
  }
  correlation::core::Trajectory
  readTrajectory(const std::string & /*filename*/,
                 correlation::readers::ProgressCallback /*progress_callback*/ = nullptr) override {
    return {};
  }
};

class MockReaderRegistry : public correlation::readers::IReaderRegistry {
public:
  MockReaderRegistry() : reader_(std::make_unique<MockCustomReader>()) {}
  [[nodiscard]] correlation::readers::BaseReader *
  getReaderForExtension(const correlation::readers::ReaderExtensionQuery &query) const override {
    if (query.extension == ".mockext") {
      return reader_.get();
    }
    return nullptr;
  }
  [[nodiscard]] correlation::readers::BaseReader *
  getReaderForExtension(std::string_view extension) const override {
    if (extension == ".mockext") {
      return reader_.get();
    }
    return nullptr;
  }
  [[nodiscard]] std::vector<std::string> getAllExtensions() const override { return {".mockext"}; }
  [[nodiscard]] const std::vector<std::unique_ptr<correlation::readers::BaseReader>> &
  getReaders() const override {
    static const std::vector<std::unique_ptr<correlation::readers::BaseReader>> EMPTY;
    return EMPTY;
  }

private:
  std::unique_ptr<MockCustomReader> reader_;
};

TEST(TrajectoryLoaderTests, InjectedReaderRegistryIsUsed) {
  const MockReaderRegistry mock_registry;
  correlation::app::TrajectoryLoader loader(mock_registry);

  const auto res = loader.loadFile("sample.mockext");
  ASSERT_TRUE(res.has_value());
  ASSERT_NE(loader.cell(), nullptr);
  EXPECT_EQ(loader.getTotalAtomCount(), 1);
  EXPECT_EQ(loader.cell()->atoms()[0].element().symbol, "Si");
}

std::filesystem::path findPlasticExamplesDir() {
  for (const auto &candidate :
       {"examples/plastic", "../examples/plastic", "../../examples/plastic"}) {
    if (std::filesystem::exists(std::filesystem::path(candidate) / "dump.compresion.0")) {
      return candidate;
    }
  }
  return {};
}

void verifyLoadedPlasticCellMetrics(const correlation::core::Trajectory &traj) {
  EXPECT_NEAR(traj.getLatticeParameters(0)[0], 70.6508, 1e-3);
  EXPECT_NEAR(traj.getLatticeParameters(1)[0], 69.9443, 1e-3);
  EXPECT_NEAR(traj.getVolume(0), 352656.30, 1.0);
}

TEST(TrajectoryLoaderTests, LoadsPlasticCompressionMultiFileSeriesInOrder) {
  const std::filesystem::path plastic_dir = findPlasticExamplesDir();
  EXPECT_FALSE(plastic_dir.empty());

  // Pass a subset in reverse/unsorted order to verify natural sort in TrajectoryLoader::loadFiles
  std::vector<std::string> const paths = {
      (plastic_dir / "dump.compresion.100000").string(),
      (plastic_dir / "dump.compresion.0").string(),
      (plastic_dir / "dump.compresion.20000").string(),
      (plastic_dir / "dump.compresion.10000").string(),
  };

  correlation::app::TrajectoryLoader loader;
  const auto res = loader.loadFiles(paths);
  EXPECT_TRUE(res.has_value());
  EXPECT_EQ(loader.getFrameCount(), 4U);
  EXPECT_EQ(loader.getTotalAtomCount(), 16000U);

  if (const auto *traj = loader.trajectory(); traj != nullptr) {
    verifyLoadedPlasticCellMetrics(*traj);
  }
}

TEST(TrajectoryLoaderTests, RenamesPlasticLammpsTypesAndPreservesOrder) {
  const std::filesystem::path plastic_dir = findPlasticExamplesDir();
  EXPECT_FALSE(plastic_dir.empty());

  std::vector<std::string> const paths = {
      (plastic_dir / "dump.compresion.0").string(),
      (plastic_dir / "dump.compresion.10000").string(),
  };

  correlation::app::TrajectoryLoader loader;
  const auto res = loader.loadFiles(paths);
  EXPECT_TRUE(res.has_value());

  const auto initial_counts = loader.getOrderedAtomCounts();
  EXPECT_EQ(initial_counts.size(), 4U);
  EXPECT_EQ(initial_counts[0].first, "1");
  EXPECT_EQ(initial_counts[0].second, 7895);
  EXPECT_EQ(initial_counts[1].first, "3");
  EXPECT_EQ(initial_counts[2].first, "4");
  EXPECT_EQ(initial_counts[3].first, "2");

  EXPECT_TRUE(loader.renameElement("1", "Cu"));
  EXPECT_TRUE(loader.renameElement("2", "Zr"));
  EXPECT_TRUE(loader.renameElement("3", "Al"));
  EXPECT_TRUE(loader.renameElement("4", "Ti"));

  // Slot order in elements() is preserved in-place so UI rows do not jump when renamed
  const auto updated_counts = loader.getOrderedAtomCounts();
  EXPECT_EQ(updated_counts.size(), 4U);
  EXPECT_EQ(updated_counts[0].first, "Cu");
  EXPECT_EQ(updated_counts[1].first, "Al");
  EXPECT_EQ(updated_counts[2].first, "Ti");
  EXPECT_EQ(updated_counts[3].first, "Zr");
}

} // namespace
