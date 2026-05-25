// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "readers/ReaderFactory.hpp"
#include "readers/BaseReader.hpp"

#include <gtest/gtest.h>
#include <memory>

namespace correlation::testing {

using namespace correlation::readers;

class MockReader : public BaseReader {
public:
  std::string getName() const override { return "MockReader"; }
  std::vector<std::string> getExtensions() const override { return {".mock", ".mck"}; }
  bool isTrajectory() const override { return false; }
  
  correlation::core::Cell readStructure(
      const std::string &filename,
      std::function<void(float, const std::string &)> progress_callback = nullptr) override {
    return correlation::core::Cell();
  }

  correlation::core::Trajectory readTrajectory(
      const std::string &filename,
      std::function<void(float, const std::string &)> progress_callback = nullptr) override {
    return correlation::core::Trajectory();
  }
};

TEST(ReaderFactoryTests, SingletonInstanceIsUnique) {
  auto &factory1 = ReaderFactory::instance();
  auto &factory2 = ReaderFactory::instance();
  EXPECT_EQ(&factory1, &factory2);
}

TEST(ReaderFactoryTests, GetRegisteredReaders) {
  auto &factory = ReaderFactory::instance();
  const auto &readers = factory.getReaders();
  EXPECT_GT(readers.size(), 0);
}

TEST(ReaderFactoryTests, GetAllExtensions) {
  auto &factory = ReaderFactory::instance();
  std::vector<std::string> extensions = factory.getAllExtensions();
  EXPECT_GT(extensions.size(), 0);
  
  // At least standard formats like .car or .cif should exist
  bool found_car = false;
  for (const auto &ext : extensions) {
    if (ext == ".car") {
      found_car = true;
      break;
    }
  }
  EXPECT_TRUE(found_car);
}

TEST(ReaderFactoryTests, GetReaderForExtension) {
  auto &factory = ReaderFactory::instance();
  
  // Try querying a standard reader by extension
  BaseReader *retrieved = factory.getReaderForExtension(".car");
  ASSERT_NE(retrieved, nullptr);
  EXPECT_FALSE(retrieved->getName().empty());
}

TEST(ReaderFactoryTests, RegisterAndLookupCustomReader) {
  auto &factory = ReaderFactory::instance();
  
  // Register a custom mock reader
  auto mock = std::make_unique<MockReader>();
  bool result = factory.registerReader(std::move(mock));
  EXPECT_TRUE(result);
  
  // Verify it resolves by extension
  BaseReader *retrieved1 = factory.getReaderForExtension(".mock");
  ASSERT_NE(retrieved1, nullptr);
  EXPECT_EQ(retrieved1->getName(), "MockReader");

  BaseReader *retrieved2 = factory.getReaderForExtension(".mck");
  ASSERT_NE(retrieved2, nullptr);
  EXPECT_EQ(retrieved2->getName(), "MockReader");
}

TEST(ReaderFactoryTests, LookupNonExistentReturnsNullptr) {
  auto &factory = ReaderFactory::instance();
  BaseReader *retrieved = factory.getReaderForExtension(".non_existent_extension");
  EXPECT_EQ(retrieved, nullptr);
}

TEST(ReaderFactoryTests, LookupEmptyExtensionReturnsNullptr) {
  auto &factory = ReaderFactory::instance();
  BaseReader *retrieved = factory.getReaderForExtension("");
  EXPECT_EQ(retrieved, nullptr);
}

} // namespace correlation::testing
