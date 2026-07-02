// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "writers/BaseWriter.hpp"
#include "writers/WriterFactory.hpp"

#include <gtest/gtest.h>
#include <memory>

namespace correlation::testing {

using namespace correlation::writers;

namespace {
class MockWriter : public BaseWriter {
public:
  [[nodiscard]] std::string getName() const override { return "MockWriter"; }
  [[nodiscard]] std::vector<std::string> getExtensions() const override { return {".mockw", ".mkw"}; }

  void write(const std::string &base_path, const correlation::analysis::DistributionFunctions &dists,
             bool smoothing) const override {
    // No-op for mock
  }
};
} // namespace

TEST(WriterFactoryTests, SingletonInstanceIsUnique) {
  auto &factory1 = WriterFactory::instance();
  auto &factory2 = WriterFactory::instance();
  EXPECT_EQ(&factory1, &factory2);
}

TEST(WriterFactoryTests, GetRegisteredWriters) {
  auto &factory = WriterFactory::instance();
  const auto &writers = factory.getWriters();
  EXPECT_GT(writers.size(), 0);
}

TEST(WriterFactoryTests, GetWriterForExtension) {
  auto &factory = WriterFactory::instance();

  // Try querying standard format by extension
  BaseWriter const *retrieved = factory.getWriterForExtension(".csv");
  ASSERT_NE(retrieved, nullptr);
  EXPECT_FALSE(retrieved->getName().empty());
}

TEST(WriterFactoryTests, GetWriterByName) {
  auto &factory = WriterFactory::instance();

  // Try querying standard format by name
  // Standard writer names could be e.g., "CSV", "HDF5", "Arrow/Parquet"
  const auto &writers = factory.getWriters();
  ASSERT_FALSE(writers.empty());
  std::string const first_writer_name = writers[0]->getName();

  BaseWriter const *retrieved = factory.getWriter(first_writer_name);
  ASSERT_NE(retrieved, nullptr);
  EXPECT_EQ(retrieved->getName(), first_writer_name);
}

TEST(WriterFactoryTests, RegisterAndLookupCustomWriter) {
  auto &factory = WriterFactory::instance();

  // Register a custom mock writer
  auto mock = std::make_unique<MockWriter>();
  bool const result = factory.registerWriter(std::move(mock));
  EXPECT_TRUE(result);

  // Verify it resolves by extension
  BaseWriter const *retrieved_ext = factory.getWriterForExtension(".mockw");
  ASSERT_NE(retrieved_ext, nullptr);
  EXPECT_EQ(retrieved_ext->getName(), "MockWriter");

  // Verify it resolves by name
  BaseWriter const *retrieved_name = factory.getWriter("MockWriter");
  ASSERT_NE(retrieved_name, nullptr);
  EXPECT_EQ(retrieved_name->getName(), "MockWriter");
}

TEST(WriterFactoryTests, LookupNonExistentReturnsNullptr) {
  auto &factory = WriterFactory::instance();
  BaseWriter const *retrieved_ext = factory.getWriterForExtension(".non_existent_writer_ext");
  EXPECT_EQ(retrieved_ext, nullptr);

  BaseWriter const *retrieved_name = factory.getWriter("NonExistentWriterName");
  EXPECT_EQ(retrieved_name, nullptr);
}

TEST(WriterFactoryTests, LookupEmptyExtensionReturnsNullptr) {
  auto &factory = WriterFactory::instance();
  BaseWriter const *retrieved = factory.getWriterForExtension("");
  EXPECT_EQ(retrieved, nullptr);
}

} // namespace correlation::testing
