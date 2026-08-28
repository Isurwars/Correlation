// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "app/UpdateChecker.hpp"

#include <gtest/gtest.h>
#include <string>

namespace {

TEST(UpdateCheckerTests, DetectsNewerVersionsCorrectly) {
  using correlation::app::UpdateChecker;

  // Newer patch version
  EXPECT_TRUE(UpdateChecker::isNewerVersion("v3.9.1", "3.9.0"));
  EXPECT_TRUE(UpdateChecker::isNewerVersion("3.9.1", "3.9.0"));

  // Newer minor version
  EXPECT_TRUE(UpdateChecker::isNewerVersion("v3.10.0", "3.9.0"));
  EXPECT_TRUE(UpdateChecker::isNewerVersion("3.10.0", "3.9.0"));

  // Newer major version
  EXPECT_TRUE(UpdateChecker::isNewerVersion("v4.0.0", "3.9.0"));
  EXPECT_TRUE(UpdateChecker::isNewerVersion("4.0.0", "3.9.0"));

  // Equal versions
  EXPECT_FALSE(UpdateChecker::isNewerVersion("v3.9.0", "3.9.0"));
  EXPECT_FALSE(UpdateChecker::isNewerVersion("3.9.0", "3.9.0"));
  EXPECT_FALSE(UpdateChecker::isNewerVersion("3.9.0", "v3.9.0"));

  // Older versions
  EXPECT_FALSE(UpdateChecker::isNewerVersion("v3.8.9", "3.9.0"));
  EXPECT_FALSE(UpdateChecker::isNewerVersion("v3.8.0", "3.9.0"));
  EXPECT_FALSE(UpdateChecker::isNewerVersion("v2.9.0", "3.9.0"));
}

TEST(UpdateCheckerTests, ParsesReleaseJsonCorrectly) {
  using correlation::app::UpdateChecker;

  const std::string mock_json = R"({
    "tag_name": "v3.9.1",
    "html_url": "https://github.com/Isurwars/Correlation/releases/tag/v3.9.1",
    "name": "Release 3.9.1"
  })";

  const auto release = UpdateChecker::parseReleaseJson(mock_json);
  ASSERT_TRUE(release.has_value());
  EXPECT_EQ(release->tag_name, "v3.9.1");
  EXPECT_EQ(release->html_url, "https://github.com/Isurwars/Correlation/releases/tag/v3.9.1");
}

TEST(UpdateCheckerTests, HandlesMissingFieldsInReleaseJson) {
  using correlation::app::UpdateChecker;

  const std::string invalid_json = R"({
    "message": "Not Found",
    "documentation_url": "https://docs.github.com/rest"
  })";

  const auto release = UpdateChecker::parseReleaseJson(invalid_json);
  EXPECT_FALSE(release.has_value());
}

TEST(UpdateCheckerTests, OpenUrlEmptyStringSafelyNoops) {
  using correlation::app::UpdateChecker;
  // Should not throw or crash on empty URL
  EXPECT_NO_THROW(UpdateChecker::openUrlInBrowser(""));
}

} // namespace
