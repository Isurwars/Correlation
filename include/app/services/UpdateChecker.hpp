/**
 * @file UpdateChecker.hpp
 * @brief Asynchronous GitHub release update checker for the application.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include <optional>
#include <string>

class AppWindow;

namespace correlation::app {

/**
 * @struct ReleaseInfo
 * @brief Holds metadata for a GitHub release.
 */
struct ReleaseInfo {
  std::string tag_name; ///< The release tag (e.g. "v3.9.1").
  std::string html_url; ///< Web page URL for the release download.
};

/**
 * @class UpdateChecker
 * @brief Checks for new releases asynchronously on GitHub without interrupting the user.
 */
class UpdateChecker {
public:
  /**
   * @brief Spawns an asynchronous background thread to check for new releases on GitHub.
   *        If a newer version is found, safely dispatches UI updates to AppWindow.
   * @param window The Slint AppWindow instance.
   * @param current_version The current application version string (e.g. "3.9.0").
   */
  static void checkForUpdatesAsync(AppWindow &window, const std::string &current_version);

  /**
   * @brief Compares two semantic version strings (e.g. "v3.9.1" vs "3.9.0").
   * @param remote_ver Remote version string.
   * @param current_ver Current version string.
   * @return true if remote_ver is strictly newer than current_ver.
   */
  static bool isNewerVersion(const std::string &remote_ver, const std::string &current_ver);

  /**
   * @brief Opens a URL in the user's default system browser.
   * @param url The URL string to open.
   */
  static void openUrlInBrowser(const std::string &url);

  /**
   * @brief Queries the GitHub API for the latest release metadata.
   * @return ReleaseInfo with tag_name and html_url if successful, std::nullopt otherwise.
   */
  static std::optional<ReleaseInfo> fetchLatestRelease();

  /**
   * @brief Parses the tag_name and html_url fields from GitHub release JSON string.
   * @param json The JSON response string.
   * @return ReleaseInfo if tag_name is found, std::nullopt otherwise.
   */
  static std::optional<ReleaseInfo> parseReleaseJson(const std::string &json);
};

} // namespace correlation::app
