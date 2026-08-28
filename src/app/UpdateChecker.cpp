/**
 * @file UpdateChecker.cpp
 * @brief Implementation of the asynchronous GitHub release update checker.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/UpdateChecker.hpp"

#if __has_include("AppWindow.h")
#include "AppWindow.h"
#endif

#include <array>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <tuple>
#include <vector>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
// clang-format off
#include <windows.h>
#include <shellapi.h>
#include <wininet.h>
// clang-format on

#ifdef _MSC_VER
#pragma comment(lib, "wininet.lib")
#pragma comment(lib, "shell32.lib")
#endif
#else
#include <fcntl.h>
#include <spawn.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#ifdef __APPLE__
#include <crt_externs.h>
#define environ (*_NSGetEnviron())
#else
extern char **environ;
#endif
#endif

namespace correlation::app {

namespace {

std::vector<int> parseSemver(std::string_view ver) {
  std::vector<int> parts;
  size_t idx = 0;
  // Skip leading non-digit characters (e.g. 'v', 'V', whitespace)
  while (idx < ver.size() && (std::isdigit(static_cast<unsigned char>(ver[idx])) == 0)) {
    ++idx;
  }

  while (idx < ver.size()) {
    if (std::isdigit(static_cast<unsigned char>(ver[idx])) != 0) {
      int val = 0;
      while (idx < ver.size() && (std::isdigit(static_cast<unsigned char>(ver[idx])) != 0)) {
        val = val * 10 + (ver[idx] - '0');
        ++idx;
      }
      parts.push_back(val);
    } else if (ver[idx] == '.') {
      ++idx;
    } else {
      // Stop at pre-release hyphen or other non-numeric separator
      break;
    }
  }

  while (parts.size() < 3) {
    parts.push_back(0);
  }
  return parts;
}

std::optional<std::string> extractJsonStringField(const std::string &json, std::string_view key) {
  const size_t key_pos = json.find(key);
  if (key_pos == std::string::npos) {
    return std::nullopt;
  }

  const size_t colon_pos = json.find(':', key_pos + key.size());
  if (colon_pos == std::string::npos) {
    return std::nullopt;
  }

  const size_t open_quote = json.find('"', colon_pos + 1);
  if (open_quote == std::string::npos) {
    return std::nullopt;
  }

  const size_t close_quote = json.find('"', open_quote + 1);
  if (close_quote == std::string::npos) {
    return std::nullopt;
  }

  return json.substr(open_quote + 1, close_quote - open_quote - 1);
}

} // namespace

bool UpdateChecker::isNewerVersion(const std::string &remote_ver, const std::string &current_ver) {
  const auto remote_parts = parseSemver(remote_ver);
  const auto current_parts = parseSemver(current_ver);

  return std::tie(remote_parts[0], remote_parts[1], remote_parts[2]) >
         std::tie(current_parts[0], current_parts[1], current_parts[2]);
}

std::optional<ReleaseInfo> UpdateChecker::parseReleaseJson(const std::string &json) {
  const auto tag = extractJsonStringField(json, "\"tag_name\"");
  if (!tag.has_value() || tag->empty()) {
    return std::nullopt;
  }

  const auto url = extractJsonStringField(json, "\"html_url\"");
  return ReleaseInfo{
      .tag_name = *tag,
      .html_url = url.value_or("https://github.com/Isurwars/Correlation/releases"),
  };
}

#ifdef _WIN32
std::optional<ReleaseInfo> UpdateChecker::fetchLatestRelease() {
  HINTERNET hInternet = InternetOpenA("Correlation-App", INTERNET_OPEN_TYPE_DIRECT, nullptr, nullptr, 0);
  if (!hInternet) {
    return std::nullopt;
  }

  HINTERNET hUrl = InternetOpenUrlA(hInternet, "https://api.github.com/repos/Isurwars/Correlation/releases/latest",
                                    "User-Agent: Correlation-App\r\nAccept: application/vnd.github.v3+json\r\n", -1,
                                    INTERNET_FLAG_RELOAD | INTERNET_FLAG_SECURE | INTERNET_FLAG_DONT_CACHE, 0);
  if (!hUrl) {
    InternetCloseHandle(hInternet);
    return std::nullopt;
  }

  std::string response;
  std::array<char, 4096> buffer{};
  DWORD bytes_read = 0;
  while (InternetReadFile(hUrl, buffer.data(), static_cast<DWORD>(buffer.size() - 1), &bytes_read) && bytes_read > 0) {
    response.append(buffer.data(), bytes_read);
  }

  InternetCloseHandle(hUrl);
  InternetCloseHandle(hInternet);

  if (response.empty()) {
    return std::nullopt;
  }
  return parseReleaseJson(response);
}

void UpdateChecker::openUrlInBrowser(const std::string &url) {
  if (url.empty()) {
    return;
  }
  ShellExecuteA(nullptr, "open", url.c_str(), nullptr, nullptr, SW_SHOWNORMAL);
}
#else
std::optional<ReleaseInfo> UpdateChecker::fetchLatestRelease() {
  std::array<int, 2> pipefd{-1, -1};
  if (pipe(pipefd.data()) != 0) {
    return std::nullopt;
  }

  posix_spawn_file_actions_t actions{};
  posix_spawn_file_actions_init(&actions);
  posix_spawn_file_actions_addclose(&actions, pipefd[0]);
  posix_spawn_file_actions_adddup2(&actions, pipefd[1], STDOUT_FILENO);
  posix_spawn_file_actions_addclose(&actions, pipefd[1]);
  posix_spawn_file_actions_addopen(&actions, STDERR_FILENO, "/dev/null", O_WRONLY, 0);

  std::array<std::string, 9> arg_strs = {"curl",
                                         "-s",
                                         "-m",
                                         "5",
                                         "-H",
                                         "User-Agent: Correlation-App",
                                         "-H",
                                         "Accept: application/vnd.github.v3+json",
                                         "https://api.github.com/repos/Isurwars/Correlation/releases/latest"};

  std::array<char *, 10> args = {
      arg_strs[0].data(), arg_strs[1].data(), arg_strs[2].data(), arg_strs[3].data(), arg_strs[4].data(),
      arg_strs[5].data(), arg_strs[6].data(), arg_strs[7].data(), arg_strs[8].data(), nullptr};

  pid_t pid = 0;
  const int spawn_res = posix_spawnp(&pid, args[0], &actions, nullptr, args.data(), environ);
  posix_spawn_file_actions_destroy(&actions);
  close(pipefd[1]);

  if (spawn_res != 0) {
    close(pipefd[0]);
    return std::nullopt;
  }

  // Read response from pipe read end
  std::string response;
  std::array<char, 4096> buffer{};
  ssize_t bytes_read = 0;
  while ((bytes_read = read(pipefd[0], buffer.data(), buffer.size())) > 0) {
    response.append(buffer.data(), static_cast<size_t>(bytes_read));
  }
  close(pipefd[0]);

  int status = 0;
  waitpid(pid, &status, 0);

  if (response.empty() || !WIFEXITED(status) || WEXITSTATUS(status) != 0) {
    return std::nullopt;
  }
  return parseReleaseJson(response);
}

void UpdateChecker::openUrlInBrowser(const std::string &url) {
  if (url.empty()) {
    return;
  }
#ifdef __APPLE__
  std::string launcher = "open";
#else
  std::string launcher = "xdg-open";
#endif
  std::vector<char> url_buf(url.begin(), url.end());
  url_buf.push_back('\0');

  std::array<char *, 3> args = {launcher.data(), url_buf.data(), nullptr};

  posix_spawn_file_actions_t actions{};
  posix_spawn_file_actions_init(&actions);
  posix_spawn_file_actions_addopen(&actions, STDIN_FILENO, "/dev/null", O_RDONLY, 0);
  posix_spawn_file_actions_addopen(&actions, STDOUT_FILENO, "/dev/null", O_WRONLY, 0);
  posix_spawn_file_actions_addopen(&actions, STDERR_FILENO, "/dev/null", O_WRONLY, 0);

  pid_t pid = 0;
  const int spawn_res = posix_spawnp(&pid, args[0], &actions, nullptr, args.data(), environ);
  posix_spawn_file_actions_destroy(&actions);
  if (spawn_res == 0) {
    waitpid(pid, nullptr, 0);
  }
}
#endif

#if __has_include("AppWindow.h")
void UpdateChecker::checkForUpdatesAsync(AppWindow &window, const std::string &current_version) {
  std::thread([&window, current_version]() {
    try {
      const auto release = fetchLatestRelease();
      if (release.has_value() && isNewerVersion(release->tag_name, current_version)) {
        const std::string tag = release->tag_name;
        const std::string url =
            release->html_url.empty() ? "https://github.com/Isurwars/Correlation/releases" : release->html_url;

        slint::invoke_from_event_loop([&window, tag, url]() {
          window.set_update_version(slint::SharedString(tag));
          window.set_update_url(slint::SharedString(url));
          window.set_update_available(true);
        });
      }
    } catch (const std::exception &e) {
      std::cerr << "Update check failed: " << e.what() << '\n';
    }
  }).detach();
}
#else
void UpdateChecker::checkForUpdatesAsync(AppWindow & /*window*/, const std::string & /*current_version*/) {
  // No-op in headless test builds without Slint AppWindow
}
#endif

} // namespace correlation::app
