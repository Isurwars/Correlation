/**
 * @file NaturalSort.hpp
 * @brief Natural alphanumeric sorting and glob pattern expansion utilities.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#ifndef _WIN32
#include <glob.h>
#endif

namespace correlation::readers {

namespace detail {

inline std::string_view extractDigits(std::string_view str_view, size_t &pos) noexcept {
  const size_t start = pos;
  while (pos < str_view.size() && std::isdigit(static_cast<unsigned char>(str_view[pos])) != 0) {
    ++pos;
  }
  return str_view.substr(start, pos - start);
}

inline std::string_view stripLeadingZeros(std::string_view num) noexcept {
  size_t z = 0;
  while (z + 1 < num.size() && num[z] == '0') {
    ++z;
  }
  return num.substr(z);
}

/**
 * @brief Checks if a string matches a wildcard pattern containing '*' and '?'.
 *
 * @param pattern Wildcard pattern.
 * @param str String to test against the pattern.
 * @return True if pattern matches str, false otherwise.
 */
inline bool wildcardMatch(std::string_view pattern, std::string_view str) noexcept {
  size_t p_idx = 0;
  size_t s_idx = 0;
  size_t star_idx = std::string_view::npos;
  size_t match_idx = 0;

  while (s_idx < str.size()) {
    if (p_idx < pattern.size() && (pattern[p_idx] == '?' || pattern[p_idx] == str[s_idx])) {
      ++p_idx;
      ++s_idx;
    } else if (p_idx < pattern.size() && pattern[p_idx] == '*') {
      star_idx = p_idx;
      match_idx = s_idx;
      ++p_idx;
    } else if (star_idx != std::string_view::npos) {
      p_idx = star_idx + 1;
      ++match_idx;
      s_idx = match_idx;
    } else {
      return false;
    }
  }

  while (p_idx < pattern.size() && pattern[p_idx] == '*') {
    ++p_idx;
  }

  return p_idx == pattern.size();
}

/**
 * @brief Matches directory entries against a component pattern for all candidate paths.
 *
 * @param component_pattern Wildcard pattern for the current path component.
 * @param is_last True if this is the terminal component in the path.
 * @param candidate_paths Current accumulated base paths.
 * @return Vector of matched valid paths.
 */
inline std::vector<std::filesystem::path>
matchComponentEntries(std::string_view component_pattern, bool is_last,
                      const std::vector<std::filesystem::path> &candidate_paths) {
  std::vector<std::filesystem::path> next_candidates;
  for (const auto &base_path : candidate_paths) {
    std::error_code error_code;
    const std::filesystem::path search_dir =
        base_path.empty() ? std::filesystem::path(".") : base_path;
    if (!std::filesystem::exists(search_dir, error_code)) {
      continue;
    }
    for (const auto &entry : std::filesystem::directory_iterator(search_dir, error_code)) {
      if (error_code) {
        break;
      }
      const std::string filename_str = entry.path().filename().string();
      if (wildcardMatch(component_pattern, filename_str)) {
        if (is_last || entry.is_directory(error_code)) {
          next_candidates.push_back(base_path.empty() ? std::filesystem::path(filename_str)
                                                      : (base_path / filename_str));
        }
      }
    }
  }
  return next_candidates;
}

/**
 * @brief Cross-platform filesystem-based wildcard expansion.
 *
 * @param pattern Wildcard pattern to expand.
 * @return Vector of matching file paths.
 * @throws std::runtime_error If the pattern matches no files.
 */
inline std::vector<std::string> expandWildcardsFilesystem(const std::string &pattern) {
  if (!pattern.contains('*') && !pattern.contains('?')) {
    return {pattern};
  }

  const std::filesystem::path pattern_path(pattern);
  std::vector<std::filesystem::path> current_candidates = {std::filesystem::path{}};

  for (auto comp_it = pattern_path.begin(); comp_it != pattern_path.end(); ++comp_it) {
    const std::string comp_str = comp_it->string();
    const bool is_last = (std::next(comp_it) == pattern_path.end());
    const bool has_wildcard = comp_str.contains('*') || comp_str.contains('?');

    if (has_wildcard) {
      current_candidates = matchComponentEntries(comp_str, is_last, current_candidates);
    } else {
      std::vector<std::filesystem::path> next_candidates;
      next_candidates.reserve(current_candidates.size());
      for (const auto &base_path : current_candidates) {
        next_candidates.push_back(base_path.empty() ? std::filesystem::path(comp_str)
                                                    : (base_path / comp_str));
      }
      current_candidates = std::move(next_candidates);
    }

    if (current_candidates.empty()) {
      break;
    }
  }

  if (current_candidates.empty()) {
    throw std::runtime_error("Wildcard pattern matched no files: " + pattern);
  }

  std::vector<std::string> results;
  results.reserve(current_candidates.size());
  for (const auto &candidate : current_candidates) {
    results.push_back(candidate.string());
  }

  return results;
}

} // namespace detail

/**
 * @brief Performs natural alphanumeric comparison of two strings.
 *
 * Ensures that numeric components are compared by value (e.g., "dump.2" < "dump.10").
 * Satisfies Strict Weak Ordering.
 */
inline bool naturalLess(std::string_view str_a, std::string_view str_b) noexcept {
  size_t i_idx = 0;
  size_t j_idx = 0;

  while (i_idx < str_a.size() && j_idx < str_b.size()) {
    const bool is_digit_a = std::isdigit(static_cast<unsigned char>(str_a[i_idx])) != 0;
    const bool is_digit_b = std::isdigit(static_cast<unsigned char>(str_b[j_idx])) != 0;

    if (is_digit_a && is_digit_b) {
      const std::string_view num_a = detail::extractDigits(str_a, i_idx);
      const std::string_view num_b = detail::extractDigits(str_b, j_idx);

      const std::string_view sig_a = detail::stripLeadingZeros(num_a);
      const std::string_view sig_b = detail::stripLeadingZeros(num_b);

      if (sig_a.size() != sig_b.size()) {
        return sig_a.size() < sig_b.size();
      }
      if (sig_a != sig_b) {
        return sig_a < sig_b;
      }
      if (num_a.size() != num_b.size()) {
        return num_a.size() > num_b.size();
      }
    } else {
      const char c_a = static_cast<char>(std::tolower(static_cast<unsigned char>(str_a[i_idx])));
      const char c_b = static_cast<char>(std::tolower(static_cast<unsigned char>(str_b[j_idx])));
      if (c_a != c_b) {
        return c_a < c_b;
      }
      ++i_idx;
      ++j_idx;
    }
  }

  return i_idx >= str_a.size() && j_idx < str_b.size();
}

/**
 * @brief Sorts a container of file paths in natural alphanumeric order.
 */
inline void naturalSort(std::vector<std::string> &paths) { std::ranges::sort(paths, naturalLess); }

/**
 * @brief Expands wildcard patterns (*, ?) in a path and sorts naturally.
 *
 * @param pattern File path or glob pattern.
 * @return Sorted vector of matching file paths.
 * @throws std::runtime_error If a wildcard pattern matches no files or globbing fails.
 */
inline std::vector<std::string> expandWildcards(const std::string &pattern) {
  if (!pattern.contains('*') && !pattern.contains('?')) {
    return {pattern};
  }

#ifndef _WIN32
  glob_t glob_result{};
  const int res = glob(pattern.c_str(), GLOB_TILDE, nullptr, &glob_result);
  if (res == GLOB_NOMATCH) {
    globfree(&glob_result);
    throw std::runtime_error("Wildcard pattern matched no files: " + pattern);
  }
  if (res != 0) {
    globfree(&glob_result);
    throw std::runtime_error("Error expanding wildcard pattern: " + pattern);
  }

  std::vector<std::string> files;
  files.reserve(glob_result.gl_pathc);
  for (size_t k = 0; k < glob_result.gl_pathc; ++k) {
    files.emplace_back(glob_result.gl_pathv[k]);
  }
  globfree(&glob_result);

  naturalSort(files);
  return files;
#else
  auto files = detail::expandWildcardsFilesystem(pattern);
  naturalSort(files);
  return files;
#endif
}

} // namespace correlation::readers
