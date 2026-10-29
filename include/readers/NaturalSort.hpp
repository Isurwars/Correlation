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
#include <glob.h>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace correlation::readers {

namespace detail {

inline std::string_view extractDigits(std::string_view s, size_t &pos) noexcept {
  const size_t start = pos;
  while (pos < s.size() && std::isdigit(static_cast<unsigned char>(s[pos])) != 0) {
    ++pos;
  }
  return s.substr(start, pos - start);
}

inline std::string_view stripLeadingZeros(std::string_view num) noexcept {
  size_t z = 0;
  while (z + 1 < num.size() && num[z] == '0') {
    ++z;
  }
  return num.substr(z);
}

} // namespace detail

/**
 * @brief Performs natural alphanumeric comparison of two strings.
 *
 * Ensures that numeric components are compared by value (e.g., "dump.2" < "dump.10").
 * Satisfies Strict Weak Ordering.
 */
inline bool naturalLess(std::string_view a, std::string_view b) noexcept {
  size_t i = 0;
  size_t j = 0;

  while (i < a.size() && j < b.size()) {
    const bool is_digit_a = std::isdigit(static_cast<unsigned char>(a[i])) != 0;
    const bool is_digit_b = std::isdigit(static_cast<unsigned char>(b[j])) != 0;

    if (is_digit_a && is_digit_b) {
      const std::string_view num_a = detail::extractDigits(a, i);
      const std::string_view num_b = detail::extractDigits(b, j);

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
      const char ca = static_cast<char>(std::tolower(static_cast<unsigned char>(a[i])));
      const char cb = static_cast<char>(std::tolower(static_cast<unsigned char>(b[j])));
      if (ca != cb) {
        return ca < cb;
      }
      ++i;
      ++j;
    }
  }

  return i >= a.size() && j < b.size();
}

/**
 * @brief Sorts a container of file paths in natural alphanumeric order.
 */
inline void naturalSort(std::vector<std::string> &paths) {
  std::sort(paths.begin(), paths.end(), naturalLess);
}

/**
 * @brief Expands wildcard patterns (*, ?) in a path using POSIX glob and sorts naturally.
 *
 * @param pattern File path or glob pattern.
 * @return Sorted vector of matching file paths.
 * @throws std::runtime_error If a wildcard pattern matches no files or globbing fails.
 */
inline std::vector<std::string> expandWildcards(const std::string &pattern) {
  if (pattern.find('*') == std::string::npos && pattern.find('?') == std::string::npos) {
    return {pattern};
  }

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
}

} // namespace correlation::readers
