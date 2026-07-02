/**
 * @file fuzz_cp2k.cpp
 * @brief libFuzzer harness for CP2KReader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/CP2KReader.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <print>
#include <string>

#include "fuzz_utils.hpp"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  if (size > static_cast<size_t>(1 * 1024 * 1024)) {
    return 0;
  }

  static thread_local correlation::fuzz::FuzzFile fuzz_file(".inp");
  fuzz_file.write(data, size);

  try {
    correlation::readers::CP2KReader reader;
    reader.readTrajectory(fuzz_file.path());
  } catch (const std::exception &e) {
    if (std::getenv("FUZZ_VERBOSE") != nullptr) {
      std::println(stderr, "Error parsing file: {}", e.what());
    }
  } catch (...) {
    if (std::getenv("FUZZ_VERBOSE") != nullptr) {
      std::println(stderr, "Unknown error parsing file");
    }
  }

  return 0;
}
