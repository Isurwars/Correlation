/**
 * @file fuzz_xyz.cpp
 * @brief libFuzzer harness for XYZReader (in-memory).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/XYZReader.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <string>

#include "fuzz_utils.hpp"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  if (size > static_cast<size_t>(1 * 1024 * 1024)) {
    return 0;
  }

  // XYZReader has an in-memory parseXYZFrame, but it's private.
  // Use the file-based entry point for full coverage of both
  // readStructure and readTrajectory.
  static thread_local correlation::fuzz::FuzzFile fuzz_file(".xyz");
  fuzz_file.write(data, size);

  try {
    correlation::readers::XYZReader reader;
    reader.readStructure(fuzz_file.path());
  } catch (const std::exception &e) {
    if (std::getenv("FUZZ_VERBOSE") != nullptr) {
      std::fprintf(stderr, "Error parsing structure: %s\n", e.what());
    }
  } catch (...) {
    if (std::getenv("FUZZ_VERBOSE") != nullptr) {
      std::fprintf(stderr, "Unknown error parsing structure\n");
    }
  }

  try {
    correlation::readers::XYZReader reader;
    reader.readTrajectory(fuzz_file.path());
  } catch (const std::exception &e) {
    if (std::getenv("FUZZ_VERBOSE") != nullptr) {
      std::fprintf(stderr, "Error parsing trajectory: %s\n", e.what());
    }
  } catch (...) {
    if (std::getenv("FUZZ_VERBOSE") != nullptr) {
      std::fprintf(stderr, "Unknown error parsing trajectory\n");
    }
  }

  return 0;
}
