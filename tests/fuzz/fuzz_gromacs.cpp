/**
 * @file fuzz_gromacs.cpp
 * @brief libFuzzer harness for GromacsReader (in-memory).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/GromacsReader.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <print>
#include <string>

#include "fuzz_utils.hpp"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  if (size > static_cast<size_t>(1 * 1024 * 1024)) {
    return 0;
  }

  try {
    const auto *char_data = std::bit_cast<const char *>(data);
    correlation::readers::GromacsReader::parseGroFrame(char_data, size);
  } catch (const std::exception &e) {
    if (std::getenv("FUZZ_VERBOSE") != nullptr) {
      std::println(stderr, "Error parsing gro frame: {}", e.what());
    }
  } catch (...) {
    if (std::getenv("FUZZ_VERBOSE") != nullptr) {
      std::println(stderr, "Unknown error parsing gro frame");
    }
  }

  static thread_local const correlation::fuzz::FuzzFile fuzz_file(".gro");
  fuzz_file.write(data, size);

  try {
    correlation::readers::GromacsReader reader;
    reader.readStructure(fuzz_file.path());
  } catch (const std::exception &e) {
    if (std::getenv("FUZZ_VERBOSE") != nullptr) {
      std::println(stderr, "Error parsing structure: {}", e.what());
    }
  } catch (...) {
    if (std::getenv("FUZZ_VERBOSE") != nullptr) {
      std::println(stderr, "Unknown error parsing structure");
    }
  }

  try {
    correlation::readers::GromacsReader reader;
    reader.readTrajectory(fuzz_file.path());
  } catch (const std::exception &e) {
    if (std::getenv("FUZZ_VERBOSE") != nullptr) {
      std::println(stderr, "Error parsing trajectory: {}", e.what());
    }
  } catch (...) {
    if (std::getenv("FUZZ_VERBOSE") != nullptr) {
      std::println(stderr, "Unknown error parsing trajectory");
    }
  }

  return 0;
}
