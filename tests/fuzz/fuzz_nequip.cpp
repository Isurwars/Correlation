/**
 * @file fuzz_nequip.cpp
 * @brief libFuzzer harness for NequipReader (in-memory).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/NequipReader.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>

#include "fuzz_utils.hpp"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  if (size > static_cast<size_t>(1 * 1024 * 1024)) {
    return 0;
  }

  static thread_local correlation::fuzz::FuzzFile fuzz_file(".nequip");
  fuzz_file.write(data, size);

  try {
    correlation::readers::NequipReader reader;
    reader.readStructure(fuzz_file.path());
  } catch (const std::exception &e) {
    if (std::getenv("FUZZ_VERBOSE") != nullptr) {
      std::cerr << "Error parsing structure: " << e.what() << '\n';
    }
  } catch (...) {
    if (std::getenv("FUZZ_VERBOSE") != nullptr) {
      std::cerr << "Unknown error parsing structure\n";
    }
  }

  try {
    correlation::readers::NequipReader reader;
    reader.readTrajectory(fuzz_file.path());
  } catch (const std::exception &e) {
    if (std::getenv("FUZZ_VERBOSE") != nullptr) {
      std::cerr << "Error parsing trajectory: " << e.what() << '\n';
    }
  } catch (...) {
    if (std::getenv("FUZZ_VERBOSE") != nullptr) {
      std::cerr << "Unknown error parsing trajectory\n";
    }
  }

  return 0;
}
