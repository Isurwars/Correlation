/**
 * @file fuzz_vasp.cpp
 * @brief libFuzzer harness for VaspReader (POSCAR/CONTCAR).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/VaspReader.hpp"

#include <cstdlib>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <string>

#include "fuzz_utils.hpp"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  if (size > static_cast<size_t>(1 * 1024 * 1024)) {
    return 0;
  }

  static thread_local correlation::fuzz::FuzzFile fuzz_file(".poscar");
  fuzz_file.write(data, size);

  try {
    correlation::readers::VaspReader reader;
    reader.readStructure(fuzz_file.path());
  } catch (const std::exception &e) {
    if (std::getenv("FUZZ_VERBOSE") != nullptr) {
      std::fprintf(stderr, "Error parsing file: %s\n", e.what());
    }
  } catch (...) {
    if (std::getenv("FUZZ_VERBOSE") != nullptr) {
      std::fprintf(stderr, "Unknown error parsing file\n");
    }
  }

  return 0;
}
