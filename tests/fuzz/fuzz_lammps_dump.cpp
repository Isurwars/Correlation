/**
 * @file fuzz_lammps_dump.cpp
 * @brief libFuzzer harness for LammpsDumpReader (in-memory).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/LammpsDumpReader.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>

#include "fuzz_utils.hpp"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  if (size > static_cast<size_t>(1 * 1024 * 1024)) {
    return 0;
}

  try {
    correlation::readers::LammpsDumpReader::parseDumpFrame(
        reinterpret_cast<const char *>(data), size);
  } catch (...) {
    // Catch-all to prevent fuzzer crashes on invalid inputs.
  }

  static thread_local correlation::fuzz::FuzzFile fuzz_file(".dump");
  fuzz_file.write(data, size);

  try {
    correlation::readers::LammpsDumpReader reader;
    reader.readStructure(fuzz_file.path());
  } catch (...) {
    // Catch-all to prevent fuzzer crashes on invalid inputs.
  }

  try {
    correlation::readers::LammpsDumpReader reader;
    reader.readTrajectory(fuzz_file.path());
  } catch (...) {
    // Catch-all to prevent fuzzer crashes on invalid inputs.
  }

  
  return 0;
}

