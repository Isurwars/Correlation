/**
 * @file fuzz_pdb.cpp
 * @brief libFuzzer harness for PdbReader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/PdbReader.hpp"

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>

#include "fuzz_utils.hpp"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  if (size > 1 * 1024 * 1024)
    return 0;

  static thread_local correlation::fuzz::FuzzFile fuzz_file(".pdb");
  fuzz_file.write(data, size);

  try {
    correlation::readers::PdbReader reader;
    reader.readStructure(fuzz_file.path());
  } catch (...) {
    // Catch-all to prevent fuzzer crashes on invalid inputs.
  }

  try {
    correlation::readers::PdbReader reader;
    reader.readTrajectory(fuzz_file.path());
  } catch (...) {
    // Catch-all to prevent fuzzer crashes on invalid inputs.
  }

  
  return 0;
}
