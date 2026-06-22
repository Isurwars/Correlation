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
#include <fstream>
#include <string>

#include "fuzz_utils.hpp"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  if (size > static_cast<size_t>(1 * 1024 * 1024)) {
    return 0;
}

  // XYZReader has an in-memory parseXYZFrame, but it's private.
  // Use the file-based entry point for full coverage of both
  // readStructure and readTrajectory.
  std::string const path = correlation::fuzz::getTempFuzzPath(".xyz");
  {
    std::ofstream f(path, std::ios::binary | std::ios::trunc);
    f.write(reinterpret_cast<const char *>(data), size);
  }

  try {
    correlation::readers::XYZReader reader;
    reader.readStructure(path);
  } catch (...) {
  }

  try {
    correlation::readers::XYZReader reader;
    reader.readTrajectory(path);
  } catch (...) {
  }

  std::remove(path.c_str());
  return 0;
}
