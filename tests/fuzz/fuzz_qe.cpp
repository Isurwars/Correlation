/**
 * @file fuzz_qe.cpp
 * @brief libFuzzer harness for QEReader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/QEReader.hpp"

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  if (size > 1 * 1024 * 1024)
    return 0;

  const std::string path = "/dev/shm/fuzz_qe.in";
  {
    std::ofstream f(path, std::ios::binary | std::ios::trunc);
    f.write(reinterpret_cast<const char *>(data), size);
  }

  try {
    correlation::readers::QEReader reader;
    reader.readTrajectory(path);
  } catch (...) {
  }

  std::remove(path.c_str());
  return 0;
}
