/**
 * @file fuzz_xdatcar.cpp
 * @brief libFuzzer harness for XdatcarReader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/XdatcarReader.hpp"

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>

#include "fuzz_utils.hpp"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  if (size > 1 * 1024 * 1024)
    return 0;

  std::string const path = correlation::fuzz::getTempFuzzPath(".xdatcar");
  {
    std::ofstream f(path, std::ios::binary | std::ios::trunc);
    f.write(reinterpret_cast<const char *>(data), size);
  }

  try {
    correlation::readers::XdatcarReader reader;
    reader.readTrajectory(path);
  } catch (...) {
  }

  std::remove(path.c_str());
  return 0;
}
