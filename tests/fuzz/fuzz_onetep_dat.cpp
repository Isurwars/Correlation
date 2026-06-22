/**
 * @file fuzz_onetep_dat.cpp
 * @brief libFuzzer harness for OnetepDatReader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/OnetepDatReader.hpp"

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>

#include "fuzz_utils.hpp"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  if (size > 1 * 1024 * 1024)
    return 0;

  std::string const path = correlation::fuzz::getTempFuzzPath(".dat");
  {
    std::ofstream f(path, std::ios::binary | std::ios::trunc);
    f.write(reinterpret_cast<const char *>(data), size);
  }

  try {
    correlation::readers::OnetepDatReader reader;
    reader.readStructure(path);
  } catch (...) {
  }

  std::remove(path.c_str());
  return 0;
}
