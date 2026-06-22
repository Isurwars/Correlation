/**
 * @file fuzz_castep_md.cpp
 * @brief libFuzzer harness for CastepMdReader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/CastepMdReader.hpp"

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

  std::string const path = correlation::fuzz::getTempFuzzPath(".md");
  {
    std::ofstream f(path, std::ios::binary | std::ios::trunc);
    f.write(reinterpret_cast<const char *>(data), size);
  }

  try {
    correlation::readers::CastepMdReader reader;
    reader.readTrajectory(path);
  } catch (...) {
  }

  std::remove(path.c_str());
  return 0;
}
