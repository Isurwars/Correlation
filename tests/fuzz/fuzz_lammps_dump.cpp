/**
 * @file fuzz_lammps_dump.cpp
 * @brief libFuzzer harness for LammpsDumpReader (in-memory).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/LammpsDumpReader.hpp"

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  if (size > 1 * 1024 * 1024)
    return 0;

  try {
    correlation::readers::LammpsDumpReader::parseDumpFrame(
        reinterpret_cast<const char *>(data), size);
  } catch (...) {
  }

  const std::string path = "/dev/shm/fuzz_lammps.dump";
  {
    std::ofstream f(path, std::ios::binary | std::ios::trunc);
    f.write(reinterpret_cast<const char *>(data), size);
  }

  try {
    correlation::readers::LammpsDumpReader reader;
    reader.readStructure(path);
  } catch (...) {
  }

  try {
    correlation::readers::LammpsDumpReader reader;
    reader.readTrajectory(path);
  } catch (...) {
  }

  std::remove(path.c_str());
  return 0;
}

