#pragma once
#include <filesystem>
#include <random>
#include <string>

namespace correlation::fuzz {

inline std::string getTempFuzzPath(const std::string &ext) {
  static thread_local std::random_device rd;
  static thread_local std::mt19937 generator(rd());
  std::uniform_int_distribution<uint64_t> distribution;
  std::string const filename = "fuzz_" + std::to_string(distribution(generator)) + ext;
  return (std::filesystem::temp_directory_path() / filename).string();
}

} // namespace correlation::fuzz
