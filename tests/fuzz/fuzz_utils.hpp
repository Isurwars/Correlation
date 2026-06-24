#pragma once
#include <filesystem>
#include <fstream>
#include <random>
#include <string>

namespace correlation::fuzz {

class FuzzFile {
public:
  explicit FuzzFile(const std::string &ext) {
    static thread_local std::random_device rd;
    static thread_local std::mt19937 generator(rd());
    std::uniform_int_distribution<uint64_t> distribution;
    std::string const filename = "fuzz_" + std::to_string(distribution(generator)) + ext;
    path_ = (std::filesystem::temp_directory_path() / filename).string();
  }

  ~FuzzFile() {
    std::remove(path_.c_str());
  }

  FuzzFile(const FuzzFile&) = delete;
  FuzzFile& operator=(const FuzzFile&) = delete;
  FuzzFile(FuzzFile&&) = delete;
  FuzzFile& operator=(FuzzFile&&) = delete;

  [[nodiscard]] const std::string& path() const noexcept { return path_; }

  void write(const uint8_t *data, size_t size) const {
    std::ofstream f(path_, std::ios::binary | std::ios::trunc);
    if (f.is_open()) {
      f.write(reinterpret_cast<const char *>(data), static_cast<std::streamsize>(size));
    }
  }

private:
  std::string path_;
};

inline std::string getTempFuzzPath(const std::string &ext) {
  static thread_local std::random_device rd;
  static thread_local std::mt19937 generator(rd());
  std::uniform_int_distribution<uint64_t> distribution;
  std::string const filename = "fuzz_" + std::to_string(distribution(generator)) + ext;
  return (std::filesystem::temp_directory_path() / filename).string();
}

} // namespace correlation::fuzz
