#pragma once
#include <filesystem>
#include <fstream>
#include <random>
#include <string>

namespace correlation::fuzz {

class FuzzFile {
public:
  explicit FuzzFile(const std::string &ext) {
    static thread_local std::random_device random_device;
    static thread_local std::mt19937 generator(random_device());
    std::uniform_int_distribution<uint64_t> distribution;
    std::string const filename = "fuzz_" + std::to_string(distribution(generator)) + ext;
    path_ = (std::filesystem::temp_directory_path() / filename).string();
  }

  ~FuzzFile() { static_cast<void>(std::remove(path_.c_str())); }

  FuzzFile(const FuzzFile &) = delete;
  FuzzFile &operator=(const FuzzFile &) = delete;
  FuzzFile(FuzzFile &&) = delete;
  FuzzFile &operator=(FuzzFile &&) = delete;

  [[nodiscard]] const std::string &path() const noexcept { return path_; }

  void write(const uint8_t *data, size_t size) const {
    std::ofstream file(path_, std::ios::binary | std::ios::trunc);
    if (file.is_open()) {
      const auto *char_data = std::bit_cast<const char *>(data);
      file.write(char_data, static_cast<std::streamsize>(size));
    }
  }

private:
  std::string path_;
};

inline std::string getTempFuzzPath(const std::string &ext) {
  static thread_local std::random_device random_device;
  static thread_local std::mt19937 generator(random_device());
  std::uniform_int_distribution<uint64_t> distribution;
  std::string const filename = "fuzz_" + std::to_string(distribution(generator)) + ext;
  return (std::filesystem::temp_directory_path() / filename).string();
}

} // namespace correlation::fuzz
