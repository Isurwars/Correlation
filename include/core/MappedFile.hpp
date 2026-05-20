/**
 * @file MappedFile.hpp
 * @brief Cross-platform memory-mapped file I/O utility.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 *
 * Provides a lightweight RAII wrapper around OS-level memory mapping:
 * - Windows: CreateFileMapping / MapViewOfFile
 * - POSIX:   mmap / munmap
 *
 * The mapped region is read-only and immutable after construction.
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace correlation::core {

/// Maximum trajectory file size (4 GiB).
inline constexpr std::uint64_t kMaxTrajectoryBytes =
    static_cast<std::uint64_t>(4) * 1024 * 1024 * 1024;

/**
 * @class MappedFile
 * @brief RAII memory-mapped read-only file view.
 *
 * The constructor opens the file and maps its entire contents into virtual
 * memory. The destructor unmaps and closes all handles. Move-only.
 *
 * Usage:
 * @code
 *   MappedFile mf("trajectory.xyz");
 *   const char* data = mf.data();
 *   size_t      len  = mf.size();
 * @endcode
 */
class MappedFile {
public:
  /**
   * @brief Opens and memory-maps a file for reading.
   * @param path  File system path.
   * @param enforce_size_limit If true, throws when file exceeds 4 GiB.
   * @throws std::runtime_error on I/O or mapping errors, or if the file
   *         exceeds the 4 GiB limit while enforce_size_limit is true.
   */
  explicit MappedFile(const std::string &path,
                      bool enforce_size_limit = true) {
#ifdef _WIN32
    file_handle_ = CreateFileA(path.c_str(), GENERIC_READ, FILE_SHARE_READ,
                               nullptr, OPEN_EXISTING,
                               FILE_ATTRIBUTE_NORMAL, nullptr);
    if (file_handle_ == INVALID_HANDLE_VALUE)
      throw std::runtime_error("MappedFile: cannot open: " + path);

    LARGE_INTEGER li;
    if (!GetFileSizeEx(file_handle_, &li)) {
      CloseHandle(file_handle_);
      throw std::runtime_error("MappedFile: cannot get size: " + path);
    }
    size_ = static_cast<std::size_t>(li.QuadPart);

    if (enforce_size_limit && static_cast<std::uint64_t>(size_) > kMaxTrajectoryBytes) {
      CloseHandle(file_handle_);
      throw std::runtime_error(
          "MappedFile: file exceeds 4 GiB trajectory limit: " + path);
    }

    mapping_ = CreateFileMappingA(file_handle_, nullptr, PAGE_READONLY,
                                  0, 0, nullptr);
    if (!mapping_) {
      CloseHandle(file_handle_);
      throw std::runtime_error("MappedFile: CreateFileMapping failed: " + path);
    }

    data_ = static_cast<const char *>(
        MapViewOfFile(mapping_, FILE_MAP_READ, 0, 0, 0));
    if (!data_) {
      CloseHandle(mapping_);
      CloseHandle(file_handle_);
      throw std::runtime_error("MappedFile: MapViewOfFile failed: " + path);
    }
#else
    fd_ = ::open(path.c_str(), O_RDONLY);
    if (fd_ < 0)
      throw std::runtime_error("MappedFile: cannot open: " + path);

    struct stat st {};
    if (::fstat(fd_, &st) != 0) {
      ::close(fd_);
      throw std::runtime_error("MappedFile: cannot stat: " + path);
    }
    size_ = static_cast<std::size_t>(st.st_size);

    if (enforce_size_limit && static_cast<std::uint64_t>(size_) > kMaxTrajectoryBytes) {
      ::close(fd_);
      throw std::runtime_error(
          "MappedFile: file exceeds 4 GiB trajectory limit: " + path);
    }

    data_ = static_cast<const char *>(
        ::mmap(nullptr, size_, PROT_READ, MAP_PRIVATE, fd_, 0));
    if (data_ == MAP_FAILED) {
      data_ = nullptr;
      ::close(fd_);
      throw std::runtime_error("MappedFile: mmap failed: " + path);
    }
#endif
  }

  ~MappedFile() { release(); }

  // Move-only.
  MappedFile(MappedFile &&other) noexcept
      : data_(other.data_), size_(other.size_)
#ifdef _WIN32
        , file_handle_(other.file_handle_), mapping_(other.mapping_)
#else
        , fd_(other.fd_)
#endif
  {
    other.data_ = nullptr;
    other.size_ = 0;
#ifdef _WIN32
    other.file_handle_ = INVALID_HANDLE_VALUE;
    other.mapping_ = nullptr;
#else
    other.fd_ = -1;
#endif
  }

  MappedFile &operator=(MappedFile &&other) noexcept {
    if (this != &other) {
      release();
      data_ = other.data_;
      size_ = other.size_;
#ifdef _WIN32
      file_handle_ = other.file_handle_;
      mapping_ = other.mapping_;
      other.file_handle_ = INVALID_HANDLE_VALUE;
      other.mapping_ = nullptr;
#else
      fd_ = other.fd_;
      other.fd_ = -1;
#endif
      other.data_ = nullptr;
      other.size_ = 0;
    }
    return *this;
  }

  MappedFile(const MappedFile &) = delete;
  MappedFile &operator=(const MappedFile &) = delete;

  /// Pointer to the start of the mapped file contents.
  [[nodiscard]] const char *data() const noexcept { return data_; }

  /// Size of the mapped region in bytes.
  [[nodiscard]] std::size_t size() const noexcept { return size_; }

  /// Convenience: pointer past the end.
  [[nodiscard]] const char *end() const noexcept { return data_ + size_; }

private:
  void release() noexcept {
#ifdef _WIN32
    if (data_)
      UnmapViewOfFile(data_);
    if (mapping_)
      CloseHandle(mapping_);
    if (file_handle_ != INVALID_HANDLE_VALUE)
      CloseHandle(file_handle_);
#else
    if (data_)
      ::munmap(const_cast<char *>(data_), size_);
    if (fd_ >= 0)
      ::close(fd_);
#endif
    data_ = nullptr;
    size_ = 0;
  }

  const char *data_{nullptr};
  std::size_t size_{0};

#ifdef _WIN32
  HANDLE file_handle_{INVALID_HANDLE_VALUE};
  HANDLE mapping_{nullptr};
#else
  int fd_{-1};
#endif
};

} // namespace correlation::core
