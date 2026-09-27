/**
 * @file ZipArchiveWriter.cpp
 * @brief Implementation of the miniz RAII wrapper.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "utils/ZipArchiveWriter.hpp"

#include <miniz.h>

#include <format>
#include <utility>

namespace correlation::utils {

struct ZipArchiveWriter::Impl {
  mz_zip_archive archive{};
  bool is_open{false};
  bool is_finalized{false};

  Impl() { mz_zip_zero_struct(&archive); }

  ~Impl() { cleanup(); }

  void cleanup() noexcept {
    if (is_open) {
      mz_zip_writer_end(&archive);
      is_open = false;
    }
  }
};

ZipArchiveWriter::ZipArchiveWriter() : impl_(std::make_unique<Impl>()) {}

ZipArchiveWriter::~ZipArchiveWriter() = default;

ZipArchiveWriter::ZipArchiveWriter(ZipArchiveWriter &&other) noexcept = default;

ZipArchiveWriter &ZipArchiveWriter::operator=(ZipArchiveWriter &&other) noexcept = default;

bool ZipArchiveWriter::isOpen() const noexcept { return impl_ && impl_->is_open; }

std::expected<void, std::string> ZipArchiveWriter::open(const std::filesystem::path &archive_path) {
  if (!impl_) {
    impl_ = std::make_unique<Impl>();
  }
  if (impl_->is_open) {
    return std::unexpected("ZipArchiveWriter is already open");
  }

  mz_bool status = mz_zip_writer_init_file(&impl_->archive, archive_path.string().c_str(), 0);
  if (!status) {
    const char *err = mz_zip_get_error_string(mz_zip_get_last_error(&impl_->archive));
    return std::unexpected(std::format("Failed to initialize zip file '{}': {}",
                                       archive_path.string(), err ? err : "unknown error"));
  }

  impl_->is_open = true;
  impl_->is_finalized = false;
  return {};
}

std::expected<void, std::string>
ZipArchiveWriter::addFileFromMemory(std::string_view entry_name, std::span<const std::byte> data) {
  if (!impl_ || !impl_->is_open) {
    return std::unexpected("Zip archive is not open");
  }
  if (impl_->is_finalized) {
    return std::unexpected("Zip archive is already finalized");
  }

  const std::string name(entry_name);
  mz_bool status = mz_zip_writer_add_mem(&impl_->archive, name.c_str(), data.data(), data.size(),
                                         MZ_DEFAULT_COMPRESSION);
  if (!status) {
    const char *err = mz_zip_get_error_string(mz_zip_get_last_error(&impl_->archive));
    return std::unexpected(
        std::format("Failed to add memory entry '{}': {}", name, err ? err : "unknown error"));
  }
  return {};
}

std::expected<void, std::string> ZipArchiveWriter::addFileFromString(std::string_view entry_name,
                                                                     std::string_view content) {
  return addFileFromMemory(
      entry_name, std::span<const std::byte>(reinterpret_cast<const std::byte *>(content.data()),
                                             content.size()));
}

std::expected<void, std::string>
ZipArchiveWriter::addFileFromDisk(std::string_view entry_name,
                                  const std::filesystem::path &source_path) {
  if (!impl_ || !impl_->is_open) {
    return std::unexpected("Zip archive is not open");
  }
  if (impl_->is_finalized) {
    return std::unexpected("Zip archive is already finalized");
  }
  if (!std::filesystem::exists(source_path)) {
    return std::unexpected(std::format("Source file does not exist: {}", source_path.string()));
  }

  const std::string name(entry_name);
  mz_bool status =
      mz_zip_writer_add_file(&impl_->archive, name.c_str(), source_path.string().c_str(), nullptr,
                             0, MZ_DEFAULT_COMPRESSION);
  if (!status) {
    const char *err = mz_zip_get_error_string(mz_zip_get_last_error(&impl_->archive));
    return std::unexpected(std::format("Failed to add file '{}' from disk '{}': {}", name,
                                       source_path.string(), err ? err : "unknown error"));
  }
  return {};
}

std::expected<void, std::string> ZipArchiveWriter::finalize() {
  if (!impl_ || !impl_->is_open) {
    return std::unexpected("Zip archive is not open");
  }
  if (impl_->is_finalized) {
    return {};
  }

  mz_bool status = mz_zip_writer_finalize_archive(&impl_->archive);
  if (!status) {
    const char *err = mz_zip_get_error_string(mz_zip_get_last_error(&impl_->archive));
    impl_->cleanup();
    return std::unexpected(
        std::format("Failed to finalize zip archive: {}", err ? err : "unknown error"));
  }

  impl_->is_finalized = true;
  impl_->cleanup();
  return {};
}

} // namespace correlation::utils
