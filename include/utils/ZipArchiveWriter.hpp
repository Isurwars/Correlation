/**
 * @file ZipArchiveWriter.hpp
 * @brief RAII wrapper for creating ZIP archives using miniz.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include <cstddef>
#include <expected>
#include <filesystem>
#include <memory>
#include <span>
#include <string>
#include <string_view>

namespace correlation::utils {

/**
 * @class ZipArchiveWriter
 * @brief Move-only RAII wrapper around miniz for creating ZIP archives.
 */
class ZipArchiveWriter {
public:
  ZipArchiveWriter();
  ~ZipArchiveWriter();

  ZipArchiveWriter(const ZipArchiveWriter &) = delete;
  ZipArchiveWriter &operator=(const ZipArchiveWriter &) = delete;

  ZipArchiveWriter(ZipArchiveWriter &&other) noexcept;
  ZipArchiveWriter &operator=(ZipArchiveWriter &&other) noexcept;

  /**
   * @brief Opens a new ZIP archive for writing on disk.
   * @param archive_path The destination path for the ZIP file.
   * @return std::expected<void, std::string> Success or error message.
   */
  std::expected<void, std::string> open(const std::filesystem::path &archive_path);

  /**
   * @brief Adds a file to the archive from an in-memory buffer.
   * @param entry_name Path inside the ZIP archive (e.g. "plots/rdf.svg").
   * @param data Binary data buffer.
   * @return std::expected<void, std::string> Success or error message.
   */
  std::expected<void, std::string> addFileFromMemory(std::string_view entry_name,
                                                     std::span<const std::byte> data);

  /**
   * @brief Adds a text file to the archive from a string or string_view.
   * @param entry_name Path inside the ZIP archive (e.g. "summary.txt").
   * @param content Text or data string.
   * @return std::expected<void, std::string> Success or error message.
   */
  std::expected<void, std::string> addFileFromString(std::string_view entry_name,
                                                     std::string_view content);

  /**
   * @brief Adds a file from the local filesystem to the archive.
   * @param entry_name Path inside the ZIP archive.
   * @param source_path Path on disk to read from.
   * @return std::expected<void, std::string> Success or error message.
   */
  std::expected<void, std::string> addFileFromDisk(std::string_view entry_name,
                                                   const std::filesystem::path &source_path);

  /**
   * @brief Finalizes and closes the archive.
   * @return std::expected<void, std::string> Success or error message.
   */
  std::expected<void, std::string> finalize();

  /**
   * @brief Checks if an archive is currently open.
   */
  [[nodiscard]] bool isOpen() const noexcept;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace correlation::utils
