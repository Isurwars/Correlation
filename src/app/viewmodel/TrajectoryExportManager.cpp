/**
 * @file TrajectoryExportManager.cpp
 * @brief Implementation of PNG-sequence ZIP and Animated GIF89a trajectory export.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/viewmodel/TrajectoryExportManager.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <format>
#include <fstream>
#include <miniz.h>
#include <unordered_map>
#include <vector>

namespace correlation::app {

namespace {

constexpr std::size_t K_PALETTE_COLORS = 256;
constexpr std::uint16_t K_CLEAR_CODE = 256;
constexpr std::uint16_t K_EOI_CODE = 257;
constexpr std::uint16_t K_FIRST_CODE = 258;
constexpr std::uint16_t K_MAX_CODE = 4095;

// 6x6x6 RGB cube (216 colors) + 40 grayscale ramp entries = 256-color global palette
[[nodiscard]] std::array<std::uint8_t, K_PALETTE_COLORS * 3> buildGlobalPalette() noexcept {
  std::array<std::uint8_t, K_PALETTE_COLORS * 3> pal{};
  std::size_t idx = 0;
  for (int red = 0; red < 6; ++red) {
    for (int green = 0; green < 6; ++green) {
      for (int blue = 0; blue < 6; ++blue) {
        pal.at(idx * 3 + 0) = static_cast<std::uint8_t>(red * 51);
        pal.at(idx * 3 + 1) = static_cast<std::uint8_t>(green * 51);
        pal.at(idx * 3 + 2) = static_cast<std::uint8_t>(blue * 51);
        ++idx;
      }
    }
  }
  for (int gray = 0; gray < 40; ++gray) {
    const auto val = static_cast<std::uint8_t>(std::clamp((gray * 255) / 39, 0, 255));
    pal.at(idx * 3 + 0) = val;
    pal.at(idx * 3 + 1) = val;
    pal.at(idx * 3 + 2) = val;
    ++idx;
  }
  return pal;
}

[[nodiscard]] std::uint8_t quantizeRgbToPalette(std::uint8_t red, std::uint8_t green,
                                                std::uint8_t blue) noexcept {
  const int max_c = std::max({red, green, blue});
  const int min_c = std::min({red, green, blue});
  if (max_c - min_c <= 8) {
    const int lum = (static_cast<int>(red) + static_cast<int>(green) + static_cast<int>(blue)) / 3;
    const int gray_idx = std::clamp((lum * 39 + 127) / 255, 0, 39);
    return static_cast<std::uint8_t>(216 + gray_idx);
  }
  const int r_bin = (static_cast<int>(red) * 5 + 127) / 255;
  const int g_bin = (static_cast<int>(green) * 5 + 127) / 255;
  const int b_bin = (static_cast<int>(blue) * 5 + 127) / 255;
  return static_cast<std::uint8_t>(r_bin * 36 + g_bin * 6 + b_bin);
}

class GifBitWriter {
public:
  void writeCode(std::uint16_t code, std::uint32_t bit_width) {
    bit_buf_ |= (static_cast<std::uint32_t>(code) << bit_count_);
    bit_count_ += bit_width;
    while (bit_count_ >= 8U) {
      flushByte(static_cast<std::uint8_t>(bit_buf_ & 0xFFU));
      bit_buf_ >>= 8U;
      bit_count_ -= 8U;
    }
  }

  void finish() {
    if (bit_count_ > 0U) {
      flushByte(static_cast<std::uint8_t>(bit_buf_ & 0xFFU));
      bit_buf_ = 0U;
      bit_count_ = 0U;
    }
    flushSubBlock();
    bytes_.push_back(0x00); // Block terminator
  }

  [[nodiscard]] const std::vector<std::uint8_t> &bytes() const noexcept { return bytes_; }

private:
  void flushByte(std::uint8_t byte_val) {
    sub_block_.push_back(byte_val);
    if (sub_block_.size() == 255) {
      flushSubBlock();
    }
  }

  void flushSubBlock() {
    if (sub_block_.empty()) {
      return;
    }
    bytes_.push_back(static_cast<std::uint8_t>(sub_block_.size()));
    bytes_.insert(bytes_.end(), sub_block_.begin(), sub_block_.end());
    sub_block_.clear();
  }

  std::uint32_t bit_buf_{0U};
  std::uint32_t bit_count_{0U};
  std::vector<std::uint8_t> sub_block_;
  std::vector<std::uint8_t> bytes_;
};

void encodeLzwFrame(std::span<const std::uint8_t> indices, std::ofstream &out) {
  const std::uint8_t min_code_size = 8;
  out.put(static_cast<char>(min_code_size));

  GifBitWriter writer;
  std::unordered_map<std::uint32_t, std::uint16_t> dict;
  dict.reserve(4096);

  std::uint32_t bit_width = 9U;
  std::uint16_t next_code = K_FIRST_CODE;
  writer.writeCode(K_CLEAR_CODE, bit_width);

  if (!indices.empty()) {
    std::uint16_t prefix = indices[0];
    for (std::size_t idx = 1; idx < indices.size(); ++idx) {
      const std::uint8_t sym = indices[idx];
      const std::uint32_t key = (static_cast<std::uint32_t>(prefix) << 8U) | sym;
      const auto iter = dict.find(key);
      if (iter != dict.end()) {
        prefix = iter->second;
      } else {
        writer.writeCode(prefix, bit_width);
        if (next_code <= K_MAX_CODE) {
          dict.emplace(key, next_code++);
          if (next_code > (1U << bit_width) && bit_width < 12U) {
            ++bit_width;
          }
        } else {
          writer.writeCode(K_CLEAR_CODE, bit_width);
          dict.clear();
          bit_width = 9U;
          next_code = K_FIRST_CODE;
        }
        prefix = sym;
      }
    }
    writer.writeCode(prefix, bit_width);
  }
  writer.writeCode(K_EOI_CODE, bit_width);
  writer.finish();

  const auto &encoded = writer.bytes();
  out.write(reinterpret_cast<const char *>(encoded.data()),
            static_cast<std::streamsize>(encoded.size()));
}

void writeU16LE(std::ofstream &out, std::uint16_t value) {
  const auto unsigned_val = static_cast<std::uint32_t>(value);
  const std::array<char, 2> buf{
      static_cast<char>(unsigned_val & 0xFFU),
      static_cast<char>((unsigned_val >> 8U) & 0xFFU),
  };
  out.write(buf.data(), 2);
}

} // namespace

bool TrajectoryExportManager::exportPngSequenceZip(
    const std::string &filepath, const std::vector<std::vector<std::uint8_t>> &frames,
    std::uint32_t width, std::uint32_t height) {
  if (filepath.empty() || frames.empty() || width == 0 || height == 0) {
    return false;
  }

  mz_zip_archive zip{};
  if (mz_zip_writer_init_file(&zip, filepath.c_str(), 0) == 0) {
    return false;
  }

  const std::size_t expected_bytes = static_cast<std::size_t>(width) * height * 4;
  bool all_ok = true;
  for (std::size_t idx = 0; idx < frames.size(); ++idx) {
    if (frames[idx].size() < expected_bytes) {
      all_ok = false;
      break;
    }
    std::size_t png_len = 0;
    void *png_mem = tdefl_write_image_to_png_file_in_memory(
        frames[idx].data(), static_cast<int>(width), static_cast<int>(height), 4, &png_len);
    if (png_mem == nullptr || png_len == 0) {
      all_ok = false;
      break;
    }
    const std::string entry_name = std::format("frame_{:04d}.png", idx);
    const mz_bool added =
        mz_zip_writer_add_mem(&zip, entry_name.c_str(), png_mem, png_len, MZ_DEFAULT_COMPRESSION);
    mz_free(png_mem);
    if (added == 0) {
      all_ok = false;
      break;
    }
  }

  const mz_bool finalized = mz_zip_writer_finalize_archive(&zip);
  mz_zip_writer_end(&zip);
  return all_ok && (finalized != 0);
}

bool TrajectoryExportManager::exportAnimatedGif(
    const std::string &filepath, const std::vector<std::vector<std::uint8_t>> &frames,
    std::uint32_t width, std::uint32_t height, std::uint16_t delay_cs) {
  if (filepath.empty() || frames.empty() || width == 0 || height == 0 || width > 65535 ||
      height > 65535) {
    return false;
  }

  std::ofstream out(filepath, std::ios::binary);
  if (!out.is_open()) {
    return false;
  }

  // GIF89a Header + Logical Screen Descriptor (global 256-color table)
  out.write("GIF89a", 6);
  writeU16LE(out, static_cast<std::uint16_t>(width));
  writeU16LE(out, static_cast<std::uint16_t>(height));
  constexpr std::array<std::uint8_t, 3> K_LSD_FLAGS{0xF7U, 0x00U, 0x00U};
  out.write(reinterpret_cast<const char *>(K_LSD_FLAGS.data()), 3);

  const auto palette = buildGlobalPalette();
  out.write(reinterpret_cast<const char *>(palette.data()),
            static_cast<std::streamsize>(palette.size()));

  // NETSCAPE2.0 Application Extension (infinite loop)
  constexpr std::array<std::uint8_t, 19> K_NETSCAPE_EXT{
      0x21U, 0xFFU, 0x0BU, 'N', 'E',   'T',   'S',   'C',   'A',   'P',
      'E',   '2',   '.',   '0', 0x03U, 0x01U, 0x00U, 0x00U, 0x00U,
  };
  out.write(reinterpret_cast<const char *>(K_NETSCAPE_EXT.data()),
            static_cast<std::streamsize>(K_NETSCAPE_EXT.size()));

  const std::size_t pixel_count = static_cast<std::size_t>(width) * height;
  std::vector<std::uint8_t> indexed(pixel_count);

  for (const auto &rgba : frames) {
    if (rgba.size() < pixel_count * 4) {
      return false;
    }
    for (std::size_t px = 0; px < pixel_count; ++px) {
      indexed[px] = quantizeRgbToPalette(rgba[px * 4 + 0], rgba[px * 4 + 1], rgba[px * 4 + 2]);
    }

    // Graphic Control Extension
    constexpr std::array<std::uint8_t, 4> K_GCE_HEAD{0x21U, 0xF9U, 0x04U, 0x04U};
    out.write(reinterpret_cast<const char *>(K_GCE_HEAD.data()), 4);
    writeU16LE(out, std::max<std::uint16_t>(2, delay_cs));
    constexpr std::array<std::uint8_t, 2> K_GCE_TAIL{0x00U, 0x00U};
    out.write(reinterpret_cast<const char *>(K_GCE_TAIL.data()), 2);

    // Image Descriptor
    out.put(static_cast<char>(0x2CU));
    writeU16LE(out, 0);
    writeU16LE(out, 0);
    writeU16LE(out, static_cast<std::uint16_t>(width));
    writeU16LE(out, static_cast<std::uint16_t>(height));
    out.put(0x00);

    encodeLzwFrame(indexed, out);
  }

  out.put(static_cast<char>(0x3BU)); // GIF Trailer
  return out.good();
}

} // namespace correlation::app
