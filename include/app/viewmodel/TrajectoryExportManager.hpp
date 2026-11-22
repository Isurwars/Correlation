/**
 * @file TrajectoryExportManager.hpp
 * @brief Lossless PNG-sequence ZIP and self-contained Animated GIF89a trajectory exporter.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace correlation::app {

/**
 * @class TrajectoryExportManager
 * @brief Encodes rendered RGBA8 trajectory frames into PNG-sequence ZIP archives or Animated GIFs.
 */
class TrajectoryExportManager {
public:
  /**
   * @brief Writes a sequence of RGBA8 frames as numbered PNGs inside a single ZIP archive.
   * @param[in] filepath Destination .zip path.
   * @param[in] frames Contiguous RGBA8 frame buffers (each size == width * height * 4).
   * @param[in] width Frame pixel width.
   * @param[in] height Frame pixel height.
   * @return True if archive was written successfully.
   */
  [[nodiscard]] static bool
  exportPngSequenceZip(const std::string &filepath,
                       const std::vector<std::vector<std::uint8_t>> &frames, std::uint32_t width,
                       std::uint32_t height);

  /**
   * @brief Encodes RGBA8 frames into a self-contained looping GIF89a animation.
   * @param[in] filepath Destination .gif path.
   * @param[in] frames Contiguous RGBA8 frame buffers (each size == width * height * 4).
   * @param[in] width Frame pixel width.
   * @param[in] height Frame pixel height.
   * @param[in] delay_cs Frame delay in centiseconds (e.g. 6 = ~16 FPS).
   * @return True if GIF file was written successfully.
   */
  [[nodiscard]] static bool exportAnimatedGif(const std::string &filepath,
                                              const std::vector<std::vector<std::uint8_t>> &frames,
                                              std::uint32_t width, std::uint32_t height,
                                              std::uint16_t delay_cs = 6);
};

} // namespace correlation::app
