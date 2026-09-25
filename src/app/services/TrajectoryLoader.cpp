/**
 * @file TrajectoryLoader.cpp
 * @brief Implementation of the TrajectoryLoader service.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/services/TrajectoryLoader.hpp"
#include "app/services/PhysicsService.hpp"
#include "readers/FileReader.hpp"
#include "readers/ReaderFactory.hpp"

#include <algorithm>
#include <filesystem>
#include <vector>

namespace correlation::app {

std::map<std::string, int> TrajectoryLoader::getAtomCounts() const {
  std::map<std::string, int> counts;
  const correlation::core::Cell *current_cell = cell();
  if (current_cell == nullptr) {
    return counts;
  }

  const auto &elements = current_cell->elements();
  std::vector<int> id_counts(elements.size(), 0);

  for (const auto &atom : current_cell->atoms()) {
    id_counts[atom.elementId()]++;
  }

  for (size_t i = 0; i < elements.size(); ++i) {
    if (id_counts[i] > 0) {
      counts[elements[i].symbol] = id_counts[i];
    }
  }

  return counts;
}

real_t TrajectoryLoader::getRecommendedTimeStep() const {
  return PhysicsService::computeRecommendedTimeStep(cell());
}
std::expected<std::string, std::string>
TrajectoryLoader::loadFile(const std::string &path, const ProgressCallback &progress) {
  try {
    std::string display_path = path;
    std::ranges::replace(display_path, '\\', '/');
    correlation::readers::FileType const type = correlation::readers::determineFileType(path);

    bool is_trajectory = false;
    std::string const ext = std::filesystem::path(path).extension().string();
    if (!ext.empty()) {
      auto *reader = correlation::readers::ReaderFactory::instance().getReaderForExtension(
          {.extension = ext, .filename = path});
      if (reader != nullptr) {
        is_trajectory = reader->isTrajectory();
      }
    } else {
      is_trajectory = (type == correlation::readers::FileType::Xdatcar);
    }

    if (is_trajectory) {
      trajectory_ = std::make_unique<correlation::core::Trajectory>(
          correlation::readers::readTrajectory(path, type, progress));
    } else {
      trajectory_ = std::make_unique<correlation::core::Trajectory>();
      trajectory_->addFrame(correlation::readers::readStructure(path, type, progress));
    }

    return "File loaded: " + display_path;
  } catch (const std::exception &e) {
    return std::unexpected(std::string("Error loading file: ") + e.what());
  }
}

} // namespace correlation::app
