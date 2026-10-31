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
#include "readers/NaturalSort.hpp"
#include "readers/ReaderFactory.hpp"

#include <algorithm>
#include <filesystem>
#include <vector>

namespace correlation::app {

TrajectoryLoader::TrajectoryLoader()
    : TrajectoryLoader(correlation::readers::ReaderFactory::instance()) {}

TrajectoryLoader::TrajectoryLoader(const correlation::readers::IReaderRegistry &registry)
    : registry_(&registry) {}

std::map<std::string, int> TrajectoryLoader::getAtomCounts() const {
  std::map<std::string, int> counts;
  for (const auto &[symbol, count] : getOrderedAtomCounts()) {
    counts[symbol] = count;
  }
  return counts;
}

std::vector<std::pair<std::string, int>> TrajectoryLoader::getOrderedAtomCounts() const {
  std::vector<std::pair<std::string, int>> counts;
  const correlation::core::Cell *current_cell = cell();
  if (current_cell == nullptr) {
    return counts;
  }

  const auto &elements = current_cell->elements();
  std::vector<int> id_counts(elements.size(), 0);

  for (const auto &atom : current_cell->atoms()) {
    id_counts[atom.elementId()]++;
  }

  counts.reserve(elements.size());
  for (size_t i = 0; i < elements.size(); ++i) {
    if (id_counts[i] > 0) {
      counts.emplace_back(elements[i].symbol, id_counts[i]);
    }
  }

  return counts;
}

bool TrajectoryLoader::renameElement(std::string_view old_symbol, std::string_view new_symbol) {
  if (!trajectory_) {
    return false;
  }
  return trajectory_->renameElement(old_symbol, new_symbol);
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
    correlation::readers::BaseReader *reader = nullptr;
    std::string const ext = std::filesystem::path(path).extension().string();
    if (!ext.empty()) {
      reader = (registry_ != nullptr)
                   ? registry_->getReaderForExtension({.extension = ext, .filename = path})
                   : correlation::readers::ReaderFactory::instance().getReaderForExtension(
                         {.extension = ext, .filename = path});
      if (reader != nullptr) {
        is_trajectory = reader->isTrajectory();
      } else {
        is_trajectory = (type == correlation::readers::FileType::Xdatcar ||
                         type == correlation::readers::FileType::LammpsDump);
      }
    } else {
      is_trajectory = (type == correlation::readers::FileType::Xdatcar ||
                       type == correlation::readers::FileType::LammpsDump);
    }

    if (reader != nullptr) {
      if (is_trajectory) {
        trajectory_ =
            std::make_shared<correlation::core::Trajectory>(reader->readTrajectory(path, progress));
      } else {
        trajectory_ = std::make_shared<correlation::core::Trajectory>();
        trajectory_->addFrame(reader->readStructure(path, progress));
      }
    } else {
      if (is_trajectory) {
        trajectory_ = std::make_shared<correlation::core::Trajectory>(
            correlation::readers::readTrajectory(path, type, progress));
      } else {
        trajectory_ = std::make_shared<correlation::core::Trajectory>();
        trajectory_->addFrame(correlation::readers::readStructure(path, type, progress));
      }
    }

    return "File loaded: " + display_path;
  } catch (const std::exception &e) {
    return std::unexpected(std::string("Error loading file: ") + e.what());
  }
}

std::expected<std::string, std::string>
TrajectoryLoader::loadFiles(const std::vector<std::string> &paths,
                            const ProgressCallback &progress) {
  try {
    if (paths.empty()) {
      return std::unexpected("No files provided to load.");
    }
    if (paths.size() == 1) {
      return loadFile(paths.front(), progress);
    }
    std::vector<std::string> sorted_paths = paths;
    correlation::readers::naturalSort(sorted_paths);
    const auto &first_path = sorted_paths.front();
    correlation::readers::FileType const type = correlation::readers::determineFileType(first_path);
    trajectory_ = std::make_shared<correlation::core::Trajectory>(
        correlation::readers::readTrajectory(sorted_paths, type, progress));
    return "Loaded " + std::to_string(sorted_paths.size()) + " files.";
  } catch (const std::exception &e) {
    return std::unexpected(std::string("Error loading files: ") + e.what());
  }
}

} // namespace correlation::app
