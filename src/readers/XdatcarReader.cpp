/**
 * @file XdatcarReader.cpp
 * @brief Implementation of the VASP XDATCAR trajectory file reader with
 *        memory-mapped lazy loading.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */
#include "readers/XdatcarReader.hpp"
#include "core/MappedFile.hpp"
#include "readers/ReaderFactory.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace correlation::readers {

namespace {

// Automatic registration
// NOLINTNEXTLINE(cert-err58-cpp, bugprone-throwing-static-initialization)
const bool registered = ReaderFactory::instance().registerReader(std::make_unique<XdatcarReader>());

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
inline size_t findLineEnd(const char *data, size_t total, size_t pos) {
  while (pos < total && data[pos] != '\n' && data[pos] != '\r') {
    ++pos;
  }
  return pos;
}

inline size_t skipLineEnding(const char *data, size_t total, size_t pos) {
  if (pos < total && data[pos] == '\r') {
    ++pos;
  }
  if (pos < total && data[pos] == '\n') {
    ++pos;
  }
  return pos;
}

inline std::string extractLine(const char *data, size_t pos, size_t lineEnd) {
  return std::string(data + pos, lineEnd - pos);
}

/**
 * @brief Holds the pre-parsed XDATCAR header information needed to parse
 *        individual frames without re-reading the header.
 */
struct XdatcarHeader {
  std::array<std::array<real_t, 3>, 3> lattice = {}; ///< Scaled lattice vectors.
  std::vector<std::string> species;                  ///< Element symbols.
  std::vector<int> atom_counts;                      ///< Count per species.
  int total_atoms{0};                                ///< Sum of atom_counts.
  std::vector<std::string> atom_species;             ///< Per-atom species assignment.
};

struct XdatcarParser {
  const char *data = nullptr;
  size_t total_size = 0;
  size_t offset = 0;
  size_t lineEnd = 0;

  XdatcarParser(const char *data, size_t total_size) : data(data), total_size(total_size) {}

  std::string nextLine() {
    lineEnd = findLineEnd(data, total_size, offset);
    std::string line = extractLine(data, offset, lineEnd);
    offset = skipLineEnding(data, total_size, lineEnd);
    return line;
  }

  std::shared_ptr<XdatcarHeader> parseHeader() {
    auto header = std::make_shared<XdatcarHeader>();

    // Line 1: Comment
    nextLine();

    // Line 2: Scaling factor
    std::string line = nextLine();
    real_t const scaling_factor = static_cast<real_t>(std::stod(line));

    // Lines 3-5: Lattice vectors
    for (int i = 0; i < 3; ++i) {
      line = nextLine();
      std::istringstream iss(line);
      if (!(iss >> header->lattice.at(i).at(0) >> header->lattice.at(i).at(1) >> header->lattice.at(i).at(2))) {
        throw std::runtime_error("XDATCAR: failed to parse lattice vector on line " + std::to_string(i + 3) + ".");
      }
    }

    // Apply scaling factor
    applyScaling(header, scaling_factor);

    // Line 6: Species names
    line = nextLine();
    {
      std::istringstream iss(line);
      std::string tok;
      while (iss >> tok) {
        header->species.push_back(tok);
      }
    }

    // Line 7: Atom counts
    line = nextLine();
    {
      std::istringstream iss(line);
      int count = 0;
      while (iss >> count) {
        header->atom_counts.push_back(count);
      }
    }

    validateAndBuildSpecies(header);

    return header;
  }

  std::vector<size_t> findFrameOffsets(const std::function<void(float, const std::string &)> &progress_callback) {
    std::vector<size_t> frame_offsets;

    while (offset < total_size) {
      // Skip empty lines
      size_t const line_end = findLineEnd(data, total_size, offset);
      std::string trimmed = extractLine(data, offset, line_end);

      // Trim leading whitespace
      size_t const start = trimmed.find_first_not_of(" \t");
      if (start == std::string::npos) {
        offset = skipLineEnding(data, total_size, line_end);
        continue;
      }
      trimmed = trimmed.substr(start);

      // Detect "Direct" keyword (case-insensitive check on first 6 chars)
      if (trimmed.size() >= 6) {
        std::string prefix = trimmed.substr(0, 6);
        std::ranges::transform(prefix, prefix.begin(), ::tolower);
        if (prefix == "direct") {
          // The frame data starts AFTER the "Direct" line — the frame_offset
          // points to the "Direct" line itself so the parser knows where it is.
          frame_offsets.push_back(offset);
        }
      }

      offset = skipLineEnding(data, total_size, line_end);

      // Report progress
      if (progress_callback && total_size > 0) {
        float const progress = static_cast<float>(offset) / static_cast<float>(total_size);
        progress_callback(progress, "Scanning XDATCAR frames...");
      }
    }

    return frame_offsets;
  }

private:
  static void applyScaling(const std::shared_ptr<XdatcarHeader> &header, real_t scaling_factor) {
    if (scaling_factor > 0.0) {
      for (auto &row : header->lattice) {
        for (int j = 0; j < 3; ++j) {
          row.at(j) *= scaling_factor;
        }
      }
    } else if (scaling_factor < 0.0) {
      real_t const target_volume = std::abs(scaling_factor);
      const auto &lattice_vecs = header->lattice;
      real_t const current_volume =
          std::abs(lattice_vecs.at(0).at(0) * (lattice_vecs.at(1).at(1) * lattice_vecs.at(2).at(2) -
                                               lattice_vecs.at(1).at(2) * lattice_vecs.at(2).at(1)) -
                   lattice_vecs.at(0).at(1) * (lattice_vecs.at(1).at(0) * lattice_vecs.at(2).at(2) -
                                               lattice_vecs.at(1).at(2) * lattice_vecs.at(2).at(0)) +
                   lattice_vecs.at(0).at(2) * (lattice_vecs.at(1).at(0) * lattice_vecs.at(2).at(1) -
                                               lattice_vecs.at(1).at(1) * lattice_vecs.at(2).at(0)));
      real_t const scale = std::cbrt(target_volume / current_volume);
      for (auto &row : header->lattice) {
        for (int j = 0; j < 3; ++j) {
          row.at(j) *= scale;
        }
      }
    }
  }

  void validateAndBuildSpecies(const std::shared_ptr<XdatcarHeader> &header) const {
    if (header->species.size() != header->atom_counts.size()) {
      throw std::runtime_error("XDATCAR: species count does not match atom count entries.");
    }

    long long total_atoms_sum = 0;
    for (int const count : header->atom_counts) {
      if (count < 0) {
        throw std::runtime_error("XDATCAR: negative atom count: " + std::to_string(count));
      }
      total_atoms_sum += count;
    }

    constexpr int kMaxAtomCount = 100'000'000;
    if (total_atoms_sum > kMaxAtomCount) {
      throw std::runtime_error("XDATCAR: total atom count exceeds limit: " + std::to_string(total_atoms_sum));
    }

    // The number of atoms cannot exceed the file size in bytes
    if (std::cmp_greater(total_atoms_sum, total_size)) {
      throw std::runtime_error("XDATCAR: total atom count (" + std::to_string(total_atoms_sum) +
                               ") exceeds file size (" + std::to_string(total_size) + " bytes)");
    }

    header->total_atoms = static_cast<int>(total_atoms_sum);

    // Build per-atom species list
    header->atom_species.reserve(header->total_atoms);
    for (size_t species_idx = 0; species_idx < header->species.size(); ++species_idx) {
      for (int atom_idx = 0; atom_idx < header->atom_counts.at(species_idx); ++atom_idx) {
        header->atom_species.push_back(header->species.at(species_idx));
      }
    }
  }
};

correlation::core::Cell parseXdatcarFrame(const std::shared_ptr<XdatcarHeader> &header, const char *frame_data,
                                          size_t frame_size) {
  size_t offset = 0;
  size_t line_end = 0;

  auto nextLn = [&]() -> std::string {
    line_end = findLineEnd(frame_data, frame_size, offset);
    std::string line = extractLine(frame_data, offset, line_end);
    offset = skipLineEnding(frame_data, frame_size, line_end);
    return line;
  };

  // Skip the "Direct configuration= N" line
  nextLn();

  const auto &lattice_vecs = header->lattice;
  correlation::core::Cell cell({lattice_vecs.at(0).at(0), lattice_vecs.at(0).at(1), lattice_vecs.at(0).at(2)},
                               {lattice_vecs.at(1).at(0), lattice_vecs.at(1).at(1), lattice_vecs.at(1).at(2)},
                               {lattice_vecs.at(2).at(0), lattice_vecs.at(2).at(1), lattice_vecs.at(2).at(2)});

  for (int atom_idx = 0; atom_idx < header->total_atoms; ++atom_idx) {
    if (offset >= frame_size) {
      break;
    }
    std::string const line = nextLn();
    std::istringstream iss(line);
    real_t pos_x = 0.0;
    real_t pos_y = 0.0;
    real_t pos_z = 0.0;
    if (!(iss >> pos_x >> pos_y >> pos_z)) {
      break;
    }
    // Convert fractional coordinates to Cartesian: pos = x*a + y*b + z*c
    correlation::math::Vector3<real_t> const pos = {
        lattice_vecs.at(0).at(0) * pos_x + lattice_vecs.at(1).at(0) * pos_y + lattice_vecs.at(2).at(0) * pos_z,
        lattice_vecs.at(0).at(1) * pos_x + lattice_vecs.at(1).at(1) * pos_y + lattice_vecs.at(2).at(1) * pos_z,
        lattice_vecs.at(0).at(2) * pos_x + lattice_vecs.at(1).at(2) * pos_y + lattice_vecs.at(2).at(2) * pos_z};
    cell.addAtom(header->atom_species.at(atom_idx), pos);
  }

  cell.wrapPositions(); // XDATCAR always uses Direct coordinates
  return cell;
}

} // namespace

// ---------------------------------------------------------------------------
// readStructure
// ---------------------------------------------------------------------------
correlation::core::Cell
XdatcarReader::readStructure(const std::string &filename,
                             std::function<void(float, const std::string &)> progress_callback) {
  auto traj = readTrajectory(filename, progress_callback);
  if (traj.getFrameCount() == 0) {
    throw std::runtime_error("XDATCAR: no frames found in file.");
  }
  return traj.getFrame(0);
}

// ---------------------------------------------------------------------------
// readTrajectory — memory-mapped lazy loading
// ---------------------------------------------------------------------------
correlation::core::Trajectory
XdatcarReader::readTrajectory(const std::string &filename,
                              std::function<void(float, const std::string &)> progress_callback) {

  if (progress_callback) {
    progress_callback(0.0F, "Reading XDATCAR file...");
  }

  auto mapped_file = std::make_shared<correlation::core::MappedFile>(filename);
  const char *data = mapped_file->data();
  const size_t total_size = mapped_file->size();

  XdatcarParser parser(data, total_size);

  auto header = parser.parseHeader();
  std::vector<size_t> frame_offsets = parser.findFrameOffsets(progress_callback);

  if (frame_offsets.empty()) {
    throw std::runtime_error("XDATCAR: no valid frames found in file.");
  }

  // Add sentinel
  frame_offsets.push_back(total_size);

  if (progress_callback) {
    progress_callback(1.0F, "XDATCAR file loaded.");
  }

  // Build the parser lambda that captures the shared header.
  auto frame_parser = [header](const char *data_begin, size_t data_size) -> correlation::core::Cell {
    return parseXdatcarFrame(header, data_begin, data_size);
  };

  return {mapped_file, std::move(frame_offsets), frame_parser, 1.0};
}

} // namespace correlation::readers
