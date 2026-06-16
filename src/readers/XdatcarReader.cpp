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
#include <cmath>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace correlation::readers {

// Automatic registration
static const bool registered = ReaderFactory::instance().registerReader(std::make_unique<XdatcarReader>()); // NOLINT(cert-err58-cpp, bugprone-throwing-static-initialization)

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
static inline size_t findLineEnd(const char *data, size_t total, size_t pos) {
  while (pos < total && data[pos] != '\n' && data[pos] != '\r') {
    ++pos;
}
  return pos;
}

static inline size_t skipLineEnding(const char *data, size_t total, size_t pos) {
  if (pos < total && data[pos] == '\r') {
    ++pos;
}
  if (pos < total && data[pos] == '\n') {
    ++pos;
}
  return pos;
}

static inline std::string extractLine(const char *data, size_t pos, size_t lineEnd) {
  return std::string(data + pos, lineEnd - pos);
}

/**
 * @brief Holds the pre-parsed XDATCAR header information needed to parse
 *        individual frames without re-reading the header.
 */
struct XdatcarHeader {
  double lattice[3][3]{};                  ///< Scaled lattice vectors.
  std::vector<std::string> species;      ///< Element symbols.
  std::vector<int> atom_counts;          ///< Count per species.
  int total_atoms{0};                    ///< Sum of atom_counts.
  std::vector<std::string> atom_species; ///< Per-atom species assignment.
};

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

  size_t offset = 0;
  size_t lineEnd = 0;

  auto nextLine = [&]() -> std::string {
    lineEnd = findLineEnd(data, total_size, offset);
    std::string line = extractLine(data, offset, lineEnd);
    offset = skipLineEnding(data, total_size, lineEnd);
    return line;
  };

  // --- Parse the header (lines 1-7) ---
  auto header = std::make_shared<XdatcarHeader>();

  // Line 1: Comment
  nextLine();

  // Line 2: Scaling factor
  std::string line = nextLine();
  double const scaling_factor = std::stod(line);

  // Lines 3-5: Lattice vectors
  for (int i = 0; i < 3; ++i) {
    line = nextLine();
    std::istringstream iss(line);
    if (!(iss >> header->lattice[i][0] >> header->lattice[i][1] >> header->lattice[i][2])) {
      throw std::runtime_error("XDATCAR: failed to parse lattice vector on line " + std::to_string(i + 3) + ".");
    }
  }

  // Apply scaling factor
  if (scaling_factor > 0.0) {
    for (auto & i : header->lattice) {
      for (int j = 0; j < 3; ++j) {
        i[j] *= scaling_factor;
}
}
  } else if (scaling_factor < 0.0) {
    double const target_volume = std::abs(scaling_factor);
    const auto &v = header->lattice;
    double const current_volume =
        std::abs(v[0][0] * (v[1][1] * v[2][2] - v[1][2] * v[2][1]) - v[0][1] * (v[1][0] * v[2][2] - v[1][2] * v[2][0]) +
                 v[0][2] * (v[1][0] * v[2][1] - v[1][1] * v[2][0]));
    double const scale = std::cbrt(target_volume / current_volume);
    for (auto & i : header->lattice) {
      for (int j = 0; j < 3; ++j) {
        i[j] *= scale;
}
}
  }

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

  if (header->species.size() != header->atom_counts.size()) {
    throw std::runtime_error("XDATCAR: species count does not match atom count entries.");
  }

  long long total_atoms_sum = 0;
  for (int const c : header->atom_counts) {
    if (c < 0) {
      throw std::runtime_error("XDATCAR: negative atom count: " + std::to_string(c));
    }
    total_atoms_sum += c;
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
  for (size_t s = 0; s < header->species.size(); ++s) {
    for (int a = 0; a < header->atom_counts[s]; ++a) {
      header->atom_species.push_back(header->species[s]);
    }
  }


  // --- Scan for "Direct" lines to find frame offsets ---
  std::vector<size_t> frame_offsets;

  while (offset < total_size) {
    // Skip empty lines
    size_t const le = findLineEnd(data, total_size, offset);
    std::string trimmed = extractLine(data, offset, le);

    // Trim leading whitespace
    size_t const start = trimmed.find_first_not_of(" \t");
    if (start == std::string::npos) {
      offset = skipLineEnding(data, total_size, le);
      continue;
    }
    trimmed = trimmed.substr(start);

    // Detect "Direct" keyword (case-insensitive check on first 6 chars)
    if (trimmed.size() >= 6) {
      std::string prefix = trimmed.substr(0, 6);
      std::transform(prefix.begin(), prefix.end(), prefix.begin(), ::tolower);
      if (prefix == "direct") {
        // The frame data starts AFTER the "Direct" line — the frame_offset
        // points to the "Direct" line itself so the parser knows where it is.
        frame_offsets.push_back(offset);
      }
    }

    offset = skipLineEnding(data, total_size, le);

    // Report progress
    if (progress_callback && total_size > 0) {
      float const progress = static_cast<float>(offset) / static_cast<float>(total_size);
      progress_callback(progress, "Scanning XDATCAR frames...");
    }
  }

  if (frame_offsets.empty()) {
    throw std::runtime_error("XDATCAR: no valid frames found in file.");
  }

  // Add sentinel
  frame_offsets.push_back(total_size);

  if (progress_callback) {
    progress_callback(1.0F, "XDATCAR file loaded.");
}

  // Build the parser lambda that captures the shared header.
  auto parser = [header](const char *d, size_t s) -> correlation::core::Cell {
    size_t off = 0;
    size_t le = 0;

    auto nextLn = [&]() -> std::string {
      le = findLineEnd(d, s, off);
      std::string ln = extractLine(d, off, le);
      off = skipLineEnding(d, s, le);
      return ln;
    };

    // Skip the "Direct configuration= N" line
    nextLn();

    const auto &v = header->lattice;
    correlation::core::Cell cell({v[0][0], v[0][1], v[0][2]}, {v[1][0], v[1][1], v[1][2]}, {v[2][0], v[2][1], v[2][2]});

    for (int i = 0; i < header->total_atoms; ++i) {
      if (off >= s) {
        break;
}
      std::string const line = nextLn();
      std::istringstream iss(line);
      double x;
      double y;
      double z;
      if (!(iss >> x >> y >> z)) {
        break;
}
      // Convert fractional coordinates to Cartesian: pos = x*a + y*b + z*c
      correlation::math::Vector3<double> const pos = {x * v[0][0] + y * v[1][0] + z * v[2][0],
                                                x * v[0][1] + y * v[1][1] + z * v[2][1],
                                                x * v[0][2] + y * v[1][2] + z * v[2][2]};
      cell.addAtom(header->atom_species[i], pos);
    }

    cell.wrapPositions(); // XDATCAR always uses Direct coordinates
    return cell;
  };

  return {mapped_file, std::move(frame_offsets), parser, 1.0};
}

} // namespace correlation::readers
