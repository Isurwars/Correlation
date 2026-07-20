/**
 * @file VaspReader.cpp
 * @brief Implementation of the VASP POSCAR/CONTCAR file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */
#include "readers/VaspReader.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "readers/ReaderFactory.hpp"

#include <array>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <fstream>
#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace correlation::readers {

namespace {

// Automatic registration
// NOLINTNEXTLINE(cert-err58-cpp, bugprone-throwing-static-initialization)
const bool registered = ReaderFactory::instance().registerReader(std::make_unique<VaspReader>());

struct VaspParser {
  std::ifstream *file = nullptr;

  explicit VaspParser(std::ifstream &file) : file(&file) {}

  correlation::core::Cell parse() const {
    std::string line;

    // Line 1: Comment
    if (!std::getline(*file, line)) {
      throw std::runtime_error("POSCAR: unexpected end of file (comment line).");
    }

    // Line 2: Scaling factor
    if (!std::getline(*file, line)) {
      throw std::runtime_error("POSCAR: unexpected end of file (scaling factor).");
    }
    real_t const scaling_factor = static_cast<real_t>(std::stod(line));

    // Lines 3-5: Lattice vectors
    auto lattice = parseLatticeVectors();
    applyScalingFactor(lattice, scaling_factor);

    correlation::core::Cell tempCell({lattice.at(0).at(0), lattice.at(0).at(1), lattice.at(0).at(2)},
                                     {lattice.at(1).at(0), lattice.at(1).at(1), lattice.at(1).at(2)},
                                     {lattice.at(2).at(0), lattice.at(2).at(1), lattice.at(2).at(2)});

    // Species names and atom counts
    auto [species, atom_counts] = parseSpeciesAndCounts();
    int const total_atoms = calculateTotalAtoms(atom_counts);

    // Read coordinate type
    bool const is_direct = parseCoordinateType();

    // Read atom positions
    parseAtomPositions(tempCell, species, atom_counts, total_atoms, is_direct);

    tempCell.wrapPositions();
    return tempCell;
  }

private:
  std::array<std::array<real_t, 3>, 3> parseLatticeVectors() const {
    std::array<std::array<real_t, 3>, 3> lattice_vectors = {};
    std::string line;
    for (int i = 0; i < 3; ++i) {
      if (!std::getline(*file, line)) {
        throw std::runtime_error("POSCAR: unexpected end of file (lattice vector).");
      }
      std::istringstream iss(line);
      if (!(iss >> lattice_vectors.at(i).at(0) >> lattice_vectors.at(i).at(1) >> lattice_vectors.at(i).at(2))) {
        throw std::runtime_error("POSCAR: failed to parse lattice vector on line " + std::to_string(i + 3) + ".");
      }
    }
    return lattice_vectors;
  }

  static void applyScalingFactor(std::array<std::array<real_t, 3>, 3> &lattice_vectors, real_t scaling_factor) {
    if (scaling_factor > 0.0) {
      for (auto &row : lattice_vectors) {
        for (real_t &val : row) {
          val *= scaling_factor;
        }
      }
    } else if (scaling_factor < 0.0) {
      real_t const target_volume = std::abs(scaling_factor);
      real_t const current_volume =
          std::abs(lattice_vectors.at(0).at(0) * (lattice_vectors.at(1).at(1) * lattice_vectors.at(2).at(2) -
                                                  lattice_vectors.at(1).at(2) * lattice_vectors.at(2).at(1)) -
                   lattice_vectors.at(0).at(1) * (lattice_vectors.at(1).at(0) * lattice_vectors.at(2).at(2) -
                                                  lattice_vectors.at(1).at(2) * lattice_vectors.at(2).at(0)) +
                   lattice_vectors.at(0).at(2) * (lattice_vectors.at(1).at(0) * lattice_vectors.at(2).at(1) -
                                                  lattice_vectors.at(1).at(1) * lattice_vectors.at(2).at(0)));
      real_t const scale = std::cbrt(target_volume / current_volume);
      for (auto &row : lattice_vectors) {
        for (real_t &val : row) {
          val *= scale;
        }
      }
    }
  }

  std::pair<std::vector<std::string>, std::vector<int>> parseSpeciesAndCounts() const {
    std::string line;
    if (!std::getline(*file, line)) {
      throw std::runtime_error("POSCAR: unexpected end of file (species/counts).");
    }

    std::vector<std::string> species;
    std::vector<int> atom_counts;
    std::istringstream iss(line);
    std::string token;
    bool is_species_line = false;

    if (iss >> token) {
      try {
        static_cast<void>(std::stoi(token));
        // VASP 4 format (no species line)
        atom_counts.push_back(std::stoi(token));
        while (iss >> token) {
          atom_counts.push_back(std::stoi(token));
        }
        for (size_t i = 0; i < atom_counts.size(); ++i) {
          species.push_back("Type" + std::to_string(i + 1));
        }
      } catch (...) {
        // VASP 5+ format (species names)
        is_species_line = true;
        species.push_back(token);
        while (iss >> token) {
          species.push_back(token);
        }
      }
    }

    if (is_species_line) {
      if (!std::getline(*file, line)) {
        throw std::runtime_error("POSCAR: unexpected end of file (atom counts).");
      }
      std::istringstream count_iss(line);
      int count = 0;
      while (count_iss >> count) {
        atom_counts.push_back(count);
      }
    }

    if (species.size() != atom_counts.size()) {
      throw std::runtime_error("POSCAR: species count (" + std::to_string(species.size()) +
                               ") does not match atom count entries (" + std::to_string(atom_counts.size()) + ").");
    }

    return {species, atom_counts};
  }

  static int calculateTotalAtoms(const std::vector<int> &atom_counts) {
    long long total_atoms_sum = 0;
    for (int const count : atom_counts) {
      if (count < 0) {
        throw std::runtime_error("POSCAR: negative atom count encountered: " + std::to_string(count));
      }
      total_atoms_sum += count;
    }
    constexpr int kMaxAtomCount = 100'000'000;
    if (total_atoms_sum > kMaxAtomCount) {
      throw std::runtime_error("POSCAR: total atom count exceeds limit: " + std::to_string(total_atoms_sum));
    }
    return static_cast<int>(total_atoms_sum);
  }

  [[nodiscard]] bool parseCoordinateType() const {
    std::string line;
    if (!std::getline(*file, line)) {
      throw std::runtime_error("POSCAR: unexpected end of file (coordinate type).");
    }

    char first_char = getFirstNonSpaceChar(line);

    if (first_char == 'S' || first_char == 's') {
      if (!std::getline(*file, line)) {
        throw std::runtime_error("POSCAR: unexpected end of file (coordinate type after selective dynamics).");
      }
      first_char = getFirstNonSpaceChar(line);
    }

    return (first_char == 'D' || first_char == 'd');
  }

  static char getFirstNonSpaceChar(const std::string &line) {
    for (char const chr : line) {
      if (std::isspace(static_cast<unsigned char>(chr)) == 0) {
        return chr;
      }
    }
    return ' ';
  }

  void parseAtomPositions(correlation::core::Cell &tempCell, const std::vector<std::string> &species,
                          const std::vector<int> &atom_counts, int total_atoms, bool is_direct) const {
    int species_idx = 0;
    int atoms_in_species = 0;
    const auto &lattice_vectors = tempCell.latticeVectors();
    std::string line;

    for (int i = 0; i < total_atoms; ++i) {
      if (!std::getline(*file, line)) {
        throw std::runtime_error("POSCAR: unexpected end of file (atom position " + std::to_string(i + 1) + " of " +
                                 std::to_string(total_atoms) + ").");
      }

      while (std::cmp_less(species_idx, atom_counts.size()) && atoms_in_species >= atom_counts.at(species_idx)) {
        atoms_in_species = 0;
        species_idx++;
      }

      std::istringstream iss(line);
      real_t pos_x = 0.0;
      real_t pos_y = 0.0;
      real_t pos_z = 0.0;
      if (!(iss >> pos_x >> pos_y >> pos_z)) {
        throw std::runtime_error("POSCAR: failed to parse atom coordinates on atom " + std::to_string(i + 1) + ".");
      }

      correlation::math::Vector3<real_t> pos;
      if (is_direct) {
        pos = {pos_x * lattice_vectors[0][0] + pos_y * lattice_vectors[1][0] + pos_z * lattice_vectors[2][0],
               pos_x * lattice_vectors[0][1] + pos_y * lattice_vectors[1][1] + pos_z * lattice_vectors[2][1],
               pos_x * lattice_vectors[0][2] + pos_y * lattice_vectors[1][2] + pos_z * lattice_vectors[2][2]};
      } else {
        pos = {pos_x, pos_y, pos_z};
      }
      tempCell.addAtom(species.at(species_idx), pos);
      atoms_in_species++;
    }
  }
};

} // namespace

correlation::core::Cell
VaspReader::readStructure(const std::string &filename,
                          std::function<void(float, const std::string &)> /*progress_callback*/) {
  return read(filename);
}

correlation::core::Trajectory
VaspReader::readTrajectory(const std::string & /*filename*/,
                           std::function<void(float, const std::string &)> /*progress_callback*/) {
  throw std::runtime_error("POSCAR/CONTCAR files are single structures, use readStructure.");
}

correlation::core::Cell VaspReader::read(const std::string &file_name) {
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" + std::strerror(errno) + ").");
  }

  VaspParser parser(myfile);
  return parser.parse();
}

} // namespace correlation::readers
