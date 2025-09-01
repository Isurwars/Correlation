// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "../include/Cell.hpp"

#include <algorithm>
#include <cmath>
#include <execution>
#include <iostream>
#include <mutex>
#include <ostream>
#include <utility>

#include "../include/Atom.hpp"
#include "../include/Constants.hpp"
#include "../include/LinearAlgebra.hpp"

//---------------------------------------------------------------------------//
//---------------------------- Cell Constructors ----------------------------//
//---------------------------------------------------------------------------//

// Default constructor
Cell::Cell()
    : lattice_parameters_{1.0, 1.0, 1.0, 90.0, 90.0, 90.0}, atoms_{},
      elements_{}, element_numbers_{}, volume_{0.0} {
  calculateLatticeVectors();
  calculateVolume();
}

// Lattice parameters constructor
Cell::Cell(const std::array<double, 6> &lat)
    : atoms_(), elements_(), element_numbers_(), volume_(0.0) {
  setLatticeParameters(lat);
  calculateLatticeVectors();
  calculateVolume();
}

// Lattice vector constructor
Cell::Cell(const linalg::Vector3<double> &v1, const linalg::Vector3<double> &v2,
           const linalg::Vector3<double> &v3)
    : atoms_{}, elements_{}, element_numbers_{}, volume_{0.0} {

  const double a = norm(v1);
  const double b = norm(v2);
  const double c = norm(v3);

  if (a <= 0.0 || b <= 0.0 || c <= 0.0) {
    throw std::invalid_argument(
        "Invalid lattice vector: zero length detected.");
  }

  // Compute angles (degrees) between lattice vectors
  const double alpha = std::acos((v2 * v3) / (b * c)) * constants::rad2deg;
  const double beta = std::acos((v1 * v3) / (a * c)) * constants::rad2deg;
  const double gamma = std::acos((v1 * v2) / (a * b)) * constants::rad2deg;

  lattice_parameters_ = {a, b, c, alpha, beta, gamma};
  v_a_ = v1;
  v_b_ = v2;
  v_c_ = v3;
  calculateVolume();
}

//---------------------------------------------------------------------------//
//--------------------------------- Methods ---------------------------------//
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
//------------------------- Set Lattice Parameters --------------------------//
//---------------------------------------------------------------------------//
void Cell::setLatticeParameters(std::array<double, 6> lat) {
  lattice_parameters_ = lat;
  calculateLatticeVectors();
  calculateVolume();
} // Cell::setLatticeparameters

//---------------------------------------------------------------------------//
//------------------------ Calculate Lattice Vectors ------------------------//
//---------------------------------------------------------------------------//
void Cell::calculateLatticeVectors() {
  // Calculate the lattice vectors from the lattice parameters.
  const std::array<double, 6> &lat = lattice_parameters();
  double A = lat[0];
  double B = lat[1];
  double C = lat[2];
  double alpha = lat[3] * constants::deg2rad;
  double beta = lat[4] * constants::deg2rad;
  double gamma = lat[5] * constants::deg2rad;

  double c_a = cos(alpha);
  double c_b = cos(beta);
  double c_g = cos(gamma);
  double s_g = sin(gamma);

  const double term_under_sqrt =
      1 - c_a * c_a - c_b * c_b - c_g * c_g + 2 * c_a * c_b * c_g;

  if (term_under_sqrt <
      1e-12) { // Use a small epsilon to handle floating-point inaccuracies
    throw std::invalid_argument(
        "Invalid lattice parameters; non-positive volume.");
  }

  v_a_ = {A, 0.0, 0.0};
  v_b_ = {B * c_g, B * s_g, 0.0};
  v_c_ = {C * c_b, C * (c_a - c_b * c_g) / s_g,
          C * std::sqrt(term_under_sqrt) / s_g};
  calculateVolume();
} // Cell::calculateLatticeVectors

//---------------------------------------------------------------------------//
//----------------------------- Calculate Volume ----------------------------//
//---------------------------------------------------------------------------//
void Cell::calculateVolume() {
  volume_ = (v_a_ * cross(v_b_, v_c_));

  if (volume_ <= 1e-8)
    throw std::logic_error("Cell volume must be positive");
} // Cell:calculateVolume

//---------------------------------------------------------------------------//
//---------------------------------  Atoms  ---------------------------------//
//---------------------------------------------------------------------------//
void Cell::addAtom(Atom at) {
  at.setID(atoms_.size());
  addElement(at.element());
  atoms_.push_back(at);
} // Cell::addAtom

void Cell::setAtoms(std::vector<Atom> ats) {
  atoms_.clear();
  for (auto &at : ats) {
    at.setID(atoms_.size());
    addAtom(at);
  }
} // Cell::setAtoms

//---------------------------------------------------------------------------//
//-------------------------------- Elements ---------------------------------//
//---------------------------------------------------------------------------//
void Cell::addElement(const std::string &ele) {
  // Guard for Empty element
  if (ele.empty())
    throw std::invalid_argument("Element name cannot be empty");

  // Find insertion position while maintaining sorted order
  auto it = std::lower_bound(elements_.begin(), elements_.end(), ele);

  // Only insert if element doesn't already exist
  if (it == elements_.end() || *it != ele)
    elements_.insert(it, ele);
} // Cell::addElement

const int Cell::element_id(const std::string &ele) {
  auto it = std::find(elements_.begin(), elements_.end(), ele);
  return (it != elements_.end()) ? std::distance(elements_.begin(), it) : -1;
} // Cell::element_id

void Cell::populateElementID() {

  // Create a lookup map for O(1) element ID access
  std::unordered_map<std::string, int> element_index_map;
  element_index_map.reserve(elements_.size());

  // Build element index map once
  for (size_t i = 0; i < elements_.size(); ++i) {
    element_index_map[elements_[i]] = i;
  }

  // Assign element IDs using the map
  for (Atom &atom : atoms_) {
    const auto it = element_index_map.find(atom.element());

    if (it == element_index_map.end()) {
      throw std::runtime_error("Missing element: " + atom.element());
    }

    atom.setElementID(it->second);
  }
} // Cell::populateElementID

void Cell::calculateElementNumbers() {
  std::vector<int> temp_num_atoms(elements().size(), 0);
  populateElementID();
  for (auto &atom : atoms_) {
    temp_num_atoms[atom.element_id()]++;
  }
  setElementsNumbers(temp_num_atoms);
} // Cell::populateElementNumbers

//---------------------------------------------------------------------------//
//---------------------------  Correct Positions ----------------------------//
//---------------------------------------------------------------------------//
void Cell::correctPositions() {

  // Create lattice matrix
  linalg::Matrix3<double> L = {{v_a_.x(), v_b_.x(), v_c_.x()},
                               {v_a_.y(), v_b_.y(), v_c_.y()},
                               {v_a_.z(), v_b_.z(), v_c_.z()}};
  linalg::Matrix3<double> L_inv = invertMatrix(L);

  for (auto &MyAtom : atoms_) {
    linalg::Vector3<double> pos = {MyAtom.position().x(), MyAtom.position().y(),
                                   MyAtom.position().z()};

    // Convert to fractional coordinates
    linalg::Vector3<double> frac = matrixVectorMultiply(L_inv, pos);

    // Wrap fractional coordinates to [0, 1)
    for (int m = 0; m < 3; ++m) {
      frac[m] -= std::floor(frac[m]);
    }

    // Convert back to Cartesian coordinates
    linalg::Vector3<double> corrected_pos = matrixVectorMultiply(L, frac);

    MyAtom.setPosition(
        {corrected_pos.x(), corrected_pos.y(), corrected_pos.z()});
  }
} // Cell::CorrectPositions

//---------------------------------------------------------------------------//
//----------------------- Correct Fractional Positions ----------------------//
//---------------------------------------------------------------------------//
void Cell::correctFracPositions() {
  linalg::Vector3<double> pos;

  for (auto &MyAtom : atoms_) {
    pos = MyAtom.position();
    double i = pos.x() - std::floor(pos.x());
    double j = pos.y() - std::floor(pos.y());
    double k = pos.z() - std::floor(pos.z());
    pos = i * v_a_ + j * v_b_ + k * v_c_;

    MyAtom.setPosition(pos);
  }
} // Cell::CorrectFracPositions

//---------------------------------------------------------------------------//
//-------------------------- Populate Bond Length ---------------------------//
//---------------------------------------------------------------------------//
void Cell::populateBondLength(double Bond_Factor) {
  // Number of elements in the Cell
  const int n = elements().size();
  // Initialize Bond length matrix (nxn) as zeros
  std::vector<std::vector<double>> temp_matrix(n, std::vector<double>(n, 0.0));

  for (int i = 0; i < n; i++) {
    double aux = covalentRadii(elements()[i]);
    for (int j = 0; j < n; j++) {
      temp_matrix[i][j] = (aux + covalentRadii(elements()[j])) * Bond_Factor;
    }
  }
  bond_length_ = temp_matrix;
} // Cell::populateBondlength

//---------------------------------------------------------------------------//
//--------------------------- Distance Population ---------------------------//
//---------------------------------------------------------------------------//
void Cell::distancePopulation(double r_cut, bool self_interaction) {
  // Validate inputs
  if (r_cut <= 0)
    throw std::invalid_argument("r_cut must be positive");
  if (v_a_.empty() || v_b_.empty() || v_c_.empty())
    throw std::logic_error("Cell vectors not initialized");
  const size_t n = elements_.size();
  if (n == 0)
    throw std::logic_error("No elements in cell");
  if (atoms_.size() == 0)
    throw std::logic_error("No atoms in cell");

  std::vector<std::vector<std::vector<double>>> temp_dist(
      n, std::vector<std::vector<double>>(n));
  const double r_cut_sq = r_cut * r_cut;

  if (bond_length_.empty())
    populateBondLength(1.2);

  calculateElementNumbers();
  correctPositions();

  // Calculate supercell dimensions
  int k_ = ceil(r_cut / v_c_.z());
  int j_ = ceil(r_cut / v_b_.y());
  int i_ = ceil(r_cut / v_a_.x());

  // Force self_interaction in case of a small supercell
  if (k_ * j_ * i_ > 8)
    self_interaction = true;

  // Precompute displacements
  std::vector<linalg::Vector3<double>> displacements;
  for (int i = -i_; i <= i_; ++i) {
    for (int j = -j_; j <= j_; ++j) {
      for (int k = -k_; k <= k_; ++k) {
        displacements.push_back(i * v_a_ + j * v_b_ + k * v_c_);
      }
    }
  }

  // Thread-safe error handling
  std::vector<std::mutex> global_mutex(n);
  std::exception_ptr async_exception = nullptr;

  auto process_atom = [&](Atom &atom_A) {
    try {
      const int id_A = atom_A.element_id();

      for (Atom &atom_B : atoms_) {
        if (&atom_A == &atom_B && !self_interaction)
          continue;

        const int id_B = atom_B.element_id();
        const auto &pos_B = atom_B.position();

        for (const auto &disp : displacements) {
          const linalg::Vector3<double> img_pos = pos_B + disp;
          const double dist = norm(atom_A.position() - img_pos);

          // Ignore self-interaction
          if (dist == 0 && (atom_A.id() == atom_B.id()))
            continue;

          // Collitions checks
          if (dist < 0.1) {
            throw std::runtime_error(
                "Atomic collision detected: " + atom_A.element() + "_" +
                std::to_string(atom_A.id()) + " and " + atom_B.element() + "_" +
                std::to_string(atom_B.id()));
          }
          if (dist < 0.5) {
            std::cerr << "\nWARNING: The atoms:\n"
                      << atom_A.element() << "_" << atom_A.id()
                      << " in position (" << atom_A.position().x() << ", "
                      << atom_A.position().y() << ", " << atom_A.position().z()
                      << "),\n"
                      << atom_B.element() << "_" << atom_B.id()
                      << " in position (" << img_pos.x() << ", " << img_pos.y()
                      << ", " << img_pos.z() << ").\n"
                      << "Have a distance less than the Bohr Radius.\n";
          }

          // Thread-safe bond addition
          if (dist <= bond_length_[id_A][id_B]) {
            atom_A.addBondedAtom(
                Atom(atom_B.element(), img_pos, atom_B.id(), id_B));
          }

          // Thread-safe distance recording
          std::lock_guard<std::mutex> lock(global_mutex[id_A]);
          temp_dist[id_A][id_B].push_back(dist);
        }
      }
    } catch (...) {
      async_exception = std::current_exception();
    }
  };

  // Parallel execution with exception propagation
  std::for_each(std::execution::par, atoms_.begin(), atoms_.end(),
                process_atom);

  if (async_exception) {
    std::rethrow_exception(async_exception);
  }

  distances_ = std::move(temp_dist);
} // Cell::distancePopulation

//---------------------------------------------------------------------------//
//------------------------- Plane Angle Population --------------------------//
//---------------------------------------------------------------------------//
void Cell::planeAnglePopulation(bool degree) {
  const double conversion_factor = degree ? constants::rad2deg : 1.0;

  const size_t num_elements = elements().size();
  if (num_elements == 0)
    throw std::logic_error("No elements in cell");
  if (atoms_.size() == 0)
    throw std::logic_error("No atoms in cell");

  // Clear previous angles and reserve space
  // NxNxN Tensor to store the planeAnglePopulation
  using AngleStorage =
      std::vector<std::vector<std::vector<std::vector<double>>>>;
  AngleStorage temp_angles(
      num_elements,
      std::vector<std::vector<std::vector<double>>>(
          num_elements, std::vector<std::vector<double>>(num_elements)));

  /*
   * This is the main loop to calculate the angles formed by every three
   * atoms. The connected atoms are calculated in Cell::distancepopulation,
   * and MUST be called first.
   *
   * The first loop iterates over every atom in the cell.
   * The second and third loop iterate over the connected atoms in every atom
   * instance. These three loops populate a 3D tensor of vectors, the indices
   * of the Tensor represent the three indices of the element in
   * Cell::elements.
   *
   * By default the angle returned is given in degrees.
   */

  for (size_t atom_idx = 0; atom_idx < atoms_.size(); ++atom_idx) {
    Atom &central_atom = atoms_[atom_idx];
    std::vector<Atom> bonded = central_atom.bonded_atoms();
    const size_t central_elem = central_atom.element_id();

    // Use indices to avoid duplicate pairs
    for (size_t i = 0; i < bonded.size(); ++i) {
      Atom &atom_a = bonded[i];
      const size_t elem_a = atom_a.element_id();

      for (size_t j = i + 1; j < bonded.size(); ++j) {
        Atom &atom_b = bonded[j];
        const size_t elem_b = atom_b.element_id();

        // Calculate angle and store in both element order permutations
        try {
          const double angle =
              central_atom.angle(atom_a, atom_b) * conversion_factor;
          temp_angles[elem_a][central_elem][elem_b].push_back(angle);
          temp_angles[elem_b][central_elem][elem_a].push_back(angle);
        } catch (const std::invalid_argument &e) {
          std::cerr << "Angle calculation error: " << e.what() << " for atoms "
                    << atom_a.id() << " and " << atom_b.id() << '\n';
        }
      }
    }
  }
  angles_ = std::move(temp_angles);
} // Cell::planeAnglePopulation
