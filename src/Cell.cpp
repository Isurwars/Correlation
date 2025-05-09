/* ----------------------------------------------------------------------------
 * Correlation: An Analysis Tool for Liquids and for Amorphous Solids
 * Copyright (c) 2013-2025 Isaías Rodríguez <isurwars@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the MIT License version as published in:
 * https://github.com/Isurwars/Correlation/blob/main/LICENSE
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 * ----------------------------------------------------------------------------
 */
#include "../include/Cell.hpp"

#include <algorithm>
#include <cmath>
#include <execution>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <numeric>
#include <ostream>
#include <regex>
#include <utility>

#include "../include/Atom.hpp"
#include "../include/Constants.hpp"
#include "../include/LinearAlgebra.hpp"
#include "../include/Smoothing.hpp"
#include "../include/StructureFactor.hpp"
#include "../include/Templates.hpp"

//---------------------------------------------------------------------------//
//---------------------------- Cell Constructors ----------------------------//
//---------------------------------------------------------------------------//

// Default constructor
Cell::Cell() {
  this->setLatticeParameters({1.0, 1.0, 1.0, 90, 90, 90});
  this->setAtoms({});
  this->setElements({});
  this->setLatticeVectors();
}

// Lattice parameters constructor
Cell::Cell(std::array<double, 6> lat) {
  this->setLatticeParameters(lat);
  this->setAtoms({});
  this->setElements({});
  this->setLatticeVectors();
}

// Lattice vector constructor
void Cell::setFromVectors(std::vector<double> v1, std::vector<double> v2,
			  std::vector<double> v3) {
  double v1_n, v2_n, v3_n, aux, a, b, c;

  v1_n = sqrt(std::inner_product(v1.begin(), v1.end(), v1.begin(), 0));
  v2_n = sqrt(std::inner_product(v2.begin(), v2.end(), v2.begin(), 0));
  v3_n = sqrt(std::inner_product(v3.begin(), v3.end(), v3.begin(), 0));
  aux = std::inner_product(v2.begin(), v2.end(), v3.begin(), 0);
  a = acos(aux / (v2_n * v3_n)) * constants::rad2deg;
  aux = std::inner_product(v1.begin(), v1.end(), v3.begin(), 0);
  b = acos(aux / (v1_n * v3_n)) * constants::rad2deg;
  aux = std::inner_product(v2.begin(), v2.end(), v1.begin(), 0);
  c = acos(aux / (v2_n * v1_n)) * constants::rad2deg;
  this->setLatticeParameters({v1_n, v2_n, v3_n, a, b, c});
  this->setAtoms({});
  this->setElements({});
  this->setLatticeVectors();
}

//---------------------------------------------------------------------------//
//------------------------------ Cell Methods -------------------------------//
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
//--------------------------- Set Lattice Vectors ---------------------------//
//---------------------------------------------------------------------------//
void Cell::setLatticeVectors() {
  // Calculate the lattice vectors from the lattice parameters.
  const std::array<double, 6> &lat = this->lattice_parameters();
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

  this->_v_a_ = {A, 0.0, 0.0};
  this->_v_b_ = {B * c_g, B * s_g, 0.0};
  this->_v_c_ = {
      C * c_b, C * (c_a - c_b * c_g) / s_g,
      C * sqrt(1 - c_a * c_a - c_b * c_b - c_g * c_g + 2 * c_a * c_b * c_g) /
	  s_g};

  this->volume =
      this->_v_a_[0] *
	  (this->_v_b_[1] * this->_v_c_[2] - this->_v_b_[2] * this->_v_c_[1]) -
      this->_v_a_[1] *
	  (this->_v_b_[0] * this->_v_c_[2] - this->_v_b_[2] * this->_v_c_[0]) +
      this->_v_a_[2] *
	  (this->_v_b_[0] * this->_v_c_[1] - this->_v_b_[1] * this->_v_c_[0]);
} // Cell::setLatticeVectors

void Cell::setLatticeParameters(std::array<double, 6> lat) {
  this->_lattice_parameters_ = lat;
  this->setLatticeVectors();
} // Cell::setLatticeparameters

//---------------------------------------------------------------------------//
//---------------------------------  Atoms  ---------------------------------//
//---------------------------------------------------------------------------//
void Cell::addAtom(Atom at) {
  at.setID(this->_atoms_.size());
  auto ele_id = findIndex(this->_elements_, at.element());
  if (!ele_id.has_value())
    this->addElement(at.element());
  this->_atoms_.push_back(at);
} // Cell::addAtom

void Cell::setAtoms(std::vector<Atom> ats) {
  this->_atoms_.clear();
  for (auto &at : ats) {
    at.setID(this->_atoms_.size());
    this->_atoms_.push_back(at);
  }
} // Cell::setAtoms

//---------------------------------------------------------------------------//
//-------------------------------- Elements ---------------------------------//
//---------------------------------------------------------------------------//
void Cell::populateElementID() {

  // Create a lookup map for O(1) element ID access
  std::unordered_map<std::string, int> element_index_map;
  element_index_map.reserve(_elements_.size());

  // Build element index map once
  for (size_t i = 0; i < this->_elements_.size(); ++i) {
    element_index_map[this->_elements_[i]] = i;
  }

  // Assign element IDs using the map
  for (Atom &atom : this->_atoms_) {
    const auto &element = atom.element();
    const auto it = element_index_map.find(element);

    if (it == element_index_map.end()) {
      throw std::runtime_error("Element '" + element +
			       "' not found in elements list");
    }

    atom.setElementID(it->second);
  }
} // Cell::populateElementID

void Cell::populateElementNumbers() {
  std::vector<int> temp_num_atoms(this->elements().size(), 0);
  this->populateElementID();
  for (auto &atom : this->_atoms_) {
    temp_num_atoms[atom.getElementID()]++;
  }
  this->setElementsNumbers(temp_num_atoms);
} // Cell::populateElementNumbers

//---------------------------------------------------------------------------//
//---------------------------  Correct Positions ----------------------------//
//---------------------------------------------------------------------------//
void Cell::correctPositions() {

  // Create lattice matrix
  Matrix3D L = {{{this->_v_a_[0], this->_v_b_[0], this->_v_c_[0]},
		 {this->_v_a_[1], this->_v_b_[1], this->_v_c_[1]},
		 {this->_v_a_[2], this->_v_b_[2], this->_v_c_[2]}}};
  Matrix3D L_inv = invertMatrix(L);

  for (auto &MyAtom : this->_atoms_) {
    Vector3D pos = {
	{MyAtom.position()[0], MyAtom.position()[1], MyAtom.position()[2]}};

    // Convert to fractional coordinates
    Vector3D frac = matrixVectorMultiply(L_inv, pos);

    // Wrap fractional coordinates to [0, 1)
    for (int m = 0; m < 3; ++m) {
      frac.data[m] -= std::floor(frac.data[m]);
    }

    // Convert back to Cartesian coordinates
    Vector3D corrected_pos = matrixVectorMultiply(L, frac);

    MyAtom.setPosition(
	{corrected_pos.data[0], corrected_pos.data[1], corrected_pos.data[2]});
  }
} // Cell::CorrectPositions

//---------------------------------------------------------------------------//
//----------------------- Correct Fractional Positions ----------------------//
//---------------------------------------------------------------------------//
void Cell::correctFracPositions() {
  std::array<double, 3> aux_pos;

  for (auto &MyAtom : this->_atoms_) {
    aux_pos = MyAtom.position();
    double i = aux_pos[0] - std::floor(aux_pos[0]);
    double j = aux_pos[1] - std::floor(aux_pos[1]);
    double k = aux_pos[2] - std::floor(aux_pos[2]);
    aux_pos[0] = i * this->_v_a_[0] + j * this->_v_b_[0] + k * this->_v_c_[0];
    aux_pos[1] = i * this->_v_a_[1] + j * this->_v_b_[1] + k * this->_v_c_[1];
    aux_pos[2] = i * this->_v_a_[2] + j * this->_v_b_[2] + k * this->_v_c_[2];

    MyAtom.setPosition(aux_pos);
  }
} // Cell::CorrectFracPositions

//---------------------------------------------------------------------------//
//-------------------------- Populate Bond Length ---------------------------//
//---------------------------------------------------------------------------//
void Cell::populateBondLength(double Bond_Factor) {
  // Number of elements in the Cell
  const int n = this->elements().size();
  // Initialize Bond length matrix (nxn) as zeros
  std::vector<std::vector<double>> temp_matrix(n, std::vector<double>(n, 0.0));

  for (int i = 0; i < n; i++) {
    double aux = covalentRadii(this->elements()[i]);
    for (int j = 0; j < n; j++) {
      temp_matrix[i][j] =
	  (aux + covalentRadii(this->elements()[j])) * Bond_Factor;
    }
  }
  this->_bond_length_ = temp_matrix;
} // Cell::populateBondlength

//---------------------------------------------------------------------------//
//--------------------------- Distance Population ---------------------------//
//---------------------------------------------------------------------------//
void Cell::distancePopulation(double r_cut, bool self_interaction) {
  // Validate inputs
  if (r_cut <= 0)
    throw std::invalid_argument("r_cut must be positive");
  if (this->_v_a_.empty() || this->_v_b_.empty() || this->_v_c_.empty())
    throw std::logic_error("Cell vectors not initialized");

  const size_t n = this->_elements_.size();
  std::vector<std::vector<std::vector<double>>> temp_dist(
      n, std::vector<std::vector<double>>(n));
  const double r_cut_sq = r_cut * r_cut;

  if (this->_bond_length_.empty())
    this->populateBondLength(1.2);

  this->populateElementNumbers();
  this->correctPositions();

  // Calculate supercell dimensions
  int k_ = ceil(r_cut / this->_v_c_[2]);
  int j_ = ceil(r_cut / this->_v_b_[1]);
  int i_ = ceil(r_cut / this->_v_a_[0]);

  // Force self_interaction in case of a small supercell
  if (k_ * j_ * i_ > 8)
    self_interaction = true;

  // Precompute displacements
  std::vector<std::array<double, 3>> displacements;
  for (int i = -i_; i <= i_; ++i) {
    for (int j = -j_; j <= j_; ++j) {
      for (int k = -k_; k <= k_; ++k) {
	displacements.push_back(
	    {i * this->_v_a_[0] + j * this->_v_b_[0] + k * this->_v_c_[0],
	     i * this->_v_a_[1] + j * this->_v_b_[1] + k * this->_v_c_[1],
	     i * this->_v_a_[2] + j * this->_v_b_[2] + k * this->_v_c_[2]});
      }
    }
  }

  // Thread-safe error handling
  std::vector<std::mutex> global_mutex(n);
  std::exception_ptr async_exception = nullptr;

  auto process_atom = [&](Atom &atom_A) {
    try {
      const int id_A = atom_A.getElementID();

      for (Atom &atom_B : this->_atoms_) {
	if (&atom_A == &atom_B && !self_interaction)
	  continue;

	const int id_B = atom_B.getElementID();
	const auto &pos_B = atom_B.position();

	for (const auto &disp : displacements) {
	  const std::array<double, 3> img_pos = {
	      pos_B[0] + disp[0], pos_B[1] + disp[1], pos_B[2] + disp[2]};

	  const double dx = atom_A.position()[0] - img_pos[0];
	  const double dy = atom_A.position()[1] - img_pos[1];
	  const double dz = atom_A.position()[2] - img_pos[2];
	  const double dist_sq = dx * dx + dy * dy + dz * dz;

	  if (dist_sq > r_cut_sq)
	    continue;

	  const double dist = std::sqrt(dist_sq);

	  // Ignore self-interaction
	  if (dist == 0 && (atom_A.getID() == atom_B.getID()))
	    continue;

	  // Collitions checks
	  if (dist < 0.1) {
	    std::cerr << "\nERROR: The atoms:\n"
		      << atom_A.element() << "_" << atom_A.getID()
		      << " in position (" << atom_A.position()[0] << ", "
		      << atom_A.position()[1] << ", " << atom_A.position()[2]
		      << "),\n"
		      << atom_B.element() << "_" << atom_B.getID()
		      << " in position (" << img_pos[0] << ", " << img_pos[1]
		      << ", " << img_pos[2] << ").\n"
		      << "Have a distance less than 10 pm.\n";
	    exit(1);
	  }
	  if (dist < 0.5) {
	    std::cerr << "\nWARNING: The atoms:\n"
		      << atom_A.element() << "_" << atom_A.getID()
		      << " in position (" << atom_A.position()[0] << ", "
		      << atom_A.position()[1] << ", " << atom_A.position()[2]
		      << "),\n"
		      << atom_B.element() << "_" << atom_B.getID()
		      << " in position (" << img_pos[0] << ", " << img_pos[1]
		      << ", " << img_pos[2] << ").\n"
		      << "Have a distance less than the Bohr Radius.\n";
	  }

	  // Thread-safe bond addition
	  if (dist <= _bond_length_[id_A][id_B]) {
	    std::lock_guard<std::mutex> lock(global_mutex[id_A]);
	    atom_A.addBondedAtom(
		Atom(atom_B.element(), img_pos, atom_B.getID(), id_B));
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
  std::for_each(std::execution::par, _atoms_.begin(), _atoms_.end(),
		process_atom);

  if (async_exception) {
    std::rethrow_exception(async_exception);
  }

  this->_distances_ = std::move(temp_dist);
} // Cell::distancePopulation

//---------------------------------------------------------------------------//
//--------------------------- Coordination Number ---------------------------//
//---------------------------------------------------------------------------//
void Cell::coordinationNumber() {
  const size_t num_elements = this->_distances_.size();
  if (num_elements == 0)
    throw std::logic_error("No elements in cell");
  const size_t num_atoms = this->_atoms_.size();
  if (num_atoms == 0)
    throw std::logic_error("No atoms in cell");

  // Step 1: Find maximum coordination number
  int max_Nc = 0;
  for (auto &atom : this->_atoms_) {
    max_Nc = std::max(max_Nc, static_cast<int>(atom.bonded_atoms().size()));
  }

  // Step 2: Initialize 3D tensor [element][element/total][coordination]
  const size_t coordination_bins = max_Nc + 2;

  // Dimensions: [element][target (elements + total)][coordination count]
  using ElementHist = std::vector<std::vector<int>>;
  std::vector<ElementHist> temp_nc(
      num_elements,
      ElementHist(num_elements + 1, std::vector<int>(coordination_bins, 0)));

  // Step 3: Build coordination histogram
  for (auto &atom : this->_atoms_) {
    const size_t element_id = atom.getElementID();
    std::vector<int> element_counts(num_elements, 0);

    // Count bonded atoms per element
    for (auto &bonded_atom : atom.bonded_atoms()) {
      const size_t bonded_element_id = bonded_atom.getElementID();
      element_counts[bonded_element_id]++;
    }

    // Update histograms
    for (size_t target_element = 0; target_element < num_elements;
	 ++target_element) {
      const int count = element_counts[target_element];
      temp_nc[element_id][target_element][count]++;
    }

    // Update total coordination count
    const int total_coord = atom.bonded_atoms().size();
    temp_nc[element_id][num_elements][total_coord]++;
  }

  // Step 4: Flatten into 2D histogram matrix
  const size_t num_columns = num_elements * num_elements + num_elements + 1;
  const size_t num_rows = coordination_bins;

  std::vector<std::vector<int>> coordination_hist(
      num_columns, std::vector<int>(num_rows, 0));

  // Fill coordination number header (first column)
  for (size_t coord = 0; coord < num_rows; ++coord) {
    coordination_hist[0][coord] = coord;
  }

  // Fill element-element interactions
  size_t hist_column = 1;
  for (size_t elem_i = 0; elem_i < num_elements; ++elem_i) {
    for (size_t elem_j = 0; elem_j < num_elements; ++elem_j) {
      coordination_hist[hist_column++] = temp_nc[elem_i][elem_j];
    }
  }

  // Fill element totals
  for (size_t elem = 0; elem < num_elements; ++elem) {
    coordination_hist[hist_column++] = temp_nc[elem][num_elements];
  }

  _Z_ = std::move(coordination_hist);
} // Cell::coorcinationNumber

//---------------------------------------------------------------------------//
//------------------------- Plane Angle Population --------------------------//
//---------------------------------------------------------------------------//
void Cell::planeAnglePopulation(bool degree) {
  const double conversion_factor = degree ? constants::rad2deg : 1.0;

  const size_t num_elements = this->elements().size();
  if (num_elements == 0)
    throw std::logic_error("No elements in cell");
  if (this->_atoms_.size() == 0)
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

  for (size_t atom_idx = 0; atom_idx < this->_atoms_.size(); ++atom_idx) {
    Atom &central_atom = this->_atoms_[atom_idx];
    std::vector<Atom> bonded = central_atom.bonded_atoms();
    const size_t central_elem = central_atom.getElementID();

    // Use indices to avoid duplicate pairs
    for (size_t i = 0; i < bonded.size(); ++i) {
      Atom &atom_a = bonded[i];
      const size_t elem_a = atom_a.getElementID();

      for (size_t j = i + 1; j < bonded.size(); ++j) {
	Atom &atom_b = bonded[j];
	const size_t elem_b = atom_b.getElementID();

	// Calculate angle and store in both element order permutations
	try {
	  const double angle =
	      central_atom.getAngle(atom_a, atom_b) * conversion_factor;
	  temp_angles[elem_a][central_elem][elem_b].push_back(angle);
	  temp_angles[elem_b][central_elem][elem_a].push_back(angle);
	} catch (const std::invalid_argument &e) {
	  std::cerr << "Angle calculation error: " << e.what() << " for atoms "
		    << atom_a.getID() << " and " << atom_b.getID() << '\n';
	}
      }
    }
  }
  this->_angles_ = std::move(temp_angles);
} // Cell::planeAnglePopulation

//---------------------------------------------------------------------------//
//--------------------- Radial Distribution Functions -----------------------//
//---------------------------------------------------------------------------//
void Cell::radialDistributionFunctions(double r_cut, double bin_width,
				       bool normalize) {
  // Validate input parameters
  if (bin_width <= 0)
    throw std::invalid_argument("Bin width must be positive");
  if (r_cut <= 0)
    throw std::invalid_argument("Cutoff radius must be positive");
  if (this->volume <= 0)
    throw std::logic_error("Cell volume must be positive");

  const size_t num_elements = this->_elements_.size();
  if (num_elements == 0)
    throw std::logic_error("No elements in cell");
  const size_t num_atoms = this->_atoms_.size();
  if (num_atoms == 0)
    throw std::logic_error("No atoms in cell");

  // Calculate basic parameters
  const size_t num_bins = static_cast<size_t>(std::floor(r_cut / bin_width));
  const size_t num_pairs = num_elements * (num_elements + 1) / 2;
  const size_t num_columns = num_pairs + 2; // +1 for r, +1 for total

  // Initialize weight factors with default value 1.0 (for normalization)
  std::vector<double> weight_factors(num_columns, 1.0);

  /*
   * _w_ij_ is the weighting factor for the partial of g_ij
   * There are to normalization commonly used:
   *     - HHS(Atlas) normalization, commonly used by
   *       experimental scientist.
   *     - Ashcroft - Waseda normalization, commonly
   *       used by theoretical scientist.
   * By default we use Ashcroft normalization:
   * _w_ij_ = c_i * c_j (Product of the numeric concentrations)
   * _w_ij_ = (#atoms_i * #atoms_j)/total_number_atoms^2
   * We offer the option to normalize to HHS by using the
   * -n, --normalize option.
   */

  // Calculate Ashcroft weights if needed
  if (!normalize) {
    const double norm_factor = 1.0 / (num_atoms * num_atoms);
    size_t weight_index = 1; // Start after r column

    for (size_t i = 0; i < num_elements; ++i) {
      const double ni = element_numbers()[i];
      for (size_t j = i; j < num_elements; ++j) {
	const double nj = element_numbers()[j];
	double weight = 2.0 * ni * nj * norm_factor;
	if (i == j)
	  weight *= 0.5;
	weight_factors[weight_index++] = weight;
      }
    }
  }
  this->_w_ij_ = weight_factors;

  // Initialize histogram with r values in first row
  std::vector<std::vector<double>> histogram(
      num_columns, std::vector<double>(num_bins, 0.0));
  for (size_t bin = 0; bin < num_bins; ++bin) {
    histogram[0][bin] = (bin + 0.5) * bin_width;
  }

  // Fill distance histogram
  size_t current_column = 1;
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      for (const auto &distance : this->_distances_[i][j]) {
	if (distance >= r_cut)
	  continue;

	const size_t bin = static_cast<size_t>(distance / bin_width);
	if (bin >= num_bins)
	  continue;

	// Count both i-j and j-i pairs unless same element
	histogram[current_column][bin] += (i == j) ? 1.0 : 2.0;
      }
      current_column++;
    }
  }

  // Calculate total Radial Distribution Function
  auto &total_dist = histogram.back();
  for (size_t col = 1; col < num_columns - 1; ++col) {
    for (size_t bin = 0; bin < num_bins; ++bin) {
      total_dist[bin] += histogram[col][bin];
    }
  }

  // Normalize histogram
  const double volume_factor = 4.0 * constants::pi * (num_atoms / this->volume);
  const double norm_factor = 1.0 / (num_atoms * bin_width);

  for (size_t col = 1; col < num_columns; ++col) {
    const double weight = normalize ? 1.0 : weight_factors[col];
    const double scaling = norm_factor / weight;

    for (size_t bin = 0; bin < num_bins; ++bin) {
      histogram[col][bin] *= scaling;
    }
  }
  this->_J_ = histogram;

  /*
   * Calculate g(r) with the inverse of the J(r) definition:
   * J(r) = 4 * pi * r^2 * rho_0 * g(r)
   * We calculate G(r) with the definition:
   * G_ij = 4 * pi * r * rho_0 * [g_ij(r) - _w_ij_]
   */

  // Calculate g(r) and G(r)
  std::vector<std::vector<double>> g_r = histogram;
  std::vector<std::vector<double>> G_r = histogram;

  for (size_t col = 1; col < num_columns; ++col) {
    for (size_t bin = 0; bin < num_bins; ++bin) {
      const double r = histogram[0][bin];
      const double r_sq = r * r;

      // Avoid division by zero for first bin
      if (r_sq < std::numeric_limits<double>::epsilon()) {
	g_r[col][bin] = 0.0;
	G_r[col][bin] = 0.0;
	continue;
      }

      // Calculate g(r)
      g_r[col][bin] /= volume_factor * r_sq;

      // Calculate G(r)
      const double reference = normalize ? 1.0 : weight_factors[col];
      G_r[col][bin] = volume_factor * r * (g_r[col][bin] - reference);
    }
  }

  this->_g_ = std::move(g_r);
  this->_G_ = std::move(G_r);

} // Cell::radialDistributionFunctions

//---------------------------------------------------------------------------//
//------------------------ Plane Angle Distribution -------------------------//
//---------------------------------------------------------------------------//
void Cell::planeAngleDistribution(double theta_cut, double bin_width) {
  // Validate input parameters
  if (bin_width <= 0.0) {
    throw std::invalid_argument("Bin width must be positive");
  }
  if (theta_cut <= 0.0) {
    throw std::invalid_argument("Theta cutoff must be positive");
  }

  /*
   * The number of columns in the output file is:
   *         n       x      (n+1)! / [2 x (n-1)!]
   * Central atom        Combination with repetition
   *                         of n in groups of 2
   * it's reduced to n x (n+1) x n /2
   * one extra column is added for Theta (angle)
   */
  const size_t num_elements = this->_elements_.size();
  const size_t num_bins = static_cast<size_t>(theta_cut / bin_width);
  const size_t num_combinations =
      num_elements * num_elements * (num_elements + 1) / 2;
  const size_t num_columns = 1 + num_combinations + 1; // [theta] + combinations

  // Initialize histogram with theta bins in first column
  std::vector<std::vector<double>> histogram(
      num_columns, std::vector<double>(num_bins, 0.0));

  // Fill theta values
  for (size_t bin = 0; bin < num_bins; ++bin) {
    histogram[0][bin] = (bin + 0.5) * bin_width;
  }

  // Fill angle counts
  size_t current_column = 1;
  for (size_t central_elem = 0; central_elem < num_elements; ++central_elem) {
    for (size_t elem_j = 0; elem_j < num_elements; ++elem_j) {
      for (size_t elem_k = elem_j; elem_k < num_elements; ++elem_k) {
	for (const double angle :
	     this->_angles_[central_elem][elem_j][elem_k]) {
	  if (angle >= theta_cut)
	    continue;

	  const size_t bin = static_cast<size_t>(angle / bin_width);
	  if (bin < num_bins) {
	    histogram[current_column][bin] += 1.0;
	  }
	}
	current_column++;
      }
    }
  }

  // Calculate normalization factor
  double normalization = 0.0;
  for (size_t col = 1; col <= num_combinations; ++col) {
    for (size_t bin = 0; bin < num_bins; ++bin) {
      normalization += histogram[col][bin];
    }
  }
  normalization *= bin_width;

  // Handle zero normalization case
  if (normalization < std::numeric_limits<double>::epsilon()) {
    this->_F_ = std::move(histogram);
    return;
  }

  // Normalize and calculate total distribution
  for (size_t col = 1; col <= num_combinations; ++col) {
    for (size_t bin = 0; bin < num_bins; ++bin) {
      histogram[col][bin] /= normalization;
      histogram.back()[bin] += histogram[col][bin];
    }
  }

  this->_F_ = std::move(histogram);
} // Cell::planeAngleDistribution

//---------------------------------------------------------------------------//
//---------------------------------- S(Q) -----------------------------------//
//---------------------------------------------------------------------------//
void Cell::SQ(double q_max, double q_bin_width, bool normalize) {
  // Validate input parameters
  if (q_max <= 0.0)
    throw std::invalid_argument("q_max must be positive");
  if (q_bin_width <= 0.0)
    throw std::invalid_argument("q_bin_width must be positive");
  if (this->_G_.empty() || this->_G_[0].empty())
    throw std::logic_error("G(r) data not initialized");

  /* n_: the number of columns is given by:
   * The combination n atoms taken in pairs,
   * plus one column for r,
   * plus one column for the total S(Q)
   */
  const size_t num_elements = this->_elements_.size();
  const size_t num_atoms = this->_atoms_.size();
  const size_t num_q_bins =
      static_cast<size_t>(std::round(q_max / q_bin_width));
  const size_t num_columns = 1 + (num_elements * (num_elements + 1)) / 2 + 1;

  // Initialize weights and S(Q) matrix
  std::vector<double> weights(num_columns, 1.0);
  std::vector<std::vector<double>> S_q(num_columns,
				       std::vector<double>(num_q_bins, 0.0));

  // Calculate weights if not normalized
  if (!normalize) {
    const double norm_factor = 1.0 / (num_atoms * num_atoms);
    size_t weight_idx = 1;

    for (size_t i = 0; i < num_elements; ++i) {
      const double ni = element_numbers()[i];
      for (size_t j = i; j < num_elements; ++j) {
	double weight = 2.0 * ni * element_numbers()[j] * norm_factor;
	if (i == j)
	  weight *= 0.5;
	weights[weight_idx++] = weight;
      }
    }
    _w_ij_ = weights;
  }

  // Initialize q values using std::transform
  std::generate(S_q[0].begin(), S_q[0].end(),
		[n = 0, q_bin_width]() mutable { return (n++) * q_bin_width; });

  // Get integration parameters
  const auto &r_values = this->_G_[0];
  const double dr = r_values[1] - r_values[0];
  const size_t num_r_points = r_values.size();

  // Lambda for sinc function computation
  auto safe_sinc = [](double x) {
    return std::abs(x) < 1e-8 ? 1.0 : std::sin(x) / x;
  };

  // Main computation loop
  for (size_t col = 1; col < num_columns; ++col) {
    for (size_t q_idx = 0; q_idx < num_q_bins; ++q_idx) {
      const double q = S_q[0][q_idx];
      double integral = 0.0;

      // Sequential integration
      for (size_t r_idx = 0; r_idx < num_r_points - 1; ++r_idx) {
	const double r = r_values[r_idx];
	/*
	 * The structure factor S(q) is calculated with:
	 * S(q) = 1 + \int{dr sinc(q * r) * r * G(r)}
	 */
	const double term = safe_sinc(q * r) * r * this->_G_[col][r_idx];
	integral += (r_idx == 0 ? 0.5 * term : term);
      }

      // Add last term with 0.5 weight
      const double last_term = safe_sinc(q * r_values.back()) *
			       r_values.back() * this->_G_[col].back();
      S_q[col][q_idx] = weights[col] + (integral + 0.5 * last_term) * dr;
    }
  }

  this->_S_ = std::move(S_q);
} // Cell::SQ

//---------------------------------------------------------------------------//
//----------------------------------- XRD -----------------------------------//
//---------------------------------------------------------------------------//
void Cell::XRD(double lambda, double theta_min, double theta_max,
	       double bin_width) {
  int n_, m_, i, j, k, col, row;
  double norm, aux;
  int n = this->elements().size();
  std::string header;

  m_ = ceil((theta_max - theta_min) / bin_width);
  n_ = n * (n + 1) / 2 + 1;

  std::vector<std::vector<double>> temp_hist(n_, std::vector<double>(m_, 0));
  // Fill the theta values of the histogram
  for (i = 0; i < m_; i++) {
    temp_hist[0][i] = theta_min + i * bin_width;
  }
  col = 0;
  // Triple loop to iterate over the distances tensor.
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      col++;
      for (const auto &it : this->_distances_[i][j]) {
	/*
	 * Bragg's Law
	 * n * lambda = 2d sin(theta)
	 * theta = asin((n * lambda) / (2 * d))
	 */
	aux = lambda / (it * 2.0);
	k = 1;
	while (k < 3) {
	  row = floor(((2 * constants::rad2deg * asin(k * aux)) - theta_min) /
		      bin_width);
	  k++;
	  if ((0 <= row) && (row < m_)) {
	    temp_hist[col][row]++;
	    if (i != j) {
	      temp_hist[col][row]++;
	    }
	  }
	}
      }
    }
  }

  // Double loop to find normalization factor
  norm = 0.0;
  for (i = 1; i < n_; i++) {
    for (j = 0; j < m_; j++) {
      norm = std::max(temp_hist[i][j], norm);
    }
  }
  norm *= 0.01;
  // Double loop to normalize PADHistogram
  for (i = 1; i < n_; i++) {
    for (j = 0; j < m_; j++) {
      temp_hist[i][j] /= norm;
    }
  }

  this->_X_ = temp_hist;
} // Cell:XRD

//--------------------------------------------------------------------------//
//------------------------------- Smoothing --------------------------------//
//--------------------------------------------------------------------------//
void Cell::Smoothing(double sigma, int _kernel_) {
  int n_, m_, col, row;

  // n_: number of columns in the histogram
  n_ = this->_g_.size();
  // m_: number of rows in the histogram
  m_ = this->_g_[0].size();
  // Array of smoothed histograms
  std::vector<std::vector<double>> temp_hist(n_, std::vector<double>(m_, 0));

  for (row = 0; row < m_; row++) {
    temp_hist[0][row] = this->_J_[0][row];
  }
  for (col = 1; col < n_; col++) {
    temp_hist[col] =
	KernelSmoothing(temp_hist[0], this->_J_[col], sigma, _kernel_);
  }
  this->_J_smoothed_ = temp_hist;

  // Smoothing g
  for (row = 0; row < m_; row++) {
    temp_hist[0][row] = this->_g_[0][row];
  }
  for (col = 1; col < n_; col++) {
    temp_hist[col] =
	KernelSmoothing(temp_hist[0], this->_g_[col], sigma, _kernel_);
  }
  this->_g_smoothed_ = temp_hist;

  // Smoothing G
  for (row = 0; row < m_; row++) {
    temp_hist[0][row] = this->_G_[0][row];
  }
  for (col = 1; col < n_; col++) {
    temp_hist[col] =
	KernelSmoothing(temp_hist[0], this->_G_[col], sigma, _kernel_);
  }
  this->_G_smoothed_ = temp_hist;

  // Smoothing F
  // n_: number of columns in the histogram
  n_ = this->_F_.size();
  // m_: number of rows in the histogram
  m_ = this->_F_[0].size();
  // Resize 2D vector of smoothed histograms
  temp_hist.resize(n_);
  for (auto &col : temp_hist) {
    col.resize(m_, 0); // resize and fill with 0
  }
  for (row = 0; row < m_; row++) {
    temp_hist[0][row] = this->_F_[0][row];
  }
  for (col = 1; col < n_; col++) {
    temp_hist[col] =
	KernelSmoothing(temp_hist[0], this->_F_[col], sigma, _kernel_);
  }
  this->_F_smoothed_ = temp_hist;

  // Smoothing S
  // n_: number of columns in the histogram
  n_ = this->_S_.size();
  // m_: number of rows in the histogram
  m_ = this->_S_[0].size();
  // Resize 2D vector of smoothed histograms
  temp_hist.resize(n_);
  for (auto &col : temp_hist) {
    col.resize(m_, 0); // resize and fill with 0
  }
  for (row = 0; row < m_; row++) {
    temp_hist[0][row] = this->_S_[0][row];
  }
  for (col = 1; col < n_; col++) {
    temp_hist[col] =
	KernelSmoothing(temp_hist[0], this->_S_[col], sigma, _kernel_);
  }
  this->_S_smoothed_ = temp_hist;

  // Smoothing X
  // n_: number of columns in the histogram
  n_ = this->_X_.size();
  // m_: number of rows in the histogram
  m_ = this->_X_[0].size();
  // Resize 2D vector of smoothed histograms
  temp_hist.resize(n_);
  for (auto &col : temp_hist) {
    col.resize(m_, 0); // resize and fill with 0
  }
  for (row = 0; row < m_; row++) {
    temp_hist[0][row] = this->_X_[0][row];
  }
  for (col = 1; col < n_; col++) {
    temp_hist[col] =
	KernelSmoothing(temp_hist[0], this->_X_[col], sigma, _kernel_);
  }
  this->_X_smoothed_ = temp_hist;

} // Cell::RDFSmoothing

//--------------------------------------------------------------------------//
//------------------------------ Voronoi Index -----------------------------//
//--------------------------------------------------------------------------//
void Cell::voronoiIndex() {
  /*
   * This function is currently and STUB and gives a rought stimate
   * to the actual voronoiindex, to compute the correct voronoi index
   * you need to compute the voronoi tesselation and then get the index
   * this is currently planned for v2.0 of the code, as the entire algorithm
   * should be generated from scratch.
   */
  for (auto &atom_A : this->_atoms_) {
    std::vector<int> atom_A_ids = atom_A.getBondedAtomsID();
    std::sort(atom_A_ids.begin(), atom_A_ids.end());
    std::cout << "A " << atom_A.getID() << ": (";
    for (auto A_id : atom_A_ids) {
      std::cout << A_id << ", ";
    }
    std::cout << " )" << std::endl;
  }
} // Cell::Voronoi
