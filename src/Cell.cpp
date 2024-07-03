/* ----------------------------------------------------------------------------
 * Correlation: An Analysis Tool for Liquids and for Amorphous Solids
 * Copyright (c) 2013-2024 Isaías Rodríguez <isurwars@gmail.com>
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
#include "Cell.hpp"

#include <algorithm>
#include <cmath>
#include <execution> // For parallel execution policies
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <ostream>
#include <regex>
#include <utility>

#include "Atom.hpp"
#include "Constants.hpp"
#include "Smoothing.hpp"
#include "Templates.hpp"

// Lattice parameters constructor
Cell::Cell(std::array<double, 6> lat) {
  this->SetLatticeParameters(lat);
  this->SetAtoms({});
  this->SetElements({});
  this->SetLatticeVectors();
}

// Default constructor
Cell::Cell() {
  this->SetLatticeParameters({1.0, 1.0, 1.0, 90, 90, 90});
  this->SetAtoms({});
  this->SetElements({});
  this->SetLatticeVectors();
}

// Lattice vector constructor
void Cell::SetFromVectors(std::vector<double> v1, std::vector<double> v2,
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
  this->SetLatticeParameters({v1_n, v2_n, v3_n, a, b, c});
  this->SetAtoms({});
  this->SetElements({});
  this->SetLatticeVectors();
}

// Calculate the lattice vectors from the lattice parameters.
void Cell::SetLatticeVectors() {
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

  this->v_a_ = {A, 0.0, 0.0};
  this->v_b_ = {B * c_g, B * s_g, 0.0};
  this->v_c_ = {
      C * c_b, C * (c_a - c_b * c_g) / s_g,
      C * sqrt(1 - c_a * c_a - c_b * c_b - c_g * c_g + 2 * c_a * c_b * c_g) /
          s_g};

  this->volume =
      this->v_a_[0] *
          (this->v_b_[1] * this->v_c_[2] - this->v_b_[2] * this->v_c_[1]) -
      this->v_a_[1] *
          (this->v_b_[0] * this->v_c_[2] - this->v_b_[2] * this->v_c_[0]) +
      this->v_a_[2] *
          (this->v_b_[0] * this->v_c_[1] - this->v_b_[1] * this->v_c_[0]);
}

// Correct the initial positions to in-cell positions
void Cell::CorrectPositions() {
  std::array<double, 3> aux_pos;
  std::vector<int> temp_num_atoms(this->elements().size(), 0);

  for (auto &MyAtom : this->_atoms_) {
    temp_num_atoms[MyAtom.element_id()]++;
    aux_pos = MyAtom.position();

    // Adjust position along the c-axis
    double k = aux_pos[2] / this->v_c_[2];
    for (int m = 0; m < 3; ++m) {
      aux_pos[m] -= k * this->v_c_[m];
    }

    // Adjust position along the b-axis
    double j = aux_pos[1] / this->v_b_[1];
    for (int m = 0; m < 3; ++m) {
      aux_pos[m] -= j * this->v_b_[m];
    }

    // Adjust position along the a-axis
    double i = aux_pos[0] / this->v_a_[0];

    // Determine the integer indices
    int i_ = (i >= -1e-15) ? static_cast<int>(std::trunc(i))
                           : static_cast<int>(std::trunc(i)) - 1;
    int j_ = (j >= -1e-15) ? static_cast<int>(std::trunc(j))
                           : static_cast<int>(std::trunc(j)) - 1;
    int k_ = (k >= -1e-15) ? static_cast<int>(std::trunc(k))
                           : static_cast<int>(std::trunc(k)) - 1;

    // Correct the position
    for (int m = 0; m < 3; ++m) {
      aux_pos[m] =
          MyAtom.position()[m] -
          (i_ * this->v_a_[m] + j_ * this->v_b_[m] + k_ * this->v_c_[m]);
    }

    // Set the new position
    MyAtom.SetPosition(aux_pos);
  }

  // Update the number of atoms per element
  this->SetElementsNumbers(temp_num_atoms);
} // Cell::CorrectPositions

// Correct the fractional positions to absolute positions
void Cell::CorrectFracPositions() {
  std::array<double, 3> aux_pos;

  for (auto &MyAtom : this->_atoms_) {
    aux_pos = MyAtom.position();
    double i = aux_pos[0];
    double j = aux_pos[1];
    double k = aux_pos[2];

    aux_pos[0] = i * this->v_a_[0] + j * this->v_b_[0] + k * this->v_c_[0];
    aux_pos[1] = i * this->v_a_[1] + j * this->v_b_[1] + k * this->v_c_[1];
    aux_pos[2] = i * this->v_a_[2] + j * this->v_b_[2] + k * this->v_c_[2];

    MyAtom.SetPosition(aux_pos);
  }
} // Cell::CorrectFracPositions

// Populate Bond_length matrix
void Cell::PopulateBondLength(double Bond_Factor) {
  std::pair<bool, int> MyId;
  int i, j;

  // Number of elements in the Cell
  const int n = this->elements().size();
  // Initialize Bond length matrix (nxn) as zeros
  std::vector<std::vector<double>> temp_matrix(n, std::vector<double>(n, 0.0));

  // Iterate in Atoms list to assign the id in the matrix to every atom.
  for (auto &MyAtom : this->_atoms_) {
    MyId = findInVector(this->elements(), MyAtom.element());
    MyAtom.SetElementId(MyId.second);
  }

  for (i = 0; i < n; i++) {
    double aux = Covalent_Radii(this->elements()[i]);
    for (j = 0; j < n; j++) {
      temp_matrix[i][j] =
          (aux + Covalent_Radii(this->elements()[j])) * Bond_Factor;
    }
  }
  this->bond_length = temp_matrix;
} // Cell::PopulateBondlength

void Cell::ReadBOND(std::string file_name) {
  /*
   * This function reads the in_bond_file to populate the bond_length Tensor.
   * The file should be in the format:
   * element_B element_A distance(in Angstroms)
   *
   * For example:
   *
   * Si Si 2.29
   * Mg Mg 2.85
   * C  C  1.55
   * C  Si 1.86
   * Si Mg 2.57
   * C  Mg 2.07
   *
   * Any missing pair of elements will use the bond_parameter as a default.
   */
  std::ifstream myfile(file_name);
  std::smatch match;
  std::pair<bool, int> MyIdA, MyIdB;

  /*
   * Every line should have two elements and a bond length separeted by spaces:
   *
   * element_A element_B bond_length
   */

  std::regex regex_bond("^([A-Z][a-z]?)"
                        "(\\s+)"
                        "([A-Z][a-z]?)"
                        "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)");

  if (myfile.is_open()) {
    /* Check if the file is open */
    std::string line;
    while (std::getline(myfile, line)) {
      /* Read line by line */
      if (std::regex_search(line, match, regex_bond)) {
        /* Bond found */
        int i, j;
        double dist;
        MyIdA =
            findInVector(this->elements(), std::string(match.str(1).data()));
        if (MyIdA.first) {
          i = MyIdA.second;
        }
        MyIdB =
            findInVector(this->elements(), std::string(match.str(3).data()));
        if (MyIdB.first) {
          j = MyIdB.second;
        }
        dist = std::stof(match.str(4).data());
        if (MyIdA.first && MyIdB.first) {
          this->bond_length[i][j] = dist;
          this->bond_length[j][i] = dist;
        }
      }
    }
  }
} // ReadBOND

// Cell::UpdateProgressBar
void Cell::UpdateProgressBar(double pos) {
  int barWidth = 50;
  int progress = round(pos * barWidth);
  std::cout << "\r[";

  for (int i = 0; i < barWidth; ++i) {
    if (i < progress) {
      std::cout << "=";
    } else if (i == progress) {
      std::cout << ">";
    } else {
      std::cout << " ";
    }
  }
  std::cout << "] " << int(pos * 100.0) << " %";
  std::cout.flush();
} // Cell::UpdateProgressBar

// DistancePopulation
void Cell::DistancePopulation(double r_cut, bool self_interaction) {
  std::array<double, 3> aux_pos;
  Atom img_atom;
  int i, j, k, i_, j_, k_;
  double aux_dist;
  double h_ = 1.0 / this->_atoms_.size();
  double progress = 0.0;

  // Number of elements in the Cell
  const int n = this->_elements_.size();
  // This matrix stores the distances between different types of elements
  std::vector<std::vector<std::vector<double>>> temp_dist(
      n, std::vector<std::vector<double>>(n));

  // Correct the atom positions to be inside the cell
  this->CorrectPositions();

  // Calculate supercell dimensions
  k_ = ceil(r_cut / this->v_c_[2]);
  j_ = ceil(r_cut / this->v_b_[1]);
  i_ = ceil(r_cut / this->v_a_[0]);

  // Force self_interaction in case of a small supercell
  if (k_ * j_ * i_ > 8)
    self_interaction = true;

  /*
   * This is the main loop to calculate the distance formed by every two atoms.
   *
   * The first loop iterates over every atom in the cell.
   * The second loop iterate over every atom in a SuperCell.
   * This is by far the most demanding cicle of the entire program and it is
   * currently in parallelization process.
   *
   * This loop can be further optimize by applying a Divide&Conquer Algorithm,
   * However the code needs a major refactor, and is one of the main objectives
   * For V2.0
   */
  auto calculate_distances = [&](Atom &atom_A) {
    int id_A = findInVector(this->_elements_, atom_A.element()).second;
    for (Atom &atom_B : this->_atoms_) {
      if (atom_A.GetID() != atom_B.GetID() || self_interaction) {
        int id_B = findInVector(this->_elements_, atom_B.element()).second;
        img_atom = atom_B;
        for (i = -i_; i <= i_; i++) {
          for (j = -j_; j <= j_; j++) {
            for (k = -k_; k <= k_; k++) {
              aux_pos = atom_B.position();
              aux_pos[0] +=
                  i * this->v_a_[0] + j * this->v_b_[0] + k * this->v_c_[0];
              aux_pos[1] += j * this->v_b_[1] + k * this->v_c_[1];
              aux_pos[2] += k * this->v_c_[2];
              img_atom.SetPosition(aux_pos);
              aux_dist = atom_A.Distance(img_atom);
              if (aux_dist <= r_cut) {
                // ignore self-interaction
                if (aux_dist == 0 && (atom_A.GetID() == atom_B.GetID())) {
                  continue;
                }
                // check for atoms collitions
                if (aux_dist < 0.1) {
                  std::cerr << "\nERROR: The atoms:\n"
                            << atom_A.element() << "_" << atom_A.GetID()
                            << " in position (" << atom_A.position()[0] << ", "
                            << atom_A.position()[1] << ", "
                            << atom_A.position()[2] << "),\n"
                            << atom_B.element() << "_" << atom_B.GetID()
                            << " in position (" << img_atom.position()[0]
                            << ", " << img_atom.position()[1] << ", "
                            << img_atom.position()[2] << ").\n"
                            << "Have a distance less than 10 pm.\n";
                  exit(1);
                }
                if (aux_dist < 0.5) {
                  std::cerr << "\nWARNING: The atoms:\n"
                            << atom_A.element() << "_" << atom_A.GetID()
                            << " in position (" << atom_A.position()[0] << ", "
                            << atom_A.position()[1] << ", "
                            << atom_A.position()[2] << "),\n"
                            << atom_B.element() << "_" << atom_B.GetID()
                            << " in position (" << img_atom.position()[0]
                            << ", " << img_atom.position()[1] << ", "
                            << img_atom.position()[2] << ").\n"
                            << "Have a distance less than the Bohr Radius.\n";
                }
                // If bonded, add it to bond vector
                if (aux_dist <= this->bond_length[atom_A.element_id()]
                                                 [img_atom.element_id()]) {
                  atom_A.AddBondedAtom(img_atom.GetImage());
                }
                // Add distance to matrix
                temp_dist[id_A][id_B].push_back(aux_dist);
              }
            }
          }
        }
      }
    }
    progress += h_;
    this->UpdateProgressBar(progress);
  };

  // Parallelize the main loop
  std::for_each(std::execution::par, std::begin(this->_atoms_),
                std::end(this->_atoms_), calculate_distances);

  std::cout << "\r[==================================================] 100 %"
            << std::endl;
  this->distances = temp_dist;
} // Cell::RDF

void Cell::CoordinationNumber() {
  int max_Nc = 0;
  int i;

  // Search for the maximum number of bonds in the cell
  for (auto &MyAtom : this->_atoms_) {
    if (max_Nc < static_cast<int>(MyAtom.bonded_atoms().size())) {
      max_Nc = static_cast<int>(MyAtom.bonded_atoms().size());
    }
  }
  const int n = this->elements().size();
  const int m = max_Nc + 2;
  // Create the Tensor nx(n+1)xm and initialize it with zeros
  std::vector<std::vector<std::vector<int>>> temp_nc(
      n, std::vector<std::vector<int>>(n + 1, std::vector<int>(m, 0)));
  // Search for the number of bonds per atom per element
  for (auto MyAtom = this->_atoms_.begin(); MyAtom != this->_atoms_.end();
       ++MyAtom) {
    std::vector<int> aux(n, 0);
    std::vector<Atom_Img> bon_aux = MyAtom->bonded_atoms();
    for (auto atom_A = bon_aux.begin(); atom_A != bon_aux.end(); ++atom_A) {
      aux[atom_A->element_id]++;
    }
    for (i = 0; i < n; i++) {
      temp_nc[MyAtom->element_id()][i][aux[i]]++;
    }
  }
  // Add total coordination per atom
  for (auto MyAtom = this->_atoms_.begin(); MyAtom != this->_atoms_.end();
       ++MyAtom) {
    temp_nc[MyAtom->element_id()][n][MyAtom->bonded_atoms().size()]++;
  }

  this->coordination = temp_nc;
} // Cell::CoordinationNumber

void Cell::PAD(bool degree) {
  double factor = 1.0;

  if (degree) {
    factor = constants::rad2deg;
  }

  // NxNxN Tensor to store the PAD
  const int n = this->elements().size();
  std::vector<std::vector<std::vector<std::vector<double>>>> temp_pad(
      n, std::vector<std::vector<std::vector<double>>>(
             n, std::vector<std::vector<double>>(n, std::vector<double>(0))));
  /*
   * This is the main loop to calculate the angles formed by every three atoms.
   * The connected atoms are calculated in Cell::RDF, and MUST be called first.
   *
   * The first loop iterates over every atom in the cell.
   * The second and third loop iterate over the connected atoms in every atom
   * instance. These three loops populate a 3D tensor of vectors, the indices of
   * the Tensor represent the three indices of the element in Cell::elements.
   *
   * By default the angle returned is given in degrees.
   */
  for (auto &MyAtom : this->_atoms_) {
    std::vector<Atom_Img> bon_aux = MyAtom.bonded_atoms();
    for (auto &atom_A : bon_aux) {
      for (auto &atom_B : bon_aux) {
        if (atom_A.atom_id != atom_B.atom_id) {
          temp_pad[atom_A.element_id][MyAtom.element_id()][atom_B.element_id]
              .push_back(MyAtom.GetAngle(atom_A, atom_B) * factor);
        }
      }
    }
  }
  this->angles = temp_pad;
} // Cell::PAD

void Cell::RDFHistogram(std::string filename, double r_cut, double bin_width,
                        bool normalize) {
  int n_, m_, i, j, col, row;
  int n = this->distances.size();
  std::string header;
  double num_atoms = this->atoms().size();

  /* m_: the number of rows is given by:
   * the cut radius/bin width
   */
  m_ = floor(r_cut / bin_width);
  /* n_: the number of columns is given by:
   * The combination n atoms taken in pairs,
   * plus one column for r,
   * plus one column for the total pair distribution
   */
  n_ = n * (n + 1) / 2 + 1 + 1;

  /*
   * w_ij is the weighting factor for the partial of G_ij
   * There are to normalization commonly used:
   *     - HHS(Atlas) normalization, commonly used by
   *       experimental scientist.
   *     - Ashcroft - Waseda normalization, commonly
   *       used by theoretical scientist.
   * By default we use Ashcroft normalization:
   * w_ij = c_i * c_j (Product of the numeric concentrations)
   * w_ij = (#atoms_i * #atoms_j)/total_number_atoms^2
   * We offer the option to normalize to HHS by using the
   * -n, --normalize option.
   */
  std::vector<double> temp_w_ij(n_, 1.0);

  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      temp_w_ij[i + j + 1] = 2.0 * this->element_numbers()[i] *
                             this->element_numbers()[j] /
                             (num_atoms * num_atoms);
      if (i == j) {
        temp_w_ij[i + j + 1] *= 0.5;
      }
    }
  }
  this->w_ij = temp_w_ij;

  std::vector<std::vector<double>> temp_hist(n_, std::vector<double>(m_, 0));
  // Fill the r values of the histogram
  for (i = 0; i < m_; i++) {
    temp_hist[0][i] = (i + 0.5) * bin_width;
  }
  col = 0;
  // Triple loop to iterate over the distances tensor.
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      col++;
      for (const auto &it : this->distances[i][j]) {
        row = floor(it / bin_width);
        if (row < m_) {
          temp_hist[col][row]++;
          if (i != j) {
            temp_hist[col][row]++;
          }
        }
      }
    }
  }

  // Calculate total Radial Distribution Function
  for (j = 0; j < m_; j++) {
    for (i = 1; i < n_ - 1; i++) {
      temp_hist[n_ - 1][j] += temp_hist[i][j];
    }
  }

  /*
   * Scale the histograms by the factor: 1 / (#atoms * bin_width)
   */

  double w_factor = num_atoms * bin_width;
  for (i = 1; i < n_; i++) {
    if (normalize) {
      w_factor = num_atoms * bin_width * this->w_ij[i];
    }
    std::transform(temp_hist[i].begin(), temp_hist[i].end(),
                   temp_hist[i].begin(),
                   [&w_factor](const auto &c) { return c / w_factor; });
  }

  this->_J_ = temp_hist;

  std::ofstream out_file(filename + "_J.csv");
  out_file << std::setprecision(6) << std::setw(13) << "r (Å),";
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      header = this->elements()[i] + "-" + this->elements()[j] + " (1/Å),";
      out_file << std::setprecision(6) << std::setw(13) << header;
    }
  }
  out_file << std::setprecision(6) << std::setw(12) << "J(r)," << std::endl
           << std::fixed;

  for (i = 0; i < m_; i++) {
    for (j = 0; j < n_; j++) {
      out_file << std::setprecision(6) << std::setw(11) << temp_hist[j][i]
               << ",";
    }
    out_file << std::endl;
  }
  out_file.close();
  out_file.clear();

  /*
   * Calculate g(r) with the inverse of the J(r) definition:
   * J(r) = 4 * pi * r^2 * rho_0 * g(r)
   */
  // numeric density
  double rho_0 = num_atoms / this->volume;
  for (col = 1; col < n_; col++) {
    for (row = 1; row < m_; row++) {
      temp_hist[col][row] /=
          4 * constants::pi * rho_0 * temp_hist[0][row] * temp_hist[0][row];
    }
  }

  this->_g_ = temp_hist;

  out_file.open(filename + "_g.csv");
  out_file << std::setprecision(6) << std::setw(13) << "r (Å),";

  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      header = this->elements()[i] + "-" + this->elements()[j] + " (1/Å),";
      out_file << std::setprecision(6) << std::setw(13) << header;
      if (normalize) {
        temp_w_ij[i + j + 1] = 1.0;
      }
    }
  }
  out_file << std::setprecision(6) << std::setw(12) << "g(r)," << std::endl
           << std::fixed;

  for (i = 0; i < m_; i++) {
    for (j = 0; j < n_; j++) {
      out_file << std::setprecision(6) << std::setw(11) << temp_hist[j][i]
               << ",";
    }
    out_file << std::endl;
  }
  out_file.close();
  out_file.clear();

  /*
   * We calculate G(r) with the definition:
   * G(r) = 4 * pi * r * rho_0 * [g(r) - 1];
   *
   * Weighted partials are calculated with:
   * G_ij = 4 * pi * r * rho_0 * [g_ij(r) - w_ij]
   */

  for (col = 1; col < n_; col++) {
    for (row = 1; row < m_; row++) {
      temp_hist[col][row] = 4 * constants::pi * rho_0 * temp_hist[0][row] *
                            (temp_hist[col][row] - temp_w_ij[col]);
    }
  }

  this->_G_ = temp_hist;

  out_file.open(filename + "_G_.csv");
  out_file << std::setprecision(6) << std::setw(13) << "r (Å),";
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      header = this->elements()[i] + "-" + this->elements()[j] + " (1/Å),";
      out_file << std::setprecision(6) << std::setw(13) << header;
    }
  }
  out_file << std::setprecision(6) << std::setw(12) << "G(r)," << std::endl
           << std::fixed;

  for (i = 0; i < m_; i++) {
    for (j = 0; j < n_; j++) {
      out_file << std::setprecision(6) << std::setw(11) << temp_hist[j][i]
               << ",";
    }
    out_file << std::endl;
  }
  out_file.close();
  out_file.clear();
} // Cell::RDFHistogram

void Cell::RDFSmoothing(std::string filename, double sigma, int _kernel_) {
  int n_, m_, i, j, col, row;
  std::string header;
  int n = this->elements().size();
  // n_: number of columns in the histogram
  n_ = this->_g_.size();
  // m_: number of rows in the histogram
  m_ = this->_g_[0].size();

  // Array of smoothed histograms
  std::vector<std::vector<double>> temp_hist(n_, std::vector<double>(m_, 0));
  /* Smoothing Loop */
  for (row = 0; row < m_; row++) {
    temp_hist[0][row] = this->_J_[0][row];
  }
  for (col = 1; col < n_; col++) {
    temp_hist[col] =
        KernelSmoothing(temp_hist[0], this->_J_[col], sigma, _kernel_);
  }

  /* printing loop */
  this->_J_smoothed_ = temp_hist;

  std::ofstream out_file(filename + "_J_smoothed.csv");
  out_file << std::setprecision(6) << std::setw(13) << "r (Å),";
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      header = this->elements()[i] + "-" + this->elements()[j] + " (1/Å),";
      out_file << std::setprecision(6) << std::setw(13) << header;
    }
  }
  out_file << std::setprecision(6) << std::setw(12) << "J(r)," << std::endl
           << std::fixed;
  for (i = 0; i < m_; i++) {
    for (j = 0; j < n_; j++) {
      out_file << std::setprecision(6) << std::setw(11) << temp_hist[j][i]
               << ",";
    }
    out_file << std::endl;
  }
  out_file.close();
  out_file.clear();

  /* Smoothing Loop */
  for (row = 0; row < m_; row++) {
    temp_hist[0][row] = this->_g_[0][row];
  }
  for (col = 1; col < n_; col++) {
    temp_hist[col] =
        KernelSmoothing(temp_hist[0], this->_g_[col], sigma, _kernel_);
  }

  /* printing loop */
  this->_g_smoothed_ = temp_hist;

  out_file.open(filename + "_g_smoothed.csv");
  out_file << std::setw(13) << "r (Å),";
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      header = this->elements()[i] + "-" + this->elements()[j] + " (1/Å),";
      out_file << std::setprecision(6) << std::setw(13) << header;
    }
  }
  out_file << std::setprecision(6) << std::setw(12) << "g(r)," << std::endl
           << std::fixed;
  for (i = 0; i < m_; i++) {
    for (j = 0; j < n_; j++) {
      out_file << std::setprecision(6) << std::setw(11) << temp_hist[j][i]
               << ",";
    }
    out_file << std::endl;
  }
  out_file.close();
  out_file.clear();

  /* Smoothing Loop */
  for (row = 0; row < m_; row++) {
    temp_hist[0][row] = this->_G_[0][row];
  }
  for (col = 1; col < n_; col++) {
    temp_hist[col] =
        KernelSmoothing(temp_hist[0], this->_G_[col], sigma, _kernel_);
  }

  /* printing loop */
  this->_g_smoothed_ = temp_hist;

  out_file.open(filename + "_G__smoothed.csv");
  out_file << std::setprecision(6) << std::setw(13) << "r (Å),";
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      header = this->elements()[i] + "-" + this->elements()[j] + " (1/Å),";
      out_file << std::setprecision(6) << std::setw(13) << header;
    }
  }
  out_file << std::setw(12) << "G(r)," << std::endl << std::fixed;
  for (i = 0; i < m_; i++) {
    for (j = 0; j < n_; j++) {
      out_file << std::setprecision(6) << std::setw(11) << temp_hist[j][i]
               << ",";
    }
    out_file << std::endl;
  }
  out_file.close();
  out_file.clear();
} // Cell::RDFSmoothing

void Cell::CoordinationNumberHistogram(std::string filename) {
  int n_, m_, i, j, col;
  int n = this->elements().size();
  std::string header;

  /* cols:
   * every element by every other element n*n
   * plus every element total n
   * plus one row for number of neighbours
   */
  n_ = n * n + n + 1;
  /* rows:
   * from 0 neighbours to n neighbours
   */
  m_ = this->coordination[0][0].size();

  std::vector<std::vector<int>> temp_hist(n_, std::vector<int>(m_, 0));
  // Fill the first n*n number of bonds values of the histogram
  for (i = 0; i < m_; i++) {
    temp_hist[0][i] = i;
  }
  col = 1;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      temp_hist[col] = this->coordination[i][j];
      col++;
    }
  }
  // fill the last n columns
  for (i = 0; i < n; i++) {
    temp_hist[col] = this->coordination[i][n];
    col++;
  }

  std::ofstream out_file(filename + "_Z.csv");
  out_file << std::setprecision(6) << std::setw(13) << "Number (#),";
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      header = this->elements()[i] + " by " + this->elements()[j] + ",";
      out_file << std::setprecision(6) << std::setw(13) << header;
    }
  }
  for (i = 0; i < n; i++) {
    header = this->elements()[i] + " by any" + ",";
    out_file << std::setprecision(6) << std::setw(13) << header;
  }
  out_file << std::endl << std::fixed;

  for (i = 0; i < m_; i++) {
    for (j = 0; j < n_; j++) {
      out_file << std::setprecision(6) << std::setw(12) << temp_hist[j][i]
               << ",";
    }
    out_file << std::endl;
  }
  out_file.close();
  out_file.clear();
} // Cell::CoordinationNumberHistogram

void Cell::VoronoiIndex() {
  for (auto &atom_A : this->_atoms_) {
    std::vector<int> atom_A_ids = atom_A.GetBondedAtomsID();
    std::sort(atom_A_ids.begin(), atom_A_ids.end());
    std::cout << "A " << atom_A.GetID() << ": (";
    for (auto A_id : atom_A_ids) {
      std::cout << A_id << ", ";
    }
    std::cout << " )" << std::endl;
  }
} // Cell::Voronoi

void Cell::SQ(std::string filename, double q_bin_width, double bin_width,
              bool normalize) {
  int n_, m_, i, j, row, col;
  int n = this->elements().size();
  std::string header;
  double Trapz;

  n_ = this->_G_.size();
  m_ = this->_G_[0].size();
  std::vector<double> temp_w_ij(n_, 1.0);
  if (!normalize) {
    temp_w_ij = this->w_ij;
  }

  /*
   * The structure factor S(q) is calculated with:
   * S(q) = 1 + 4*pi*rho_0*(q^-1)*\int{ dr r*sin(qr)*[g(r) - 1]
   * or S(q) = 1 + (q^{-1})*\int{dr sin(q * r) * G(r)}
   */
  std::vector<std::vector<double>> temp_S(n_, std::vector<double>(m_, 0));

  // Fill the q values of the histogram
  for (i = 0; i < m_; i++) {
    temp_S[0][i] = 0.0 + (i + 0.0) * q_bin_width;
  }
  double rho_0 = this->atoms().size() / this->volume;

  /*
   * Double loop on q (rows) and G_ij (cols)
   */
  for (col = 1; col < n_; col++) {
    for (row = 0; row < m_; row++) {
      /*
       * Integration with Trapezoidal_rule
       */
      Trapz = 0.0;
      if (temp_S[0][row] < 0.2) {
        for (i = 1; i < m_; i++) {
          Trapz += (1 - 0.5 * pow(temp_S[0][row] * this->_g_[0][i], 2)) *
                   this->_g_[col][i];
        }
        temp_S[col][row] = temp_w_ij[col] + Trapz * bin_width * rho_0;
      } else {
        for (i = 1; i < (m_ - 1); i++) {
          Trapz +=
              std::sin(temp_S[0][row] * this->_G_[0][i]) * this->_G_[col][i];
        }
        Trapz += 0.5 * std::sin(temp_S[0][row] * this->_G_[0][m_ - 1]) *
                 this->_G_[col][m_ - 1];
        temp_S[col][row] =
            temp_w_ij[col] + Trapz * (bin_width / temp_S[0][row]);
      }
    }
  }

  this->_S_ = temp_S;

  std::ofstream out_file(filename + "_S_q.csv");
  out_file << std::setprecision(5) << std::setw(13) << "q (1/Å),";
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      header = this->elements()[i] + "-" + this->elements()[j] + " (Å),";
      out_file << std::setprecision(5) << std::setw(13) << header;
    }
  }
  out_file << std::endl << std::fixed;

  for (i = 0; i < m_; i++) {
    for (j = 0; j < n_; j++) {
      out_file << std::setprecision(5) << std::setw(11) << temp_S[j][i] << ",";
    }
    out_file << std::endl;
  }
  out_file.close();
} // Cell::SQ

void Cell::XRD(std::string filename, double lambda, double theta_min,
               double theta_max, double bin_width) {
  int n_, m_, i, j, k, col, row;
  double norm, aux;
  int n = this->elements().size();
  std::string header;

  m_ = ceil((theta_max - theta_min) / bin_width);
  n_ = n * (n + 1) / 2 + 1;

  std::vector<std::vector<double>> temp_hist(n_, std::vector<double>(m_, 0));
  // Fill the r values of the histogram
  for (i = 0; i < m_; i++) {
    temp_hist[0][i] = theta_min + i * bin_width;
  }
  col = 0;
  // Triple loop to iterate over the distances tensor.
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      col++;
      for (const auto &it : this->distances[i][j]) {
        /*
         * Bragg's Law
         * n * lambda = 2d sin(theta)
         * theta = asin((n * lambda) / (2 * d))
         */
        aux = lambda / (it * 2.0);
        k = 1;
        std::cout << "d:" << it << std::endl;
        while ((k * aux) <= 1.0) {
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
  std::ofstream out_file(filename + "_XRD.csv");
  out_file << std::setprecision(5) << std::setw(13) << "2-theta (°),";
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      header = this->elements()[i] + "-" + this->elements()[j] + " (%),";
      out_file << std::setprecision(5) << std::setw(13) << header;
    }
  }
  out_file << std::endl << std::fixed;

  for (i = 0; i < m_; i++) {
    for (j = 0; j < n_; j++) {
      out_file << std::setprecision(5) << std::setw(12) << temp_hist[j][i]
               << ",";
    }
    out_file << std::endl;
  }

  out_file.close();
} // Cell:XRD

void Cell::PADHistogram(std::string filename, double theta_cut,
                        double bin_width) {
  int n_, m_, i, j, k, h, col, row;
  double norm;
  int n = this->elements().size();
  std::string header;

  /*
   * The number of columns in the output file is:
   *         n       x      (n+1)! / [2 x (n-1)!]
   * Central atom        Combination with repetition
   *                         of n in groups of 2
   * it's reduced to n x (n+1) x n /2
   * one extra column is added for Theta (angle)
   */

  n_ = 1 + (n * n * (n + 1) / 2);
  // from 0 to theta_cut degrees rows
  m_ = std::round(theta_cut / bin_width);
  // Matrix to store the Histograms n_ + 1 columns, m_ rows
  std::vector<std::vector<double>> temp_hist(n_ + 1,
                                             std::vector<double>(m_, 0));
  // Fill the theta values of the histogram
  for (i = 0; i < m_; i++) {
    temp_hist[0][i] = (i + 0.5) * bin_width;
  }
  col = 0;
  // Quadruple loop to iterate over the 3D angle tensor.
  for (i = 0; i < n; i++) {
    // i iterates over all central atoms
    for (j = 0; j < n; j++) {
      // j iterates over all initial atoms
      for (k = j; k < n; k++) {
        // k iterates only over half + 1 of the spectrum
        col++;
        for (const auto &it : this->angles[j][i][k]) {
          row = floor(it / bin_width);
          if (row < m_) {
            temp_hist[col][row]++;
          }
        }
        // Remove double count when j == k
        if (j == k) {
          for (h = 0; h < m_; h++) {
            temp_hist[col][h] /= 2.0;
          }
        }
      }
    }
  }
  // Double loop to find normalization factor
  norm = 0.0;
  for (i = 1; i < n_; i++) {
    for (j = 0; j < m_; j++) {
      norm += temp_hist[i][j];
    }
  }
  norm *= bin_width;

  // Double loop to normalize PADHistogram
  for (i = 1; i < n_; i++) {
    for (j = 0; j < m_; j++) {
      temp_hist[i][j] /= norm;
    }
  }

  // Double loop to calculate total PADHistogram
  for (i = 1; i < n_; i++) {
    for (j = 0; j < m_; j++) {
      temp_hist[n_][j] += temp_hist[i][j];
    }
  }

  this->_f_theta_ = temp_hist;
  std::ofstream out_file(filename + "_PAD.csv");
  out_file << std::setprecision(6) << std::setw(13) << "theta (°),";
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      for (k = j; k < n; k++) {
        header = this->elements()[j] + "-" + this->elements()[i] + "-" +
                 this->elements()[k] + ",";
        out_file << std::setprecision(6) << std::setw(12) << header;
      }
    }
  }
  out_file << std::setprecision(6) << std::setw(12) << "f(theta)," << std::endl
           << std::fixed;

  for (i = 0; i < m_; i++) {
    for (j = 0; j < n_ + 1; j++) {
      out_file << std::setprecision(6) << std::setw(11) << temp_hist[j][i]
               << ",";
    }
    out_file << std::endl;
  }

  out_file.close();
} // Cell::PADHistogram
