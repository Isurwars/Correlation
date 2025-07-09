// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "../include/DistributionFunctions.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "../include/Constants.hpp"
#include "../include/Smoothing.hpp"

//---------------------------------------------------------------------------//
//------------------------ Cell Action Constructors -------------------------//
//---------------------------------------------------------------------------//

DistributionFunctions::DistributionFunctions(const Cell &cell) : cell_(cell) {}

//--------------------------------------------------------------------------//
//--------------------------------- Methods --------------------------------//
//--------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
//--------------------------- Coordination Number ---------------------------//
//---------------------------------------------------------------------------//
void DistributionFunctions::coordinationNumber() {
  const size_t num_elements = cell_.distances().size();
  if (num_elements == 0)
    throw std::logic_error("No elements in cell");
  const size_t num_atoms = cell_.atoms().size();
  if (num_atoms == 0)
    throw std::logic_error("No atoms in cell");

  // Step 1: Find maximum coordination number
  int max_Nc = 0;
  for (auto &atom : cell_.atoms()) {
    max_Nc = std::max(max_Nc, static_cast<int>(atom.bonded_atoms().size()));
  }

  // Step 2: Initialize 3D tensor [element][element/total][coordination]
  const size_t coordination_bins = max_Nc + 2;

  // Dimensions: [element][target (elements + total)][coordination count]
  using ElementHist = std::vector<std::vector<double>>;
  std::vector<ElementHist> temp_nc(
      num_elements,
      ElementHist(num_elements + 1, std::vector<double>(coordination_bins, 0)));

  // Step 3: Build coordination histogram
  for (auto &atom : cell_.atoms()) {
    const size_t element_id = atom.element_id();
    std::vector<int> element_counts(num_elements, 0);

    // Count bonded atoms per element
    for (auto &bonded_atom : atom.bonded_atoms()) {
      const size_t bonded_element_id = bonded_atom.element_id();
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

  std::vector<std::vector<double>> coordination_hist(
      num_columns, std::vector<double>(num_rows, 0));

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

  Z_ = std::move(coordination_hist);
} // DistributionFunctions::coorcinationNumber

//---------------------------------------------------------------------------//
//--------------------- Radial Distribution Functions -----------------------//
//---------------------------------------------------------------------------//
void DistributionFunctions::calculateRDF(double r_cut, double bin_width,
                                         bool normalize) {
  // Validate input parameters
  if (bin_width <= 0)
    throw std::invalid_argument("Bin width must be positive");
  if (r_cut <= 0)
    throw std::invalid_argument("Cutoff radius must be positive");
  if (cell_.volume() <= 0)
    throw std::logic_error("Cell volume must be positive");

  const size_t num_elements = cell_.elements().size();
  if (num_elements == 0)
    throw std::logic_error("No elements in cell");
  const size_t num_atoms = cell_.atoms().size();
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
      const double ni = cell_.element_numbers()[i];
      for (size_t j = i; j < num_elements; ++j) {
        const double nj = cell_.element_numbers()[j];
        double weight = 2.0 * ni * nj * norm_factor;
        if (i == j)
          weight *= 0.5;
        weight_factors[weight_index++] = weight;
      }
    }
  }
  cell_.setWeightFactors(weight_factors);

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
      for (const auto &distance : cell_.distances()[i][j]) {
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
  const double volume_factor =
      4.0 * constants::pi * (num_atoms / cell_.volume());
  const double norm_factor = 1.0 / (num_atoms * bin_width);

  for (size_t col = 1; col < num_columns; ++col) {
    const double weight = normalize ? 1.0 : weight_factors[col];
    const double scaling = norm_factor / weight;

    for (size_t bin = 0; bin < num_bins; ++bin) {
      histogram[col][bin] *= scaling;
    }
  }
  J_ = histogram;

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

  g_ = std::move(g_r);
  G_ = std::move(G_r);

} // DistributionFunctions::radialDistributionFunctions

//---------------------------------------------------------------------------//
//------------------------ Plane Angle Distribution -------------------------//
//---------------------------------------------------------------------------//
void DistributionFunctions::calculatePAD(double theta_cut, double bin_width) {
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
  const size_t num_elements = cell_.elements().size();
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
             cell_.angles()[elem_j][central_elem][elem_k]) {
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
    F_ = std::move(histogram);
    return;
  }

  // Normalize and calculate total distribution
  for (size_t col = 1; col <= num_combinations; ++col) {
    for (size_t bin = 0; bin < num_bins; ++bin) {
      histogram[col][bin] /= normalization;
      histogram.back()[bin] += histogram[col][bin];
    }
  }

  F_ = std::move(histogram);
} // DistributionFunctions::calculatePAD

//---------------------------------------------------------------------------//
//---------------------------------- S(Q) -----------------------------------//
//---------------------------------------------------------------------------//
void DistributionFunctions::calculateSQ(double q_max, double q_bin_width,
                                        double r_max, bool normalize) {
  // Validate input parameters
  if (q_max <= 0.0)
    throw std::invalid_argument("q_max must be positive");
  if (q_bin_width <= 0.0)
    throw std::invalid_argument("q_bin_width must be positive");
  if (G_.empty() || G_[0].empty())
    throw std::logic_error("G(r) data not initialized");

  const size_t num_elements = cell_.elements().size();
  const size_t num_atoms = cell_.atoms().size();
  const size_t num_q_bins =
      static_cast<size_t>(std::round(q_max / q_bin_width));
  const size_t num_columns = 1 + (num_elements * (num_elements + 1)) / 2 + 1;

  // Initialize weights and S(Q) matrix
  std::vector<double> weights(num_columns, 1.0);
  std::vector<std::vector<double>> S_q(num_columns,
                                       std::vector<double>(num_q_bins, 0.0));

  // Set q values in first row
  std::generate(S_q[0].begin(), S_q[0].end(),
                [n = 0, q_bin_width]() mutable { return (n++) * q_bin_width; });

  // Calculate weights if needed
  if (!normalize) {
    const double norm_factor = 1.0 / (num_atoms * num_atoms);
    size_t weight_idx = 1;

    for (size_t i = 0; i < num_elements; ++i) {
      const double ni = cell_.element_numbers()[i];
      for (size_t j = i; j < num_elements; ++j) {
        double weight = 2.0 * ni * cell_.element_numbers()[j] * norm_factor;
        if (i == j)
          weight *= 0.5;
        weights[weight_idx++] = weight;
      }
    }
    cell_.setWeightFactors(weights);
  }

  // Get integration parameters
  const auto &r_values = G_[0];
  const size_t num_r_points = r_values.size();
  if (num_r_points < 2)
    throw std::logic_error("Insufficient r points for integration");

  const double dr = r_values[1] - r_values[0];

  // Precompute trapezoidal weights
  std::vector<double> trapz_weights(num_r_points, dr);
  trapz_weights.front() = 0.5 * dr;
  trapz_weights.back() = 0.5 * dr;

  // Precompute weighted G(r) for each column (skip r-values column)
  std::vector<std::vector<double>> weighted_rG(
      num_columns - 1, std::vector<double>(num_r_points));
  for (size_t col = 1; col < num_columns; ++col) {
    for (size_t r_idx = 0; r_idx < num_r_points; ++r_idx) {
      if (r_values[r_idx] <= r_max) {
        weighted_rG[col - 1][r_idx] =
            trapz_weights[r_idx] * r_values[r_idx] * G_[col][r_idx];
      }
    }
  }

  // Precompute sinc arguments table (memory optimization)
  std::vector<std::vector<double>> qr_table(num_q_bins,
                                            std::vector<double>(num_r_points));
  for (size_t q_idx = 0; q_idx < num_q_bins; ++q_idx) {
    const double q = S_q[0][q_idx];
    for (size_t r_idx = 0; r_idx < num_r_points; ++r_idx) {
      qr_table[q_idx][r_idx] = q * r_values[r_idx];
    }
  }

// Main computation loop - parallelized
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (size_t q_idx = 0; q_idx < num_q_bins; ++q_idx) {
    for (size_t col = 1; col < num_columns; ++col) {
      double integral = 0.0;
      const auto &wG_col = weighted_rG[col - 1];

      for (size_t r_idx = 0; r_idx < num_r_points; ++r_idx) {
        const double x = qr_table[q_idx][r_idx];
        double sinc_val = (std::abs(x) < 1e-8) ? 1.0 : std::sin(x) / x;
        integral += sinc_val * wG_col[r_idx];
      }
      S_q[col][q_idx] = weights[col] + integral;
    }
  }

  S_ = std::move(S_q);
} // DistributionFunctions::calculateSQ

//---------------------------------------------------------------------------//
//------------------------------- S(Q) Dirac --------------------------------//
//---------------------------------------------------------------------------//
void DistributionFunctions::calculateSQDirac(double q_max, double q_bin_width,
                                             bool normalize) {
  // Validate input parameters
  if (q_max <= 0.0)
    throw std::invalid_argument("q_max must be positive");
  if (q_bin_width <= 0.0)
    throw std::invalid_argument("q_bin_width must be positive");

  const size_t num_elements = cell_.elements().size();
  const size_t num_atoms = cell_.atoms().size();
  const size_t num_q_bins =
      static_cast<size_t>(std::round(q_max / q_bin_width));
  const size_t num_columns = 1 + (num_elements * (num_elements + 1)) / 2 + 1;

  // Initialize weights and S(Q) matrix
  std::vector<double> weights(num_columns, 1.0);
  std::vector<std::vector<double>> S_q(num_columns,
                                       std::vector<double>(num_q_bins, 0.0));

  // Set q values in first row
  std::generate(S_q[0].begin(), S_q[0].end(),
                [n = 0, q_bin_width]() mutable { return (n++) * q_bin_width; });

  // Calculate weights if needed
  if (!normalize) {
    const double norm_factor = 1.0 / (num_atoms * num_atoms);
    size_t weight_idx = 1;

    for (size_t i = 0; i < num_elements; ++i) {
      const double ni = cell_.element_numbers()[i];
      for (size_t j = i; j < num_elements; ++j) {
        double weight = 2.0 * ni * cell_.element_numbers()[j] * norm_factor;
        if (i == j)
          weight *= 0.5;
        weights[weight_idx++] = weight;
      }
    }
    cell_.setWeightFactors(weights);
  }

  /*
   * Dirac Formula is a direct calculation over atom positions
   * S(q) = (1/N) <Sum_i>
   */

} // calculateSQDirac

//---------------------------------------------------------------------------//
//----------------------------------- XRD -----------------------------------//
//---------------------------------------------------------------------------//
void DistributionFunctions::calculateXRD(double lambda, double theta_min,
                                         double theta_max, double bin_width) {
  int n_, m_, i, j, k, col, row;
  double norm, aux;
  int n = cell_.elements().size();

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
      for (const auto &it : cell_.distances()[i][j]) {
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

  X_ = temp_hist;
} // DistributionFunctions:calculateXRD

//--------------------------------------------------------------------------//
//------------------------------- Smoothing --------------------------------//
//--------------------------------------------------------------------------//
void DistributionFunctions::Smoothing(double sigma, int _kernel_) {
  size_t row, col;
  // n_: number of columns in the histogram
  size_t n_ = g_.size();
  // m_: number of rows in the histogram
  size_t m_ = g_[0].size();
  // Array of smoothed histograms
  std::vector<std::vector<double>> temp_hist(n_, std::vector<double>(m_, 0));

  for (row = 0; row < m_; row++) {
    temp_hist[0][row] = J_[0][row];
  }
  for (col = 1; col < n_; col++) {
    temp_hist[col] = KernelSmoothing(temp_hist[0], J_[col], sigma, _kernel_);
  }
  J_smoothed_ = temp_hist;

  // Smoothing g
  for (row = 0; row < m_; row++) {
    temp_hist[0][row] = g_[0][row];
  }
  for (col = 1; col < n_; col++) {
    temp_hist[col] = KernelSmoothing(temp_hist[0], g_[col], sigma, _kernel_);
  }
  g_smoothed_ = temp_hist;

  // Smoothing G
  for (row = 0; row < m_; row++) {
    temp_hist[0][row] = G_[0][row];
  }
  for (col = 1; col < n_; col++) {
    temp_hist[col] = KernelSmoothing(temp_hist[0], G_[col], sigma, _kernel_);
  }
  G_smoothed_ = temp_hist;

  // Smoothing F
  // n_: number of columns in the histogram
  n_ = F_.size();
  // m_: number of rows in the histogram
  m_ = F_[0].size();
  // Resize 2D vector of smoothed histograms
  temp_hist.resize(n_);
  for (auto &col_hist : temp_hist) {
    col_hist.resize(m_, 0); // resize and fill with 0
  }
  for (row = 0; row < m_; row++) {
    temp_hist[0][row] = F_[0][row];
  }
  for (col = 1; col < n_; col++) {
    temp_hist[col] = KernelSmoothing(temp_hist[0], F_[col], sigma, _kernel_);
  }
  F_smoothed_ = temp_hist;

  // Smoothing S
  // n_: number of columns in the histogram
  n_ = S_.size();
  // m_: number of rows in the histogram
  m_ = S_[0].size();
  // Resize 2D vector of smoothed histograms
  temp_hist.resize(n_);
  for (auto &col_hist : temp_hist) {
    col_hist.resize(m_, 0); // resize and fill with 0
  }
  for (row = 0; row < m_; row++) {
    temp_hist[0][row] = S_[0][row];
  }
  for (col = 1; col < n_; col++) {
    temp_hist[col] = KernelSmoothing(temp_hist[0], S_[col], sigma, _kernel_);
  }
  S_smoothed_ = temp_hist;

  // Smoothing X
  // n_: number of columns in the histogram
  n_ = X_.size();
  // m_: number of rows in the histogram
  m_ = X_[0].size();
  // Resize 2D vector of smoothed histograms
  temp_hist.resize(n_);
  for (auto &col_hist : temp_hist) {
    col_hist.resize(m_, 0); // resize and fill with 0
  }
  for (row = 0; row < m_; row++) {
    temp_hist[0][row] = X_[0][row];
  }
  for (col = 1; col < n_; col++) {
    temp_hist[col] = KernelSmoothing(temp_hist[0], X_[col], sigma, _kernel_);
  }
  X_smoothed_ = temp_hist;

} // DistributionFunctions::RDFSmoothing

//--------------------------------------------------------------------------//
//------------------------------ Voronoi Index -----------------------------//
//--------------------------------------------------------------------------//
void DistributionFunctions::voronoiIndex() {
  /*
   * This function is currently and STUB and gives a rought stimate
   * to the actual voronoiindex, to compute the correct voronoi index
   * you need to compute the voronoi tesselation and then get the index
   * this is currently planned for v2.0 of the code, as the entire algorithm
   * should be generated from scratch.
   */

} // DistributionFunctions::Voronoi
