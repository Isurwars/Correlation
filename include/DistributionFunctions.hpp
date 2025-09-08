#ifndef INCLUDE_DISTRIBUTION_FUNCTIONS_HPP_
#define INCLUDE_DISTRIBUTION_FUNCTIONS_HPP_
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <map>
#include <string>
#include <vector>

#include "Cell.hpp"
#include "NeighborList.hpp"
#include "Smoothing.hpp"

// A structure to hold all data related to a single histogram.
// This includes the x-axis values (bins) and the y-axis values for
// both the raw and smoothed data.
struct Histogram {
  std::vector<double> bins;
  // Maps a partial key (e.g., "Si-O" or "Total") to its histogram values
  std::map<std::string, std::vector<double>> partials;
  std::map<std::string, std::vector<double>> smoothed_partials;
};

class DistributionFunctions {
public:
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//

  explicit DistributionFunctions(const Cell &cell,
                                 const NeighborList &neighbors);

  //-------------------------------------------------------------------------//
  //------------------------------- Accessors -------------------------------//
  //-------------------------------------------------------------------------//
  const Cell &cell() const { return cell_; }

  /**
   * @brief Access a specific calculated histogram by name.
   * @param name The name of the distribution function (e.g., "g(r)", "S(Q)").
   * @return A const reference to the Histogram data.
   * @throws std::out_of_range if the name is not found.
   */
  const Histogram &getHistogram(const std::string &name) const {
    if (histograms_.count(name) == 0) {
      throw std::out_of_range("Histogram not found: " + name);
    }
    return histograms_.at(name);
  }

  const std::map<std::string, Histogram> &getAllHistograms() const {
    return histograms_;
  }

  //-------------------------------------------------------------------------//
  //--------------------------- Calculation Methods -------------------------//
  //-------------------------------------------------------------------------//

  void calculateCoordinationNumber();

  void calculateRDF(double r_cut = 20.0, double bin_width = 0.05,
                    bool normalize = false);

  void calculatePAD(double theta_cut = 180.0, double bin_width = 1.0);

  void calculateSQ(double q_max = 25.0, double q_bin_width = 0.05,
                   double r_integration_max = 8.0);

  void calculateXRD(double lambda = 1.5406, double theta_min = 5.0,
                    double theta_max = 90.0, double bin_width = 1.0);

  /**
   * @brief Applies kernel smoothing to a specified distribution function.
   * @param name The name of the histogram to smooth (e.g., "g(r)").
   * @param sigma The standard deviation (width) of the Gaussian kernel.
   */
  void smooth(const std::string &, double, KernelType);
  void smoothAll(double, KernelType = KernelType::Gaussian);

  // Provides a list of keys for all computed histograms.
  std::vector<std::string> getAvailableHistograms() const;

private:
  // Helper to create a key for partial maps (e.g., "Si-O")
  std::string getPartialKey(int type1, int type2) const;

  const Cell &cell_;
  const NeighborList &neighbors_;

  // A map to store all calculated histograms by name.
  // This is more flexible and scalable than individual member variables.
  std::map<std::string, Histogram> histograms_;

  // Stores coordination number data separately as it's not a typical histogram
  std::vector<std::vector<double>> coordination_data_;
};
#endif // INCLUDE_DISTRIBUTION_FUNCTIONS_HPP_
