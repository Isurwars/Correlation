// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Cell.hpp"
#include "Smoothing.hpp"
#include "StructureAnalyzer.hpp"

class Trajectory;
class TrajectoryAnalyzer;

struct AnalysisSettings {
  double r_max = 20.0;
  double r_bin_width = 0.02;
  double q_max = 20.0;
  double q_bin_width = 0.02;
  double r_int_max = 10.0;
  double angle_bin_width = 1.0;
  bool smoothing = true;
  double smoothing_sigma = 0.1;
  KernelType smoothing_kernel = KernelType::Gaussian;
};

// A structure to hold all data related to a single histogram.
struct Histogram {
  std::vector<double> bins;
  std::string bin_label;
  // Maps a partial key (e.g., "Si-O" or "Total") to its histogram values
  std::map<std::string, std::vector<double>> partials;
  std::map<std::string, std::vector<double>> smoothed_partials;
};

class DistributionFunctions {
public:
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//

  explicit DistributionFunctions(
      Cell &, double = 0.0,
      const std::vector<std::vector<double>> & = {});

  // Move constructor
  DistributionFunctions(DistributionFunctions &&other) noexcept;

  // Move assignment operator
  DistributionFunctions &operator=(DistributionFunctions &&other) noexcept;

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
  const Histogram &getHistogram(const std::string &name) const;

  const std::map<std::string, Histogram> &getAllHistograms() const {
    return histograms_;
  }

  std::vector<std::string> getAvailableHistograms() const;
  //-------------------------------------------------------------------------//
  //--------------------------- Calculation Methods -------------------------//
  //-------------------------------------------------------------------------//

  void calculateCoordinationNumber();
  /**
   * @brief Calculates the radial distribution function g(r)
   * @param r_max Maximum radius to calculate up to
   * @param r_bin_width Width of each bin in Angstroms
   * @throws std::invalid_argument if parameters are invalid
   * @throws std::logic_error if cell volume is invalid
   */
  void calculateRDF(double r_max = 20.0, double bin_width = 0.05);

  void calculatePAD(double bin_width = 1.0);

  /**
   * @brief Calculates the Velocity Autocorrelation Function (VACF).
   * @param max_correlation_frames Maximum lag frames. -1 uses default (half trajectory).
   */
  void calculateVACF(const Trajectory &traj, int max_correlation_frames = -1);

  void calculateVDOS();

  void calculateSQ(double q_max = 25.0, double q_bin_width = 0.05,
                   double r_integration_max = 8.0);

  void calculateXRD(double lambda = 1.5406, double theta_min = 5.0,
                    double theta_max = 90.0, double bin_width = 1.0);

  void smooth(const std::string &name, double sigma,
              KernelType kernel = KernelType::Gaussian);
  void smoothAll(double sigma, KernelType kernel = KernelType::Gaussian);

  void setStructureAnalyzer(const StructureAnalyzer *analyzer);

  //-------------------------------------------------------------------------//
  //----------------------------- Accumulation ------------------------------//
  //-------------------------------------------------------------------------//
  
  /**
   * @brief Adds the histograms from another DistributionFunctions object to this one.
   *        Used for averaging over multiple frames.
   * @param other The other DistributionFunctions object to add.
   */
  void add(const DistributionFunctions &other);

  /**
   * @brief Scales all histograms in this object by a factor.
   *        Used for normalizing averaged distributions.
   * @param factor The scaling factor.
   */
  void scale(double factor);

  /**
   * @brief Computes the mean distribution functions over a trajectory.
   *        Uses parallel execution to speed up calculation.
   */
  static std::unique_ptr<DistributionFunctions> computeMean(
      Trajectory &trajectory, const TrajectoryAnalyzer &analyzer,
      size_t start_frame, const AnalysisSettings &settings,
      std::function<void(float)> progress_callback = nullptr);

private:
  const StructureAnalyzer *neighbors() const;
  void ensureNeighborsComputed(double r_cut);
  std::string getPartialKey(int type1, int type2) const;
  std::string getInversePartialKey(int type1, int type2) const;
  void calculateAshcroftWeights();

  Cell &cell_;
  const StructureAnalyzer *neighbors_ref_{nullptr};
  std::unique_ptr<StructureAnalyzer> neighbors_owned_;
  double current_cutoff_{-1.0};
  std::vector<std::vector<double>> bond_cutoffs_sq_;

  std::map<std::string, Histogram> histograms_;
  std::map<std::string, double> ashcroft_weights_;
};
