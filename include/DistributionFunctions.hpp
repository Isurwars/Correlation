/**
 * @file DistributionFunctions.hpp
 * @brief Manager for distribution function calculations (RDF, PAD, S(Q), etc.).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Cell.hpp"
#include "math/Smoothing.hpp"
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
  double dihedral_bin_width = 1.0;
  size_t max_ring_size = 8;
  // Maps calculator ID (e.g., "RDF", "SQ") to whether it is enabled.
  // An empty map means all calculators are enabled by default.
  std::map<std::string, bool> active_calculators;
  bool isActive(const std::string& id) const {
    if (active_calculators.empty()) return true; // default: all enabled
    auto it = active_calculators.find(id);
    return it != active_calculators.end() && it->second;
  }
  bool smoothing = true;
  double smoothing_sigma = 0.1;
  correlation::math::KernelType smoothing_kernel =
      correlation::math::KernelType::Gaussian;
};

// A structure to hold all data related to a single histogram.
struct Histogram {
  std::vector<double> bins;
  std::string title;
  std::string x_label;
  std::string y_label;
  std::string x_unit;
  std::string y_unit;
  std::string description;
  std::string file_suffix;
  // Maps a partial key (e.g., "Si-O" or "Total") to its histogram values
  std::map<std::string, std::vector<double>> partials;
  std::map<std::string, std::vector<double>> smoothed_partials;
};

/**
 * @brief Manages the calculation and storage of various distribution functions.
 *
 * This class includes methods for calculating:
 * - Radial Distribution Functions (RDF) g(r), J(r), G(r)
 * - Plane and Angle Distributions (PAD)
 * - Validation Autocorrelation Function (VACF)
 * - Vibrational Density of States (VDOS)
 * - Structure Factor S(Q)
 * - X-Ray Diffraction (XRD) patterns
 *
 * It acts as a container for all analysis results associated with a Cell or
 * Trajectory.
 */
class DistributionFunctions {
public:
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//

  /**
   * @brief Constructs a DistributionFunctions object.
   * @param cell Reference to the cell (structure) being analyzed.
   * @param cutoff Optional cutoff radius for neighbor calculations (if not pre
   * computed).
   * @param bond_cutoffs Optional bond cutoffs for neighbor calculations.
   */
  explicit DistributionFunctions(
      Cell &cell, double cutoff = 0.0,
      const std::vector<std::vector<double>> &bond_cutoffs = {});

  /**
   * @brief Move constructor.
   */
  DistributionFunctions(DistributionFunctions &&other) noexcept;

  /**
   * @brief Move assignment operator.
   */
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

  const std::map<std::string, double> &getAshcroftWeights() const {
    return ashcroft_weights_;
  }

  void addHistogram(const std::string &name, Histogram &&histogram);

  const StructureAnalyzer *neighbors() const;

  std::vector<std::string> getAvailableHistograms() const;
  //-------------------------------------------------------------------------//
  //--------------------------- Calculation Methods -------------------------//
  //-------------------------------------------------------------------------//

  /**
   * @brief Calculates the Coordination Number (CN) distribution.
   * Requires that neighbors have been computed.
   */
  void calculateCoordinationNumber();
  /**
   * @brief Calculates the radial distribution function g(r)
   * @param r_max Maximum radius to calculate up to
   * @param bin_width Width of each bin in Angstroms
   * @throws std::invalid_argument if parameters are invalid
   * @throws std::logic_error if cell volume is invalid
   */
  void calculateRDF(double r_max = 20.0, double bin_width = 0.05);

  /**
   * @brief Calculates the Plane Angle Distribution (PAD).
   * @param bin_width Width of angular bins in degrees.
   */
  void calculatePAD(double bin_width = 1.0);

  /**
   * @brief Calculates the Dihedral Angle Distribution (DAD).
   * @param bin_width Width of angular bins in degrees.
   */
  void calculateDAD(double bin_width = 1.0);

  /**
   * @brief Calculates the Velocity Autocorrelation Function (VACF).
   * @param traj The trajectory containing pre-calculated velocities.
   * @param max_correlation_frames Maximum lag frames. -1 uses default (half
   * trajectory).
   * @param start_frame Starting frame index.
   * @param end_frame Ending frame index (exclusive).
   */
  void calculateVACF(const Trajectory &traj, int max_correlation_frames = -1,
                     size_t start_frame = 0,
                     size_t end_frame = static_cast<size_t>(-1));

  /**
   * @brief Calculates the Vibrational Density of States (VDOS) from the VACF.
   * Requires calculateVACF to be called first.
   */
  void calculateVDOS();

  /**
   * @brief Calculates the XRD pattern.
   * @param lambda X-ray wavelength in Angstroms.
   * @param theta_min Minimum 2-theta angle.
   * @param theta_max Maximum 2-theta angle.
   * @param bin_width Bin width for 2-theta.
   */
  void calculateXRD(double lambda = 1.5406, double theta_min = 5.0,
                    double theta_max = 90.0, double bin_width = 1.0);

  /**
   * @brief Smooths a specific histogram.
   * @param name Name of the histogram (e.g., "g(r)").
   * @param sigma Smoothing width.
   * @param kernel The smoothing kernel type.
   */
  void smooth(const std::string &name, double sigma,
              correlation::math::KernelType kernel =
                  correlation::math::KernelType::Gaussian);

  /**
   * @brief Smooths all available histograms.
   * @param sigma Smoothing width.
   * @param kernel The smoothing kernel type.
   */
  void smoothAll(double sigma,
                 correlation::math::KernelType kernel =
                     correlation::math::KernelType::Gaussian);

  /**
   * @brief Uses an external StructureAnalyzer for neighborhood/bond info.
   * @param analyzer A pre-computed external analyzer instance.
   */
  void setStructureAnalyzer(const StructureAnalyzer *analyzer);

  /**
   * @brief Takes ownership of a StructureAnalyzer for neighborhood/bond info.
   * @param analyzer A pre-computed analyzer instance.
   */
  void setStructureAnalyzerOwned(std::unique_ptr<StructureAnalyzer> analyzer);

  //-------------------------------------------------------------------------//
  //----------------------------- Accumulation ------------------------------//
  //-------------------------------------------------------------------------//

  /**
   * @brief Adds the histograms from another DistributionFunctions object to
   * this one. Used for averaging over multiple frames.
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
      std::function<void(float, const std::string &)> progress_callback =
          nullptr);

private:
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
