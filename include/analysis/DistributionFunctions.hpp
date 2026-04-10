/**
 * @file DistributionFunctions.hpp
 * @brief Manager for distribution function calculations (RDF, PAD, S(Q), etc.).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "analysis/StructureAnalyzer.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "math/Smoothing.hpp"

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace correlation::analysis {

class TrajectoryAnalyzer;

/**
 * @brief Configuration settings for distribution function analysis.
 */
struct AnalysisSettings {
  double r_max = 20.0;           ///< Maximum radius for RDF calculations (Angstroms).
  double r_bin_width = 0.02;     ///< Bin width for radial distributions (Angstroms).
  double q_max = 20.0;           ///< Maximum momentum transfer for S(Q) (Angstroms^-1).
  double q_bin_width = 0.02;     ///< Bin width for S(Q) (Angstroms^-1).
  double r_int_max = 10.0;       ///< Cutoff for integration-based properties.
  double angle_bin_width = 1.0;  ///< Bin width for bond angle distributions (degrees).
  double dihedral_bin_width = 1.0; ///< Bin width for dihedral distributions (degrees).
  size_t max_ring_size = 8;      ///< Maximum size of rings to search for.

  /// Maps calculator ID (e.g., "RDF", "SQ") to whether it is enabled.
  /// An empty map means all calculators are enabled by default.
  std::map<std::string, bool> active_calculators;

  /**
   * @brief Helper to check if a specific calculator is enabled.
   * @param id The identifier of the calculator (e.g. "RDF").
   * @return True if active or if no active calculators are specified.
   */
  bool isActive(const std::string &id) const {
    if (active_calculators.empty())
      return true; // default: all enabled
    auto it = active_calculators.find(id);
    return it != active_calculators.end() && it->second;
  }

  bool smoothing = true;         ///< Whether to apply post-processing smoothing.
  double smoothing_sigma = 0.1;  ///< Gaussian smoothing standard deviation.
  correlation::math::KernelType smoothing_kernel =
      correlation::math::KernelType::Gaussian; ///< The kernel to use for smoothing.
};

// A structure to hold all data related to a single histogram.
/**
 * @brief Container for a single calculated distribution function.
 */
struct Histogram {
  std::vector<double> bins;      ///< The x-axis values (radii, angles, etc.).
  std::string title;             ///< Descriptive title for the plot.
  std::string x_label;           ///< Label for the x-axis.
  std::string y_label;           ///< Label for the y-axis.
  std::string x_unit;            ///< Physical unit for the x-axis (e.g. "A").
  std::string y_unit;            ///< Physical unit for the y-axis (e.g. "A^-3").
  std::string description;       ///< Internal description of what this data represents.
  std::string file_suffix;       ///< Default suffix for saving this histogram to disk.

  /// Maps a partial key (e.g., "Si-O" or "Total") to its histogram values.
  std::map<std::string, std::vector<double>> partials;

  /// Maps a partial key to its smoothed histogram values.
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
 * It acts as a container for all analysis results associated with a
 * correlation::core::Cell or correlation::core::Trajectory.
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
      correlation::core::Cell &cell, double cutoff = 0.0,
      const std::vector<std::vector<double>> &bond_cutoffs = {});

  /**
   * @brief Move constructor.
   * @param other The DistributionFunctions object to move from.
   */
  DistributionFunctions(DistributionFunctions &&other) noexcept;

  /**
   * @brief Move assignment operator.
   * @param other The DistributionFunctions object to move from.
   * @return Reference to this object.
   */
  DistributionFunctions &operator=(DistributionFunctions &&other) noexcept;

  //-------------------------------------------------------------------------//
  //------------------------------- Accessors -------------------------------//
  //-------------------------------------------------------------------------//

  /**
   * @brief Access the underlying simulation cell.
   * @return Constant reference to the Cell object.
   */
  const correlation::core::Cell &cell() const { return cell_; }

  /**
   * @brief Access a specific calculated histogram by name.
   * @param name The name of the distribution function (e.g., "g(r)", "S(Q)").
   * @return A const reference to the Histogram data.
   * @throws std::out_of_range if the name is not found.
   */
  const Histogram &getHistogram(const std::string &name) const;

  /**
   * @brief Access all calculated histograms.
   * @return A map of histogram names to Histogram objects.
   */
  const std::map<std::string, Histogram> &getAllHistograms() const {
    return histograms_;
  }

  /**
   * @brief Gets the Ashcroft-Langreth weights used for S(Q) partials.
   * @return A map of element pair strings to weight values.
   */
  const std::map<std::string, double> &getAshcroftWeights() const {
    return ashcroft_weights_;
  }

  /**
   * @brief Manually add a histogram to the collection.
   * @param name Name of the histogram.
   * @param histogram The Histogram data to move into the collection.
   */
  void addHistogram(const std::string &name, Histogram &&histogram);

  /**
   * @brief Returns the structure analyzer providing neighbor information.
   * @return Pointer to the current StructureAnalyzer.
   */
  const StructureAnalyzer *neighbors() const;

  /**
   * @brief Gets a list of names for all currently available histograms.
   * @return Vector of histogram name strings.
   */
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
  void calculateVACF(const correlation::core::Trajectory &traj,
                     int max_correlation_frames = -1, size_t start_frame = 0,
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
  void smoothAll(double sigma, correlation::math::KernelType kernel =
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
   * @param trajectory The trajectory to analyze.
   * @param analyzer A pre-configured TrajectoryAnalyzer for neighbor info.
   * @param start_frame Index of the frame to start analysis from.
   * @param settings The analysis parameters (bin widths, cutoffs, etc.).
   * @param progress_callback Optional callback to report completion progress.
   * @return A unique_ptr to the newly created and populated DistributionFunctions object.
   */
  static std::unique_ptr<DistributionFunctions> computeMean(
      correlation::core::Trajectory &trajectory,
      const TrajectoryAnalyzer &analyzer, size_t start_frame,
      const AnalysisSettings &settings,
      std::function<void(float, const std::string &)> progress_callback =
          nullptr);

private:
  /**
   * @brief Ensures neighbors are computed for the given cutoff.
   * @param r_cut Cutoff radius.
   */
  void ensureNeighborsComputed(double r_cut);

  /**
   * @brief Generates a unique string key for a pair of element types.
   * @param type1 First element ID.
   * @param type2 Second element ID.
   * @return Canonical string key (e.g. "Si-O").
   */
  std::string getPartialKey(int type1, int type2) const;

  /**
   * @brief Generates the inverse string key for a pair of element types.
   * @param type1 First element ID.
   * @param type2 Second element ID.
   * @return Inverse string key (e.g. "O-Si").
   */
  std::string getInversePartialKey(int type1, int type2) const;

  /**
   * @brief Calculates weights for partial distributions based on concentrations.
   */
  void calculateAshcroftWeights();

  correlation::core::Cell &cell_; ///< Reference to the cell being analyzed.
  const StructureAnalyzer *neighbors_ref_{nullptr}; ///< Pointer to external neighbor info.
  std::unique_ptr<StructureAnalyzer> neighbors_owned_; ///< Owned neighbor info.
  double current_cutoff_{-1.0}; ///< Last cutoff used for neighbor searching.
  std::vector<std::vector<double>> bond_cutoffs_sq_; ///< Cached squared bond cutoffs.

  std::map<std::string, Histogram> histograms_; ///< Storage for all analysis results.
  std::map<std::string, double> ashcroft_weights_; ///< Scalar weights for S(Q) calculations.
};

} // namespace correlation::analysis
