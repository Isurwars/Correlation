// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "DistributionFunctions.hpp"

#include <algorithm>
#include <stdexcept>

#include "Smoothing.hpp"
#include "Trajectory.hpp"
#include "TrajectoryAnalyzer.hpp"
#include "calculators/CNCalculator.hpp"
#include "calculators/DADCalculator.hpp"
#include "calculators/MDCalculator.hpp"
#include "calculators/PADCalculator.hpp"
#include "calculators/RDFCalculator.hpp"
#include "calculators/SQCalculator.hpp"
#include "calculators/VACFCalculator.hpp"
#include "calculators/VDOSCalculator.hpp"
#include "calculators/XRDCalculator.hpp"

#include <atomic>
#include <future>
#include <mutex>

//---------------------------------------------------------------------------//
//------------------------------- Constructor -------------------------------//
//---------------------------------------------------------------------------//

DistributionFunctions::DistributionFunctions(
    Cell &cell, double cutoff,
    const std::vector<std::vector<double>> &bond_cutoffs)
    : cell_(cell), neighbors_owned_(nullptr), current_cutoff_(0.0),
      bond_cutoffs_sq_(bond_cutoffs) {
  if (cutoff > 0.0) {
    if (bond_cutoffs.empty()) {
      throw std::invalid_argument(
          "Bond cutoffs must be provided if cutoff > 0.0");
    }
    ensureNeighborsComputed(cutoff);
  }
  calculateAshcroftWeights();
}

// Move Constructor
DistributionFunctions::DistributionFunctions(
    DistributionFunctions &&other) noexcept
    : cell_(other.cell_), neighbors_ref_(other.neighbors_ref_),
      neighbors_owned_(std::move(other.neighbors_owned_)),
      current_cutoff_(other.current_cutoff_),
      histograms_(std::move(other.histograms_)),
      ashcroft_weights_(std::move(other.ashcroft_weights_)) {

  other.current_cutoff_ = -1.0;
  other.neighbors_ref_ = nullptr;
}

// Move Assignment Operator
DistributionFunctions &
DistributionFunctions::operator=(DistributionFunctions &&other) noexcept {
  if (this != &other) {
    neighbors_ref_ = other.neighbors_ref_;
    neighbors_owned_ = std::move(other.neighbors_owned_);
    current_cutoff_ = other.current_cutoff_;
    histograms_ = std::move(other.histograms_);
    ashcroft_weights_ = std::move(other.ashcroft_weights_);

    other.current_cutoff_ = -1.0;
    other.neighbors_ref_ = nullptr;
  }
  return *this;
}

//--------------------------------------------------------------------------//
//---------------------------- Helper Functions ----------------------------//
//--------------------------------------------------------------------------//

void DistributionFunctions::setStructureAnalyzer(
    const StructureAnalyzer *analyzer) {
  neighbors_ref_ = analyzer;
  if (neighbors_ref_) {
    neighbors_owned_.reset();
  }
}

const StructureAnalyzer *DistributionFunctions::neighbors() const {
  if (neighbors_ref_) {
    return neighbors_ref_;
  }
  return neighbors_owned_.get();
}

const Histogram &
DistributionFunctions::getHistogram(const std::string &name) const {
  auto it = histograms_.find(name);
  if (it == histograms_.end()) {
    throw std::out_of_range("Histogram '" + name + "' not found.");
  }
  return it->second;
}

void DistributionFunctions::ensureNeighborsComputed(double r_max) {
  // If we have an external reference, we assume it's valid and sufficient.
  if (neighbors_ref_) {
    return;
  }

  if (!neighbors_owned_ || r_max > current_cutoff_) {
    neighbors_owned_ = std::make_unique<StructureAnalyzer>(
        cell_, r_max, bond_cutoffs_sq_, true);
    current_cutoff_ = r_max;
  }
}

std::vector<std::string> DistributionFunctions::getAvailableHistograms() const {
  std::vector<std::string> keys;
  keys.reserve(histograms_.size());
  // Iterate through the map of histograms and extract the key for each entry.
  for (const auto &[key, val] : histograms_) {
    keys.push_back(key);
  }
  return keys;
}

std::string DistributionFunctions::getPartialKey(int type1, int type2) const {
  const auto &elements = cell_.elements();
  // Ensure consistent ordering for pairs (e.g., Si-O is same as O-Si)
  if (type1 > type2)
    std::swap(type1, type2);
  return elements[type1].symbol + "-" + elements[type2].symbol;
}

std::string DistributionFunctions::getInversePartialKey(int type1,
                                                        int type2) const {
  const auto &elements = cell_.elements();
  // Ensure consistent ordering for pairs (e.g., Si-O is same as O-Si)
  if (type1 < type2)
    std::swap(type1, type2);
  return elements[type1].symbol + "-" + elements[type2].symbol;
}

void DistributionFunctions::calculateAshcroftWeights() {
  const auto &atoms = cell_.atoms();
  if (atoms.empty()) {
    throw std::invalid_argument("Cell must contain atoms");
  }
  const size_t num_atoms = atoms.size();

  // 1. Count the number of atoms for each element symbol.
  std::map<std::string, int> element_counts;
  for (const auto &atom : atoms) {
    element_counts[atom.element().symbol]++;
  }

  // 2. Get the list of unique elements present.
  const auto &elements = cell_.elements();

  // 3. Calculate w_ij = (N_i * N_j) / N_total^2 for all pairs.
  // This is the concentration-based weighting factor.
  for (size_t i = 0; i < elements.size(); ++i) {
    for (size_t j = i; j < elements.size(); ++j) {
      const auto &element_i = elements[i];
      const auto &element_j = elements[j];

      const double count_i =
          static_cast<double>(element_counts.at(element_i.symbol));
      const double count_j =
          static_cast<double>(element_counts.at(element_j.symbol));

      double weight = (count_i * count_j) / (num_atoms * num_atoms);
      if (i != j) {
        weight *= 2.0;
      }

      // Get the canonical key (e.g., "Si-O", not "O-Si")
      std::string key = getPartialKey(element_i.id.value, element_j.id.value);
      ashcroft_weights_[key] = weight;
    }
  }
  // These weights are crucial for the Faber-Ziman formalism in S(Q)
  // calculation.
}

//---------------------------------------------------------------------------//
//---------------------------- Smoothing Methods ----------------------------//
//---------------------------------------------------------------------------//
void DistributionFunctions::smooth(const std::string &name, double sigma,
                                   KernelType kernel) {
  if (histograms_.find(name) == histograms_.end()) {
    throw std::runtime_error("Histogram '" + name +
                             "' not found for smoothing.");
  }

  auto &hist = histograms_.at(name);
  hist.smoothed_partials.clear();

  for (const auto &[key, partial_values] : hist.partials) {
    const double dx = hist.bins[1] - hist.bins[0];
    const double min_sigma = std::max(dx, sigma);
    hist.smoothed_partials[key] =
        KernelSmoothing(hist.bins, partial_values, min_sigma, kernel);
  }
}

void DistributionFunctions::smoothAll(double sigma, KernelType kernel) {
  for (const auto &[name, histogram] : histograms_) {
    smooth(name, sigma, kernel);
  }
}

//---------------------------------------------------------------------------//
//----------------------------- Calculation CN ------------------------------//
//---------------------------------------------------------------------------//

void DistributionFunctions::calculateCoordinationNumber() {
  histograms_["CN"] = CNCalculator::calculate(cell_, neighbors());
}

//---------------------------------------------------------------------------//
//----------------------------- Calculation RDF -----------------------------//
//---------------------------------------------------------------------------//

void DistributionFunctions::calculateRDF(double r_max, double r_bin_width) {
  auto results = RDFCalculator::calculate(cell_, neighbors(), ashcroft_weights_,
                                          r_max, r_bin_width);
  for (auto &[name, histogram] : results) {
    histograms_[name] = std::move(histogram);
  }
}

//---------------------------------------------------------------------------//
//----------------------------- Calculation PAD -----------------------------//
//---------------------------------------------------------------------------//

void DistributionFunctions::calculatePAD(double bin_width) {
  histograms_["f(theta)"] =
      PADCalculator::calculate(cell_, neighbors(), bin_width);
}

//---------------------------------------------------------------------------//
//----------------------------- Calculation DAD -----------------------------//
//---------------------------------------------------------------------------//

void DistributionFunctions::calculateDAD(double bin_width) {
  histograms_["DAD"] = DADCalculator::calculate(cell_, neighbors(), bin_width);
}

//---------------------------------------------------------------------------//
//----------------------------- Calculation MD ------------------------------//
//---------------------------------------------------------------------------//

void DistributionFunctions::calculateMD(size_t max_ring_size) {
  histograms_["MD"] =
      MDCalculator::calculate(neighbors()->neighborGraph(), max_ring_size);
}

//---------------------------------------------------------------------------//
//----------------------------- Calculation VACF ----------------------------//
//---------------------------------------------------------------------------//

void DistributionFunctions::calculateVACF(const Trajectory &traj,
                                          int max_correlation_frames,
                                          size_t start_frame,
                                          size_t end_frame) {
  auto results = VACFCalculator::calculate(traj, max_correlation_frames,
                                           start_frame, end_frame);
  for (auto &[name, histogram] : results) {
    histograms_[name] = std::move(histogram);
  }
}

//---------------------------------------------------------------------------//
//----------------------------- Calculation VDOS ----------------------------//
//---------------------------------------------------------------------------//

void DistributionFunctions::calculateVDOS() {
  if (histograms_.find("VACF") == histograms_.end()) {
    throw std::logic_error(
        "Cannot calculate VDOS. Please calculate VACF first.");
  }
  histograms_["VDOS"] = VDOSCalculator::calculate(histograms_.at("VACF"));
}

//---------------------------------------------------------------------------//
//---------------------------- Calculation S(Q) -----------------------------//
//---------------------------------------------------------------------------//

void DistributionFunctions::calculateSQ(double q_max, double q_bin_width,
                                        double r_integration_max) {
  if (histograms_.find("g(r)") == histograms_.end()) {
    throw std::logic_error(
        "Cannot calculate S(Q). Please calculate g(r) first by calling "
        "calculateRDF().");
  }
  histograms_["S(Q)"] =
      SQCalculator::calculate(histograms_.at("g(r)"), cell_, ashcroft_weights_,
                              q_max, q_bin_width, r_integration_max);
}

void DistributionFunctions::calculateXRD(double lambda, double theta_min,
                                         double theta_max, double bin_width) {
  if (histograms_.find("g(r)") == histograms_.end()) {
    throw std::logic_error(
        "Cannot calculate XRD. Please calculate g(r) first by calling "
        "calculateRDF().");
  }
  histograms_["XRD"] =
      XRDCalculator::calculate(histograms_.at("g(r)"), cell_, ashcroft_weights_,
                               lambda, theta_min, theta_max, bin_width);
}

//---------------------------------------------------------------------------//
//----------------------------- Accumulation --------------------------------//
//---------------------------------------------------------------------------//

void DistributionFunctions::add(const DistributionFunctions &other) {
  // Iterate over all histograms in the other object
  for (const auto &[name, other_hist] : other.histograms_) {
    // Check if this object has the corresponding histogram
    auto it = histograms_.find(name);
    if (it == histograms_.end()) {
      // If not found, copy it entirely (e.g. if this is the first accumulation
      // but it shouldn't happen if we initialize with one frame)
      histograms_[name] = other_hist;
      continue;
    }

    auto &this_hist = it->second;

    // Consistency check (optional but recommended)
    if (this_hist.bins.size() != other_hist.bins.size()) {
      continue; // Cannot add incompatible histograms
    }

    // Accumulate Partials
    for (const auto &[key, other_partial] : other_hist.partials) {
      if (this_hist.partials.find(key) == this_hist.partials.end()) {
        this_hist.partials[key] = other_partial;
      } else {
        auto &this_partial = this_hist.partials[key];
        for (size_t i = 0; i < this_partial.size(); ++i) {
          this_partial[i] += other_partial[i];
        }
      }
    }

    // Clear smoothed partials as they are no longer valid after accumulation
    this_hist.smoothed_partials.clear();
  }
}

void DistributionFunctions::scale(double factor) {
  for (auto &[name, hist] : histograms_) {
    for (auto &[key, partial] : hist.partials) {
      for (size_t i = 0; i < partial.size(); ++i) {
        partial[i] *= factor;
      }
    }
    // Also scale smoothed partials if they exist (though they shouldn't if we
    // just accumulated)
    for (auto &[key, smoothed_partial] : hist.smoothed_partials) {
      for (size_t i = 0; i < smoothed_partial.size(); ++i) {
        smoothed_partial[i] *= factor;
      }
    }
  }
}
//---------------------------------------------------------------------------//
//----------------------------- Mean Calculation ----------------------------//
//---------------------------------------------------------------------------//

std::unique_ptr<DistributionFunctions> DistributionFunctions::computeMean(
    Trajectory &trajectory, const TrajectoryAnalyzer &analyzer,
    size_t start_frame, const AnalysisSettings &settings,
    std::function<void(float)> progress_callback) {

  const auto &analyzers = analyzer.getAnalyzers();
  if (analyzers.empty()) {
    return nullptr;
  }

  const size_t num_frames = analyzers.size();
  std::vector<std::unique_ptr<DistributionFunctions>> results(num_frames);
  std::vector<std::future<void>> futures;

  // Retrieve bond cutoffs from the analyzer.
  // Note: TrajectoryAnalyzer uses same cutoffs for all frames usually.
  auto bond_cutoffs = analyzer.getBondCutoffsSQ();

  // Mutex for progress callback (if needed, though atomic is better for
  // count)
  std::mutex callback_mutex;
  std::atomic<size_t> completed_frames{0};

  for (size_t i = 0; i < num_frames; ++i) {
    futures.push_back(std::async(std::launch::async, [&, i]() {
      size_t frame_idx = start_frame + i;
      if (frame_idx >= trajectory.getFrames().size()) {
        // Should not happen if TrajectoryAnalyzer was constructed
        // correctly
        return;
      }

      Cell &frame = trajectory.getFrames()[frame_idx];

      // Create DF for this frame
      auto frame_df =
          std::make_unique<DistributionFunctions>(frame, 0.0, bond_cutoffs);
      frame_df->setStructureAnalyzer(analyzers[i].get());

      frame_df->calculateCoordinationNumber();
      frame_df->calculateRDF(settings.r_max, settings.r_bin_width);
      if (settings.angle_bin_width > 0) {
        frame_df->calculatePAD(settings.angle_bin_width);
        frame_df->calculateDAD(settings.angle_bin_width);
      }
      frame_df->calculateMD(settings.max_ring_size);
      if (settings.q_max > 0) {
        frame_df->calculateSQ(settings.q_max, settings.q_bin_width,
                              settings.r_int_max);
      }

      results[i] = std::move(frame_df);

      size_t current_completed = ++completed_frames;
      if (progress_callback) {
        // Avoid flooding the callback, update maybe every 1% or just
        // simple throttle? Or lock.
        std::lock_guard<std::mutex> lock(callback_mutex);
        progress_callback(static_cast<float>(current_completed) /
                          static_cast<float>(num_frames));
      }
    }));
  }

  // Wait for all tasks
  for (auto &f : futures) {
    if (f.valid()) {
      f.get();
    }
  }

  // Accumulate results
  if (results.empty() || !results[0]) {
    return nullptr;
  }

  auto &final_df = results[0];
  for (size_t i = 1; i < num_frames; ++i) {
    if (results[i]) {
      final_df->add(*results[i]);
    }
  }

  if (num_frames > 1) {
    final_df->scale(1.0 / static_cast<double>(num_frames));
  }

  if (settings.smoothing) {
    final_df->smoothAll(settings.smoothing_sigma, settings.smoothing_kernel);
  }

  return std::move(final_df);
}
