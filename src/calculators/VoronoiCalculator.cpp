/**
 * @file VoronoiCalculator.cpp
 * @brief Implementation of the Voronoi Tessellation calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/VoronoiCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Constants.hpp"
#include "c_loops.hh"
#include "cell.hh"
#include "container_prd.hh"

#include <algorithm>
#include <cmath>
#include <format>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
// NOLINTNEXTLINE(cert-err58-cpp)
const bool registered = CalculatorFactory::registerTypeSafe<VoronoiCalculator>("VoronoiCalculator");
} // namespace

void VoronoiCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                       const correlation::analysis::AnalysisSettings & /*settings*/) const {
  auto hists = calculate(dists.cell(), dists.neighbors());
  for (auto &[name, hist] : hists) {
    dists.addHistogram(name, std::move(hist));
  }
}

// ---------------------------------------------------------------------------
// computeVoronoiCells — validation, coordinate mapping, voro++ computation
// ---------------------------------------------------------------------------
VoronoiCalculator::CellData VoronoiCalculator::computeVoronoiCells(const correlation::core::Cell &cell) {
  const auto &atoms = cell.atoms();
  size_t const num_atoms = atoms.size();
  if (num_atoms == 0) {
    throw std::invalid_argument("Cannot calculate Voronoi Tessellation for an empty cell.");
  }

  const auto &lattice = cell.latticeVectors();
  real_t const volume = cell.volume();
  if (volume <= 1e-9) {
    throw std::logic_error("Cell volume must be positive and finite.");
  }

  // Get aligned box parameters from lattice vectors.
  // In the aligned frame of reference:
  // a = (bx, 0, 0)
  // b = (bxy, by, 0)
  // c = (bxz, byz, bz)
  real_t const bx_ = lattice[0].x();
  real_t const bxy = lattice[1].x();
  real_t const by_ = lattice[1].y();
  real_t const bxz = lattice[2].x();
  real_t const byz = lattice[2].y();
  real_t const bz_ = lattice[2].z();

  if (bx_ <= 1e-9 || by_ <= 1e-9 || bz_ <= 1e-9) {
    throw std::runtime_error("Invalid or non-orthogonal/skewed cell dimensions for Voronoi calculation.");
  }

  // Map atom positions to aligned Cartesian coordinates via fractional coordinates
  // Perform the wrapping and aligned position calculations in real_t precision
  math::Matrix3<real_t> const lattice_d(lattice);
  math::Matrix3<real_t> const inv_lattice_d = math::invert(lattice_d);
  std::vector<std::array<real_t, 3>> aligned_positions(num_atoms);
  for (size_t i = 0; i < num_atoms; ++i) {
    math::Vector3<real_t> const pos_d(atoms[i].position());
    auto frac = inv_lattice_d * pos_d;
    // Wrap to [0, 1) fundamental domain
    frac.x() -= std::floor(frac.x());
    frac.y() -= std::floor(frac.y());
    frac.z() -= std::floor(frac.z());

    real_t const ax_ = frac.x() * bx_ + frac.y() * bxy + frac.z() * bxz;
    real_t const ay_ = frac.y() * by_ + frac.z() * byz;
    real_t const az_ = frac.z() * bz_;
    aligned_positions[i] = {ax_, ay_, az_};
  }

  // Dynamic grid estimator (aim for ~6 particles per grid block)
  real_t const optimal_block_vol = 6.0 / (static_cast<real_t>(num_atoms) / volume);
  real_t const block_side = std::max(static_cast<real_t>(1.0), std::cbrt(optimal_block_vol));
  int const nx_ = std::max(1, static_cast<int>(std::round(bx_ / block_side)));
  int const ny_ = std::max(1, static_cast<int>(std::round(by_ / block_side)));
  int const nz_ = std::max(1, static_cast<int>(std::round(bz_ / block_side)));

  // Setup periodic container
  voro::container_periodic con(bx_, bxy, by_, bxz, byz, bz_, nx_, ny_, nz_, 8);

  // Put atoms into container
  for (size_t i = 0; i < num_atoms; ++i) {
    con.put(static_cast<int>(i), aligned_positions[i][0], aligned_positions[i][1], aligned_positions[i][2]);
  }

  // Compute Voronoi cells
  CellData data;
  data.volumes.assign(num_atoms, 0.0);
  data.sphericities.assign(num_atoms, 0.0);
  data.coordination_numbers.assign(num_atoms, 0);
  data.signatures.assign(num_atoms, "");

  voro::voronoicell voro_cell;
  voro::c_loop_all_periodic voro_loop(con);
  if (!voro_loop.start()) {
    return data;
  }

  for (bool first = true; first || voro_loop.inc(); first = false) {
    if (!con.compute_cell(voro_cell, voro_loop)) {
      continue;
    }

    int const pid = voro_loop.pid();
    if (pid < 0 || pid >= static_cast<int>(num_atoms)) {
      continue;
    }

    real_t const vol = voro_cell.volume();
    real_t const area = voro_cell.surface_area();
    real_t const sphericity =
        (area > 1e-9) ? (std::pow(correlation::math::pi, 1.0 / 3.0) * std::pow(6.0 * vol, 2.0 / 3.0)) / area : 0.0;

    std::vector<int> orders;
    voro_cell.face_orders(orders);
    std::vector<double> areas;
    voro_cell.face_areas(areas);

    int n_3 = 0;
    int n_4 = 0;
    int n_5 = 0;
    int n_6 = 0;
    int significant_faces = 0;

    for (size_t i = 0; i < orders.size(); ++i) {
      if (i < areas.size() && areas[i] < 1e-4) {
        continue;
      }
      significant_faces++;
      int const order = orders[i];
      if (order == 3) {
        ++n_3;
      } else if (order == 4) {
        ++n_4;
      } else if (order == 5) {
        ++n_5;
      } else if (order == 6) {
        ++n_6;
      }
    }

    data.volumes[pid] = vol;
    data.sphericities[pid] = sphericity;
    data.coordination_numbers[pid] = significant_faces;
    data.signatures[pid] = std::format("({}, {}, {}, {})", n_3, n_4, n_5, n_6);
  }

  return data;
}

// ---------------------------------------------------------------------------
// buildSignatureMap — count, sort, and describe polyhedral signatures
// ---------------------------------------------------------------------------
std::pair<std::vector<std::pair<std::string, int>>, std::string>
VoronoiCalculator::buildSignatureMap(const std::vector<std::string> &signatures) {
  std::map<std::string, int> counts;
  for (const auto &sig : signatures) {
    if (!sig.empty()) {
      counts[sig]++;
    }
  }

  std::vector<std::pair<std::string, int>> sorted(counts.begin(), counts.end());
  std::ranges::sort(sorted, [](const auto &frst, const auto &secnd) { return frst.second > secnd.second; });

  std::string desc = "Polyhedral signature frequencies:\n";
  for (size_t i = 0; i < sorted.size(); ++i) {
    desc += std::format("  {}: {} (count: {})\n", i, sorted[i].first, sorted[i].second);
  }

  return {sorted, desc};
}

// ---------------------------------------------------------------------------
// makeHistogram — factory for initialised Histogram objects
// ---------------------------------------------------------------------------
correlation::analysis::Histogram
VoronoiCalculator::makeHistogram(const std::string &title, const std::string &x_label, const std::string &y_label,
                                 const std::string &x_unit, const std::string &y_unit, const std::string &description,
                                 const std::string &file_suffix, const std::vector<real_t> &bins,
                                 const std::vector<std::string> &element_symbols) {
  correlation::analysis::Histogram hist;
  hist.title = title;
  hist.x_label = x_label;
  hist.y_label = y_label;
  hist.x_unit = x_unit;
  hist.y_unit = y_unit;
  hist.description = description;
  hist.file_suffix = file_suffix;
  hist.bins = bins;

  size_t const num_bins = bins.size();
  for (const auto &sym : element_symbols) {
    hist.partials[sym].assign(num_bins, 0.0);
  }
  hist.partials["Total"].assign(num_bins, 0.0);

  return hist;
}

// ---------------------------------------------------------------------------
// populateHistogram — bin numeric per-atom values into a histogram
// ---------------------------------------------------------------------------
void VoronoiCalculator::populateHistogram(correlation::analysis::Histogram &hist, const BinRange &range,
                                          const std::vector<real_t> &values,
                                          const std::vector<correlation::core::Atom> &atoms) {
  for (size_t i = 0; i < values.size(); ++i) {
    real_t const val = values[i];
    if (val < range.range_min || val >= range.range_max) {
      continue;
    }
    auto const idx = static_cast<size_t>(val / range.bin_width);
    if (idx >= range.num_bins) {
      continue;
    }
    const std::string &sym = atoms[i].element().symbol;
    hist.partials[sym][idx] += 1.0;
    hist.partials["Total"][idx] += 1.0;
  }
}

// ---------------------------------------------------------------------------
// calculate — orchestrator
// ---------------------------------------------------------------------------
std::map<std::string, correlation::analysis::Histogram>
VoronoiCalculator::calculate(const correlation::core::Cell &cell,
                             const correlation::analysis::StructureAnalyzer * /*neighbors*/) {
  // 1. Compute Voronoi tessellation
  CellData data = computeVoronoiCells(cell);
  size_t const num_atoms = cell.atoms().size();

  // Unique chemical elements
  std::vector<std::string> element_symbols;
  for (const auto &elem : cell.elements()) {
    element_symbols.push_back(elem.symbol);
  }

  // 2. Build signature map
  auto [sorted_sigs, sig_desc] = buildSignatureMap(data.signatures);
  std::map<std::string, size_t> signature_indices;
  for (size_t i = 0; i < sorted_sigs.size(); ++i) {
    signature_indices[sorted_sigs[i].first] = i;
  }

  // 3. Initialize histograms
  std::map<std::string, correlation::analysis::Histogram> results;

  // 3a. Volume histogram
  real_t max_vol = *std::max_element(data.volumes.begin(), data.volumes.end());
  if (max_vol <= 0.0) {
    max_vol = 100.0;
  }
  real_t const max_vol_range = max_vol * 1.2;
  size_t const vol_bins = 100;
  real_t const vol_d = max_vol_range / static_cast<real_t>(vol_bins);

  std::vector<real_t> vol_bin_centers(vol_bins);
  for (size_t i = 0; i < vol_bins; ++i) {
    vol_bin_centers[i] = (static_cast<real_t>(i) + 0.5) * vol_d;
  }
  results["Voronoi Volume"] =
      makeHistogram("Voronoi Cell Volume Distribution", "Volume", "Probability Density", "Å³", "probability",
                    "Probability distribution of Voronoi cell volumes.", "_vvol", vol_bin_centers, element_symbols);

  // 3b. Sphericity histogram
  size_t const sph_bins = 100;
  real_t const sph_d = 1.0 / static_cast<real_t>(sph_bins);

  std::vector<real_t> sph_bin_centers(sph_bins);
  for (size_t i = 0; i < sph_bins; ++i) {
    sph_bin_centers[i] = (static_cast<real_t>(i) + 0.5) * sph_d;
  }
  results["Voronoi Sphericity"] = makeHistogram(
      "Voronoi Cell Sphericity Distribution", "Sphericity", "Probability Density", "dimensionless", "probability",
      "Probability distribution of Voronoi cell sphericities (Psi).", "_vsph", sph_bin_centers, element_symbols);

  // 3c. Coordination number histogram
  int max_cn = *std::max_element(data.coordination_numbers.begin(), data.coordination_numbers.end());
  size_t const cn_bins = std::max(25, max_cn + 2);

  std::vector<real_t> cn_bin_values(cn_bins);
  std::iota(cn_bin_values.begin(), cn_bin_values.end(), 0.0);
  results["Voronoi Coordination Number"] = makeHistogram(
      "Voronoi Coordination Number Distribution", "Coordination Number", "Probability", "faces", "probability",
      "Probability distribution of Voronoi cell face counts.", "_vcn", cn_bin_values, element_symbols);

  // 3d. Signatures histogram
  size_t const sig_bins = std::max(size_t{1}, sorted_sigs.size());

  std::vector<real_t> sig_bin_values(sig_bins);
  std::iota(sig_bin_values.begin(), sig_bin_values.end(), 0.0);
  results["Voronoi Signatures"] =
      makeHistogram("Voronoi Polyhedral Signatures", "Signature Index", "Probability", "index", "probability", sig_desc,
                    "_vsig", sig_bin_values, element_symbols);

  // 4. Populate histograms
  populateHistogram(results["Voronoi Volume"],
                    {.bin_width = vol_d, .num_bins = vol_bins, .range_min = 0.0, .range_max = max_vol_range},
                    data.volumes, cell.atoms());
  populateHistogram(results["Voronoi Sphericity"],
                    {.bin_width = sph_d, .num_bins = sph_bins, .range_min = 0.0, .range_max = 1.0}, data.sphericities,
                    cell.atoms());

  // Coordination numbers — convert to real_t for the shared helper
  std::vector<real_t> cn_as_double(data.coordination_numbers.begin(), data.coordination_numbers.end());
  populateHistogram(
      results["Voronoi Coordination Number"],
      {.bin_width = 1.0, .num_bins = cn_bins, .range_min = 0.0, .range_max = static_cast<real_t>(cn_bins)},
      cn_as_double, cell.atoms());

  // Signatures — unique lookup logic, kept inline
  auto &sig_hist = results["Voronoi Signatures"];
  for (size_t i = 0; i < num_atoms; ++i) {
    const std::string &sig = data.signatures[i];
    if (sig.empty() || !signature_indices.contains(sig)) {
      continue;
    }
    size_t const idx = signature_indices[sig];
    if (idx >= sig_bins) {
      continue;
    }
    const std::string &sym = cell.atoms()[i].element().symbol;
    sig_hist.partials[sym][idx] += 1.0;
    sig_hist.partials["Total"][idx] += 1.0;
  }

  // 5. Normalize histograms
  real_t const factor = 1.0 / static_cast<real_t>(num_atoms);
  for (auto &[name, hist] : results) {
    for (auto &[key, vec] : hist.partials) {
      for (auto &val : vec) {
        val *= factor;
      }
    }
  }

  return results;
}

} // namespace correlation::calculators
