/**
 * @file SDFCalculator.cpp
 * @brief Implementation of the Spatial Distribution Function calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/SDFCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Precision.hpp"

#include <cmath>
#include <sys/stat.h>
#include <vector>

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
const bool registered = CalculatorFactory::registerTypeSafe<SDFCalculator>("SDFCalculator");
} // namespace

void SDFCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                   const correlation::analysis::AnalysisSettings &settings) const {

  const auto &cell = dists.cell();
  // SDF is a 3D grid and requires a coarser resolution than 1D radial distributions.
  // We enforce a minimum grid spacing of 0.5 Å to prevent memory explosion.
  const real_t d_x = static_cast<real_t>(std::max(settings.r_bin_width > 0.0 ? settings.r_bin_width : 0.5, 0.5));

  // For a general implementation, we build a 3D grid based on the cell
  // dimensions
  real_t const l_x = cell.lattice_parameters()[0];
  real_t const l_y = cell.lattice_parameters()[1];
  real_t const l_z = cell.lattice_parameters()[2];

  if (l_x <= 0 || l_y <= 0 || l_z <= 0) {
    return;
  }

  auto const n_x = static_cast<size_t>(std::ceil(l_x / d_x));
  auto const n_y = static_cast<size_t>(std::ceil(l_y / d_x));
  auto const n_z = static_cast<size_t>(std::ceil(l_z / d_x));
  size_t const total_bins = n_x * n_y * n_z;

  if (total_bins == 0 || total_bins > 1000000) {
    // Prevent memory explosion if grid is too fine
    return;
  }

  correlation::analysis::Histogram sdf_hist;
  sdf_hist.x_label = "Flattened 3D Index";
  sdf_hist.title = "Spatial Distribution Function (3D)";
  sdf_hist.y_label = "Density";
  sdf_hist.x_unit = "-";
  sdf_hist.y_unit = "atoms/A^3";
  sdf_hist.description = "Flattened 3D Spatial Distribution Grid (" + std::to_string(n_x) + "x" + std::to_string(n_y) +
                         "x" + std::to_string(n_z) + ")";
  sdf_hist.file_suffix = "_sdf";

  // Dummy bins to satisfy dimensions
  sdf_hist.bins.resize(total_bins, 0.0);
  // Optional: store nx, ny, nz in the first 3 bins just to preserve them across
  // frames Since we accumulate/scale, this is tricky. We'll rely on the
  // description string.

  // Voxel volume = total cell volume / number of voxels (correct for triclinic).
  const real_t d_V = cell.volume() / static_cast<real_t>(total_bins);
  const auto &inv_lv = cell.inverseLatticeVectors();

  for (const auto &atom : cell.atoms()) {
    std::string const sym = atom.element().symbol;
    if (!sdf_hist.partials.contains(sym)) {
      sdf_hist.partials[sym].assign(total_bins, 0.0);
    }

    // Convert Cartesian → fractional coordinates (correct for non-orthogonal cells)
    auto frac = inv_lv * atom.position();
    // Wrap to [0, 1) fractional
    real_t const f_x = frac.x() - std::floor(frac.x());
    real_t const f_y = frac.y() - std::floor(frac.y());
    real_t const f_z = frac.z() - std::floor(frac.z());

    size_t const i_x = static_cast<size_t>(f_x * static_cast<real_t>(n_x)) % n_x;
    size_t const i_y = static_cast<size_t>(f_y * static_cast<real_t>(n_y)) % n_y;
    size_t const i_z = static_cast<size_t>(f_z * static_cast<real_t>(n_z)) % n_z;

    size_t const idx = i_x * (n_y * n_z) + i_y * n_z + i_z;
    sdf_hist.partials[sym][idx] += static_cast<real_t>(1.0) / d_V; // Density contribution per frame

    // Also accumulate to total
    if (!sdf_hist.partials.contains("Total")) {
      sdf_hist.partials["Total"].assign(total_bins, 0.0);
    }
    sdf_hist.partials["Total"][idx] += static_cast<real_t>(1.0) / d_V;
  }

  if (!sdf_hist.partials.contains("Total")) {
    sdf_hist.partials["Total"].assign(total_bins, 0.0);
  }

  dists.addHistogram("SDF", std::move(sdf_hist));
}

} // namespace correlation::calculators
