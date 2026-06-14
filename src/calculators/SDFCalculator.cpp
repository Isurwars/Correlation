/**
 * @file SDFCalculator.cpp
 * @brief Implementation of the Spatial Distribution Function calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/SDFCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"

#include <cmath>
#include <vector>

namespace correlation::calculators {

namespace {
bool registered = CalculatorFactory::instance().registerCalculator(std::make_unique<SDFCalculator>());
} // namespace

void SDFCalculator::calculateFrame(correlation::analysis::DistributionFunctions &df,
                                   const correlation::analysis::AnalysisSettings &settings) const {

  const auto &cell = df.cell();
  // SDF is a 3D grid and requires a coarser resolution than 1D radial distributions.
  // We enforce a minimum grid spacing of 0.5 Å to prevent memory explosion.
  const double dx = std::max(settings.r_bin_width > 0.0 ? settings.r_bin_width : 0.5, 0.5);

  // For a general implementation, we build a 3D grid based on the cell
  // dimensions
  double const lx = cell.lattice_parameters()[0];
  double const ly = cell.lattice_parameters()[1];
  double const lz = cell.lattice_parameters()[2];

  if (lx <= 0 || ly <= 0 || lz <= 0) {
    return;
}

  auto const nx = static_cast<size_t>(std::ceil(lx / dx));
  auto const ny = static_cast<size_t>(std::ceil(ly / dx));
  auto const nz = static_cast<size_t>(std::ceil(lz / dx));
  size_t const total_bins = nx * ny * nz;

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
  sdf_hist.description = "Flattened 3D Spatial Distribution Grid (" + std::to_string(nx) + "x" + std::to_string(ny) +
                         "x" + std::to_string(nz) + ")";
  sdf_hist.file_suffix = "_sdf";

  // Dummy bins to satisfy dimensions
  sdf_hist.bins.resize(total_bins, 0.0);
  // Optional: store nx, ny, nz in the first 3 bins just to preserve them across
  // frames Since we accumulate/scale, this is tricky. We'll rely on the
  // description string.

  // Voxel volume = total cell volume / number of voxels (correct for triclinic).
  const double dV = cell.volume() / static_cast<double>(total_bins);
  const auto &inv_lv = cell.inverseLatticeVectors();

  for (const auto &atom : cell.atoms()) {
    std::string const sym = atom.element().symbol;
    if (!sdf_hist.partials.contains(sym)) {
      sdf_hist.partials[sym].assign(total_bins, 0.0);
    }

    // Convert Cartesian → fractional coordinates (correct for non-orthogonal cells)
    auto frac = inv_lv * atom.position();
    // Wrap to [0, 1) fractional
    double const fx = frac.x() - std::floor(frac.x());
    double const fy = frac.y() - std::floor(frac.y());
    double const fz = frac.z() - std::floor(frac.z());

    size_t const ix = static_cast<size_t>(fx * nx) % nx; // NOLINT(bugprone-narrowing-conversions)
    size_t const iy = static_cast<size_t>(fy * ny) % ny; // NOLINT(bugprone-narrowing-conversions)
    size_t const iz = static_cast<size_t>(fz * nz) % nz; // NOLINT(bugprone-narrowing-conversions)

    size_t const idx = ix * (ny * nz) + iy * nz + iz;
    sdf_hist.partials[sym][idx] += 1.0 / dV; // Density contribution per frame

    // Also accumulate to total
    if (!sdf_hist.partials.contains("Total")) {
      sdf_hist.partials["Total"].assign(total_bins, 0.0);
    }
    sdf_hist.partials["Total"][idx] += 1.0 / dV;
  }

  if (!sdf_hist.partials.contains("Total")) {
    sdf_hist.partials["Total"].assign(total_bins, 0.0);
  }

  df.addHistogram("SDF", std::move(sdf_hist));
}

} // namespace correlation::calculators
