/**
 * @file SDFCalculator.cpp
 * @brief Implementation of the Spatial Distribution Function calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "calculators/SDFCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"

#include <cmath>
#include <vector>

namespace correlation::calculators {

namespace {
bool registered = CalculatorFactory::instance().registerCalculator(
    std::make_unique<SDFCalculator>());
} // namespace

void SDFCalculator::calculateFrame(
    correlation::analysis::DistributionFunctions &df,
    const correlation::analysis::AnalysisSettings &settings) const {

  const auto &cell = df.cell();
  const double dx = settings.r_bin_width > 0.0 ? settings.r_bin_width : 0.5;

  // For a general implementation, we build a 3D grid based on the cell
  // dimensions
  double lx = cell.lattice_parameters()[0];
  double ly = cell.lattice_parameters()[1];
  double lz = cell.lattice_parameters()[2];

  if (lx <= 0 || ly <= 0 || lz <= 0)
    return;

  size_t nx = static_cast<size_t>(std::ceil(lx / dx));
  size_t ny = static_cast<size_t>(std::ceil(ly / dx));
  size_t nz = static_cast<size_t>(std::ceil(lz / dx));
  size_t total_bins = nx * ny * nz;

  if (total_bins == 0 || total_bins > 100000000) {
    // Prevent memory explosion if grid is too fine
    return;
  }

  correlation::analysis::Histogram sdf_hist;
  sdf_hist.x_label = "Flattened 3D Index";
  sdf_hist.title = "Spatial Distribution Function (3D)";
  sdf_hist.y_label = "Density";
  sdf_hist.x_unit = "-";
  sdf_hist.y_unit = "atoms/A^3";
  sdf_hist.description = "Flattened 3D Spatial Distribution Grid (" +
                         std::to_string(nx) + "x" + std::to_string(ny) + "x" +
                         std::to_string(nz) + ")";
  sdf_hist.file_suffix = "_sdf";

  // Dummy bins to satisfy dimensions
  sdf_hist.bins.resize(total_bins, 0.0);
  // Optional: store nx, ny, nz in the first 3 bins just to preserve them across
  // frames Since we accumulate/scale, this is tricky. We'll rely on the
  // description string.

  double dV = (lx / nx) * (ly / ny) * (lz / nz);

  for (const auto &atom : cell.atoms()) {
    std::string sym = atom.element().symbol;
    if (sdf_hist.partials.find(sym) == sdf_hist.partials.end()) {
      sdf_hist.partials[sym].assign(total_bins, 0.0);
    }

    auto pos = atom.position();
    // Apply periodic boundary conditions to wrap into [0, L)
    double rx = pos.x() - lx * std::floor(pos.x() / lx);
    double ry = pos.y() - ly * std::floor(pos.y() / ly);
    double rz = pos.z() - lz * std::floor(pos.z() / lz);

    size_t ix = static_cast<size_t>(rx / lx * nx) % nx;
    size_t iy = static_cast<size_t>(ry / ly * ny) % ny;
    size_t iz = static_cast<size_t>(rz / lz * nz) % nz;

    size_t idx = ix * (ny * nz) + iy * nz + iz;
    sdf_hist.partials[sym][idx] += 1.0 / dV; // Density contribution per frame

    // Also accumulate to total
    if (sdf_hist.partials.find("Total") == sdf_hist.partials.end()) {
      sdf_hist.partials["Total"].assign(total_bins, 0.0);
    }
    sdf_hist.partials["Total"][idx] += 1.0 / dV;
  }

  df.addHistogram("SDF", std::move(sdf_hist));
}

} // namespace correlation::calculators
