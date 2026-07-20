/**
 * @file wasm_bindings.cpp
 * @brief Emscripten bindings for running Correlation in the browser.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 *
 * Compile with the Emscripten toolchain:
 *   emcmake cmake .. -DBUILD_WASM=ON -DBUILD_TESTING=OFF
 *   emmake make correlation_wasm
 *
 * Exposes the following to JavaScript via embind:
 *   - Cell, Trajectory, DistributionFunctions (construction + calculation)
 *   - read()  — parse a file buffer into a Trajectory
 *   - Histogram data extraction for plotting
 */

#ifdef __EMSCRIPTEN__

#include "analysis/DistributionFunctions.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "readers/FileReader.hpp"

#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <fstream>
#include <string>
#include <vector>

using namespace emscripten;
using namespace correlation::core;
using namespace correlation::analysis;
using namespace correlation::readers;

// ---------------------------------------------------------------------------
// Helper: write a string buffer to the Emscripten virtual filesystem, read it
// back via the normal reader pipeline.
// ---------------------------------------------------------------------------
static Trajectory readFromBuffer(const std::string &data, const std::string &filename) {
  // Write to virtual FS.
  {
    std::ofstream f("/" + filename, std::ios::binary);
    f.write(data.c_str(), data.size());
  }
  FileType ft = determineFileType(filename);
  return readTrajectory("/" + filename, ft);
}

// ---------------------------------------------------------------------------
// Helper: extract histogram bins as a JS Float64Array.
// ---------------------------------------------------------------------------
static val getBinsJS(const Histogram &h) { return val(typed_memory_view(h.bins.size(), h.bins.data())); }

// ---------------------------------------------------------------------------
// Helper: extract a partial as a JS Float64Array.
// ---------------------------------------------------------------------------
static val getPartialJS(const Histogram &h, const std::string &key) {
  auto it = h.partials.find(key);
  if (it == h.partials.end())
    return val::null();
  return val(typed_memory_view(it->second.size(), it->second.data()));
}

// ---------------------------------------------------------------------------
// Helper: list partial keys.
// ---------------------------------------------------------------------------
static val getPartialKeysJS(const Histogram &h) {
  val keys = val::array();
  unsigned idx = 0;
  for (const auto &[k, _] : h.partials) {
    keys.set(idx++, k);
  }
  return keys;
}

// ============================================================================
// Emscripten bindings
// ============================================================================
EMSCRIPTEN_BINDINGS(correlation_wasm) {

  // ---- Cell ----
  class_<Cell>("Cell").constructor<>().function("atomCount", &Cell::atomCount).function("getVolume", &Cell::volume);

  // ---- Trajectory ----
  class_<Trajectory>("Trajectory")
      .constructor<>()
      .function("numFrames", &Trajectory::getFrameCount)
      .property("timeStep", &Trajectory::getTimeStep, &Trajectory::setTimeStep);

  // ---- Histogram ----
  class_<Histogram>("Histogram")
      .property("title", &Histogram::title)
      .property("xLabel", &Histogram::x_label)
      .property("yLabel", &Histogram::y_label)
      .function("getBins", &getBinsJS)
      .function("getPartial", &getPartialJS)
      .function("getPartialKeys", &getPartialKeysJS);

  // ---- AnalysisSettings ----
  class_<AnalysisSettings>("AnalysisSettings")
      .constructor<>()
      .property("rMax", &AnalysisSettings::r_max)
      .property("rBinWidth", &AnalysisSettings::r_bin_width)
      .property("qMax", &AnalysisSettings::q_max)
      .property("qBinWidth", &AnalysisSettings::q_bin_width)
      .property("angleBinWidth", &AnalysisSettings::angle_bin_width)
      .property("lefCutoff", &AnalysisSettings::lef_cutoff)
      .property("lefSigma", &AnalysisSettings::lef_sigma);

  // ---- DistributionFunctions ----
  class_<DistributionFunctions>("DistributionFunctions")
      .constructor<Cell &, real_t, const std::vector<std::vector<real_t>> &>()
      .function("calculateRDF", &DistributionFunctions::calculateRDF)
      .function("calculatePAD", &DistributionFunctions::calculatePAD)
      .function("getHistogram",
                select_overload<const Histogram &(const std::string &) const>(&DistributionFunctions::getHistogram))
      .function("getAvailableHistograms", &DistributionFunctions::getAvailableHistograms);

  // ---- Free functions ----
  function("readFromBuffer", &readFromBuffer);
}

#endif // __EMSCRIPTEN__
