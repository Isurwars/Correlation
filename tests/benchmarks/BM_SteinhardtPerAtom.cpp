/**
 * @file BM_SteinhardtPerAtom.cpp
 * @brief Microbenchmarks for Steinhardt per-atom order parameters and halo rasterization.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "CrystalTestHelper.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "calculators/order/SteinhardtCalculator.hpp"
#include "render3d/Camera3D.hpp"
#include "render3d/Scene3D.hpp"
#include "render3d/cpu/SoftwareFramebuffer.hpp"
#include "render3d/cpu/SoftwareRasterizer.hpp"

#include <benchmark/benchmark.h>

namespace correlation::benchmarks {

namespace {

constexpr uint32_t FRAMEBUFFER_WIDTH = 800;
constexpr uint32_t FRAMEBUFFER_HEIGHT = 600;

void benchSteinhardtPerAtom10k(benchmark::State &state) {
  // 14x14x14 FCC supercell = 4 * 14^3 = 10,976 atoms
  const auto cell = testing::crystals::createFCCCell(3.61, "Cu", 14, 14, 14);
  const real_t cutoff = 3.0;
  const analysis::BondCutoffMatrix cutoff_matrix = {
      {{.min_sq = 0.0, .max_sq = cutoff * cutoff}},
  };
  const analysis::StructureAnalyzer analyzer(cell, cutoff, cutoff_matrix, false);

  for ([[maybe_unused]] const auto iter : state) {
    const auto params = calculators::SteinhardtCalculator::calculatePerAtom(cell, &analyzer);
    benchmark::DoNotOptimize(params.Q6_bar.data());
    benchmark::ClobberMemory();
  }
  state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(cell.atoms().size()));
}

void benchRasterizerHaloGlow(benchmark::State &state) {
  render3d::SoftwareFramebuffer framebuffer(FRAMEBUFFER_WIDTH, FRAMEBUFFER_HEIGHT);
  const render3d::SoftwareRasterizer rasterizer;

  render3d::Camera3D camera;
  camera.setViewport(static_cast<float>(FRAMEBUFFER_WIDTH),
                     static_cast<float>(FRAMEBUFFER_HEIGHT));

  // 6x6x6 FCC supercell = 864 atoms
  const auto cell = testing::crystals::createFCCCell(3.61, "Cu", 6, 6, 6);
  const auto scene = render3d::Scene3D::fromCell(cell, 1.15, false, false);
  if (scene == nullptr) {
    state.SkipWithError("Failed to generate scene");
    return;
  }
  camera.fitToBox(scene->boundsMin(), scene->boundsMax());
  camera.orbit(45.0F, 30.0F);

  // Apply halo properties to all atoms
  for (auto &atom : scene->atoms()) {
    atom.halo_color = render3d::RgbaColor{.r = 68, .g = 1, .b = 84, .a = 200};
    atom.halo_radius_ratio = 1.35F;
    atom.dim_factor = 1.0F;
  }

  for ([[maybe_unused]] const auto iter : state) {
    framebuffer.clear(render3d::RgbaColor{.r = 20, .g = 20, .b = 25, .a = 255}, 1.0F);
    rasterizer.render(*scene, camera, framebuffer);
    benchmark::DoNotOptimize(framebuffer.colorBuffer().data());
    benchmark::ClobberMemory();
  }
  state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(scene->atoms().size()));
}

} // namespace

BENCHMARK(benchSteinhardtPerAtom10k)->Unit(benchmark::kMillisecond);
BENCHMARK(benchRasterizerHaloGlow)->Unit(benchmark::kMillisecond);

} // namespace correlation::benchmarks
