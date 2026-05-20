/**
 * @file GPUSQCalculator.hpp
 * @brief GPU-accelerated structure factor S(Q) calculator with CPU fallback.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 *
 * When compiled with CORRELATION_USE_CUDA, this calculator performs the
 * plane-wave summation on an Nvidia GPU. If no compatible device is found
 * at runtime, it transparently falls back to the CPU-based
 * StructureFactorCalculator.
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"

namespace correlation::calculators {

/**
 * @class GPUSQCalculator
 * @brief Computes S(Q) on the GPU (CUDA) with automatic CPU fallback.
 *
 * Build requirements:
 *   - Compile with `-DBUILD_WITH_CUDA=ON` and a CUDA toolkit.
 *   - The .cu translation unit (GPUSQCalculator.cu) contains the kernel.
 *
 * Runtime behaviour:
 *   - On construction the class probes `cudaGetDeviceCount()`.
 *   - If no device is found, `calculateFrame()` forwards to the CPU
 *     StructureFactorCalculator.
 */
class GPUSQCalculator : public BaseCalculator {
public:
  GPUSQCalculator();

  std::string getName() const override { return "S(Q) — GPU Accelerated"; }
  std::string getShortName() const override { return "S_Q_GPU"; }
  std::string getGroup() const override { return "Scattering"; }
  std::string getDescription() const override {
    return "GPU-accelerated static structure factor S(Q). Falls back to CPU "
           "when no compatible Nvidia GPU is detected.";
  }

  bool isFrameCalculator() const override { return true; }
  bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(
      correlation::analysis::DistributionFunctions &df,
      const correlation::analysis::AnalysisSettings &settings) const override;

  /// Returns true when a usable CUDA device was found at construction.
  [[nodiscard]] bool hasGPU() const noexcept { return has_gpu_; }

private:
  bool has_gpu_{false};
};

} // namespace correlation::calculators
