/**
 * @file GPUSteinhardtCalculator.hpp
 * @brief GPU-accelerated Steinhardt bond-order parameter calculator (Q4, Q6, W6) with CPU fallback.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 *
 * When compiled with CORRELATION_USE_CUDA, this calculator performs spherical
 * harmonic evaluations and bond-order accumulations on the GPU. If no
 * compatible GPU is found at runtime, it falls back to the CPU SteinhardtCalculator.
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"

namespace correlation::calculators {

/**
 * @class GPUSteinhardtCalculator
 * @brief Computes Steinhardt parameters (Q4, Q6, W6) on GPU (CUDA/HIP) with automatic CPU fallback.
 */
class GPUSteinhardtCalculator : public BaseCalculator {
public:
  GPUSteinhardtCalculator();

  [[nodiscard]] std::string getName() const override { return "Steinhardt Parameter — GPU Accelerated"; }
  [[nodiscard]] std::string getShortName() const override { return "Steinhardt_GPU"; }
  [[nodiscard]] std::string getGroup() const override { return "Structural"; }
  [[nodiscard]] std::string getDescription() const override {
    return "GPU-accelerated Steinhardt bond-orientational parameters (Q4, Q6, W6). "
           "Falls back to CPU when no compatible GPU is detected.";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Checks if a functional CUDA/HIP device was detected during instantiation.
   * @return True if GPU backend is available.
   */
  [[nodiscard]] bool hasGPU() const noexcept { return has_gpu_; }

private:
  bool has_gpu_{false};
};

} // namespace correlation::calculators
