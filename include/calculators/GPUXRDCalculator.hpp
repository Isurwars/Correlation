/**
 * @file GPUXRDCalculator.hpp
 * @brief GPU-accelerated X-Ray Diffraction (XRD) calculator with automatic CPU fallback.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "core/Cell.hpp"

namespace correlation::calculators {

/**
 * @struct GPUXRDParams
 * @brief Configuration parameters for GPU-accelerated XRD calculation.
 */
struct GPUXRDParams {
  real_t lambda{1.5406};  ///< X-ray source wavelength in Angstroms (default: Cu-Kalpha = 1.5406 A).
  real_t theta_min{10.0}; ///< Minimum 2theta diffraction angle in degrees.
  real_t theta_max{80.0}; ///< Maximum 2theta diffraction angle in degrees.
  real_t bin_width{0.1};  ///< Angular resolution delta(2theta) in degrees.
};

namespace gpu {

/**
 * @brief Computes XRD diffraction pattern using direct Debye scattering summation on GPU with CPU fallback.
 * @param[in] cell The atomic unit cell or cluster.
 * @param[in] params Diffraction parameters (wavelength, theta bounds, resolution).
 * @return Histogram containing 2theta bins and diffraction intensity.
 */
correlation::analysis::Histogram compute_xrd_gpu(const correlation::core::Cell &cell, const GPUXRDParams &params = {});

} // namespace gpu

/**
 * @class GPUXRDCalculator
 * @brief Computes XRD patterns on GPU (CUDA/HIP) using direct Debye scattering with CPU fallback.
 *
 * Build requirements:
 *   - Compile with `-DBUILD_WITH_CUDA=ON` or `-DBUILD_WITH_HIP=ON`.
 *
 * Runtime behaviour:
 *   - Probes for a compatible GPU device on construction and execution.
 *   - If no device is detected, falls back seamlessly to the CPU XRDCalculator.
 */
class GPUXRDCalculator : public BaseCalculator {
public:
  GPUXRDCalculator();

  [[nodiscard]] std::string_view getName() const override { return "XRD — GPU Accelerated"; }
  [[nodiscard]] std::string_view getShortName() const override { return "XRD_GPU"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Diffraction"; }
  [[nodiscard]] std::string_view getDescription() const override {
    return "GPU-accelerated X-Ray Diffraction (XRD) pattern calculation using direct Debye scattering.";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  [[nodiscard]] bool hasGPU() const noexcept { return has_gpu_; }

private:
  bool has_gpu_{false};
};

} // namespace correlation::calculators
