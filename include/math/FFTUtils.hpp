/**
 * @file FFTUtils.hpp
 * @brief Fast Fourier Transform (FFT) and autocorrelation utilities.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include <complex>
#include <stdexcept>
#include <vector>

#if defined(CORRELATION_USE_FFTW3)
#include <fftw3.h>
#include <mutex>
#include <unordered_map>
#elif defined(CORRELATION_USE_MKL)
#include <mkl_dfti.h>
#include <unordered_map>
#endif

namespace correlation::math {

/**
 * @brief Helper to check if a size is a product of small prime factors (2, 3, 5).
 * Highly optimized for FFTW3 and MKL.
 */
inline bool isGoodFFTSize(size_t n) {
  if (n == 0) {
    return false;
  }
  while (n % 2 == 0) {
    n /= 2;
  }
  while (n % 3 == 0) {
    n /= 3;
  }
  while (n % 5 == 0) {
    n /= 5;
  }
  return n == 1;
}

/**
 * @brief Finds the next "good" FFT size >= target.
 */
inline size_t findNextGoodFFTSize(size_t target) {
#if defined(CORRELATION_USE_FFTW3) || defined(CORRELATION_USE_MKL)
  while (!isGoodFFTSize(target)) {
    target++;
  }
  return target;
#else
  // For hand-coded Cooley-Tukey, we must pad to the next power of two
  size_t len = 1;
  while (len < target) {
    len <<= 1;
  }
  return len;
#endif
}

#if defined(CORRELATION_USE_FFTW3)

struct FFTWPlanCache {
  std::unordered_map<size_t, fftw_plan> forward_plans;
  std::unordered_map<size_t, fftw_plan> backward_plans;

  // Rule of Five: Delete copy and move operations to prevent double-destruction of fftw_plans
  FFTWPlanCache() = default;
  FFTWPlanCache(const FFTWPlanCache &) = delete;
  FFTWPlanCache &operator=(const FFTWPlanCache &) = delete;
  FFTWPlanCache(FFTWPlanCache &&) = delete;
  FFTWPlanCache &operator=(FFTWPlanCache &&) = delete;

  ~FFTWPlanCache() {
    for (auto &[key, plan] : forward_plans) {
      fftw_destroy_plan(plan);
    }
    for (auto &[key, plan] : backward_plans) {
      fftw_destroy_plan(plan);
    }
  }
};

inline void computeFFT(std::vector<std::complex<double>> &data, bool invert) {
  size_t size = data.size();
  if (size == 0) {
    return;
  }

  thread_local FFTWPlanCache cache;
  auto &plans = invert ? cache.backward_plans : cache.forward_plans;
  auto iterator = plans.find(size);
  fftw_plan plan = nullptr;

  if (iterator != plans.end()) {
    plan = iterator->second;
  } else {
    static std::mutex planner_mutex;
    std::lock_guard<std::mutex> lock(planner_mutex);

    plan = fftw_plan_dft_1d(
        static_cast<int>(size),
        reinterpret_cast<fftw_complex *>(data.data()), // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
        reinterpret_cast<fftw_complex *>(data.data()), // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
        invert ? FFTW_BACKWARD : FFTW_FORWARD, FFTW_ESTIMATE);
    if (plan == nullptr) {
      throw std::runtime_error("Failed to create FFTW plan for size " + std::to_string(size));
    }
    plans[size] = plan;
  }

  fftw_execute_dft(
      plan,
      reinterpret_cast<fftw_complex *>(data.data()),  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
      reinterpret_cast<fftw_complex *>(data.data())); // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)

  if (invert) {
    double scale = 1.0 / static_cast<double>(size);
    for (auto &val : data) {
      val *= scale;
    }
  }
}

#elif defined(CORRELATION_USE_MKL)

struct MKLDescriptorCache {
  std::unordered_map<size_t, DFTI_DESCRIPTOR_HANDLE> forward_handles;
  std::unordered_map<size_t, DFTI_DESCRIPTOR_HANDLE> backward_handles;

  ~MKLDescriptorCache() {
    for (auto &p : forward_handles) {
      DftiFreeDescriptor(&p.second);
    }
    for (auto &p : backward_handles) {
      DftiFreeDescriptor(&p.second);
    }
  }
};

inline void computeFFT(std::vector<std::complex<double>> &data, bool invert) {
  size_t size = data.size();
  if (size == 0) {
    return;
  }

  thread_local MKLDescriptorCache cache;
  auto &handles = invert ? cache.backward_handles : cache.forward_handles;
  auto it = handles.find(size);
  DFTI_DESCRIPTOR_HANDLE handle = nullptr;

  if (it != handles.end()) {
    handle = it->second;
  } else {
    MKL_LONG status = DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, static_cast<MKL_LONG>(size));
    if (status != DFTI_NO_ERROR) {
      throw std::runtime_error("MKL DftiCreateDescriptor failed");
    }
    if (invert) {
      double scale = 1.0 / static_cast<double>(size);
      DftiSetValue(handle, DFTI_BACKWARD_SCALE, scale);
    }
    status = DftiCommitDescriptor(handle);
    if (status != DFTI_NO_ERROR) {
      DftiFreeDescriptor(&handle);
      throw std::runtime_error("MKL DftiCommitDescriptor failed");
    }
    handles[size] = handle;
  }

  MKL_LONG status;
  if (invert) {
    status = DftiComputeBackward(handle, data.data());
  } else {
    status = DftiComputeForward(handle, data.data());
  }

  if (status != DFTI_NO_ERROR) {
    throw std::runtime_error("MKL DFTI computation failed");
  }
}

#else

/**
 * @brief Simple Radix-2 Cooley-Tukey FFT implementation with bit-reversal
 * optimization.
 *
 * @param data The input/output vector of complex numbers. Modifies in-place.
 * @param invert If true, performs an inverse FFT and scales the result by 1/N.
 */
inline void computeFFT(std::vector<std::complex<double>> &data, bool invert) {
  size_t size = data.size();
  if (size == 0) {
    return;
  }

  if ((size & (size - 1)) != 0) {
    throw std::invalid_argument("FFT length must be a power of 2");
  }

  // Bit-reversal permutation (Optimized)
  for (size_t i = 1, j = 0; i < size; i++) {
    size_t bit = size >> 1;
    for (; (j & bit) != 0; bit >>= 1) {
      j ^= bit;
    }
    j ^= bit;
    if (i < j) {
      std::swap(data[i], data[j]);
    }
  }

  // Cooley-Tukey with precomputed twiddle basics
  for (size_t len = 2; len <= size; len <<= 1) {
    double angle = correlation::math::two_pi / static_cast<double>(len) * (invert ? -1 : 1);
    std::complex<double> wlen(std::cos(angle), std::sin(angle));
    for (size_t i = 0; i < size; i += len) {
      std::complex<double> twiddle(1.0, 0.0);
      for (size_t j = 0; j < len / 2; j++) {
        std::complex<double> even_val = data[i + j];
        std::complex<double> odd_val = data[i + j + len / 2] * twiddle;
        data[i + j] = even_val + odd_val;
        data[i + j + len / 2] = even_val - odd_val;
        twiddle *= wlen;
      }
    }
  }

  if (invert) {
    for (auto &val : data) {
      val /= static_cast<double>(size);
    }
  }
}

#endif

/**
 * @brief Computes the autocorrelation of a 1D real sequence using FFT.
 *
 * @param signal    Input signal.
 * @param workspace Reusable scratch buffer. Resized automatically when needed.
 *                  Pass the same buffer across repeated calls (e.g. per-atom
 *                  loops) to avoid repeated heap allocation.
 * @return Autocorrelation array of length n.
 */
inline std::vector<double> autocorrelate(const std::vector<double> &signal,
                                         std::vector<std::complex<double>> &workspace) {
  const size_t size = signal.size();
  if (size == 0) {
    return {};
  }

  size_t len = findNextGoodFFTSize(2 * size - 1);

  // Reuse caller-supplied workspace; only reallocates when length grows.
  workspace.assign(len, {0.0, 0.0});
  for (size_t i = 0; i < size; ++i) {
    workspace[i] = {signal[i], 0.0};
  }

  computeFFT(workspace, false);

  // In-place |X[k]|² — no intermediate vector needed.
  for (size_t i = 0; i < len; ++i) {
    const double real_part = workspace[i].real();
    const double imag_part = workspace[i].imag();
    workspace[i] = {real_part * real_part + imag_part * imag_part, 0.0};
  }

  computeFFT(workspace, true);

  std::vector<double> result(size);
  for (size_t i = 0; i < size; ++i) {
    result[i] = workspace[i].real();
  }
  return result;
}

/**
 * @brief Convenience overload — allocates its own workspace.
 *
 * Prefer the two-argument overload in tight loops to recycle the allocation.
 *
 * @param signal Input signal vector.
 * @return The autocorrelation vector.
 */
inline std::vector<double> autocorrelate(const std::vector<double> &signal) {
  std::vector<std::complex<double>> workspace;
  return autocorrelate(signal, workspace);
}

} // namespace correlation::math
