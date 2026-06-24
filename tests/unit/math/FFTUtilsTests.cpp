// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "math/FFTUtils.hpp"

#include <complex>
#include <gtest/gtest.h>
#include <vector>

namespace correlation::testing {

using namespace correlation::math;

namespace {

class FFTUtilsTests : public ::testing::Test {};

TEST_F(FFTUtilsTests, ComputeFFTHandlesPowerOfTwoAndInvert) {
  // Test forward and inverse FFT on size 8
  std::vector<std::complex<double>> signal = {{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0},
                                              {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
  auto original = signal;

  // Forward FFT
  computeFFT(signal, false);

  // Inverse FFT
  computeFFT(signal, true);

  for (size_t i = 0; i < signal.size(); ++i) {
    EXPECT_NEAR(signal[i].real(), original[i].real(), 1e-9);
    EXPECT_NEAR(signal[i].imag(), original[i].imag(), 1e-9);
  }
}

#if defined(CORRELATION_USE_FFTW3) || defined(CORRELATION_USE_MKL)
TEST_F(FFTUtilsTests, ComputeFFTHandlesNonPowerOfTwo) {
  std::vector<std::complex<double>> signal = {{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}, {5.0, 0.0}, {6.0, 0.0}};
  auto original = signal;
  computeFFT(signal, false);
  computeFFT(signal, true);
  for (size_t i = 0; i < signal.size(); ++i) {
    EXPECT_NEAR(signal[i].real(), original[i].real(), 1e-9);
    EXPECT_NEAR(signal[i].imag(), original[i].imag(), 1e-9);
  }
}
#else
TEST_F(FFTUtilsTests, ComputeFFTThrowsOnNonPowerOfTwo) {
  std::vector<std::complex<double>> bad_signal(3, {1.0, 0.0});
  EXPECT_THROW(computeFFT(bad_signal, false), std::invalid_argument);

  std::vector<std::complex<double>> bad_signal_zero_power(6, {1.0, 0.0});
  EXPECT_THROW(computeFFT(bad_signal_zero_power, false), std::invalid_argument);
}
#endif

TEST_F(FFTUtilsTests, ComputeFFTHandlesEmptyInput) {
  std::vector<std::complex<double>> empty_signal;
  EXPECT_NO_THROW(computeFFT(empty_signal, false));
}

TEST_F(FFTUtilsTests, AutocorrelateEmptyReturnsEmpty) {
  std::vector<double> empty_x;
  auto res = autocorrelate(empty_x);
  EXPECT_TRUE(res.empty());
}

TEST_F(FFTUtilsTests, AutocorrelateMatchesMathematicalDefinition) {
  std::vector<double> positive_positions = {1.0, 2.0, 3.0};

  // Linear autocorrelation:
  // R[0] = 1*1 + 2*2 + 3*3 = 14
  // R[1] = 1*2 + 2*3 = 8
  // R[2] = 1*3 = 3
  auto result = autocorrelate(positive_positions);

  ASSERT_EQ(result.size(), positive_positions.size());
  EXPECT_NEAR(result[0], 14.0, 1e-9);
  EXPECT_NEAR(result[1], 8.0, 1e-9);
  EXPECT_NEAR(result[2], 3.0, 1e-9);
}

TEST_F(FFTUtilsTests, AutocorrelateReusesWorkspaceCorrectly) {
  std::vector<double> pos_1 = {1.0, 2.0};
  std::vector<double> pos_2 = {1.0, 2.0, 3.0, 4.0};
  std::vector<std::complex<double>> workspace;

  auto result_1 = autocorrelate(pos_1, workspace);
  EXPECT_EQ(result_1.size(), 2);
#if defined(CORRELATION_USE_FFTW3) || defined(CORRELATION_USE_MKL)
  EXPECT_GE(workspace.size(), 3);
#else
  EXPECT_GE(workspace.size(), 4); // padded to 2^2
#endif

  auto result_2 = autocorrelate(pos_2, workspace);
  EXPECT_EQ(result_2.size(), 4);
#if defined(CORRELATION_USE_FFTW3) || defined(CORRELATION_USE_MKL)
  EXPECT_GE(workspace.size(), 8);
#else
  EXPECT_GE(workspace.size(), 8); // padded to 2^3
#endif
}

// --- Extreme / Edge-Case Tests ---

TEST_F(FFTUtilsTests, ComputeFFTSizeOne) {
  // Size 1 is a valid power of two (2^0)
  std::vector<std::complex<double>> signal = {{7.5, -2.3}};
  auto original = signal;

  // Forward FFT of a single element is itself
  computeFFT(signal, false);
  EXPECT_NEAR(signal[0].real(), original[0].real(), 1e-9);
  EXPECT_NEAR(signal[0].imag(), original[0].imag(), 1e-9);

  // Inverse FFT should return to original
  computeFFT(signal, true);
  EXPECT_NEAR(signal[0].real(), original[0].real(), 1e-9);
  EXPECT_NEAR(signal[0].imag(), original[0].imag(), 1e-9);
}

} // namespace

} // namespace correlation::testing
