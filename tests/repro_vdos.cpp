
TEST(DynamicsAnalyzerTest, VDOSIsZeroAtZeroFrequency) {
  // Create a simple constant VACF (DC signal)
  // Fourier transform of a constant is a delta at f=0.
  // However, VDOS definition usually implies a density of states weighted by
  // frequency squared or similar, or just the Fourier transform. The issue
  // report says "Multiply by frequency", which implies VDOS(f) = f * FT(VACF).
  // If so, at f=0, VDOS should be 0 regardless of the VACF content.

  std::vector<double> vacf(100, 1.0); // Constant VACF
  double dt = 1.0;

  auto [frequencies, intensities, _] =
      DynamicsAnalyzer::calculateVDOS(vacf, dt);

  ASSERT_FALSE(frequencies.empty());
  ASSERT_EQ(frequencies[0], 0.0);

  // Before the fix, this might be non-zero because it's just the integral of
  // VACF (which is sum of 1.0s)
  EXPECT_NEAR(intensities[0], 0.0, 1e-9)
      << "VDOS at 0 THz should be zero due to frequency prefactor";
}
