// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <gtest/gtest.h>

#include "../include/AppBackend.hpp"
#include <stdexcept>

// Test fixture for the AppBackend class
class Test15_AppBackend : public ::testing::Test {};

TEST_F(Test15_AppBackend, DefaultConstructorInitializesCorrectly) {
  // Arrange & Act
  AppBackend backend;
  ProgramOptions opts = backend.options();

  // Assert
  EXPECT_DOUBLE_EQ(opts.r_max, AppDefaults::R_MAX);
  EXPECT_DOUBLE_EQ(opts.r_bin_width, AppDefaults::R_BIN_WIDTH);
  EXPECT_DOUBLE_EQ(opts.q_max, AppDefaults::Q_MAX);
  EXPECT_DOUBLE_EQ(opts.q_bin_width, AppDefaults::Q_BIN_WIDTH);
  EXPECT_DOUBLE_EQ(opts.r_int_max, AppDefaults::R_INT_MAX);
  EXPECT_DOUBLE_EQ(opts.angle_bin_width, AppDefaults::ANGLE_BIN_WIDTH);
  EXPECT_DOUBLE_EQ(opts.smoothing_sigma, AppDefaults::SMOOTHING_SIGMA);
  EXPECT_EQ(opts.smoothing_kernel, AppDefaults::SMOOTHING_KERNEL);
  EXPECT_DOUBLE_EQ(opts.time_step, AppDefaults::TIME_STEP);

  EXPECT_EQ(backend.cell(), nullptr);
  EXPECT_EQ(backend.getFrameCount(), 0);
  EXPECT_EQ(backend.getTotalAtomCount(), 0);
}

TEST_F(Test15_AppBackend, SetOptionsModifiesState) {
  // Arrange
  AppBackend backend;
  ProgramOptions opts;
  opts.r_max = 50.0;
  opts.smoothing_sigma = 0.5;
  opts.min_frame = 10;
  opts.max_frame = 20;

  // Act
  backend.setOptions(opts);

  // Assert
  ProgramOptions new_opts = backend.options();
  EXPECT_DOUBLE_EQ(new_opts.r_max, 50.0);
  EXPECT_DOUBLE_EQ(new_opts.smoothing_sigma, 0.5);
  EXPECT_EQ(new_opts.min_frame, 10);
  EXPECT_EQ(new_opts.max_frame, 20);
}

TEST_F(Test15_AppBackend, LoadInvalidFileThrowsException) {
  // Arrange
  AppBackend backend;

  // Act & Assert
  EXPECT_THROW(backend.load_file("nonexistent_file_that_should_not_exist.xyz"),
               std::runtime_error);
}

TEST_F(Test15_AppBackend, RecommendedTimeStepWithNoCellReturnsDefault) {
  // Arrange
  AppBackend backend;

  // Act
  double time_step = backend.getRecommendedTimeStep();

  // Assert
  EXPECT_DOUBLE_EQ(time_step, AppDefaults::TIME_STEP);
}

// Additional rigorous test for recommended step logic should be added by
// loading a mock cell, but since AppBackend deals with reading from file,
// creating an artificial trajectory internally is complex here. We can assume
// the math is straightforward.
