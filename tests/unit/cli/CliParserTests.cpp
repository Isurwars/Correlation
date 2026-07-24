// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "cli/CliParser.hpp"

#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace {

// ---------------------------------------------------------------------------
// Helper: build a fake argv from a list of strings.
// Returns a vector<char*> whose .data() can be passed to parseArgs().
// The caller keeps `storage` alive for the duration of the test.
// ---------------------------------------------------------------------------
struct ArgBuilder {
  std::vector<std::string> storage;
  std::vector<char *> argv;

  // Construct from initializer list: {"correlation-cli", "--help", ...}
  explicit ArgBuilder(std::initializer_list<std::string> args) : storage(args) {
    for (auto &store : storage) {
      argv.push_back(store.data());
    }
  }

  [[nodiscard]] int argc() const { return static_cast<int>(argv.size()); }
  char **data() { return argv.data(); }
};
class CliParserTests : public ::testing::Test {};
} // namespace

// ===== Default values =====

TEST_F(CliParserTests, DefaultsWithInputFileOnly) {
  ArgBuilder args{"correlation-cli", "input.poscar"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));

  EXPECT_EQ(opts.input_file, "input.poscar");
  EXPECT_EQ(opts.output_base, "input"); // stem of input.poscar
  EXPECT_THAT(opts.r_max, correlation::testing::IsRealEq(20.0));
  EXPECT_THAT(opts.r_bin_width, correlation::testing::IsRealEq(0.02));
  EXPECT_THAT(opts.q_max, correlation::testing::IsRealEq(20.0));
  EXPECT_THAT(opts.q_bin_width, correlation::testing::IsRealEq(0.02));
  EXPECT_THAT(opts.angle_bin_width, correlation::testing::IsRealEq(1.0));
  EXPECT_THAT(opts.dihedral_bin_width, correlation::testing::IsRealEq(1.0));
  EXPECT_FALSE(opts.has_dihedral_bin);
  EXPECT_EQ(opts.min_frame, 0);
  EXPECT_EQ(opts.max_frame, -1);
  EXPECT_TRUE(opts.csv);
  EXPECT_FALSE(opts.hdf5);
  EXPECT_FALSE(opts.parquet);
  EXPECT_TRUE(opts.smoothing);
  EXPECT_FALSE(opts.quiet);
  EXPECT_TRUE(opts.disable_groups.empty());
  EXPECT_THAT(opts.time_step, correlation::testing::IsRealEq(1.0));
  EXPECT_THAT(opts.r_int_max, correlation::testing::IsRealEq(10.0));
  EXPECT_EQ(opts.max_ring_size, 8);
  EXPECT_THAT(opts.smoothing_sigma, correlation::testing::IsRealEq(0.1));
  EXPECT_EQ(opts.smoothing_kernel, correlation::math::KernelType::Gaussian);
  EXPECT_FALSE(opts.show_version);
  EXPECT_FALSE(opts.show_help);
}

// ===== Help and version flags =====

TEST_F(CliParserTests, HelpLongFlag) {
  ArgBuilder args{"correlation-cli", "--help"};
  correlation::cli::CliOptions opts;

  EXPECT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_TRUE(opts.show_help);
}

TEST_F(CliParserTests, HelpShortFlag) {
  ArgBuilder args{"correlation-cli", "-h"};
  correlation::cli::CliOptions opts;

  EXPECT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_TRUE(opts.show_help);
}

TEST_F(CliParserTests, VersionLongFlag) {
  ArgBuilder args{"correlation-cli", "--version"};
  correlation::cli::CliOptions opts;

  EXPECT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_TRUE(opts.show_version);
}

TEST_F(CliParserTests, VersionShortFlag) {
  ArgBuilder args{"correlation-cli", "-v"};
  correlation::cli::CliOptions opts;

  EXPECT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_TRUE(opts.show_version);
}

// ===== Error cases =====

TEST_F(CliParserTests, NoArgsReturnsFalse) {
  ArgBuilder args{"correlation-cli"};
  correlation::cli::CliOptions opts;

  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
}

TEST_F(CliParserTests, UnknownOptionReturnsFalse) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--unknown-flag"};
  correlation::cli::CliOptions opts;

  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
}

TEST_F(CliParserTests, OnlyDashFlagNoFileReturnsFalse) {
  // A single unknown dash-prefixed arg with no input file
  ArgBuilder args{"correlation-cli", "--no-smoothing"};
  correlation::cli::CliOptions opts;

  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
}

// ===== Numeric options =====

TEST_F(CliParserTests, RMaxOption) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--r-max", "35.5"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_THAT(opts.r_max, correlation::testing::IsRealEq(35.5));
}

TEST_F(CliParserTests, RBinOption) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--r-bin", "0.05"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_THAT(opts.r_bin_width, correlation::testing::IsRealEq(0.05));
}

TEST_F(CliParserTests, QMaxOption) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--q-max", "30.0"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_THAT(opts.q_max, correlation::testing::IsRealEq(30.0));
}

TEST_F(CliParserTests, QBinOption) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--q-bin", "0.1"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_THAT(opts.q_bin_width, correlation::testing::IsRealEq(0.1));
}

TEST_F(CliParserTests, AngleBinOption) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--angle-bin", "2.0"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_THAT(opts.angle_bin_width, correlation::testing::IsRealEq(2.0));
}

TEST_F(CliParserTests, DihedralBinOption) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--dihedral-bin", "3.0"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_THAT(opts.dihedral_bin_width, correlation::testing::IsRealEq(3.0));
  EXPECT_TRUE(opts.has_dihedral_bin);
}

TEST_F(CliParserTests, DihedralBinDefaultsToCopyOfAngleBin) {
  ArgBuilder args{"correlation-cli", "f.poscar"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_FALSE(opts.has_dihedral_bin);
  // Default dihedral_bin_width equals the default angle_bin_width
  EXPECT_THAT(opts.dihedral_bin_width, correlation::testing::IsRealEq(opts.angle_bin_width));
}

TEST_F(CliParserTests, TimeStepOption) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--time-step", "0.5"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_THAT(opts.time_step, correlation::testing::IsRealEq(0.5));
}

TEST_F(CliParserTests, RIntMaxOption) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--r-int-max", "15.0"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_THAT(opts.r_int_max, correlation::testing::IsRealEq(15.0));
}

TEST_F(CliParserTests, MaxRingSizeOption) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--max-ring-size", "12"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_EQ(opts.max_ring_size, 12);
}

TEST_F(CliParserTests, SmoothingSigmaOption) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--smoothing-sigma", "0.5"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_THAT(opts.smoothing_sigma, correlation::testing::IsRealEq(0.5));
}

// ===== Frame options =====

TEST_F(CliParserTests, MinFrameConvertedToZeroBased) {
  // User passes 1-based frame 5 → should store 4
  ArgBuilder args{"correlation-cli", "f.poscar", "--min-frame", "5"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_EQ(opts.min_frame, 4);
}

TEST_F(CliParserTests, MinFrameZeroClamped) {
  // User passes 0 (invalid 1-based) → 0-1 = -1, clamped to 0
  ArgBuilder args{"correlation-cli", "f.poscar", "--min-frame", "0"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_EQ(opts.min_frame, 0);
}

TEST_F(CliParserTests, MaxFrameOption) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--max-frame", "100"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_EQ(opts.max_frame, 100);
}

// ===== Boolean toggle options =====

TEST_F(CliParserTests, CsvToggle) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--no-csv", "--hdf5"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_FALSE(opts.csv);
  EXPECT_TRUE(opts.hdf5);
}

TEST_F(CliParserTests, CsvEnableExplicit) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--no-csv", "--csv"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_TRUE(opts.csv);
}

TEST_F(CliParserTests, Hdf5Toggle) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--hdf5"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_TRUE(opts.hdf5);
}

TEST_F(CliParserTests, Hdf5DisableToggle) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--hdf5", "--no-hdf5"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_FALSE(opts.hdf5);
}

TEST_F(CliParserTests, ParquetToggle) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--parquet"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_TRUE(opts.parquet);
}

TEST_F(CliParserTests, ParquetDisableToggle) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--parquet", "--no-parquet"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_FALSE(opts.parquet);
}

TEST_F(CliParserTests, NoSmoothingOption) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--no-smoothing"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_FALSE(opts.smoothing);
}

// ===== Quiet flag =====

TEST_F(CliParserTests, QuietLongFlag) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--quiet"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_TRUE(opts.quiet);
}

TEST_F(CliParserTests, QuietShortFlag) {
  ArgBuilder args{"correlation-cli", "f.poscar", "-q"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_TRUE(opts.quiet);
}

// ===== Output option =====

TEST_F(CliParserTests, OutputLongFlag) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--output", "/tmp/results"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_EQ(opts.output_base, "/tmp/results");
}

TEST_F(CliParserTests, OutputShortFlag) {
  ArgBuilder args{"correlation-cli", "f.poscar", "-o", "out_base"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_EQ(opts.output_base, "out_base");
}

TEST_F(CliParserTests, DefaultOutputBaseIsStemOfInput) {
  ArgBuilder args{"correlation-cli", "/data/experiment/sample.poscar"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_EQ(opts.output_base, (std::filesystem::path("/data/experiment") / "sample").lexically_normal().string());
}

TEST_F(CliParserTests, DisableGroupsOption) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--disable-groups", "radial,scattering"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_EQ(opts.disable_groups, "radial,scattering");
}

TEST_F(CliParserTests, DisableGroupsShortOption) {
  ArgBuilder args{"correlation-cli", "f.poscar", "-g", "structural"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_EQ(opts.disable_groups, "structural");
}

// ===== Smoothing kernel parsing =====

TEST_F(CliParserTests, KernelGaussian) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--smoothing-kernel", "gaussian"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_EQ(opts.smoothing_kernel, correlation::math::KernelType::Gaussian);
}

TEST_F(CliParserTests, KernelGauss) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--smoothing-kernel", "gauss"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_EQ(opts.smoothing_kernel, correlation::math::KernelType::Gaussian);
}

TEST_F(CliParserTests, KernelGaussianCaseInsensitive) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--smoothing-kernel", "GAUSSIAN"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_EQ(opts.smoothing_kernel, correlation::math::KernelType::Gaussian);
}

TEST_F(CliParserTests, KernelBump) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--smoothing-kernel", "bump"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_EQ(opts.smoothing_kernel, correlation::math::KernelType::Bump);
}

TEST_F(CliParserTests, KernelTriweight) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--smoothing-kernel", "triweight"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_EQ(opts.smoothing_kernel, correlation::math::KernelType::Triweight);
}

TEST_F(CliParserTests, KernelUnknownDefaultsToGaussian) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--smoothing-kernel", "invalid_kernel"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_EQ(opts.smoothing_kernel, correlation::math::KernelType::Gaussian);
}

// ===== Combined options =====

TEST_F(CliParserTests, MultipleOptionsCombined) {
  ArgBuilder args{"correlation-cli",
                  "sim.xyz",
                  "--r-max",
                  "30.0",
                  "--r-bin",
                  "0.05",
                  "--q-max",
                  "25.0",
                  "--no-csv",
                  "--hdf5",
                  "--quiet",
                  "--smoothing-kernel",
                  "bump",
                  "--disable-groups",
                  "scattering,rings",
                  "-o",
                  "/results/out"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_EQ(opts.input_file, "sim.xyz");
  EXPECT_THAT(opts.r_max, correlation::testing::IsRealEq(30.0));
  EXPECT_THAT(opts.r_bin_width, correlation::testing::IsRealEq(0.05));
  EXPECT_THAT(opts.q_max, correlation::testing::IsRealEq(25.0));
  EXPECT_FALSE(opts.csv);
  EXPECT_TRUE(opts.hdf5);
  EXPECT_TRUE(opts.quiet);
  EXPECT_EQ(opts.smoothing_kernel, correlation::math::KernelType::Bump);
  EXPECT_EQ(opts.disable_groups, "scattering,rings");
  EXPECT_EQ(opts.output_base, "/results/out");
}

// ===== Help with input file (help takes precedence) =====

TEST_F(CliParserTests, HelpAfterInputFileStillSetsHelp) {
  ArgBuilder args{"correlation-cli", "f.poscar", "--help"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_TRUE(opts.show_help);
  EXPECT_EQ(opts.input_file, "f.poscar");
}

TEST_F(CliParserTests, VersionAfterInputFileStillSetsVersion) {
  ArgBuilder args{"correlation-cli", "f.poscar", "-v"};
  correlation::cli::CliOptions opts;

  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_TRUE(opts.show_version);
}

// ===== Validation failure tests =====

TEST_F(CliParserTests, RMaxMustBePositive) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--r-max", "0.0"};
  correlation::cli::CliOptions opts;
  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));

  ArgBuilder args2{"correlation-cli", "input.poscar", "--r-max", "-5.0"};
  EXPECT_FALSE(correlation::cli::parseArgs(args2.argc(), args2.data(), opts));
}

TEST_F(CliParserTests, RBinMustBePositive) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--r-bin", "0.0"};
  correlation::cli::CliOptions opts;
  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));

  ArgBuilder args2{"correlation-cli", "input.poscar", "--r-bin", "-0.1"};
  EXPECT_FALSE(correlation::cli::parseArgs(args2.argc(), args2.data(), opts));
}

TEST_F(CliParserTests, RBinCannotBeGreaterThanOrEqualRMax) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--r-max", "10.0", "--r-bin", "10.0"};
  correlation::cli::CliOptions opts;
  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));

  ArgBuilder args2{"correlation-cli", "input.poscar", "--r-max", "10.0", "--r-bin", "12.5"};
  EXPECT_FALSE(correlation::cli::parseArgs(args2.argc(), args2.data(), opts));
}

TEST_F(CliParserTests, QMaxMustBePositive) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--q-max", "0.0"};
  correlation::cli::CliOptions opts;
  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));

  ArgBuilder args2{"correlation-cli", "input.poscar", "--q-max", "-1.0"};
  EXPECT_FALSE(correlation::cli::parseArgs(args2.argc(), args2.data(), opts));
}

TEST_F(CliParserTests, QBinMustBePositive) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--q-bin", "0.0"};
  correlation::cli::CliOptions opts;
  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));

  ArgBuilder args2{"correlation-cli", "input.poscar", "--q-bin", "-0.01"};
  EXPECT_FALSE(correlation::cli::parseArgs(args2.argc(), args2.data(), opts));
}

TEST_F(CliParserTests, QBinCannotBeGreaterThanOrEqualQMax) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--q-max", "15.0", "--q-bin", "15.0"};
  correlation::cli::CliOptions opts;
  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));

  ArgBuilder args2{"correlation-cli", "input.poscar", "--q-max", "15.0", "--q-bin", "16.0"};
  EXPECT_FALSE(correlation::cli::parseArgs(args2.argc(), args2.data(), opts));
}

TEST_F(CliParserTests, AngleBinMustBePositive) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--angle-bin", "0.0"};
  correlation::cli::CliOptions opts;
  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
}

TEST_F(CliParserTests, AngleBinCannotExceed180) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--angle-bin", "180.1"};
  correlation::cli::CliOptions opts;
  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
}

TEST_F(CliParserTests, DihedralBinMustBePositive) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--dihedral-bin", "0.0"};
  correlation::cli::CliOptions opts;
  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
}

TEST_F(CliParserTests, DihedralBinCannotExceed360) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--dihedral-bin", "360.5"};
  correlation::cli::CliOptions opts;
  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
}

TEST_F(CliParserTests, TimeStepMustBePositive) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--time-step", "0.0"};
  correlation::cli::CliOptions opts;
  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
}

TEST_F(CliParserTests, RIntMaxMustBePositive) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--r-int-max", "0.0"};
  correlation::cli::CliOptions opts;
  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
}

TEST_F(CliParserTests, MaxRingSizeMustBePositive) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--max-ring-size", "0"};
  correlation::cli::CliOptions opts;
  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
}

TEST_F(CliParserTests, MaxFrameCannotBeLessThanMinusOne) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--max-frame", "-2"};
  correlation::cli::CliOptions opts;
  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
}

TEST_F(CliParserTests, MinFrameCannotBeGreaterThanMaxFrame) {
  // Use a min-frame of 5 (maps to index 4) and max-frame of 3 (maps to index 3)
  ArgBuilder args{"correlation-cli", "input.poscar", "--min-frame", "5", "--max-frame", "3"};
  correlation::cli::CliOptions opts;
  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
}

TEST_F(CliParserTests, SmoothingSigmaCannotBeNegative) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--smoothing-sigma", "-0.1"};
  correlation::cli::CliOptions opts;
  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
}

TEST_F(CliParserTests, SmoothingSigmaMustBePositiveWhenSmoothingEnabled) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--smoothing-sigma", "0.0"};
  correlation::cli::CliOptions opts;
  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
}

TEST_F(CliParserTests, AtLeastOneOutputFormatMustBeEnabled) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--no-csv", "--no-hdf5", "--no-parquet"};
  correlation::cli::CliOptions opts;
  EXPECT_FALSE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
}

TEST_F(CliParserTests, LiquidMaterialDefaults) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--material", "liquid"};
  correlation::cli::CliOptions opts;
  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_EQ(opts.material_type, 1);
  EXPECT_THAT(opts.r_bin_width, correlation::testing::IsRealEq(0.05));
  EXPECT_THAT(opts.q_bin_width, correlation::testing::IsRealEq(0.05));
  EXPECT_THAT(opts.angle_bin_width, correlation::testing::IsRealEq(2.0));
  EXPECT_THAT(opts.dihedral_bin_width, correlation::testing::IsRealEq(2.0));
  EXPECT_THAT(opts.smoothing_sigma, correlation::testing::IsRealEq(0.15));
}

TEST_F(CliParserTests, CrystallineMaterialDefaults) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--material", "crystalline"};
  correlation::cli::CliOptions opts;
  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_EQ(opts.material_type, 2);
  EXPECT_THAT(opts.r_bin_width, correlation::testing::IsRealEq(0.002));
  EXPECT_THAT(opts.q_bin_width, correlation::testing::IsRealEq(0.002));
  EXPECT_THAT(opts.angle_bin_width, correlation::testing::IsRealEq(0.1));
  EXPECT_THAT(opts.dihedral_bin_width, correlation::testing::IsRealEq(0.1));
  EXPECT_THAT(opts.smoothing_sigma, correlation::testing::IsRealEq(0.01));
}

TEST_F(CliParserTests, MaterialDefaultsWithOverrides) {
  ArgBuilder args{"correlation-cli", "input.poscar", "--material",        "crystalline",
                  "--r-bin",         "0.01",         "--smoothing-sigma", "0.05"};
  correlation::cli::CliOptions opts;
  ASSERT_TRUE(correlation::cli::parseArgs(args.argc(), args.data(), opts));
  EXPECT_EQ(opts.material_type, 2);
  EXPECT_THAT(opts.r_bin_width, correlation::testing::IsRealEq(0.01));       // Overridden
  EXPECT_THAT(opts.q_bin_width, correlation::testing::IsRealEq(0.002));      // Crystal default
  EXPECT_THAT(opts.angle_bin_width, correlation::testing::IsRealEq(0.1));    // Crystal default
  EXPECT_THAT(opts.dihedral_bin_width, correlation::testing::IsRealEq(0.1)); // Crystal default
  EXPECT_THAT(opts.smoothing_sigma, correlation::testing::IsRealEq(0.05));   // Overridden
}
