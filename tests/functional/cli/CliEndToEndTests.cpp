// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <array>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <string_view>

#ifdef _WIN32
#define popen _popen
#define pclose _pclose
#endif

namespace {

// ---------------------------------------------------------------------------
// The path to the correlation-cli binary is injected at compile time via
// -DCORRELATION_CLI_PATH="..." in CMakeLists.txt.
// ---------------------------------------------------------------------------
#ifndef CORRELATION_CLI_PATH
#error "CORRELATION_CLI_PATH must be defined — check tests/CMakeLists.txt"
#endif

constexpr std::string_view CLI_BIN = CORRELATION_CLI_PATH;

// Helper to find the test data directory relative to the build dir.
std::string getTestDataDir() {
  std::vector<std::string> const candidates = {
      "../../tests/data/", // build/tests -> tests/data
      "../tests/data/",    // build -> tests/data
      "tests/data/",       // project root
      "data/",             // if cwd is tests
  };
  for (const auto &dir : candidates) {
    if (std::filesystem::exists(dir + "vasp/Si.poscar")) {
      return dir + "vasp/";
    }
  }
  return "../../tests/data/vasp/";
}

// Run a shell command and return its exit code.
int runCli(const std::string &args) {
#ifdef _WIN32
  std::string cmd = std::string(CLI_BIN) + " " + args + " 2>nul";
#else
  std::string const cmd = std::string(CLI_BIN) + " " + args + " 2>/dev/null";
#endif
  // NOLINTNEXTLINE(cert-env33-c)
  int const status = std::system(cmd.c_str());

#ifdef _WIN32
  return status;
#else
  return WEXITSTATUS(status);
#endif
}

// Run a shell command and capture its stdout.
std::string runCliCapture(const std::string &args) {
#ifdef _WIN32
  std::string cmd = std::string(CLI_BIN) + " " + args + " 2>nul";
#else
  std::string const cmd = std::string(CLI_BIN) + " " + args + " 2>/dev/null";
#endif
  std::array<char, 512> buffer{};
  std::string result;
  // NOLINTNEXTLINE(cert-env33-c)
  FILE *pipe = popen(cmd.c_str(), "r");
  if (pipe == nullptr) {
    return "";
  }
  while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe) != nullptr) {
    result += buffer.data();
  }
  pclose(pipe);
  return result;
}

class CliEndToEndTests : public ::testing::Test {
protected:
  [[nodiscard]] const std::string &dataDir() const { return data_dir_; }

  void SetUp() override { data_dir_ = getTestDataDir(); }

private:
  std::string data_dir_;
};

// ===== Basic flag tests =====

TEST_F(CliEndToEndTests, HelpReturnsZero) { EXPECT_EQ(runCli("--help"), 0); }

TEST_F(CliEndToEndTests, VersionReturnsZero) { EXPECT_EQ(runCli("--version"), 0); }

TEST_F(CliEndToEndTests, VersionPrintsVersionString) {
  std::string const output = runCliCapture("--version");
  EXPECT_NE(output.find("Correlation version"), std::string::npos) << "Actual output: " << output;
}

TEST_F(CliEndToEndTests, NoArgsReturnsNonZero) { EXPECT_NE(runCli(""), 0); }

TEST_F(CliEndToEndTests, UnknownOptionReturnsNonZero) { EXPECT_NE(runCli("--this-does-not-exist"), 0); }

// ===== File loading tests =====

TEST_F(CliEndToEndTests, NonexistentFileReturnsNonZero) { EXPECT_NE(runCli("nonexistent_file_xyz.poscar --quiet"), 0); }

TEST_F(CliEndToEndTests, ValidFileRunsSuccessfully) {
  // Use a temp directory for output so we don't pollute the test data dir
  auto tmp_dir = std::filesystem::temp_directory_path() / "correlation_e2e_test";
  std::filesystem::create_directories(tmp_dir);
  std::string const out_base = (tmp_dir / "result").string();

  std::string const input = dataDir() + "Si.poscar";
  int const run_cli_status = runCli(input + " --quiet -o " + out_base + " --r-max 10 --r-bin 0.1");

  // Clean up
  std::filesystem::remove_all(tmp_dir);

  EXPECT_EQ(run_cli_status, 0) << "CLI should succeed on valid Si.poscar input";
}

TEST_F(CliEndToEndTests, ShortHelpFlag) { EXPECT_EQ(runCli("-h"), 0); }

TEST_F(CliEndToEndTests, ShortVersionFlag) { EXPECT_EQ(runCli("-v"), 0); }

TEST_F(CliEndToEndTests, DefaultExecutesAllCalculators) {
  auto tmp_dir = std::filesystem::temp_directory_path() / "correlation_e2e_default_groups";
  std::filesystem::create_directories(tmp_dir);
  std::string const out_base = (tmp_dir / "result").string();

  std::string const input = dataDir() + "Si.poscar";
  // Run with no disable-groups flags -> should default to running all groups
  int const run_cli_status = runCli(input + " --quiet -o " + out_base + " --r-max 10 --r-bin 0.1");

  EXPECT_EQ(run_cli_status, 0);

  // Radial files should exist
  EXPECT_TRUE(std::filesystem::exists(out_base + "_g.csv"));
  EXPECT_TRUE(std::filesystem::exists(out_base + "_J.csv"));
  EXPECT_TRUE(std::filesystem::exists(out_base + "_G.csv"));

  // Structural/rings files should also exist because they are enabled by default
  EXPECT_TRUE(std::filesystem::exists(out_base + "_CN.csv"));

  std::filesystem::remove_all(tmp_dir);
}

TEST_F(CliEndToEndTests, DisableRadialAndScatteringGroups) {
  auto tmp_dir = std::filesystem::temp_directory_path() / "correlation_e2e_disable_radial";
  std::filesystem::create_directories(tmp_dir);
  std::string const out_base = (tmp_dir / "result").string();

  std::string const input = dataDir() + "Si.poscar";
  // Run with radial and scattering groups disabled
  int const run_cli_status = runCli(input + " --quiet -o " + out_base + " --disable-groups radial,scattering");

  EXPECT_EQ(run_cli_status, 0);

  // CN (structural) should exist
  EXPECT_TRUE(std::filesystem::exists(out_base + "_CN.csv"));

  // Radial files should NOT exist
  EXPECT_FALSE(std::filesystem::exists(out_base + "_g.csv"));
  EXPECT_FALSE(std::filesystem::exists(out_base + "_J.csv"));
  EXPECT_FALSE(std::filesystem::exists(out_base + "_G_reduced.csv"));

  std::filesystem::remove_all(tmp_dir);
}

TEST_F(CliEndToEndTests, DisableStructuralGroup) {
  auto tmp_dir = std::filesystem::temp_directory_path() / "correlation_e2e_disable_structural";
  std::filesystem::create_directories(tmp_dir);
  std::string const out_base = (tmp_dir / "result").string();

  std::string const input = dataDir() + "Si.poscar";
  // Run disabling structural group
  int const run_cli_status =
      runCli(input + " --quiet -o " + out_base + " --disable-groups structural --r-max 10 --r-bin 0.1");

  EXPECT_EQ(run_cli_status, 0);

  // Radial should exist
  EXPECT_TRUE(std::filesystem::exists(out_base + "_g.csv"));

  // CN (structural) should NOT exist
  EXPECT_FALSE(std::filesystem::exists(out_base + "_CN.csv"));

  std::filesystem::remove_all(tmp_dir);
}

} // namespace
