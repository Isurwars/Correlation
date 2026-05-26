// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <array>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <gtest/gtest.h>
#include <string>

namespace {

// ---------------------------------------------------------------------------
// The path to the correlation-cli binary is injected at compile time via
// -DCORRELATION_CLI_PATH="..." in CMakeLists.txt.
// ---------------------------------------------------------------------------
#ifndef CORRELATION_CLI_PATH
#error "CORRELATION_CLI_PATH must be defined — check tests/CMakeLists.txt"
#endif

const std::string CLI_BIN = CORRELATION_CLI_PATH;

// Helper to find the test data directory relative to the build dir.
std::string getTestDataDir() {
  std::vector<std::string> candidates = {
      "../../tests/data/", // build/tests -> tests/data
      "../tests/data/",    // build -> tests/data
      "tests/data/",       // project root
      "data/",             // if cwd is tests
  };
  for (const auto &dir : candidates) {
    if (std::filesystem::exists(dir + "Si.poscar")) {
      return dir;
    }
  }
  return "../../tests/data/";
}

// Run a shell command and return its exit code.
int runCli(const std::string &args) {
  std::string cmd = CLI_BIN + " " + args + " 2>/dev/null";
  int status = std::system(cmd.c_str());
#ifdef _WIN32
  return status;
#else
  return WEXITSTATUS(status);
#endif
}

// Run a shell command and capture its stdout.
std::string runCliCapture(const std::string &args) {
  std::string cmd = CLI_BIN + " " + args + " 2>/dev/null";
  std::array<char, 512> buffer{};
  std::string result;
  FILE *pipe = popen(cmd.c_str(), "r");
  if (!pipe)
    return "";
  while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe) !=
         nullptr) {
    result += buffer.data();
  }
  pclose(pipe);
  return result;
}

} // namespace

class CliEndToEndTests : public ::testing::Test {
protected:
  std::string data_dir_;
  void SetUp() override { data_dir_ = getTestDataDir(); }
};

// ===== Basic flag tests =====

TEST_F(CliEndToEndTests, HelpReturnsZero) { EXPECT_EQ(runCli("--help"), 0); }

TEST_F(CliEndToEndTests, VersionReturnsZero) {
  EXPECT_EQ(runCli("--version"), 0);
}

TEST_F(CliEndToEndTests, VersionPrintsVersionString) {
  std::string output = runCliCapture("--version");
  EXPECT_NE(output.find("Correlation version"), std::string::npos)
      << "Actual output: " << output;
}

TEST_F(CliEndToEndTests, NoArgsReturnsNonZero) { EXPECT_NE(runCli(""), 0); }

TEST_F(CliEndToEndTests, UnknownOptionReturnsNonZero) {
  EXPECT_NE(runCli("--this-does-not-exist"), 0);
}

// ===== File loading tests =====

TEST_F(CliEndToEndTests, NonexistentFileReturnsNonZero) {
  EXPECT_NE(runCli("nonexistent_file_xyz.poscar --quiet"), 0);
}

TEST_F(CliEndToEndTests, ValidFileRunsSuccessfully) {
  // Use a temp directory for output so we don't pollute the test data dir
  auto tmp_dir =
      std::filesystem::temp_directory_path() / "correlation_e2e_test";
  std::filesystem::create_directories(tmp_dir);
  std::string out_base = (tmp_dir / "result").string();

  std::string input = data_dir_ + "Si.poscar";
  int rc = runCli(input + " --quiet -o " + out_base +
                  " --calculators RDF --r-max 10 --r-bin 0.1");

  // Clean up
  std::filesystem::remove_all(tmp_dir);

  EXPECT_EQ(rc, 0) << "CLI should succeed on valid Si.poscar input";
}

TEST_F(CliEndToEndTests, ShortHelpFlag) { EXPECT_EQ(runCli("-h"), 0); }

TEST_F(CliEndToEndTests, ShortVersionFlag) { EXPECT_EQ(runCli("-v"), 0); }
