// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <array>
#include <cstdio>
#include <filesystem>
#include <gtest/gtest.h>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>
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
  std::string cmd = std::string(CLI_BIN) + " " + args;

  STARTUPINFOA si = {sizeof(STARTUPINFOA)};
  PROCESS_INFORMATION pi = {};

  si.dwFlags = STARTF_USESTDHANDLES;
  si.hStdError = nullptr;

  if (CreateProcessA(nullptr, cmd.data(), nullptr, nullptr, FALSE, CREATE_NO_WINDOW, nullptr, nullptr, &si, &pi) == 0) {
    return -1;
  }

  WaitForSingleObject(pi.hProcess, INFINITE);
  DWORD exitCode = 0;
  GetExitCodeProcess(pi.hProcess, &exitCode);

  CloseHandle(pi.hProcess);
  CloseHandle(pi.hThread);
  return static_cast<int>(exitCode);

#else // POSIX (Linux / macOS)
  pid_t const pid = fork();
  if (pid == 0) {
    // Child Process: Redirect stderr to /dev/null safely without shell invocation
    int const dev_null = open("/dev/null", O_WRONLY, 0);
    if (dev_null != -1) {
      dup2(dev_null, STDERR_FILENO);
      close(dev_null);
    }

    std::vector<std::string> tokens;
    tokens.emplace_back(CLI_BIN);

    std::istringstream iss(args);
    std::string token;
    while (iss >> token) {
      tokens.push_back(token);
    }

    std::vector<char *> argv;
    argv.reserve(tokens.size() + 1);
    for (auto &tokenRef : tokens) {
      argv.push_back(tokenRef.data());
    }
    argv.push_back(nullptr);

    execv(std::string(CLI_BIN).c_str(), argv.data());
    _exit(127); // Exit child if execv fails
  }

  if (pid > 0) {
    int status = 0;
    waitpid(pid, &status, 0);
    return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
  }

  return -1;
#endif
}

// Run a shell command and capture its stdout.
std::string runCliCapture(const std::string &args) {
#ifdef _WIN32
  std::string cmd = std::string(CLI_BIN) + " " + args;

  SECURITY_ATTRIBUTES saAttr{};
  saAttr.nLength = sizeof(SECURITY_ATTRIBUTES);
  saAttr.bInheritHandle = TRUE;

  HANDLE hReadPipe = nullptr;
  HANDLE hWritePipe = nullptr;
  if (CreatePipe(&hReadPipe, &hWritePipe, &saAttr, 0) == 0) {
    return "";
  }
  SetHandleInformation(hReadPipe, HANDLE_INHERIT, 0);

  STARTUPINFOA si{};
  si.cb = sizeof(STARTUPINFOA);
  si.dwFlags = STARTF_USESTDHANDLES | STARTF_USESHOWWINDOW;
  si.hStdOutput = hWritePipe;
  si.hStdError = nullptr;
  si.wShowWindow = SW_HIDE;

  PROCESS_INFORMATION pi{};
  if (CreateProcessA(nullptr, cmd.data(), nullptr, nullptr, TRUE, CREATE_NO_WINDOW, nullptr, nullptr, &si, &pi) == 0) {
    CloseHandle(hReadPipe);
    CloseHandle(hWritePipe);
    return "";
  }

  CloseHandle(hWritePipe);

  std::string result;
  std::array<char, 512> buffer{};
  DWORD bytesRead = 0;

  while (ReadFile(hReadPipe, buffer.data(), static_cast<DWORD>(buffer.size()), &bytesRead, nullptr) != 0 &&
         bytesRead > 0) {
    result.append(buffer.data(), bytesRead);
  }

  CloseHandle(hReadPipe);
  CloseHandle(pi.hProcess);
  CloseHandle(pi.hThread);
  return result;

#else // POSIX (Linux / macOS)
  std::array<int, 2> pipefd{};
  if (pipe(pipefd.data()) == -1) {
    return "";
  }

  pid_t const pid = fork();
  if (pid == -1) {
    close(pipefd[0]);
    close(pipefd[1]);
    return "";
  }

  if (pid == 0) {
    dup2(pipefd[1], STDOUT_FILENO);
    close(pipefd[0]);
    close(pipefd[1]);

    int const dev_null = open("/dev/null", O_WRONLY, 0);
    if (dev_null != -1) {
      dup2(dev_null, STDERR_FILENO);
      close(dev_null);
    }

    std::vector<std::string> tokens;
    tokens.emplace_back(CLI_BIN);

    std::istringstream iss(args);
    std::string token;
    while (iss >> token) {
      tokens.push_back(token);
    }

    std::vector<char *> argv;
    argv.reserve(tokens.size() + 1);
    for (auto &tokenRef : tokens) {
      argv.push_back(tokenRef.data());
    }
    argv.push_back(nullptr);

    execv(std::string(CLI_BIN).c_str(), argv.data());
    _exit(127);
  }

  close(pipefd[1]);

  std::string result;
  std::array<char, 512> buffer{};
  ssize_t bytesRead = 0;

  while ((bytesRead = read(pipefd[0], buffer.data(), buffer.size())) > 0) {
    result.append(buffer.data(), static_cast<size_t>(bytesRead));
  }

  close(pipefd[0]);
  int status = 0;
  waitpid(pid, &status, 0);
  return result;
#endif
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
