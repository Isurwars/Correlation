#!/usr/bin/env python3
"""Unit tests for Correlation lifecycle hooks. Run: `python3 -m unittest test_hooks` from this dir."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from hook_runner import (
    audit_cmake_file,
    audit_cpp_file,
    audit_slint_file,
    estimate_cognitive_complexity,
    evaluate_command_safety,
    handle_pre_tool_use,
    parse_porcelain_z,
)


def _write_temp(suffix: str, content: str) -> Path:
    with tempfile.NamedTemporaryFile(suffix=suffix, mode="w", delete=False) as handle:
        handle.write(content)
        return Path(handle.name)


class TestHookSafety(unittest.TestCase):
    """PreToolUse command evaluation."""

    def assert_unsafe(self, commands: list[str], fragment: str = "") -> None:
        for cmd in commands:
            safe, reason = evaluate_command_safety(cmd)
            self.assertFalse(safe, f"Expected unsafe for: {cmd}")
            if fragment:
                self.assertIn(fragment, reason, cmd)

    def test_destructive_rm_patterns(self) -> None:
        self.assert_unsafe(
            [
                "rm -rf build",
                "rm -fr build/",
                "rm -r -f target/",
                "rm -f -r some_dir",
                "rm -Rf /tmp/foo",
                "rm --recursive --force .",
                "rm --force --recursive node_modules",
                'rm -rf "path with spaces"',
            ],
            "recursive forced deletion",
        )

    def test_protected_paths_deletion(self) -> None:
        self.assert_unsafe(
            [
                "rm -r build",
                "rmdir build",
                "unlink build",
                "rm CMakeCache.txt",
                "rm -r graphify-out",
                "rmdir graphify-out/",
                "rm build/CMakeCache.txt",
            ]
        )

    def test_bypass_attempts_are_blocked(self) -> None:
        self.assert_unsafe(
            [
                "/usr/bin/rm -rf foo",
                "rm -r ./build",
                "rm -r build/../build",
                "sudo rm -rf foo",
                "env FOO=1 rm -rf foo",
                "ls | xargs rm -rf",
                "echo hi && rm -rf foo",
                "true;rm -rf foo",
                "$(rm -rf foo)",
                'bash -c "rm -rf foo"',
                "find build -delete",
                "git clean -fdx",
                "git clean -fX",
                "cmake -E rm -rf build",
                "echo x > build/CMakeCache.txt",
                "echo {} >graphify-out/graph.json",
            ]
        )

    def test_allowed_commands(self) -> None:
        allowed = [
            "cmake --build build",
            "ctest --test-dir build",
            "ninja -C build",
            "git rm src/old_file.cpp",
            "python3 format_terms.py",
            "ls -la build/",
            "git status",
            "ctest --output-on-failure",
            "git clean -n",
            'grep -rn "rm -rf" docs/',
            "echo 'a; rm -rf b' > notes.txt",
            "rm scratch.txt",
            "ls .git",
        ]
        for cmd in allowed:
            safe, reason = evaluate_command_safety(cmd)
            self.assertTrue(safe, f"Expected safe for: {cmd}, got: {reason}")

    def test_pre_tool_use_payload_contract(self) -> None:
        deny = handle_pre_tool_use(
            {"toolCall": {"name": "run_command", "args": {"CommandLine": "rm -rf build"}}}
        )
        self.assertEqual(deny["decision"], "deny")
        allow = handle_pre_tool_use({"toolCall": {"name": "view_file", "args": {}}})
        self.assertEqual(allow["decision"], "allow")


class TestGitDiscovery(unittest.TestCase):
    """PostToolUse modified-file discovery."""

    def test_parse_porcelain_z(self) -> None:
        output = " M src/a.cpp\0?? new dir/b.hpp\0R  ui/new.slint\0ui/old.slint\0"
        self.assertEqual(parse_porcelain_z(output), ["src/a.cpp", "new dir/b.hpp", "ui/new.slint"])


class TestCodeQualityAudits(unittest.TestCase):
    """PostToolUse source linters."""

    def audit(self, suffix: str, content: str, auditor) -> list[str]:
        path = _write_temp(suffix, content)
        try:
            return auditor(path)
        finally:
            path.unlink(missing_ok=True)

    def test_cpp_zero_nolint(self) -> None:
        for snippet in (
            "// NOLINTNEXTLINE(readability-function-cognitive-complexity)\nvoid f() {}\n",
            "int x = 0; /* NOLINT */\n",
            "// NOLINTEND\n",
        ):
            violations = self.audit(".cpp", snippet, audit_cpp_file)
            self.assertTrue(any("Zero-NOLINT" in v for v in violations), snippet)

    def test_nolint_in_string_is_not_a_comment_false_positive(self) -> None:
        violations = self.audit(".cpp", 'const char* s = "NOLINT";\n', audit_cpp_file)
        self.assertFalse(any("Zero-NOLINT" in v for v in violations))

    def test_cuda_files_audited(self) -> None:
        violations = self.audit(".cu", "// NOLINT\n__global__ void k() {}\n", audit_cpp_file)
        self.assertTrue(any("Zero-NOLINT" in v for v in violations))

    def test_cpp_cognitive_complexity(self) -> None:
        nested_code = """
        void complex_func(int a, int b) {
            if (a > 0) {
                for (int i = 0; i < a; ++i) {
                    while (b > 0) {
                        if (a && b) {
                            if (a || b) {
                                switch (i) {
                                    case 1: if (a) {} break;
                                    case 2: if (b) {} break;
                                    default: break;
                                }
                            }
                        }
                    }
                }
            }
        }
        """
        violations = self.audit(".cpp", nested_code, audit_cpp_file)
        self.assertTrue(any("Cognitive Complexity" in v for v in violations))

    def test_braceless_ifs_do_not_inflate_nesting(self) -> None:
        guards = "\n".join(f"    if (x == {i}) return {i};" for i in range(12))
        code = f"int f(int x) {{\n{guards}\n    for (int i = 0; i < x; ++i) {{ x += i; }}\n    return x;\n}}\n"
        (result,) = estimate_cognitive_complexity(code)
        self.assertEqual(result[1], "f")
        self.assertEqual(result[2], 13)  # 12 flat ifs + 1 flat for

    def test_lambda_and_init_braces_are_not_control_nesting(self) -> None:
        code = """
        void g(std::vector<int>& v) {
            auto cfg = Config{1, 2};
            std::ranges::for_each(v, [](int& e) { e += 1; });
            if (v.empty()) { return; }
        }
        """
        (result,) = estimate_cognitive_complexity(code)
        self.assertEqual(result[2], 1)

    def test_else_if_scoring(self) -> None:
        code = "void h(int a) { if (a) { } else if (a > 1) { } else { } }"
        (result,) = estimate_cognitive_complexity(code)
        self.assertEqual(result[2], 3)  # if +1, else-if +1, else +1

    def test_cmake_file_glob(self) -> None:
        violations = self.audit("CMakeLists.txt", "file(GLOB SOURCES src/*.cpp)\n", audit_cmake_file)
        self.assertTrue(any("CMake Anti-Pattern" in v for v in violations))

    def test_slint_export_block(self) -> None:
        violations = self.audit(".slint", "export {}\n", audit_slint_file)
        self.assertTrue(any("Slint Anti-Pattern" in v for v in violations))


if __name__ == "__main__":
    unittest.main()
