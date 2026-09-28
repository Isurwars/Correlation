#!/usr/bin/env python3
"""Automated unit test suite for Correlation lifecycle hooks."""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path
import unittest

from hook_runner import (
    evaluate_command_safety,
    audit_cpp_file,
    audit_slint_file,
    audit_cmake_file,
)


class TestHookSafety(unittest.TestCase):
    """Verifies command evaluation logic in PreToolUse."""

    def test_destructive_rm_patterns(self):
        destructive = [
            "rm -rf build",
            "rm -fr build/",
            "rm -r -f target/",
            "rm -f -r some_dir",
            "rm -Rf /tmp/foo",
            "rm --recursive --force .",
            "rm --force --recursive node_modules",
            "rm -rf \"path with spaces\"",
        ]
        for cmd in destructive:
            safe, reason = evaluate_command_safety(cmd)
            self.assertFalse(safe, f"Expected unsafe for: {cmd}")
            self.assertIn("recursive forced deletion", reason)

    def test_protected_paths_deletion(self):
        protected = [
            "rm -r build",
            "rmdir build",
            "unlink build",
            "rm CMakeCache.txt",
            "rm -r graphify-out",
            "rmdir graphify-out/",
        ]
        for cmd in protected:
            safe, reason = evaluate_command_safety(cmd)
            self.assertFalse(safe, f"Expected unsafe for: {cmd}")

    def test_allowed_commands(self):
        allowed = [
            "cmake --build build",
            "ctest --test-dir build",
            "ninja -C build",
            "git rm src/old_file.cpp",
            "python3 format_terms.py",
            "ls -la build/",
            "git status",
            "ctest --output-on-failure",
        ]
        for cmd in allowed:
            safe, reason = evaluate_command_safety(cmd)
            self.assertTrue(safe, f"Expected safe for: {cmd}, got: {reason}")


class TestCodeQualityAudits(unittest.TestCase):
    """Verifies code quality linters in PostToolUse."""

    def test_cpp_zero_nolint(self):
        with tempfile.NamedTemporaryFile(suffix=".cpp", mode="w", delete=False) as f:
            f.write("// NOLINTNEXTLINE(readability-function-cognitive-complexity)\nvoid test() {}\n")
            path = Path(f.name)

        try:
            violations = audit_cpp_file(path)
            self.assertTrue(any("Zero-NOLINT" in v for v in violations))
        finally:
            path.unlink(missing_ok=True)

    def test_cpp_cognitive_complexity(self):
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
        with tempfile.NamedTemporaryFile(suffix=".cpp", mode="w", delete=False) as f:
            f.write(nested_code)
            path = Path(f.name)

        try:
            violations = audit_cpp_file(path)
            self.assertTrue(any("Cognitive Complexity" in v for v in violations))
        finally:
            path.unlink(missing_ok=True)

    def test_cmake_file_glob(self):
        with tempfile.NamedTemporaryFile(suffix="CMakeLists.txt", mode="w", delete=False) as f:
            f.write("file(GLOB SOURCES src/*.cpp)\n")
            path = Path(f.name)

        try:
            violations = audit_cmake_file(path)
            self.assertTrue(any("CMake Anti-Pattern" in v for v in violations))
        finally:
            path.unlink(missing_ok=True)

    def test_slint_export_block(self):
        with tempfile.NamedTemporaryFile(suffix=".slint", mode="w", delete=False) as f:
            f.write("export {}\n")
            path = Path(f.name)

        try:
            violations = audit_slint_file(path)
            self.assertTrue(any("Slint Anti-Pattern" in v for v in violations))
        finally:
            path.unlink(missing_ok=True)


if __name__ == "__main__":
    unittest.main()
