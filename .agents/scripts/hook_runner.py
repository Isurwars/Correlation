#!/usr/bin/env python3
"""Antigravity Lifecycle Hook Runner for Correlation Workspace.

Enforces workspace safety invariants, zero-NOLINT policies, cognitive
complexity gates, and repository integrity via deterministic hooks.
"""

from __future__ import annotations

import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any

WORKSPACE_ROOT = Path(__file__).resolve().parent.parent.parent
STATE_FILE = WORKSPACE_ROOT / ".agents" / ".hook_state.json"

PROTECTED_PATHS = {
    "build",
    "build/",
    "cmakecache.txt",
    "graphify-out",
    "graphify-out/",
    ".git",
    ".git/",
}

SAFE_BUILD_COMMANDS = {
    "cmake",
    "ctest",
    "ninja",
    "make",
    "meson",
}

CONTROL_FLOW_KEYWORDS = {"if", "for", "while", "switch", "catch"}
COGNITIVE_COMPLEXITY_THRESHOLD = 25


def _tokenize_command(cmd: str) -> list[str]:
    """Safely tokenizes shell command line handling quotes."""
    try:
        return shlex.split(cmd)
    except Exception:
        return cmd.split()


def _is_destructive_rm(tokens: list[str]) -> tuple[bool, str]:
    """Detects rm commands equipped with both recursive and force flags."""
    for i, token in enumerate(tokens):
        if token == "rm" and (i == 0 or tokens[i - 1] not in ("git",)):
            has_recursive = False
            has_force = False
            for arg in tokens[i + 1 :]:
                if arg in (";", "&&", "||", "|"):
                    break
                if arg == "--recursive":
                    has_recursive = True
                elif arg in ("--force", "-f"):
                    has_force = True
                elif arg.startswith("-") and not arg.startswith("--"):
                    flags = arg[1:]
                    if "r" in flags or "R" in flags:
                        has_recursive = True
                    if "f" in flags:
                        has_force = True

            if has_recursive and has_force:
                return (
                    True,
                    f"Blocked destructive command: '{' '.join(tokens[i:])}' matches recursive forced deletion (rm -rf).",
                )
    return False, ""


def _targets_protected_path(tokens: list[str]) -> tuple[bool, str]:
    """Detects commands attempting to wipe or delete protected build/graphify artifacts."""
    if not tokens:
        return False, ""

    primary_cmd = os.path.basename(tokens[0])
    if primary_cmd in SAFE_BUILD_COMMANDS:
        return False, ""

    destructive_cmds = {"rm", "rmdir", "unlink", "shred"}
    if primary_cmd in destructive_cmds:
        for arg in tokens[1:]:
            normalized = arg.strip("\"'").rstrip("/").lower()
            if normalized in ("build", "cmakecache.txt", "graphify-out", ".git"):
                return (
                    True,
                    f"Blocked deletion of protected workspace path: '{arg}'. Direct deletion is prohibited by AGENTS.md.",
                )

    full_str = " ".join(tokens).lower()
    for protected in ("cmakecache.txt", "graphify-out/graph.json", "graphify-out/graph_report.md"):
        if f"> {protected}" in full_str or f">{protected}" in full_str:
            return (
                True,
                f"Blocked file truncation/overwrite of protected file: '{protected}'.",
            )

    return False, ""


def evaluate_command_safety(command_line: str) -> tuple[bool, str]:
    """Evaluates whether a shell command is safe to execute."""
    if not command_line or not command_line.strip():
        return True, ""

    subcommands = re.split(r"[;&|]+", command_line)
    for subcmd in subcommands:
        tokens = _tokenize_command(subcmd.strip())
        if not tokens:
            continue

        is_dest_rm, rm_reason = _is_destructive_rm(tokens)
        if is_dest_rm:
            return False, rm_reason

        targets_prot, prot_reason = _targets_protected_path(tokens)
        if targets_prot:
            return False, prot_reason

    return True, ""


def audit_cpp_file(file_path: Path) -> list[str]:
    """Scans C++ source/header file for zero-NOLINT and cognitive complexity invariants."""
    violations: list[str] = []
    if not file_path.is_file():
        return violations

    try:
        content = file_path.read_text(encoding="utf-8", errors="replace")
    except Exception as exc:
        return [f"{file_path}: Unable to read file ({exc})"]

    lines = content.splitlines()
    rel_path = file_path.relative_to(WORKSPACE_ROOT) if file_path.is_relative_to(WORKSPACE_ROOT) else file_path

    # 1. Zero-NOLINT check
    nolint_pattern = re.compile(r"//\s*NOLINT(NEXTLINE|BEGIN)?(\([^)]*\))?")
    for idx, line in enumerate(lines, start=1):
        if nolint_pattern.search(line):
            violations.append(
                f"[{rel_path}:{idx}] -> Zero-NOLINT Violation -> Remove inline suppression comment."
            )

    # 2. Cognitive Complexity Heuristic
    current_func: str | None = None
    func_start_line = 0
    func_complexity = 0
    nesting_depth = 0
    brace_depth = 0
    func_base_depth = 0

    func_sig_re = re.compile(
        r"^(?:(?:auto|void|int|bool|double|float|[a-zA-Z_]\w*(?:::[a-zA-Z_]\w*)*(?:<[^>]+>)?)\s+)?([a-zA-Z_]\w*::[a-zA-Z_~]\w*|[a-zA-Z_]\w+)\s*\([^)]*\)\s*(?:const)?\s*(?:noexcept)?\s*\{?"
    )

    for idx, line in enumerate(lines, start=1):
        stripped = line.strip()
        if stripped.startswith("//") or stripped.startswith("/*") or stripped.startswith("*"):
            continue

        # Check for function start when at top level or namespace level
        if current_func is None and brace_depth <= 2 and not stripped.startswith(":"):
            match = func_sig_re.search(stripped)
            if match and not stripped.endswith(";"):
                current_func = match.group(1)
                func_start_line = idx
                func_complexity = 0
                nesting_depth = 0
                func_base_depth = brace_depth

        open_braces = line.count("{")
        close_braces = line.count("}")
        brace_depth += open_braces - close_braces

        if current_func is not None:
            # Control flow keywords add 1 + nesting depth
            for kw in CONTROL_FLOW_KEYWORDS:
                if re.search(rf"\b{kw}\b", stripped):
                    func_complexity += 1 + nesting_depth
                    nesting_depth += 1

            # Logical operators add 1
            func_complexity += len(re.findall(r"(&&|\|\|)", stripped))

            if close_braces > 0:
                nesting_depth = max(0, nesting_depth - close_braces)

            if brace_depth <= func_base_depth and (open_braces > 0 or close_braces > 0):
                if func_complexity > COGNITIVE_COMPLEXITY_THRESHOLD:
                    violations.append(
                        f"[{rel_path}:{func_start_line}] -> Cognitive Complexity ({func_complexity} > {COGNITIVE_COMPLEXITY_THRESHOLD}) in function '{current_func}' -> Decompose into subroutines."
                    )
                current_func = None
                brace_depth = max(0, brace_depth)

    return violations


def audit_slint_file(file_path: Path) -> list[str]:
    """Audits Slint UI file for deprecated patterns or anti-patterns."""
    violations: list[str] = []
    if not file_path.is_file():
        return violations

    try:
        content = file_path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return violations

    lines = content.splitlines()
    rel_path = file_path.relative_to(WORKSPACE_ROOT) if file_path.is_relative_to(WORKSPACE_ROOT) else file_path

    # Check for legacy export blocks outside component scopes
    for idx, line in enumerate(lines, start=1):
        stripped = line.strip()
        if re.match(r"^export\s*\{\s*\}", stripped):
            violations.append(
                f"[{rel_path}:{idx}] -> Slint Anti-Pattern -> Deprecated empty export block."
            )

    return violations


def audit_cmake_file(file_path: Path) -> list[str]:
    """Audits CMake build file for globbing anti-patterns."""
    violations: list[str] = []
    if not file_path.is_file():
        return violations

    try:
        content = file_path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return violations

    lines = content.splitlines()
    rel_path = file_path.relative_to(WORKSPACE_ROOT) if file_path.is_relative_to(WORKSPACE_ROOT) else file_path

    for idx, line in enumerate(lines, start=1):
        stripped = line.strip()
        if re.search(r"\bfile\s*\(\s*GLOB(_RECURSE)?\b", stripped, re.IGNORECASE):
            violations.append(
                f"[{rel_path}:{idx}] -> CMake Anti-Pattern -> file(GLOB) used. Explicitly list source files."
            )

    return violations


def get_modified_files() -> list[Path]:
    """Discovers modified and untracked files in the workspace via git."""
    try:
        res = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=str(WORKSPACE_ROOT),
            capture_output=True,
            text=True,
            timeout=5,
        )
        if res.returncode != 0:
            return []
        files = []
        for line in res.stdout.splitlines():
            if len(line) >= 4:
                file_rel = line[3:].strip()
                # Handle renamed files
                if " -> " in file_rel:
                    file_rel = file_rel.split(" -> ")[1]
                path = WORKSPACE_ROOT / file_rel
                if path.is_file():
                    files.append(path)
        return files
    except Exception:
        return []


def audit_workspace() -> list[str]:
    """Runs quality audits across all modified C++, Slint, and CMake files."""
    violations: list[str] = []
    modified_files = get_modified_files()

    for path in modified_files:
        ext = path.suffix.lower()
        name = path.name.lower()
        if ext in (".cpp", ".hpp", ".h", ".cc", ".cxx"):
            violations.extend(audit_cpp_file(path))
        elif ext == ".slint":
            violations.extend(audit_slint_file(path))
        elif name == "cmakelists.txt" or ext == ".cmake":
            violations.extend(audit_cmake_file(path))

    return violations


def handle_pre_tool_use(payload: dict[str, Any]) -> dict[str, Any]:
    """Validates tool execution before execution."""
    tool_call = payload.get("toolCall", {})
    name = tool_call.get("name", "")
    args = tool_call.get("args", {})

    if name == "run_command":
        cmd = args.get("CommandLine", "")
        safe, reason = evaluate_command_safety(cmd)
        if not safe:
            return {
                "decision": "deny",
                "reason": f"[SAFETY INTERCEPTOR] {reason}",
            }

    return {"decision": "allow"}


def handle_post_tool_use(payload: dict[str, Any]) -> dict[str, Any]:
    """Audits workspace state post tool execution and records violations."""
    violations = audit_workspace()

    state = {
        "violations": violations,
        "stepIdx": payload.get("stepIdx", 0),
    }

    try:
        STATE_FILE.parent.mkdir(parents=True, exist_ok=True)
        STATE_FILE.write_text(json.dumps(state, indent=2), encoding="utf-8")
    except Exception:
        pass

    return {}


def handle_pre_invocation(payload: dict[str, Any]) -> dict[str, Any]:
    """Injects notifications before model turn if violations were detected."""
    if not STATE_FILE.is_file():
        return {}

    try:
        state = json.loads(STATE_FILE.read_text(encoding="utf-8"))
        violations = state.get("violations", [])
        if violations:
            msg = "[HOOK AUDIT ALERT] Active quality gate violations detected:\n" + "\n".join(
                f"- {v}" for v in violations[:10]
            )
            return {
                "injectSteps": [
                    {
                        "ephemeralMessage": msg,
                    }
                ]
            }
    except Exception:
        pass

    return {}


def main() -> None:
    mode = sys.argv[1] if len(sys.argv) > 1 else "pre_tool_use"

    try:
        raw_input = sys.stdin.read()
        payload = json.loads(raw_input) if raw_input.strip() else {}
    except Exception:
        payload = {}

    if mode == "pre_tool_use":
        result = handle_pre_tool_use(payload)
    elif mode == "post_tool_use":
        result = handle_post_tool_use(payload)
    elif mode == "pre_invocation":
        result = handle_pre_invocation(payload)
    else:
        result = {}

    sys.stdout.write(json.dumps(result))
    sys.stdout.flush()


if __name__ == "__main__":
    main()
