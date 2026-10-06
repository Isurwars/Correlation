#!/usr/bin/env python3
"""Antigravity Lifecycle Hook Runner for the Correlation workspace.

Modes (argv[1]):
    pre_tool_use    -- deny destructive shell commands (PreToolUse, run_command).
    post_tool_use   -- audit modified C++/Slint/CMake files; cache results (PostToolUse).
    pre_invocation  -- inject an ephemeral alert when the violation set changes (PreInvocation).

The cognitive-complexity check is a fast *heuristic*. clang-tidy
(`readability-function-cognitive-complexity`) remains the authoritative gate.
"""

from __future__ import annotations

import bisect
import hashlib
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

COGNITIVE_COMPLEXITY_THRESHOLD = 25
MAX_ALERT_LINES = 10

CPP_EXTENSIONS = {".cpp", ".hpp", ".h", ".cc", ".cxx", ".cu", ".cuh"}

# --------------------------------------------------------------------------------------
# PreToolUse: command safety
# --------------------------------------------------------------------------------------

PROTECTED_DIRS = {"build", "graphify-out", ".git"}
PROTECTED_FILES = {"cmakecache.txt", "graphify-out/graph.json", "graphify-out/graph_report.md"}
SAFE_BUILD_COMMANDS = {"cmake", "ctest", "ninja", "make", "meson"}
DESTRUCTIVE_COMMANDS = {"rm", "rmdir", "unlink", "shred"}
WRAPPER_COMMANDS = {"sudo", "env", "nohup", "time", "xargs", "command", "exec", "nice", "doas"}
NESTED_SHELLS = {"sh", "bash", "zsh", "dash"}
OPERATOR_TOKENS = {";", "&&", "||", "|", "&", "|&", ";;", "(", ")", "\n"}
REDIRECT_TOKENS = {">", ">>", ">|", "&>", "&>>"}


def _tokenize_command(cmd: str) -> list[str]:
    """Tokenizes a shell line, emitting operators (`;`, `&&`, `|`, `>`) as separate tokens."""
    try:
        lexer = shlex.shlex(cmd, posix=True, punctuation_chars=True)
        lexer.whitespace_split = True
        return list(lexer)
    except ValueError:
        return cmd.split()


def _split_subcommands(tokens: list[str]) -> list[list[str]]:
    """Splits a token stream on shell control operators."""
    groups: list[list[str]] = [[]]
    for token in tokens:
        if token in OPERATOR_TOKENS:
            groups.append([])
        else:
            groups[-1].append(token)
    return [g for g in groups if g]


def _strip_wrappers(tokens: list[str]) -> list[str]:
    """Removes prefixes such as `sudo`, `env VAR=1`, `xargs -0` that wrap the real command."""
    idx = 0
    while idx < len(tokens):
        tok = tokens[idx]
        base = os.path.basename(tok)
        if base in WRAPPER_COMMANDS or tok.startswith("-") or ("=" in tok and idx > 0):
            idx += 1
            continue
        if "=" in tok and idx == 0:  # leading VAR=value assignment
            idx += 1
            continue
        break
    return tokens[idx:]


def _normalize_path(arg: str) -> str:
    """Lower-cased, workspace-relative, normalized path (`./build/` -> `build`)."""
    raw = os.path.expanduser(arg.strip("\"'"))
    path = Path(os.path.normpath(raw))
    if path.is_absolute():
        try:
            path = path.relative_to(WORKSPACE_ROOT)
        except ValueError:
            return str(path).lower()
    return str(path).lower()


def _is_protected(arg: str) -> bool:
    norm = _normalize_path(arg)
    return norm in PROTECTED_FILES or norm in PROTECTED_DIRS or norm.endswith("/cmakecache.txt")


def _rm_flags(args: list[str]) -> tuple[bool, bool]:
    recursive = force = False
    for arg in args:
        if arg == "--":
            break
        if arg in ("--recursive",):
            recursive = True
        elif arg in ("--force",):
            force = True
        elif arg.startswith("-") and not arg.startswith("--"):
            recursive |= "r" in arg or "R" in arg
            force |= "f" in arg
    return recursive, force


def _check_rm(cmd: str, args: list[str], tokens: list[str]) -> str:
    if cmd == "rm":
        recursive, force = _rm_flags(args)
        if recursive and force:
            return (
                f"Blocked destructive command: '{' '.join(tokens)}' matches recursive "
                "forced deletion (rm -rf)."
            )
    if cmd in DESTRUCTIVE_COMMANDS:
        for arg in args:
            if not arg.startswith("-") and _is_protected(arg):
                return (
                    f"Blocked deletion of protected workspace path: '{arg}'. "
                    "Direct deletion is prohibited by AGENTS.md."
                )
    return ""


def _check_special_deleters(cmd: str, args: list[str]) -> str:
    if cmd == "find" and "-delete" in args and any(_is_protected(a) for a in args):
        return "Blocked 'find -delete' targeting a protected workspace path."
    if cmd == "git" and args[:1] == ["clean"]:
        flags = "".join(a for a in args[1:] if a.startswith("-") and not a.startswith("--"))
        if "x" in flags or "X" in flags:
            return "Blocked 'git clean -x/-X': it wipes ignored build/ and graphify-out/ trees."
    if cmd == "cmake" and args[:1] == ["-E"] and args[1:2] and args[1] in (
        "rm",
        "remove",
        "remove_directory",
    ):
        if any(_is_protected(a) for a in args[2:]):
            return "Blocked 'cmake -E rm' targeting a protected workspace path."
    return ""


def _check_redirects(tokens: list[str]) -> str:
    for idx, tok in enumerate(tokens[:-1]):
        if tok in REDIRECT_TOKENS:
            target = _normalize_path(tokens[idx + 1])
            if any(target == f or target.endswith("/" + f) for f in PROTECTED_FILES):
                return f"Blocked file truncation/overwrite of protected file: '{tokens[idx + 1]}'."
    return ""


def _evaluate_subcommand(tokens: list[str], depth: int) -> str:
    redirect_reason = _check_redirects(tokens)
    if redirect_reason:
        return redirect_reason

    tokens = _strip_wrappers(tokens)
    if not tokens:
        return ""
    cmd = os.path.basename(tokens[0])
    args = tokens[1:]

    if cmd in NESTED_SHELLS and "-c" in args and depth < 3:
        idx = args.index("-c")
        if idx + 1 < len(args):
            safe, reason = _evaluate(args[idx + 1], depth + 1)
            return "" if safe else reason

    special = _check_special_deleters(cmd, args)
    if special:
        return special
    if cmd in SAFE_BUILD_COMMANDS:
        return ""
    return _check_rm(cmd, args, tokens)


def _evaluate(command_line: str, depth: int) -> tuple[bool, str]:
    if not command_line or not command_line.strip():
        return True, ""
    for sub in _split_subcommands(_tokenize_command(command_line)):
        reason = _evaluate_subcommand(sub, depth)
        if reason:
            return False, reason
    return True, ""


def evaluate_command_safety(command_line: str) -> tuple[bool, str]:
    """Evaluates whether a shell command is safe to execute."""
    return _evaluate(command_line, 0)


# --------------------------------------------------------------------------------------
# PostToolUse: source audits
# --------------------------------------------------------------------------------------

_NOLINT_RE = re.compile(r"(//|/\*).*?\bNOLINT(NEXTLINE|BEGIN|END)?\b")
_STRIP_RE = re.compile(
    r"//[^\n]*|/\*.*?\*/|\"(?:\\.|[^\"\\\n])*\"|'(?:\\.|[^'\\\n])*'|^[ \t]*#[^\n]*",
    re.DOTALL | re.MULTILINE,
)
_TOKEN_RE = re.compile(r"\b(?:if|else|for|while|switch|catch)\b|&&|\|\||[?{}();]")
_FUNC_HEADER_RE = re.compile(
    r"\)\s*(?:const|noexcept(?:\([^)]*\))?|override|final|mutable|&|&&|\s|->\s*[\w:<>,\s*&]+)*$"
)
_NON_FUNC_HEADER_RE = re.compile(r"^\s*(?:namespace|class|struct|enum|union|extern)\b|=\s*$")
_FUNC_NAME_RE = re.compile(r"([A-Za-z_~][\w:~]*)\s*\(")
_CONTROL_KEYWORDS = {"if", "for", "while", "switch", "catch"}


def _blank(match: re.Match[str]) -> str:
    """Replaces a comment/literal with spaces while preserving line breaks."""
    return re.sub(r"[^\n]", " ", match.group(0))


def _rel(path: Path) -> Path:
    return path.relative_to(WORKSPACE_ROOT) if path.is_relative_to(WORKSPACE_ROOT) else path


def _read(path: Path) -> str | None:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None


class _ComplexityScanner:
    """Brace-stack approximation of clang-tidy cognitive complexity.

    Nesting only grows for braces opened by control statements, so brace-less
    `if (x) return;` and lambda/initializer braces no longer skew the score.
    """

    def __init__(self, code: str) -> None:
        self.code = code
        self.line_starts = [0] + [m.end() for m in re.finditer(r"\n", code)]
        self.stack: list[str] = []  # 'func' | 'ctrl' | 'block'
        self.paren_depth = 0
        self.pending_ctrl = False
        self.after_else = False
        self.seg_start = 0
        self.func: dict[str, Any] | None = None
        self.last_logical = ""
        self.results: list[tuple[int, str, int]] = []

    def _line(self, pos: int) -> int:
        return bisect.bisect_right(self.line_starts, pos)

    def _nesting(self) -> int:
        if self.func is None:
            return 0
        return self.stack[self.func["base"] + 1 :].count("ctrl")

    def _open_brace(self, pos: int) -> None:
        header = self.code[self.seg_start : pos]
        if self.func is None:
            if _FUNC_HEADER_RE.search(header) and not _NON_FUNC_HEADER_RE.search(header):
                name_match = _FUNC_NAME_RE.search(header)
                self.func = {
                    "base": len(self.stack),
                    "name": name_match.group(1) if name_match else "<anonymous>",
                    "line": self._line(pos),
                    "score": 0,
                }
                self.stack.append("func")
            else:
                self.stack.append("block")
        else:
            self.stack.append("ctrl" if self.pending_ctrl else "block")
        self.pending_ctrl = False
        self.after_else = False

    def _close_brace(self) -> None:
        if not self.stack:
            return
        kind = self.stack.pop()
        if kind == "func" and self.func is not None:
            self.results.append((self.func["line"], self.func["name"], self.func["score"]))
            self.func = None

    def _keyword(self, tok: str) -> None:
        if self.func is None:
            return
        if tok == "else":
            self.func["score"] += 1
            self.after_else = True
            self.pending_ctrl = True
            return
        if tok == "if" and self.after_else:
            self.after_else = False  # `else if` already scored
        else:
            self.func["score"] += 1 + self._nesting()
        self.pending_ctrl = True

    def _paren_or_semicolon(self, tok: str) -> bool:
        """Handles `(`, `)` and `;`; returns True when the segment boundary must reset."""
        if tok == "(":
            self.paren_depth += 1
            return False
        if tok == ")":
            self.paren_depth = max(0, self.paren_depth - 1)
            return False
        if self.paren_depth == 0:  # tok == ';' outside a for-header
            self.pending_ctrl = False
            self.after_else = False
            return True
        return False

    def _logical(self, tok: str) -> None:
        """Scores `?` with nesting; scores a run of identical `&&`/`||` operators once."""
        if self.func is None:
            return
        if tok == "?":
            self.func["score"] += 1 + self._nesting()
            self.last_logical = ""
        elif tok != self.last_logical:
            self.func["score"] += 1
            self.last_logical = tok

    def _dispatch(self, tok: str, pos: int) -> bool:
        """Processes one token; returns True when the header segment must reset."""
        if tok in ("&&", "||", "?"):
            self._logical(tok)
            return False
        if tok not in "()":
            self.last_logical = ""
        if tok in _CONTROL_KEYWORDS or tok == "else":
            self._keyword(tok)
            return False
        if tok in "();":
            return self._paren_or_semicolon(tok)
        if tok == "{":
            self._open_brace(pos)
        else:
            self._close_brace()
        return True

    def scan(self) -> list[tuple[int, str, int]]:
        for match in _TOKEN_RE.finditer(self.code):
            if self._dispatch(match.group(0), match.start()):
                self.seg_start = match.end()
        return self.results


def estimate_cognitive_complexity(content: str) -> list[tuple[int, str, int]]:
    """Returns (line, function, approx_score) for every detected function body."""
    return _ComplexityScanner(_STRIP_RE.sub(_blank, content)).scan()


def audit_cpp_file(file_path: Path) -> list[str]:
    """Scans a C++/CUDA file for NOLINT suppressions and likely complexity violations."""
    content = _read(file_path) if file_path.is_file() else None
    if content is None:
        return []
    rel_path = _rel(file_path)
    violations = [
        f"[{rel_path}:{idx}] -> Zero-NOLINT Violation -> Remove inline suppression comment."
        for idx, line in enumerate(content.splitlines(), start=1)
        if _NOLINT_RE.search(line)
    ]
    for line, name, score in estimate_cognitive_complexity(content):
        if score > COGNITIVE_COMPLEXITY_THRESHOLD:
            violations.append(
                f"[{rel_path}:{line}] -> Cognitive Complexity ~{score} (heuristic) > "
                f"{COGNITIVE_COMPLEXITY_THRESHOLD} in '{name}' -> Confirm with "
                f"`clang-tidy -p build {rel_path}`; decompose if confirmed."
            )
    return violations


def audit_slint_file(file_path: Path) -> list[str]:
    """Audits a Slint UI file for deprecated patterns."""
    content = _read(file_path) if file_path.is_file() else None
    if content is None:
        return []
    rel_path = _rel(file_path)
    return [
        f"[{rel_path}:{idx}] -> Slint Anti-Pattern -> Deprecated empty export block."
        for idx, line in enumerate(content.splitlines(), start=1)
        if re.match(r"^export\s*\{\s*\}", line.strip())
    ]


def audit_cmake_file(file_path: Path) -> list[str]:
    """Audits a CMake file for source-globbing anti-patterns."""
    content = _read(file_path) if file_path.is_file() else None
    if content is None:
        return []
    rel_path = _rel(file_path)
    return [
        f"[{rel_path}:{idx}] -> CMake Anti-Pattern -> file(GLOB) used. Explicitly list source files."
        for idx, line in enumerate(content.splitlines(), start=1)
        if re.search(r"\bfile\s*\(\s*GLOB(_RECURSE)?\b", line, re.IGNORECASE)
    ]


def parse_porcelain_z(output: str) -> list[str]:
    """Parses `git status --porcelain=v1 -z` output into workspace-relative paths."""
    entries = output.split("\0")
    paths: list[str] = []
    idx = 0
    while idx < len(entries):
        entry = entries[idx]
        idx += 1
        if len(entry) < 4:
            continue
        status, rel = entry[:2], entry[3:]
        if "R" in status or "C" in status:
            idx += 1  # skip the original path of a rename/copy
        paths.append(rel)
    return paths


def get_modified_files() -> list[Path]:
    """Discovers modified and untracked files (including files inside new directories)."""
    try:
        res = subprocess.run(
            ["git", "status", "--porcelain=v1", "-z", "-uall"],
            cwd=str(WORKSPACE_ROOT),
            capture_output=True,
            text=True,
            timeout=4,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return []
    if res.returncode != 0:
        return []
    return [p for p in (WORKSPACE_ROOT / rel for rel in parse_porcelain_z(res.stdout)) if p.is_file()]


def _is_auditable(path: Path) -> bool:
    ext = path.suffix.lower()
    return (
        ext in CPP_EXTENSIONS
        or ext in (".slint", ".cmake")
        or path.name.lower() == "cmakelists.txt"
    )


def _audit_path(path: Path) -> list[str]:
    ext = path.suffix.lower()
    if ext in CPP_EXTENSIONS:
        return audit_cpp_file(path)
    if ext == ".slint":
        return audit_slint_file(path)
    if path.name.lower() == "cmakelists.txt" or ext == ".cmake":
        return audit_cmake_file(path)
    return []


def _fingerprint(paths: list[Path]) -> str:
    digest = hashlib.sha1()
    for path in sorted(paths):
        try:
            stat = path.stat()
        except OSError:
            continue
        digest.update(f"{path}:{stat.st_mtime_ns}:{stat.st_size};".encode())
    return digest.hexdigest()


def _digest(violations: list[str]) -> str:
    return hashlib.sha1("\n".join(violations).encode()).hexdigest()


def _load_state() -> dict[str, Any]:
    try:
        return json.loads(STATE_FILE.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def _save_state(state: dict[str, Any]) -> None:
    try:
        STATE_FILE.parent.mkdir(parents=True, exist_ok=True)
        STATE_FILE.write_text(json.dumps(state, indent=2), encoding="utf-8")
    except OSError:
        pass


def audit_workspace(modified: list[Path] | None = None) -> list[str]:
    """Runs quality audits across all modified C++, Slint, and CMake files."""
    files = get_modified_files() if modified is None else modified
    violations: list[str] = []
    for path in files:
        violations.extend(_audit_path(path))
    return violations


# --------------------------------------------------------------------------------------
# Hook handlers
# --------------------------------------------------------------------------------------


def handle_pre_tool_use(payload: dict[str, Any]) -> dict[str, Any]:
    """Denies destructive shell commands before execution."""
    tool_call = payload.get("toolCall", {})
    if tool_call.get("name", "") == "run_command":
        safe, reason = evaluate_command_safety(tool_call.get("args", {}).get("CommandLine", ""))
        if not safe:
            return {"decision": "deny", "reason": f"[SAFETY INTERCEPTOR] {reason}"}
    return {"decision": "allow"}


def handle_post_tool_use(payload: dict[str, Any]) -> dict[str, Any]:
    """Re-audits only when the set/mtimes of modified files changed (cheap no-op otherwise)."""
    modified = [p for p in get_modified_files() if _is_auditable(p)]
    fingerprint = _fingerprint(modified)
    state = _load_state()
    if state.get("fingerprint") == fingerprint:
        return {}
    state.update(
        {
            "fingerprint": fingerprint,
            "violations": audit_workspace(modified),
            "stepIdx": payload.get("stepIdx", 0),
        }
    )
    _save_state(state)
    return {}


def handle_pre_invocation(_payload: dict[str, Any]) -> dict[str, Any]:
    """Injects an alert only when the violation set differs from the last one reported."""
    state = _load_state()
    violations = state.get("violations", [])
    digest = _digest(violations)
    if not violations or state.get("notifiedDigest") == digest:
        return {}
    state["notifiedDigest"] = digest
    _save_state(state)
    lines = "\n".join(f"- {v}" for v in violations[:MAX_ALERT_LINES])
    extra = len(violations) - MAX_ALERT_LINES
    if extra > 0:
        lines += f"\n- ... {extra} more"
    return {
        "injectSteps": [
            {"ephemeralMessage": f"[HOOK AUDIT ALERT] Quality gate findings:\n{lines}"}
        ]
    }


HANDLERS = {
    "pre_tool_use": handle_pre_tool_use,
    "post_tool_use": handle_post_tool_use,
    "pre_invocation": handle_pre_invocation,
}


def main() -> None:
    mode = sys.argv[1] if len(sys.argv) > 1 else "pre_tool_use"
    try:
        raw_input = sys.stdin.read()
        payload = json.loads(raw_input) if raw_input.strip() else {}
    except ValueError:
        payload = {}
    handler = HANDLERS.get(mode)
    result = handler(payload) if handler else {}
    sys.stdout.write(json.dumps(result))
    sys.stdout.flush()


if __name__ == "__main__":
    main()
