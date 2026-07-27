---
name: git-commit-standards
description: Enforces Conventional Commits format, validates scope tags against project modules, and generates consistent commit messages from staged diffs.
---

# Git Commit Standards Skill

This skill generates and validates commit messages following Conventional Commits for the Correlation project.

## 1. Commit Message Generation Protocol

When asked to generate a commit message:

1. **Inspect Staged Changes:** Run `git diff --cached --stat` to identify modified files.
2. **Classify Change Type:** Map the dominant change category to a Conventional Commits type:

| Type | Trigger Pattern |
| :--- | :--- |
| `feat` | New functions, classes, UI components, or CLI options |
| `fix` | Bug corrections, error handling fixes, crash prevention |
| `refactor` | Code restructuring without behavior change |
| `perf` | Algorithm optimization, cache alignment, SIMD changes |
| `test` | New or modified test cases |
| `docs` | Doxygen comments, README, SKILL.md changes |
| `build` | CMakeLists.txt, FetchContent, compiler flag changes |
| `chore` | Formatting, linting, dependency bumps |

3. **Select Scope:** Match modified file paths to project scopes:

| Scope | Path Pattern |
| :--- | :--- |
| `rdf`, `pdf`, `sq`, `pad`, `rings`, `msd`, `vacf`, `vdos`, `lef` | `src/calculators/<name>*` |
| `ui` | `ui/**/*.slint` |
| `slint` | `ui/material/**` |
| `bindings` | `src/bindings/**` |
| `cmake` | `**/CMakeLists.txt`, `cmake/**` |
| `cli` | `src/cli/**` |
| `io` | `src/io/**`, `src/parsers/**` |
| `core` | `src/core/**`, `include/**` |
| `tests` | `tests/**` |

4. **Format Output:**
```
<type>(<scope>): <imperative description, max 72 chars>

<body: explain WHY, not WHAT — 1-3 sentences max>
```

## 2. Validation Rules

- Subject line must be ≤ 72 characters.
- Subject must use imperative mood ("add", "fix", "remove" — not "added", "fixes", "removed").
- Body must explain rationale, not repeat the diff.
- Breaking changes must include `BREAKING CHANGE:` footer or `!` suffix.

## 3. Anti-Patterns

- ❌ `"Updated files"` — no description of what changed.
- ❌ `"Fixed bug"` — no scope or specificity.
- ❌ `"feat: Added new feature to do the thing"` — past tense, vague.
- ✅ `"feat(rdf): add Lorch modification function for g(r) smoothing"` — correct.
