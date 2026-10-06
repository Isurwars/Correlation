---
name: commit
description: Generate and validate Conventional Commits from git staged diffs. Use when crafting commits, preparing pull requests, or validating commit message format.
---

# Conventional Commit Standards

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
| `docs` | Doxygen comments, README, documentation changes |
| `build` | CMakeLists.txt, dependencies, compiler flag changes |
| `ci` | CI/CD GitHub workflows and deployment pipelines |
| `chore` | Formatting, linting, agent hooks, dependency bumps |

3. **Select Scope:** Match modified file paths to project scopes:

| Scope | Path Pattern |
| :--- | :--- |
| `rdf`, `pdf`, `sq`, `pad`, `rings`, `msd`, `vacf`, `vdos`, `lef` | `src/calculators/<name>*`, `include/calculators/<name>*` |
| `analysis`, `calculators`, `math`, `physics`, `mlip` | `src/analysis/**`, `src/math/**`, `src/physics/**`, `src/mlip/**` |
| `readers`, `writers`, `io` | `src/readers/**`, `src/writers/**`, `include/readers/**`, `include/writers/**` |
| `ui`, `slint`, `plotters` | `ui/**`, `src/app/**`, `src/plotters/**` |
| `bindings` | `src/bindings/**`, `python/**` |
| `cli` | `src/cli/**`, `include/cli/**` |
| `core`, `utils` | `src/core/**`, `src/utils/**`, `include/core/**` |
| `cmake` | `**/CMakeLists.txt`, `cmake/**`, `CMakePresets.json` |
| `tests` | `tests/**` |
| `docs` | `docs/**`, `*.md` |

4. **Format Output:**
```
<type>(<scope>): <imperative description, max 72 chars>

<body: explain WHY, not WHAT — 1-3 sentences max>
```

## 2. Validation Rules

- Subject line must be ≤ 72 characters.
- Subject must use imperative mood ("add", "fix", "remove" — not "added", "fixes", "removed").
- Body must explain rationale, not repeat the diff.
- Breaking changes must include `BREAKING CHANGE:` footer or `!` suffix: `feat(core)!: description`.

## 3. Anti-Patterns

- ❌ `"Updated files"` — no description of what changed.
- ❌ `"Fixed bug"` — no scope or specificity.
- ❌ `"feat: Added new feature to do the thing"` — past tense, vague.
- ✅ `"feat(rdf): add Lorch modification function for g(r) smoothing"` — correct.
