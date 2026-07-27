---
name: changelog-generator
description: Generates structured CHANGELOG.md files and release notes from Conventional Commits history and git tags.
---

# Changelog Generator Skill

This skill parses git commit history using Conventional Commits tags to produce user-facing release notes and update `CHANGELOG.md`.

## 1. Execution Workflow

### Step 1: Detect Last Tag
Identify the latest semantic release tag:
```bash
git describe --tags --abbrev=0 2>/dev/null || git rev-list --max-parents=0 HEAD
```

### Step 2: Extract Commit Logs
Extract all commits since the last tag up to `HEAD`:
```bash
git log <last-tag>..HEAD --pretty=format:"%h %s"
```

### Step 3: Categorize Changes
Group commit entries into standard changelog sections:

| Category | Commit Types | Section Header |
| :--- | :--- | :--- |
| **Features** | `feat` | `### 🚀 Features` |
| **Bug Fixes** | `fix` | `### 🐛 Bug Fixes` |
| **Performance** | `perf` | `### ⚡ Performance Improvements` |
| **Refactoring** | `refactor` | `### ♻️ Code Refactoring` |
| **Build & Dependencies** | `build`, `ci` | `### 🛠️ Build & CI` |
| **Documentation** | `docs` | `### 📚 Documentation` |
| **Breaking Changes** | `BREAKING CHANGE:`, `type!:` | `### ⚠️ Breaking Changes` |

---

## 2. Formatting Specification

```markdown
## [vX.Y.Z] - YYYY-MM-DD

### 🚀 Features
- **rdf:** add Lorch windowing function for Fourier transform smoothing ([a1b2c3d](file:///commit/a1b2c3d))
- **ui:** implement responsive DropDownMenu scrollbar policy ([e5f6g7h](file:///commit/e5f6g7h))

### 🐛 Bug Fixes
- **ui:** fix absolute position coordinate calculation in AppWindow ([1a2b3c4](file:///commit/1a2b3c4))

### ⚡ Performance Improvements
- **openmp:** align thread-local histogram accumulators to 64 bytes ([5d6e7f8](file:///commit/5d6e7f8))
```

---

## 3. Automation Protocol

1. Read latest commits.
2. Format markdown release notes.
3. Prepend the new release section below the header of `CHANGELOG.md`.
4. Validate that no breaking changes are missing from the summary.
