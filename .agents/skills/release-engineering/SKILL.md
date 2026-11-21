---
name: release-engineering
description: End-to-end release workflow for Correlation: 8-file version bump, docs/CHANGELOG generation from Conventional Commits, CMake release verification, annotated git tagging, and Linux/Windows/macOS packaging. Use when preparing releases or updating changelogs.
---

# Release Engineering & Packaging

## 1. Exact 8-File Version Bump Matrix
Every release MUST update all 8 tracked version locations in lockstep:
1. `CMakeLists.txt` (`project(Correlation VERSION X.Y.Z LANGUAGES CXX C)`)
2. `ui/AppWindow.slint` (`in-out property <string> app_version: "vX.Y.Z";`)
3. `packaging/PKGBUILD` (`pkgver=X.Y.Z`)
4. `python/correlation/__init__.py` (`__version__ = "X.Y.Z"`)
5. `python/matcorr/__init__.py` (`__version__ = "X.Y.Z"`)
6. `README.md` (version shield badge)
7. `docs/CONTRIBUTING.md` (version shield badge)
8. `docs/CHANGELOG` *(NOT `CHANGELOG.md`!)* — prepend `## [X.Y.Z] - YYYY-MM-DD` grouped by Features, Performance Improvements, Bug Fixes, Code Refactoring, and Testing with short commit hashes.

## 2. Changelog Extraction
```bash
LAST_TAG=$(git describe --tags --abbrev=0)
git log "${LAST_TAG}..HEAD" --oneline
```

## 3. Pre-Release Verification Gate
1. Confirm zero stale version strings outside `docs/CHANGELOG`:
   ```bash
   git grep "<PREV_VERSION>"
   ```
2. Reconfigure release preset to regenerate `build-release/packaging/*` templates (`PKGBUILD`, `correlation.iss`, `resource.rc`, `meta.yaml`):
   ```bash
   cmake --preset release
   ```
3. Commit & create annotated tag:
   ```bash
   git add CMakeLists.txt README.md docs/CHANGELOG docs/CONTRIBUTING.md packaging/PKGBUILD python/correlation/__init__.py python/matcorr/__init__.py ui/AppWindow.slint
   git commit -m "chore(release): bump version to X.Y.Z and update changelog"
   git tag -a vX.Y.Z -m "Release vX.Y.Z"
   ```
4. Push from an authenticated terminal (`git push origin main vX.Y.Z`) to trigger `.github/workflows/release.yml` (Linux AppImage/DEB/RPM/TGZ + AUR, Windows NSIS `.exe`, macOS `.dmg`).

## 4. Local Packaging (CPack)
```bash
cmake --build --preset release
cpack --config build-release/CPackConfig.cmake -G "DEB;RPM;TGZ"
```
