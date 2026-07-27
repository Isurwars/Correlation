---
name: release-packaging
description: Guides creation of distributable Linux packages (AppImage, DEB, RPM) and release artifacts for the Correlation desktop application.
---

# Release Packaging Skill

This skill automates the creation of distributable release artifacts for the Correlation application on Linux platforms.

## 1. Pre-Packaging Checklist

Before creating a release package:

1. **Version Bump:** Ensure `CMakeLists.txt` contains the correct version: `project(Correlation VERSION x.y.z)`.
2. **Release Build:** Compile with optimizations and stripped debug symbols:
   ```bash
   cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr
   cmake --build build -j$(nproc)
   ```
3. **Test Verification:** All test targets must pass:
   ```bash
   cmake --build build --target correlation_unit_tests correlation_functional_tests correlation_gui_tests -j$(nproc)
   ./build/tests/correlation_unit_tests
   ./build/tests/correlation_functional_tests
   ./build/tests/correlation_gui_tests
   ```

---

## 2. AppImage Creation

AppImage is the preferred format for portable Linux distribution.

### Steps
1. **Install `linuxdeploy`:**
   ```bash
   wget -q https://github.com/linuxdeploy/linuxdeploy/releases/download/continuous/linuxdeploy-x86_64.AppImage
   chmod +x linuxdeploy-x86_64.AppImage
   ```

2. **Install to AppDir:**
   ```bash
   DESTDIR=AppDir cmake --install build
   ```

3. **Create Desktop Entry:**
   ```ini
   [Desktop Entry]
   Name=Correlation
   Exec=correlation
   Icon=correlation
   Type=Application
   Categories=Science;Physics;Education;
   Comment=Liquid and Amorphous Solid Structural Analysis
   ```

4. **Bundle with linuxdeploy:**
   ```bash
   ./linuxdeploy-x86_64.AppImage \
     --appdir AppDir \
     --desktop-file AppDir/usr/share/applications/correlation.desktop \
     --icon-file AppDir/usr/share/icons/correlation.png \
     --output appimage
   ```

---

## 3. DEB Package Creation

For Debian/Ubuntu distributions.

### Using CPack (CMake-Integrated)
Add to `CMakeLists.txt`:
```cmake
set(CPACK_GENERATOR "DEB")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Isaías Rodríguez <isurwars@gmail.com>")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6 (>= 2.31), libstdc++6 (>= 11)")
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "Structural analysis tool for liquids and amorphous solids")
include(CPack)
```

Build the package:
```bash
cd build && cpack -G DEB
```

---

## 4. RPM Package Creation

For Fedora/RHEL/CentOS distributions.

### Using CPack
```cmake
set(CPACK_GENERATOR "RPM")
set(CPACK_RPM_PACKAGE_LICENSE "AGPL-3.0-only")
set(CPACK_RPM_PACKAGE_GROUP "Applications/Science")
set(CPACK_RPM_PACKAGE_REQUIRES "glibc >= 2.31, libstdc++ >= 11")
include(CPack)
```

Build the package:
```bash
cd build && cpack -G RPM
```

---

## 5. Release Artifact Naming Convention

```
correlation-<version>-<arch>.<format>
```

Examples:
- `correlation-2.1.0-x86_64.AppImage`
- `correlation_2.1.0_amd64.deb`
- `correlation-2.1.0-1.x86_64.rpm`

## 6. Post-Packaging Verification

1. **AppImage:** Run on a clean system or container to verify all shared libraries are bundled.
2. **DEB:** Install with `sudo dpkg -i <package>.deb` and verify with `which correlation`.
3. **RPM:** Install with `sudo rpm -i <package>.rpm` and verify execution.
