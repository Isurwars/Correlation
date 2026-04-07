#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status

# Check if a version number was provided
if [ -z "$1" ]; then
  echo "Usage: ./release.sh <version> (e.g., 2.5.0)"
  exit 1
fi

VERSION=$1
echo "🚀 Prepping Correlation Release v$VERSION..."

# 1. Update CMakeLists.txt (Matches 'project(Correlation VERSION x.y.z)')
sed -i -E "s/project\(Correlation VERSION [0-9]+\.[0-9]+\.[0-9]+\)/project(Correlation VERSION $VERSION)/g" CMakeLists.txt

# 2. Re-generate packaging/correlation.iss from a template
# (Assuming you have a correlation.iss.in template where you use @PROJECT_VERSION@ as a placeholder)
if [ -f "packaging/correlation.iss.in" ]; then
    sed "s/@PROJECT_VERSION@/$VERSION/g" packaging/correlation.iss.in > packaging/correlation.iss
    echo "Generated correlation.iss from template."
else
    # Fallback: Just replace the version string directly in the .iss file
    sed -i -E "s/AppVersion=[0-9]+\.[0-9]+\.[0-9]+/AppVersion=$VERSION/g" packaging/correlation.iss
fi

# 3. Update README.md badge
# (Assuming a shields.io badge like: ![Version](https://img.shields.io/badge/version-2.4.0-green))
sed -i -E "s/badge\/version-v?[0-9]+\.[0-9]+\.[0-9]+-[A-Za-z0-9]+/badge\/version-$VERSION-green/g" README.md

# 4. Update Arch Linux packaging files
if [ -f "packaging/PKGBUILD" ]; then
    sed -i -E "s/^pkgver=[0-9]+\.[0-9]+\.[0-9]+/pkgver=$VERSION/g" packaging/PKGBUILD
    echo "Updated PKGBUILD version."
fi

if [ -f "packaging/.SRCINFO" ]; then
    sed -i -E "s/pkgver = [0-9]+\.[0-9]+\.[0-9]+/pkgver = $VERSION/g" packaging/.SRCINFO
    echo "Updated .SRCINFO version."
fi

# 5. Update Linux desktop file
if [ -f "packaging/correlation.desktop" ]; then
    sed -i -E "s/^Version=[0-9]+\.[0-9]+\.[0-9]+/Version=$VERSION/g" packaging/correlation.desktop
    echo "Updated correlation.desktop version."
fi

echo "✅ Version updated to $VERSION across all relevant files!"