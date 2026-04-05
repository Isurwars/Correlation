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
# (Assuming you have a correlation.iss.in template where you use @VERSION@ as a placeholder)
if [ -f "packaging/correlation.iss.in" ]; then
    sed "s/@VERSION@/$VERSION/g" packaging/correlation.iss.in > packaging/correlation.iss
    echo "Generated correlation.iss from template."
else
    # Fallback: Just replace the version string directly in the .iss file
    sed -i -E "s/AppVersion=[0-9]+\.[0-9]+\.[0-9]+/AppVersion=$VERSION/g" packaging/correlation.iss
fi

# 3. Update README.md badge
# (Assuming a shields.io badge like: ![Version](https://img.shields.io/badge/version-v2.4.0-blue))
sed -i -E "s/badge\/version-v[0-9]+\.[0-9]+\.[0-9]+-[A-Za-z0-9]+/badge\/version-v$VERSION-blue/g" README.md

# 4. Stage changes
#git add CMakeLists.txt packaging/correlation.iss README.md

# 5. Show diff and wait for confirmation
#echo "---------------------------------------------------"
#git diff --cached
#echo "---------------------------------------------------"
#read -p "Does the diff look good? Press [Enter] to commit and push, or [Ctrl+C] to abort..."

# 6. Commit, tag, and push
#git commit -m "Release v$VERSION"
git tag "v$VERSION"
#git push && 
git push --tags

echo "✅ Version $VERSION pushed successfully! Please trigger the release action manually from the GitHub Actions tab."