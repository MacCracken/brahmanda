#!/usr/bin/env bash
set -euo pipefail

# version-bump.sh — Update version across all files
#
# Usage: ./scripts/version-bump.sh 1.1.0

if [ $# -ne 1 ]; then
    echo "Usage: $0 <version>"
    exit 1
fi

VERSION="$1"

echo "$VERSION" > VERSION
sed -i "s/^version = \".*\"/version = \"$VERSION\"/" Cargo.toml
cargo generate-lockfile 2>/dev/null || true

echo "Version bumped to $VERSION"
echo "  - VERSION file"
echo "  - Cargo.toml"
echo "  - Cargo.lock"
