#!/bin/bash
# ==========================================================================
# Nabbit macOS Installer Build Script
# ==========================================================================
# Creates a self-contained Nabbit.app bundle and packages it into a .dmg.
#
# Prerequisites:
#   - macOS 12+ (Monterey or later)
#   - Python 3.12+ (brew install python@3.12)
#   - Xcode Command Line Tools (xcode-select --install)
#   - ~5 GB free disk space
#
# Usage:
#   bash installer/build_macos.sh
#
# Environment variables:
#   NABBIT_LITE=1        — Skip ablang/PyTorch (smaller ~200MB bundle)
#   CODESIGN_IDENTITY=X  — Sign with a Developer ID instead of ad-hoc
#   PYTHON=python3.12    — Override the Python interpreter path
#   VERSION=1.0.0        — Override the version string
# ==========================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="$PROJECT_DIR/build"
DIST_DIR="$PROJECT_DIR/dist"
VENV_DIR="$BUILD_DIR/venv_build"
PYTHON="${PYTHON:-python3}"
VERSION="${VERSION:-1.0.0}"
APP_NAME="Nabbit"
DMG_NAME="${APP_NAME}-${VERSION}-macOS.dmg"
LITE="${NABBIT_LITE:-0}"
CODESIGN_IDENTITY="${CODESIGN_IDENTITY:--}"

echo "============================================================"
echo "  Nabbit macOS Installer Build"
echo "============================================================"
echo "  Project:  $PROJECT_DIR"
echo "  Python:   $PYTHON"
echo "  Version:  $VERSION"
echo "  Lite:     $LITE"
echo "  Signing:  $CODESIGN_IDENTITY"
echo "============================================================"
echo ""

# ---------------------------------------------------------------------------
# Step 1: Check prerequisites
# ---------------------------------------------------------------------------
echo ">>> Step 1: Checking prerequisites..."

# Python 3.12+
if ! command -v "$PYTHON" &>/dev/null; then
    echo "ERROR: $PYTHON not found. Install Python 3.12+:"
    echo "  brew install python@3.12"
    exit 1
fi

PY_VERSION=$("$PYTHON" -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
PY_MAJOR=$("$PYTHON" -c "import sys; print(sys.version_info.major)")
PY_MINOR=$("$PYTHON" -c "import sys; print(sys.version_info.minor)")

if [ "$PY_MAJOR" -lt 3 ] || { [ "$PY_MAJOR" -eq 3 ] && [ "$PY_MINOR" -lt 12 ]; }; then
    echo "ERROR: Python 3.12+ required, found $PY_VERSION"
    exit 1
fi
echo "  Python $PY_VERSION ✓"

# Xcode Command Line Tools
if ! xcode-select -p &>/dev/null; then
    echo "ERROR: Xcode Command Line Tools not installed."
    echo "  Run: xcode-select --install"
    exit 1
fi
echo "  Xcode CLT ✓"

# hdiutil (should always be available on macOS)
if ! command -v hdiutil &>/dev/null; then
    echo "ERROR: hdiutil not found (required for DMG creation)"
    exit 1
fi
echo "  hdiutil ✓"

echo ""

# ---------------------------------------------------------------------------
# Step 2: Create build virtual environment
# ---------------------------------------------------------------------------
echo ">>> Step 2: Creating build virtual environment..."

if [ -d "$VENV_DIR" ]; then
    echo "  Removing existing build venv..."
    rm -rf "$VENV_DIR"
fi

"$PYTHON" -m venv "$VENV_DIR"
# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"

pip install --upgrade pip setuptools wheel --quiet

echo "  venv created at $VENV_DIR ✓"
echo ""

# ---------------------------------------------------------------------------
# Step 3: Install dependencies
# ---------------------------------------------------------------------------
echo ">>> Step 3: Installing dependencies..."

# PyInstaller
pip install pyinstaller --quiet
echo "  PyInstaller ✓"

# Core dependencies
pip install numpy pandas scipy matplotlib plotly --quiet
echo "  Core deps (numpy, pandas, scipy, matplotlib, plotly) ✓"

# Optional but bundled
pip install biopython scikit-learn --quiet
echo "  biopython, scikit-learn ✓"

# ablang + PyTorch (unless lite mode)
if [ "$LITE" != "1" ]; then
    echo "  Installing CPU-only PyTorch (this may take a moment)..."
    pip install torch --index-url https://download.pytorch.org/whl/cpu --quiet
    echo "  PyTorch (CPU-only) ✓"

    pip install ablang --quiet
    echo "  ablang ✓"
else
    echo "  LITE mode: skipping ablang/PyTorch"
fi

echo ""

# ---------------------------------------------------------------------------
# Step 4: Generate .icns icon
# ---------------------------------------------------------------------------
echo ">>> Step 4: Generating app icon..."

ICONSET_DIR="$SCRIPT_DIR/Nabbit.iconset"
ICNS_FILE="$SCRIPT_DIR/Nabbit.icns"
SVG_FILE="$PROJECT_DIR/nabbit.svg"

if [ -f "$SVG_FILE" ]; then
    rm -rf "$ICONSET_DIR"
    mkdir -p "$ICONSET_DIR"

    # Try to render SVG to PNG using Python
    "$VENV_DIR/bin/python" -c "
import subprocess, os, sys

svg = '$SVG_FILE'
iconset = '$ICONSET_DIR'
sizes = [16, 32, 64, 128, 256, 512]

# Try cairosvg first
try:
    import cairosvg
    for s in sizes:
        cairosvg.svg2png(url=svg, write_to=os.path.join(iconset, f'icon_{s}x{s}.png'),
                         output_width=s, output_height=s)
        cairosvg.svg2png(url=svg, write_to=os.path.join(iconset, f'icon_{s}x{s}@2x.png'),
                         output_width=s*2, output_height=s*2)
    print('  Generated icons via cairosvg')
    sys.exit(0)
except ImportError:
    pass

# Fallback: use qlmanage (macOS Quick Look) to render SVG to PNG
import tempfile
tmp = tempfile.mkdtemp()
subprocess.run(['qlmanage', '-t', '-s', '1024', '-o', tmp, svg],
               capture_output=True)

# Find the rendered PNG
rendered = None
for f in os.listdir(tmp):
    if f.endswith('.png'):
        rendered = os.path.join(tmp, f)
        break

if not rendered:
    # Last resort: generate a simple placeholder icon with matplotlib
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(10.24, 10.24), dpi=100)
    ax.set_facecolor('#1a1a2e')
    fig.patch.set_facecolor('#1a1a2e')
    ax.text(0.5, 0.5, 'N', transform=ax.transAxes,
            fontsize=400, fontweight='bold', color='#a855f7',
            ha='center', va='center', family='sans-serif')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axis('off')
    rendered = os.path.join(tmp, 'nabbit_icon.png')
    fig.savefig(rendered, dpi=100, bbox_inches='tight', pad_inches=0)
    plt.close()
    print('  Generated placeholder icon via matplotlib')

# Resize to all required sizes using sips
import shutil
for s in sizes:
    for suffix, actual in [(f'{s}x{s}', s), (f'{s}x{s}@2x', s*2)]:
        out = os.path.join(iconset, f'icon_{suffix}.png')
        shutil.copy2(rendered, out)
        subprocess.run(['sips', '-z', str(actual), str(actual), out],
                       capture_output=True)

print('  Generated icons via qlmanage/sips')
" 2>/dev/null || echo "  Warning: icon generation had issues, will use default"

    # Convert iconset to icns
    if [ -d "$ICONSET_DIR" ] && ls "$ICONSET_DIR"/*.png &>/dev/null; then
        iconutil -c icns "$ICONSET_DIR" -o "$ICNS_FILE" 2>/dev/null && \
            echo "  Icon: $ICNS_FILE ✓" || \
            echo "  Warning: iconutil failed, building without custom icon"
        rm -rf "$ICONSET_DIR"
    fi
else
    echo "  No nabbit.svg found, building without custom icon"
fi

echo ""

# ---------------------------------------------------------------------------
# Step 5: Verify nabbit.py compiles
# ---------------------------------------------------------------------------
echo ">>> Step 5: Verifying nabbit.py syntax..."

"$VENV_DIR/bin/python" -c "import py_compile; py_compile.compile('$PROJECT_DIR/nabbit.py', doraise=True)"
echo "  nabbit.py compiles ✓"
echo ""

# ---------------------------------------------------------------------------
# Step 6: Run PyInstaller
# ---------------------------------------------------------------------------
echo ">>> Step 6: Running PyInstaller (this may take several minutes)..."

cd "$PROJECT_DIR"

# Export NABBIT_LITE so the spec file can read it
export NABBIT_LITE="$LITE"

"$VENV_DIR/bin/pyinstaller" --clean --noconfirm "$SCRIPT_DIR/Nabbit.spec"

if [ ! -d "$DIST_DIR/Nabbit.app" ]; then
    echo "ERROR: PyInstaller did not produce Nabbit.app"
    exit 1
fi
echo "  Nabbit.app built ✓"
echo ""

# ---------------------------------------------------------------------------
# Step 7: Smoke test
# ---------------------------------------------------------------------------
echo ">>> Step 7: Smoke testing the bundle..."

# Test 1: Binary loads without crashing (run with --help to trigger imports)
if "$DIST_DIR/Nabbit.app/Contents/MacOS/Nabbit" --help &>/dev/null; then
    echo "  --help test passed ✓"
else
    echo "  WARNING: --help test failed (may still work, continuing...)"
fi

echo ""

# ---------------------------------------------------------------------------
# Step 8: Code signing
# ---------------------------------------------------------------------------
echo ">>> Step 8: Code signing..."

codesign --sign "$CODESIGN_IDENTITY" --force --deep "$DIST_DIR/Nabbit.app" 2>/dev/null && \
    echo "  Signed with identity: $CODESIGN_IDENTITY ✓" || \
    echo "  Warning: code signing failed (app may still run with right-click > Open)"

echo ""

# ---------------------------------------------------------------------------
# Step 9: Create DMG
# ---------------------------------------------------------------------------
echo ">>> Step 9: Creating DMG installer..."

DMG_STAGING="$BUILD_DIR/dmg_staging"
rm -rf "$DMG_STAGING"
mkdir -p "$DMG_STAGING"

# Copy the app
cp -R "$DIST_DIR/Nabbit.app" "$DMG_STAGING/"

# Create Applications symlink for drag-to-install
ln -s /Applications "$DMG_STAGING/Applications"

# Remove any existing DMG
rm -f "$DIST_DIR/$DMG_NAME"

# Create compressed read-only DMG
hdiutil create \
    -volname "$APP_NAME" \
    -srcfolder "$DMG_STAGING" \
    -ov \
    -format UDBZ \
    "$DIST_DIR/$DMG_NAME"

# Clean up staging
rm -rf "$DMG_STAGING"

echo ""
echo "============================================================"
echo "  BUILD COMPLETE"
echo "============================================================"
echo ""
echo "  App:  $DIST_DIR/Nabbit.app"
echo "  DMG:  $DIST_DIR/$DMG_NAME"
echo "  Size: $(du -sh "$DIST_DIR/$DMG_NAME" | cut -f1)"
echo ""
echo "  To install:"
echo "    1. Double-click $DMG_NAME"
echo "    2. Drag Nabbit to Applications"
echo "    3. Right-click Nabbit.app > Open (first launch only)"
echo ""
if [ "$CODESIGN_IDENTITY" = "-" ]; then
    echo "  NOTE: App is ad-hoc signed. Users will need to right-click > Open"
    echo "  on first launch, or run: xattr -cr /Applications/Nabbit.app"
    echo ""
    echo "  For Gatekeeper-trusted distribution, rebuild with:"
    echo "    CODESIGN_IDENTITY='Developer ID Application: ...' bash $0"
fi
echo "============================================================"
