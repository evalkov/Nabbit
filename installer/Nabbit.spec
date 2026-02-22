# -*- mode: python ; coding: utf-8 -*-
"""
Nabbit PyInstaller Spec File
=============================
Builds a macOS .app bundle containing Python 3.12+, all dependencies,
and nabbit.py. The resulting app launches the web UI on double-click.

Usage:
    pyinstaller --clean --noconfirm installer/Nabbit.spec

Or use the build script:
    bash installer/build_macos.sh
"""
import os
import sys
from PyInstaller.utils.hooks import collect_data_files, collect_submodules

block_cipher = None

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
NABBIT_DIR = os.path.abspath(os.path.join(SPECPATH, '..'))
NABBIT_PY = os.path.join(NABBIT_DIR, 'nabbit.py')
LAUNCHER = os.path.join(NABBIT_DIR, 'installer', 'nabbit_launcher.py')
ICON_FILE = os.path.join(NABBIT_DIR, 'installer', 'Nabbit.icns')

# ---------------------------------------------------------------------------
# Whether to include ablang + PyTorch (set NABBIT_LITE=1 to skip)
# ---------------------------------------------------------------------------
INCLUDE_ABLANG = os.environ.get('NABBIT_LITE', '0') != '1'

# ---------------------------------------------------------------------------
# Hidden imports
# ---------------------------------------------------------------------------
# nabbit.py uses try/except for optional imports, so PyInstaller's static
# analysis can't detect them. List everything explicitly.
hidden_imports = [
    # --- Core (normally auto-detected, but be explicit) ---
    'numpy',
    'numpy.core._methods',
    'numpy.lib.format',
    'numpy.random',
    'pandas',
    'pandas._libs.tslibs.base',
    'pandas._libs.tslibs.np_datetime',
    'scipy',
    'scipy.stats',
    'scipy.special',
    'scipy.sparse',
    'scipy.sparse.csgraph',
    'scipy.sparse.csgraph._validation',
    'scipy.spatial',
    'scipy.spatial.distance',
    'scipy.cluster',
    'scipy.cluster.hierarchy',
    'scipy._lib.messagestream',
    'matplotlib',
    'matplotlib.pyplot',
    'matplotlib.colors',
    'matplotlib.cm',
    'matplotlib.backends',
    'matplotlib.backends.backend_agg',

    # --- Plotly (dashboard) ---
    'plotly',
    'plotly.graph_objects',
    'plotly.subplots',
    'plotly.io',

    # --- Biopython (optional but bundled) ---
    'Bio',
    'Bio.Phylo',
    'Bio.Phylo.TreeConstruction',
    'Bio.Phylo.BaseTree',
    'Bio.Align',
    'Bio.Align.substitution_matrices',
    'Bio.Seq',
    'Bio.SeqRecord',
    'Bio.SeqIO',

    # --- scikit-learn (optional but bundled) ---
    'sklearn',
    'sklearn.decomposition',
    'sklearn.decomposition._pca',
    'sklearn.utils',
    'sklearn.utils._cython_blas',
    'sklearn.utils._typedefs',
    'sklearn.neighbors._typedefs',
    'sklearn.neighbors._partition_nodes',

    # --- stdlib that PyInstaller sometimes misses ---
    'http.server',
    'webbrowser',
    'multiprocessing',
    'multiprocessing.pool',
    'concurrent.futures',
    'gzip',
    'tempfile',
    'uuid',
    'logging',
    'argparse',
    'json',
    'csv',
    'collections',
    'pathlib',
    'typing',
    'itertools',
    're',
]

# ablang + PyTorch hidden imports (large, optional)
if INCLUDE_ABLANG:
    hidden_imports += [
        'ablang',
        'torch',
        'torch.nn',
        'torch.nn.functional',
        'torch.utils',
        'torch.jit',
    ]
    # Collect all torch submodules since it uses heavy lazy-loading
    hidden_imports += collect_submodules('torch')

# ---------------------------------------------------------------------------
# Data files
# ---------------------------------------------------------------------------
datas = []

# nabbit.py itself (needed so the frozen binary can locate it as a module)
datas += [(NABBIT_PY, '.')]

# Include the SVG logo if it exists
nabbit_svg = os.path.join(NABBIT_DIR, 'nabbit.svg')
if os.path.exists(nabbit_svg):
    datas += [(nabbit_svg, '.')]

# Package data files that PyInstaller doesn't auto-collect
datas += collect_data_files('plotly')
datas += collect_data_files('scipy')
datas += collect_data_files('sklearn')
datas += collect_data_files('Bio')

if INCLUDE_ABLANG:
    datas += collect_data_files('ablang')

# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------
a = Analysis(
    [LAUNCHER],
    pathex=[NABBIT_DIR],
    binaries=[],
    datas=datas,
    hiddenimports=hidden_imports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        # GUI toolkits not needed (Nabbit uses a web browser)
        'tkinter', '_tkinter',
        'PyQt5', 'PyQt6', 'PySide2', 'PySide6',
        # Dev/notebook tools
        'IPython', 'notebook', 'jupyter', 'jupyter_client', 'jupyter_core',
        # Test frameworks (keep unittest — numpy.testing needs it)
        'test', 'tests', 'pytest',
        # Unused torch backends (if included)
        'torch.distributed', 'torch.testing',
        'caffe2',
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

# ---------------------------------------------------------------------------
# PYZ (Python bytecode archive)
# ---------------------------------------------------------------------------
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

# ---------------------------------------------------------------------------
# EXE
# ---------------------------------------------------------------------------
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,   # onedir mode (fast startup, no extraction)
    name='Nabbit',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,               # UPX corrupts numpy/scipy shared libs
    console=False,            # No Terminal window on macOS
    icon=ICON_FILE if os.path.exists(ICON_FILE) else None,
)

# ---------------------------------------------------------------------------
# COLLECT (onedir bundle)
# ---------------------------------------------------------------------------
coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    name='Nabbit',
)

# ---------------------------------------------------------------------------
# macOS APP bundle
# ---------------------------------------------------------------------------
app = BUNDLE(
    coll,
    name='Nabbit.app',
    icon=ICON_FILE if os.path.exists(ICON_FILE) else None,
    bundle_identifier='com.nabbit.app',
    info_plist={
        'CFBundleName': 'Nabbit',
        'CFBundleDisplayName': 'Nabbit',
        'CFBundleVersion': '1.0.0',
        'CFBundleShortVersionString': '1.0.0',
        'CFBundleExecutable': 'Nabbit',
        'CFBundleIdentifier': 'com.nabbit.app',
        'CFBundleInfoDictionaryVersion': '6.0',
        'CFBundlePackageType': 'APPL',
        'LSMinimumSystemVersion': '12.0',
        'NSHighResolutionCapable': True,
        'LSBackgroundOnly': False,
        'NSAppleEventsUsageDescription':
            'Nabbit needs to open file-picker dialogs.',
        'LSApplicationCategoryType':
            'public.app-category.medical',
    },
)
