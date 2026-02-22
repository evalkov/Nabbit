#!/usr/bin/env python3
"""
Nabbit macOS App Launcher
=========================
Entry point for the PyInstaller-frozen Nabbit application.

When double-clicked (no CLI args): starts the web launcher UI.
When re-invoked by the web server subprocess: passes args through to the pipeline.
"""
import multiprocessing
import sys
import os


def main():
    # freeze_support() MUST be called before anything else in a frozen app,
    # otherwise spawned worker processes (ProcessPoolExecutor) will crash.
    multiprocessing.freeze_support()

    # When launched as a macOS .app with no arguments, default to --serve mode
    # with port 0 (auto-select a free port).
    if getattr(sys, 'frozen', False) and len(sys.argv) == 1:
        sys.argv = [sys.argv[0], '--serve', '--port', '0']

    # Add the bundle's resource directory to the path so nabbit.py can be found
    if getattr(sys, 'frozen', False):
        bundle_dir = getattr(sys, '_MEIPASS', os.path.dirname(sys.executable))
        if bundle_dir not in sys.path:
            sys.path.insert(0, bundle_dir)

    from nabbit import main as nabbit_main
    nabbit_main()


if __name__ == '__main__':
    main()
