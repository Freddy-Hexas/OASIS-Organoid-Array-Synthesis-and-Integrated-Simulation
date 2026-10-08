"""Portable workspace paths shared by the web server and local tools.

All defaults are relative to this package. Relative environment paths are
also relative to the package, so a different working directory is harmless.
"""
from pathlib import Path
import os


PACKAGE_DIR = Path(__file__).resolve().parent.parent


def configured_path(name, default):
    value = os.environ.get(name, '').strip()
    path = Path(value).expanduser() if value else default
    if not path.is_absolute():
        path = PACKAGE_DIR / path
    return path.resolve()


DATA_DIR = configured_path('GDS_WORKBENCH_DATA_DIR', PACKAGE_DIR / 'data')
RUNS_DIR = configured_path('GDS_WORKBENCH_RUNS_DIR', PACKAGE_DIR / 'outputs' / 'runs')
GEOMETRY_DIR = PACKAGE_DIR / 'outputs' / 'geometry'
