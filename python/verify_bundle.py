"""Compatibility entry point for the package-level source integrity check."""
from pathlib import Path
import runpy
runpy.run_path(str(Path(__file__).resolve().parents[1]/'verify_release.py'),run_name='__main__')
