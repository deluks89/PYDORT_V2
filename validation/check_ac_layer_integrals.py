"""Compatibility entry point for the current curved AC integral checks."""
from pathlib import Path
import runpy

if __name__ == '__main__':
    runpy.run_path(str(Path(__file__).with_name('check_ac_chord_integrals.py')), run_name='__main__')
