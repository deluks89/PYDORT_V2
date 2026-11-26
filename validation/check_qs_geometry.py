"""Compile and run the QS geometry regression against the built native objects."""
from pathlib import Path
import os
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    (ROOT / 'tmp').mkdir(exist_ok=True)
    compiler = shutil.which('gfortran') or 'C:/mingw64/bin/gfortran.exe'
    backend = ROOT / 'Fortran/vlidort_qs'
    env = dict(os.environ, TEMP=str(ROOT / 'tmp'), TMP=str(ROOT / 'tmp'))
    with tempfile.TemporaryDirectory(prefix='geometry_', dir=ROOT / 'tmp') as temp:
        executable = Path(temp) / 'geometry_check.exe'
        objects = [backend / 'obj' / (name+'.o') for name in (
            'MSSeg_pars', 'FOMS_pars', 'FOMS_Geometry_Routines', 'FOMS_Geometry_Masters')]
        command = [compiler, '-static', '-mcmodel=large', '-I'+str(backend / 'mod'),
                   str(ROOT / 'validation/check_qs_limb_geometry.f90'),
                   *map(str, objects), '-o', str(executable)]
        subprocess.run(command, check=True, cwd=temp, env=env)
        subprocess.run([str(executable)], check=True, cwd=temp, env=env)


if __name__ == '__main__':
    main()
