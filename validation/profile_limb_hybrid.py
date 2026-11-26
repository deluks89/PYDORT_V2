"""Profile isolated, instrumented native copies; never modify production binaries.

Times distinguish input export, process creation/entry, driver execution,
FO/geometry/MS kernels, and complete Standard source solves. QS MS records
homogeneous and particular solution calls separately.
"""
from pathlib import Path
from collections import defaultdict
import argparse
import ctypes
from ctypes import wintypes
import json
import os
import re
import shutil
import struct
import subprocess
import time

import numpy as np
import validate_limb_hybrid as h

ROOT = h.ROOT
TIMER_DECL = '''
   integer(kind=8) :: prof_t0, prof_t1, prof_rate
   real(kind=8) :: prof_c0, prof_c1
'''


def begin():
    return 'call cpu_time(prof_c0)\ncall system_clock(prof_t0,prof_rate)\n'


def end(tag, extra='0'):
    return ("call cpu_time(prof_c1)\ncall system_clock(prof_t1)\n"
            f"write(*,*) 'PROFILE {tag}',prof_c1-prof_c0, &\n"
            f" real(prof_t1-prof_t0,kind=8)/real(prof_rate,kind=8), {extra}\n")


def instrument_calls(source, pattern):
    source = re.sub(r'(?im)^(\s*implicit none\s*)$', lambda m: m[0]+TIMER_DECL, source, count=1)
    lines, result, i = source.splitlines(keepends=True), [], 0
    while i < len(lines):
        match = re.match(r'\s*call\s+(\w+)', lines[i], re.I)
        if match and re.search(pattern, match[1], re.I):
            name = match[1].upper()
            result.append(begin())
            while True:
                line = lines[i]
                result.append(line)
                i += 1
                code = line.split('!')[0].rstrip()
                if code and not code.endswith('&'):
                    break
            result.append(end(name, 'ngeoms' if ('PARTIC' in name or name == 'MSSEG_LP_FOURIER') else '0'))
        else:
            result.append(lines[i])
            i += 1
    return ''.join(result)


def instrument_driver(source):
    source = re.sub(r'(?im)^(\s*implicit none\s*)$', lambda m: m[0]+TIMER_DECL, source, count=1)
    marker = re.search(r'(?im)^\s*GMASK\s*=', source).start()
    source = source[:marker]+"\nwrite(*,*) 'PROFILE_ENTRY'\nflush(6)\n"+begin()+source[marker:]
    final = re.search(r'(?im)^\s*(?:contains|end program DRIVER)\s*$', source).start()
    # QS exits through an explicit STOP before END PROGRAM.
    stops = list(re.finditer(r'(?im)^\s*stop\s*$', source[:final]))
    if stops:
        final = stops[-1].start()
    return source[:final]+end('DRIVER_TOTAL')+source[final:]


def build(work, compiler):
    env = dict(os.environ, TEMP=str(ROOT/'tmp'), TMP=str(ROOT/'tmp'))
    specifications = {
        'vlidort_qs': [
            ('DRIVER.f90', None),
            ('vqs_main/masters/vqs_masters_V1p4.f90', r'^FOMS_'),
            ('vqs_main/masters/vqs_lp_masters_V1p4.f90', r'^FOMS_'),
            ('vqs_main/foms/FOMS_MS_VectorRT_Mk3_Exp.f90', r'^MSSEG_.*FOURIER$'),
            ('vqs_main/foms/FOMS_MS_LP_VectorRT_Mk3_Exp.f90', r'^MSSEG_.*FOURIER$'),
        ],
        'vlidort_standard': [('DRIVER.f90', None), ('sparse_aircraft.f90', r'^VLIDORT_(LPS_)?MASTER$')],
    }
    for backend, entries in specifications.items():
        original = ROOT/'Fortran'/backend
        folder = work/backend
        folder.mkdir(parents=True, exist_ok=True)
        for sub in ('INPUT', 'INPUT/xs', 'INPUT/aerosol_moments',
                    'INPUT/aerosol_model_jacobians', 'OUTPUT', 'mod'):
            (folder/sub).mkdir(parents=True, exist_ok=True)
        replacements = {}
        flags = ['-c', '-std=legacy', '-mcmodel=large', '-O3'] if backend == 'vlidort_qs' else [
            '-c', '-O2', '-fno-automatic']
        for filename, pattern in entries:
            path = original/filename
            source = path.read_text()
            source = instrument_driver(source) if pattern is None else instrument_calls(source, pattern)
            copied = folder/path.name
            copied.write_text(source)
            obj = copied.with_suffix('.o')
            command = [compiler, *flags, '-I'+str(folder/'mod'), '-I'+str(original/'mod'),
                       '-J'+str(folder/'mod'), str(copied), '-o', str(obj)]
            result = subprocess.run(command, capture_output=True, text=True, cwd=folder, env=env)
            (folder/(path.stem+'_build.log')).write_text(result.stdout+result.stderr)
            if result.returncode:
                raise RuntimeError(result.stdout+result.stderr)
            replacements[obj.name] = obj
        objects = [replacements.get(p.name, p) for p in sorted((original/'obj').glob('*.o'))]
        command = [compiler, *map(str, objects), '-static', '-Wl,--stack,536870912',
                   '-o', str(folder/'DRIVER.exe')]
        subprocess.run(command, check=True, cwd=folder, env=env)
        print('Built profiled', backend, flush=True)


def process_cpu(process):
    stamps = [wintypes.FILETIME() for _ in range(4)]
    function = ctypes.windll.kernel32.GetProcessTimes
    function.argtypes = [wintypes.HANDLE]+[ctypes.POINTER(wintypes.FILETIME)]*4
    if not function(int(process._handle), *[ctypes.byref(t) for t in stamps]):
        raise ctypes.WinError()
    return sum((t.dwHighDateTime << 32)+t.dwLowDateTime for t in stamps[2:])*1e-7


def image_sizes(path):
    data = path.read_bytes()
    pe = struct.unpack_from('<I', data, 0x3c)[0]
    return dict(exe_file_MB=len(data)/1e6,
                pe_image_MB=struct.unpack_from('<I', data, pe+24+56)[0]/1e6)


def run(folder, log):
    start = time.perf_counter()
    proc = subprocess.Popen([str(folder/'DRIVER.exe')], cwd=folder,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    created = time.perf_counter()
    text, entry = [], None
    for line in proc.stdout:
        if 'PROFILE_ENTRY' in line:
            entry = time.perf_counter()-start
        text.append(line)
    proc.wait(timeout=900)
    elapsed = time.perf_counter()-start
    cpu = process_cpu(proc)
    io_start = time.perf_counter()
    log.write_text(''.join(text))
    log_time = time.perf_counter()-io_start
    if proc.returncode:
        raise RuntimeError(''.join(text))
    groups = defaultdict(lambda: dict(calls=0, cpu_s=0.0, wall_s=0.0, geometries=0))
    for line in text:
        parts = line.split()
        if len(parts) == 5 and parts[0] == 'PROFILE':
            group = groups[parts[1]]
            group['calls'] += 1
            group['cpu_s'] += float(parts[2])
            group['wall_s'] += float(parts[3])
            group['geometries'] += int(parts[4])
    result = dict(wall_s=elapsed, process_cpu_s=cpu, popen_s=created-start,
                  entry_s=entry, stages=dict(groups), **image_sizes(folder/'DRIVER.exe'))
    path = folder/'OUTPUT/radiances.out'
    io_start = time.perf_counter()
    if path.exists() and path.stat().st_size:
        result['radiance'] = np.loadtxt(path).reshape(-1).tolist()
    result['log_write_s'] = log_time
    result['result_read_s'] = time.perf_counter()-io_start
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reuse-build', action='store_true')
    parser.add_argument('--vector', action='store_true')
    parser.add_argument('--wavelengths', type=int, default=1)
    parser.add_argument('--repeats', type=int, default=2)
    parser.add_argument('--startup-control', action='store_true')
    args = parser.parse_args()
    (ROOT/'tmp').mkdir(exist_ok=True)
    if args.startup_control:
        startup_control()
        return
    work = ROOT/'tmp/limb_profile_build'
    work.mkdir(exist_ok=True)
    output = ROOT/'Examples/testing/limb_hybrid/profile'/f"{'vector' if args.vector else 'scalar'}_{args.wavelengths}wl"
    output.mkdir(parents=True, exist_ok=True)
    if not args.reuse_build:
        build(work, shutil.which('gfortran') or 'C:/mingw64/bin/gfortran.exe')
    for backend in ('vlidort_qs', 'vlidort_standard'):
        for sub in ('INPUT/xs', 'INPUT/aerosol_moments', 'INPUT/aerosol_model_jacobians'):
            (work/backend/sub).mkdir(parents=True, exist_ok=True)
    original = h.rt.QS_BACKEND_PATH, h.rt.STANDARD_BACKEND_PATH
    h.rt.QS_BACKEND_PATH = (work/'vlidort_qs').as_posix()+'/'
    h.rt.STANDARD_BACKEND_PATH = (work/'vlidort_standard').as_posix()+'/'
    results = []
    try:
        p = h.scene(15000, 91, 90, args.vector)
        p.wavelengths = np.linspace(477, 478, args.wavelengths)
        p.surface_albedo.change_wavelengths(p.wavelengths)
        p.aerosols.run_models(p.wavelengths)
        for name, mode, ss, points, radiance in [
            ('qs_empty', 'qs', False, 5, False),
            ('standard_empty', 'aircraft_limb_ms', False, 5, False),
            ('qs_full', 'qs', False, 5, True),
            ('qs_fo', 'qs', True, 5, True),
            ('standard_5', 'aircraft_limb_ms', False, 5, True),
            ('standard_9', 'aircraft_limb_ms', False, 9, True),
        ]:
            p.settings_rt.do_radiance = radiance
            start = time.perf_counter()
            folder = h.generate(p, work, mode, points, ss)
            export_time = time.perf_counter()-start
            for repeat in range(args.repeats):
                row = dict(case=name, repeat=repeat, vector=args.vector, wavelengths=args.wavelengths,
                           export_s=export_time, **run(folder, output/f'{name}_{repeat}.log'))
                results.append(row)
                (output/'profile.json').write_text(json.dumps(results, indent=2))
                print(json.dumps(row), flush=True)
    finally:
        h.rt.QS_BACKEND_PATH, h.rt.STANDARD_BACKEND_PATH = original


def startup_control():
    import sys
    folder = ROOT/'tmp/startup_control'
    folder.mkdir(exist_ok=True)
    source = folder/'empty.f90'
    text = "program empty\nwrite(*,*) 'PROFILE_ENTRY'\nflush(6)\nend program empty\n"
    compiler = shutil.which('gfortran') or 'C:/mingw64/bin/gfortran.exe'
    env = dict(os.environ, TEMP=str(ROOT/'tmp'), TMP=str(ROOT/'tmp'))
    if not source.exists() or source.read_text() != text or not (folder/'DRIVER.exe').exists():
        source.write_text(text)
        subprocess.run([compiler, '-static', str(source), '-o', str(folder/'DRIVER.exe')], check=True, env=env)
    rows = [dict(case='empty_fortran', **run(folder, folder/'empty.log'))]
    for name, command in [
        ('python', [sys._base_executable, '-c', "print('PROFILE_ENTRY',flush=True)"]),
    ]:
        start = time.perf_counter()
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        created = time.perf_counter()
        entry = None
        lines = []
        for line in process.stdout:
            lines.append(line)
            if line.strip() == 'PROFILE_ENTRY':
                entry = time.perf_counter()-start
        process.wait(timeout=60)
        rows.append(dict(case=name, popen_s=created-start, entry_s=entry,
                         wall_s=time.perf_counter()-start, process_cpu_s=process_cpu(process),
                         returncode=process.returncode, output=''.join(lines)))
    output = ROOT/'Examples/testing/limb_hybrid/profile/startup_control.json'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(rows, indent=2))
    print(json.dumps(rows, indent=2), flush=True)


if __name__ == '__main__':
    main()
