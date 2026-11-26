"""Profile Standard polarized AC internals without changing the production executable."""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import numpy as np
import validate_limb_hybrid as h
import validate_combined_ac as c
import profile_limb_hybrid as prof

ROOT=h.ROOT
folder=ROOT/'tmp/combined_ac_profile'
folder.mkdir(parents=True,exist_ok=True)
(folder/'mod').mkdir(exist_ok=True)
standard=ROOT/'Fortran/vlidort_standard'
source=(standard/'vlidort_main/linearized/vlidort_lps_masters_V2p8p3.f90').read_text()
pattern=(r'^VLIDORT_(L_)?[QU]HOM_SOLUTION$|^(LP_)?BVP_(MATRIXSETUP|SOLUTION)_MASTER$|'
         r'^VLIDORT_(LP_|L_)?[CG](BEAM|USER)_SOLUTION$|^(DN|UP)USER_PROFILEWF$|'
         r'^VLIDORT_(DN|UP)USER_INTENSITY$')
source=prof.instrument_calls(source,pattern)
# The selected routines also occur in the separate Fourier subroutine.
source=source.replace(prof.TIMER_DECL,'')
source=re.sub(r'(?im)^([ \t]*implicit none)[ \t]*$',lambda m:m[0]+prof.TIMER_DECL,source)
copied=folder/'vlidort_lps_masters_V2p8p3.f90'
copied.write_text(source)
compiler='C:/mingw64/bin/gfortran.exe'
env=dict(os.environ,TEMP=str(ROOT/'tmp'),TMP=str(ROOT/'tmp'))
result=subprocess.run([compiler,'-c','-O2','-fno-automatic',
    '-I'+str(standard/'mod'),'-J'+str(folder/'mod'),str(copied),'-o',str(copied.with_suffix('.o'))],
    capture_output=True,text=True,env=env)
(folder/'compile.log').write_text(result.stdout+result.stderr)
assert result.returncode==0,result.stderr
obsolete={'MSSeg_pars.o','FOMS_pars.o','vqs_outputs_def.o','vqs_lin_outputs_def.o'}
objects=[copied.with_suffix('.o') if p.name==copied.with_suffix('.o').name else p
         for p in sorted((standard/'obj').glob('*.o')) if p.name not in obsolete]
subprocess.run([compiler,*map(str,objects),'-static','-Wl,--stack,536870912','-o',str(folder/'DRIVER.exe')],
               check=True,env=env)
hashes={name:hashlib.sha256((ROOT/'Fortran'/name/'DRIVER.exe').read_bytes()).hexdigest()
        for name in ['vlidort_qs','vlidort_standard']}
output=ROOT/'Examples/testing/limb_hybrid/combined_profile'
output.mkdir(parents=True,exist_ok=True)
with h.isolated(hashes) as work:
    import shutil
    shutil.copy2(folder/'DRIVER.exe',work/'vlidort_standard/DRIVER.exe')
    p=h.scene(15000,91,90,True)
    p.settings_rt.do_deltam_phase_function_scaling=True
    p.settings_rt.jacobian_list=['absorption','haze','haze_asym0']
    native=c.generate(p,work,'aircraft_combined',3)
    result=c.run(native,output/'native.log')
    reference=np.load(ROOT/'Examples/testing/limb_hybrid/combined_vector/h15000_v91_a90_ac3.npz')
    for key in ['radiance','absorption','haze','haze_asym0']:
        np.testing.assert_allclose(result[key],reference['model_'+key],rtol=1e-12,atol=1e-14)
    totals={}
    for line in (output/'native.log').read_text().splitlines():
        m=re.search(r'PROFILE (\w+)\s+([\d.Ee+\-]+)\s+([\d.Ee+\-]+)',line)
        if m:
            item=totals.setdefault(m[1],dict(calls=0,wall_seconds=0.0))
            item['calls']+=1;item['wall_seconds']+=float(m[3])
    (output/'profile.json').write_text(json.dumps(dict(stages=totals,
        native_fo_seconds=result['fo_seconds'],native_ms_seconds=result['ms_seconds']),indent=2))
    print(json.dumps(totals,indent=2),flush=True)
