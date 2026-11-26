"""Build a complete isolated binary64 Standard copy, with dependency-safe concurrency."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import os
import re
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'tmp/standard_binary64/Fortran'

def main():
    for name in ('vlidort_standard','vlidort_qs','driver_support'):
        shutil.copytree(ROOT/'Fortran'/name, BASE/name,dirs_exist_ok=True,
            ignore=shutil.ignore_patterns('obj','mod','INPUT','OUTPUT','*.exe','*.o','*.mod','*.log'))
    folder=BASE/'vlidort_standard'
    text=(folder/'makefile').read_text().replace(' -freal-8-real-16','')
    (folder/'makefile').write_text(text)
    flat=text.replace('\\\n',' ')
    variables=dict(re.findall(r'^([A-Z_]+) = (.*)$',flat,re.M))
    def expand(value):
        while '$(' in value:
            value=re.sub(r'\$\((\w+)\)',lambda m:variables[m[1]],value)
        return value
    sources=[folder/p for key in ('REGULAR_SOURCES','LINEARIZED_SOURCES','FO_SHARED_SOURCES','MS_SHARED_SOURCES','DRIVER_SOURCES')
             for p in expand(variables[key]).split()]
    modules={}; used={}
    for source in sources:
        code=source.read_text()
        for name in re.findall(r'^\s*module\s+(?!procedure\b)(\w+)',code,re.M|re.I):
            modules[name.lower()]=source
        used[source]=set(x.lower() for x in re.findall(r'^\s*use(?:\s*,[^:]+::)?\s+(\w+)',code,re.M|re.I))
    deps={s:{modules[m] for m in used[s] if m in modules and modules[m]!=s} for s in sources}
    for sub in ('obj','mod','build_logs'): (folder/sub).mkdir(exist_ok=True)
    env=dict(os.environ);env['TEMP']=env['TMP']=str(BASE.parent)
    compiler='C:/mingw64/bin/gfortran.exe'
    flags=['-c','-O2','-Wall','-fno-automatic','-Imod','-Jmod']
    def compile_one(source):
        result=subprocess.run([compiler,*flags,str(source),'-o',str(folder/'obj'/source.with_suffix('.o').name)],
                              cwd=folder,env=env,capture_output=True,text=True)
        (folder/'build_logs'/source.with_suffix('.log').name).write_text(result.stdout+result.stderr)
        if result.returncode: raise RuntimeError(str(source)+'\n'+result.stderr)
        return source
    remaining=set(sources);done=set();running={}
    with ThreadPoolExecutor(max_workers=8) as pool:
        while remaining or running:
            ready=[s for s in sources if s in remaining and deps[s]<=done]
            for source in ready[:8-len(running)]:
                remaining.remove(source);running[pool.submit(compile_one,source)]=source
            if not running: raise RuntimeError('Unresolved dependencies: '+str(remaining))
            finished,_=wait(running,return_when=FIRST_COMPLETED)
            for future in finished:
                done.add(future.result());running.pop(future)
                print(f'Compiled {len(done)}/{len(sources)}: {future.result().name}',flush=True)
    objects=[str(folder/'obj'/s.with_suffix('.o').name) for s in sources]
    subprocess.run([compiler,*objects,'-static','-Wl,--stack,536870912','-o',str(folder/'DRIVER.exe')],
                   cwd=folder,env=env,check=True)
    print('Isolated binary64 build complete.',flush=True)

if __name__=='__main__':main()
