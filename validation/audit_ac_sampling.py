"""Audit AC angular sampling using an isolated, logging-only native build."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,os,re,shutil,subprocess
import numpy as np
from compare_limb_altitudes import make_scene
from compare_limb_angles import OUT as ORIGINAL,ANGLES
import validate_limb_hybrid as h
import validate_combined_ac as c
ROOT=h.ROOT
OUT=ORIGINAL/'sampling_audit'
def main():
 OUT.mkdir(exist_ok=True)
 build=ROOT/'tmp/ac_sampling_audit_build';build.mkdir(exist_ok=True)
 native=ROOT/'Fortran/vlidort_standard'
 manifest=json.loads((ORIGINAL/'manifest.json').read_text())
 hashes={n:hashlib.sha256((ROOT/'Fortran'/n/'DRIVER.exe').read_bytes()).hexdigest() for n in ['vlidort_qs','vlidort_standard']}
 assert hashes==manifest['hashes'],'Production binary differs from original sweep'
 shutil.copytree(native/'mod',build/'mod',dirs_exist_ok=True)
 code=(native/'sparse_aircraft.f90').read_text()
 anchor='      allocate(all_ms(ng,4,nl),all_jac(nq,nl,ng,4,nl))'
 assert code.count(anchor)==1
 code=code.replace(anchor,"      write(*,*) 'AC AUDIT points/legs/total:',points,nlegs,ng\n"+
 "      do sample=1,ng\n"+
 "         write(*,*) 'AC AUDIT node:',sample,all_count(sample),all_angles(sample,:)\n"+
 "      enddo\n"+anchor)
 (build/'sparse_aircraft.f90').write_text(code)
 env=dict(os.environ);env['TEMP']=env['TMP']=str(ROOT/'tmp')
 compiler='C:/mingw64/bin/gfortran.exe'
 subprocess.run([compiler,'-c','-O2','-Wall','-fno-automatic','-Imod','-Jmod','sparse_aircraft.f90','-o','sparse_aircraft.o'],cwd=build,env=env,check=True)
 flat=(native/'makefile').read_text().replace('\\\n',' ')
 variables=dict(re.findall(r'^([A-Z_]+) = (.*)$',flat,re.M))
 def expand(v):
  while '$(' in v:v=re.sub(r'\$\((\w+)\)',lambda m:variables[m[1]],v)
  return v
 sources=[p for k in ['REGULAR_SOURCES','LINEARIZED_SOURCES','FO_SHARED_SOURCES','MS_SHARED_SOURCES','DRIVER_SOURCES'] for p in expand(variables[k]).split()]
 objects=[str(build/'sparse_aircraft.o') if Path(p).name=='sparse_aircraft.f90' else str(native/'obj'/Path(p).with_suffix('.o').name) for p in sources]
 subprocess.run([compiler,*objects,'-static','-Wl,--stack,536870912','-o',str(build/'DRIVER.exe')],cwd=build,env=env,check=True)
 symbols=subprocess.run(['C:/mingw64/bin/nm.exe',str(native/'DRIVER.exe')],capture_output=True,text=True,check=True).stdout
 (OUT/'production_symbols.txt').write_text(symbols)
 forbidden=[l for l in symbols.splitlines() if re.search(r'foms_ms_.*vectorrt|vqs_.*master',l,re.I)]
 assert not forbidden,forbidden
 print('Isolated logging-only build ready; no full QS master symbols.',flush=True)
 p=make_scene();p.settings_rt.jacobian_list=['absorption','haze']
 jobs=[]
 with h.isolated(hashes) as work,ThreadPoolExecutor(max_workers=4) as pool:
  def execute(folder,angle,height,points,tag):
   result=c.run(folder,OUT/(tag+'.log'))
   log=(OUT/(tag+'.log')).read_text()
   matches=re.findall(r'AC AUDIT points/legs/total:\s*(\d+)\s+(\d+)\s+(\d+)',log)
   assert len(matches)==1,(tag,matches)
   requested,legs,total=map(int,matches[0]);assert requested==points and total<=points*legs
   row=dict(vza=angle,height_m=height,requested=points,legs=legs,total_angular_samples=total)
   with np.load(ORIGINAL/f'v{angle}_h{height}_aircraft_combined.npz') as old:
    for key in ['radiance','absorption','haze']:
     error=float(np.max(abs(result[key]-old[key]))/np.max(abs(old[key])))
     row[key+'_relative_difference_from_saved']=error
     if points==3:np.testing.assert_allclose(result[key],old[key],rtol=1e-12,atol=1e-14)
   (OUT/(tag+'.json')).write_text(json.dumps(row,indent=2))
   print(json.dumps(row),flush=True)
   return row
  for angle in ANGLES:
   for height in manifest['heights_m']:
    height=int(height)
    p.viewing_geometry.vza=float(angle);p.viewing_geometry.observer_height_agl=float(height)
    tag=f'v{angle}_h{height}_p3';folder=work/tag
    stage=c.generate(p,work,'aircraft_combined',3)
    shutil.copytree(stage/'INPUT',folder/'INPUT');(folder/'OUTPUT').mkdir()
    shutil.copy2(build/'DRIVER.exe',folder/'DRIVER.exe')
    jobs.append(pool.submit(execute,folder,angle,height,3,tag))
  rows=[j.result() for j in jobs]
  # Changing the sample count must affect a nondegenerate case.
  for angle in [89,91]:
   for points in [2,5]:
    height=15000;p.viewing_geometry.vza=float(angle);p.viewing_geometry.observer_height_agl=float(height)
    tag=f'v{angle}_h{height}_p{points}';folder=work/tag
    stage=c.generate(p,work,'aircraft_combined',points)
    shutil.copytree(stage/'INPUT',folder/'INPUT');(folder/'OUTPUT').mkdir()
    shutil.copy2(build/'DRIVER.exe',folder/'DRIVER.exe')
    row=execute(folder,angle,height,points,tag)
    assert row['haze_relative_difference_from_saved']>1e-12
    rows.append(row)
 (OUT/'results.json').write_text(json.dumps(rows,indent=2))
 (OUT/'provenance.json').write_text(json.dumps(dict(production_hashes=hashes,original_manifest_hashes=manifest['hashes'],audit_sha256=hashlib.sha256((build/'DRIVER.exe').read_bytes()).hexdigest(),instrumentation='Only added writes immediately before the one aircraft_ms_sources call.',full_qs_master_symbols=forbidden),indent=2))
 print('AC SAMPLING AUDIT COMPLETE',flush=True)
if __name__=='__main__':main()
