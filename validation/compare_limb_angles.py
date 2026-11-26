"""Sweep limb viewing angle/aircraft height: radiance, O4 and two Jacobian families."""
from concurrent.futures import ThreadPoolExecutor
import argparse,hashlib,json,shutil,time
import numpy as np
from compare_limb_altitudes import make_scene,MODES,OUT as PREVIOUS
import validate_limb_hybrid as h
import validate_combined_ac as c

OUT=h.ROOT/'Examples/testing/limb_angle_comparison'
ANGLES=[80,85,88,89,91,92,95,100]
def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--angles',type=float,nargs='+',default=ANGLES)
 args=parser.parse_args()
 OUT.mkdir(parents=True,exist_ok=True)
 t=time.perf_counter();p=make_scene();p.settings_rt.jacobian_list=['absorption','haze']
 base=json.loads((PREVIOUS/'manifest.json').read_text());heights=np.array(base['heights_m'])
 grid=p.vertical_grid.copy();np.testing.assert_array_equal(grid,base['vertical_grid_m'])
 hashes={n:hashlib.sha256((h.ROOT/'Fortran'/n/'DRIVER.exe').read_bytes()).hexdigest() for n in ['vlidort_qs','vlidort_standard']}
 manifest=base|dict(angles=ANGLES,hashes=hashes,jacobians=['absorption','haze'],
  initialization_seconds=time.perf_counter()-t,
  parallelism='Up to three isolated native processes for ONE geometry; never multiple QS cases. Input generation serial.',
  timing='CPU is per native process; wall times measured concurrently, not serial benchmark times. Export and staging logged separately.',
  unsupported_standard='VZA>90 and impact radius>=surface radius: public Standard would fall back to QS; no substitute result is counted.',
  threshold_definition='1% and 5% maximum across relative I/O4 errors and per-profile-peak-normalized absorption/haze Jacobian errors.')
 mp=OUT/'manifest.json'
 if mp.exists():
  old=json.loads(mp.read_text())
  assert old['hashes']==hashes and old['vertical_grid_m']==grid.tolist(),'Stale results: executable or grid changed.'
 else:mp.write_text(json.dumps(manifest,indent=2))
 rows=[]
 with h.isolated(hashes) as work,ThreadPoolExecutor(max_workers=3) as pool:
  folders={}
  for mode,_ in MODES:
   folder=work/('evaluation_'+mode);folders[mode]=folder
   (folder/'INPUT').mkdir(parents=True);(folder/'OUTPUT').mkdir()
   shutil.copy2(work/('vlidort_qs' if mode=='qs' else 'vlidort_standard')/'DRIVER.exe',folder/'DRIVER.exe')
  for angle in args.angles:
   for i,height in enumerate(heights):
    p.viewing_geometry.vza=float(angle);p.viewing_geometry.observer_height_agl=float(height)
    impact=(p.settings_rt.planet_radius*1000+height)*np.sin(np.deg2rad(angle))
    turning=bool(angle>90 and impact>=p.settings_rt.planet_radius*1000+grid[0])
    path_type='upward' if angle<90 else ('turning' if turning else 'ground_intersecting')
    case=f'v{angle:g}_h{height:g}';jobs={};pending={};started=time.perf_counter()
    for mode,label in MODES[i%3:]+MODES[:i%3]:
     tag=case+'_'+mode;checkpoint=OUT/(tag+'.json')
     if checkpoint.exists():
      rows.append(json.loads(checkpoint.read_text()));continue
     common=dict(vza=float(angle),height_m=float(height),mode=mode,path_type=path_type,impact_altitude_m=float(impact-p.settings_rt.planet_radius*1000))
     if mode=='standard' and turning:
      row=common|dict(status='unsupported',reason='Downward limb ray misses surface; Standard public API falls back to QS.')
      checkpoint.write_text(json.dumps(row,indent=2));rows.append(row);continue
     export_start=time.perf_counter()
     stage=c.generate(p,work,mode,3)
     folder=folders[mode]
     shutil.copytree(stage/'INPUT',folder/'INPUT',dirs_exist_ok=True)
     for path in (folder/'OUTPUT').iterdir():
      if path.is_file():path.unlink()
     export=time.perf_counter()-export_start
     pending[mode]=(tag,common,export)
    # All input folders are independent before any native call starts.
    for mode,(tag,common,export) in pending.items():
     jobs[mode]=pool.submit(c.run,folders[mode],OUT/(tag+'.log'))
    for mode,future in jobs.items():
     tag,common,export=pending[mode]
     try:
      result=future.result()
      for key in ['absorption','haze']:
       assert result[key].shape==(len(grid)-1,) and np.all(np.isfinite(result[key]))
      rad=float(result['radiance'][0])
      scd=float(np.dot(p.get_o4_profile()*np.diff(grid)*100,-result['absorption'][::-1]/rad))
      assert rad>0 and np.isfinite(scd) and scd>0
      np.savez(OUT/(tag+'.npz'),grid=grid,o4_scd=scd,**result)
      row=common|dict(status='ok',radiance=rad,o4_scd=scd,cpu_seconds=result['cpu_seconds'],
       native_wall_seconds=result['seconds'],export_staging_seconds=export,
       wall_seconds=export+result['seconds'],fo_seconds=result['fo_seconds'],ms_seconds=result['ms_seconds'])
     except Exception as exc:
      row=common|dict(status='failed',reason=str(exc))
     (OUT/(tag+'.json')).write_text(json.dumps(row,indent=2));rows.append(row)
    statuses={r['mode']:r['status'] for r in rows if r['vza']==angle and r['height_m']==height}
    print(json.dumps(dict(case=case,status=statuses,batch_wall_seconds=time.perf_counter()-started)),flush=True)
    (OUT/'results.json').write_text(json.dumps(rows,indent=2))
 print('ANGLE SWEEP COMPLETE',flush=True)
if __name__=='__main__':main()
