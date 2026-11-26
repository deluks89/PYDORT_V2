"""Compare full absorption and aerosol extinction Jacobian profiles for 20 limb heights."""
import hashlib,json,time
import numpy as np
from compare_limb_altitudes import make_scene,MODES,OUT as PREVIOUS
import validate_limb_hybrid as h
import validate_combined_ac as c

OUT=h.ROOT/'Examples/testing/limb_altitude_jacobians'
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 start=time.perf_counter()
 p=make_scene();p.settings_rt.jacobian_list=['absorption','haze']
 grid=p.vertical_grid.copy()
 manifest=json.loads((PREVIOUS/'manifest.json').read_text())
 heights=np.array(manifest['heights_m'])
 np.testing.assert_array_equal(grid,manifest['vertical_grid_m'])
 assert len(heights)==20 and np.all(np.isin(heights,grid))
 hashes={n:hashlib.sha256((h.ROOT/'Fortran'/n/'DRIVER.exe').read_bytes()).hexdigest() for n in ['vlidort_qs','vlidort_standard']}
 manifest.update(hashes=hashes,jacobians=['absorption','haze'],
  initialization_seconds=time.perf_counter()-start,
  absorption_definition='d(I/E_sun)/d(layer absorption optical depth), sr^-1',
  haze_definition='d(I/E_sun)/d(layer aerosol extinction coefficient at 477 nm), km sr^-1',
  array_order='NPZ native arrays top-down; all_profiles.csv and figures bottom-up',
  error_definition='100*(J_model-J_QS)/max_layer(abs(J_QS)), separately per aircraft height and Jacobian family')
 (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2))
 rows=[]
 with h.isolated(hashes) as work:
  for i,height in enumerate(heights):
   p.viewing_geometry.observer_height_agl=float(height)
   for mode,label in MODES[i%3:]+MODES[:i%3]:
    t=time.perf_counter();folder=c.generate(p,work,mode,3);generation=time.perf_counter()-t
    tag=f'h{height:g}_{mode}'
    result=c.run(folder,OUT/(tag+'.log'))
    for key in ['absorption','haze']:
     assert result[key].shape==(len(grid)-1,)
     assert np.all(np.isfinite(result[key]))
    old=np.load(PREVIOUS/(tag+'.npz'))
    np.testing.assert_allclose(result['radiance'],old['radiance'],rtol=1e-10,atol=1e-13)
    assert np.max(np.abs(result['absorption']-old['absorption']))/np.max(np.abs(old['absorption']))<1e-8
    np.savez(OUT/(tag+'.npz'),grid=grid,**result)
    row=dict(height_m=float(height),mode=mode,cpu_seconds=result['cpu_seconds'],
     native_wall_seconds=result['seconds'],input_seconds=generation,
     wall_seconds=generation+result['seconds'],fo_seconds=result['fo_seconds'],ms_seconds=result['ms_seconds'])
    rows.append(row);h.write_csv(OUT/'timings.csv',rows)
    print(json.dumps(row),flush=True)
 assert len(rows)==60
 print('PASS: all 60 Jacobian runs finite; radiance and absorption reproduce prior results.',flush=True)
if __name__=='__main__':main()
