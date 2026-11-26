"""Refine exported layers without resampling optical or thermodynamic profiles."""
from pathlib import Path
import hashlib,json
import numpy as np
import validate_limb_hybrid as h
import validate_combined_ac as c
OUT=h.ROOT/'Examples/testing/limb_hybrid/horizon_accuracy_audit'
OUT.mkdir(parents=True,exist_ok=True)
hashes={n:hashlib.sha256((h.ROOT/'Fortran'/n/'DRIVER.exe').read_bytes()).hexdigest() for n in ['vlidort_standard','vlidort_qs']}
def subdivide(folder,count):
 path=folder/'INPUT/vertical_grid.inp'
 grid=np.loadtxt(path,skiprows=2)
 newgrid=np.r_[grid[:-1],np.linspace(grid[-2],grid[-1],count+1)[1:]]
 parent=np.r_[np.arange(len(grid)-2),np.repeat(len(grid)-2,count)]
 np.savetxt(path,newgrid,header=f'{len(parent)} ! n_layers\nLayer boundary altitudes [m]',comments='')
 for name in ['pressure_profile.inp','temperature_profile.inp','vapor_pressure_profile.inp','gas_concentration_profiles.inp','aerosol_extinction_profiles.inp','jacobian_layer_flags.inp']:
  path=folder/'INPUT'/name
  lines=path.read_text().splitlines()
  data=lines[3:]
  assert len(data)==len(grid)-1,(name,len(data))
  path.write_text('\n'.join([f'{len(parent)} ! n_rows',*lines[1:3],*[data[k] for k in parent]])+'\n')
 return parent
rows=[]
with h.isolated(hashes) as work:
 p=h.scene(550,90,90,False)
 for count in [1,4,16,64]:
  result={}
  for mode in ['qs','standard_native_fo','aircraft_combined']:
   folder=c.generate(p,work,mode,3)
   parent=subdivide(folder,count)
   r=c.run(folder,OUT/f'refine{count}_{mode}.log')
   r['haze_parent']=np.bincount(parent,weights=r['haze'],minlength=15)
   np.savez(OUT/f'refine{count}_{mode}.npz',parent=parent,**r)
   result[mode]=r
  row=dict(subdivisions=count,layer_spacing_m=2000/count)
  for mode,r in result.items():
   row[mode+'_I']=float(r['radiance'][0])
   row[mode+'_I_error_pct']=100*float(r['radiance'][0]/result['qs']['radiance'][0]-1)
   row[mode+'_haze_error_pct']=100*float(np.max(np.abs(r['haze_parent']-result['qs']['haze_parent']))/np.max(np.abs(result['qs']['haze_parent'])))
  rows.append(row);h.write_csv(OUT/'vertical_refinement.csv',rows);print(json.dumps(row),flush=True)
 # Separate angular sampling from vertical resolution.
 for angle in [90.,90.00001]:
  p.viewing_geometry.vza=angle
  for points in [9,33]:
   r=c.run(c.generate(p,work,'aircraft_combined',points),OUT/f'angle{angle:.8f}_points{points}.log')
   np.savez(OUT/f'angle{angle:.8f}_points{points}.npz',**r)
   print('POINTS',angle,points,float(r['radiance'][0]),flush=True)
print('REFINEMENT COMPLETE',flush=True)
