"""Diagnose AC versus local Standard MS around a horizontal viewing ray."""
from pathlib import Path
import sys, os, hashlib, json, re, argparse, csv
import numpy as np
import validate_limb_hybrid as h
import validate_combined_ac as c
OUT=h.ROOT/'Examples/testing/limb_hybrid/horizon_accuracy_audit'
OUT.mkdir(parents=True,exist_ok=True)
hashes={n:hashlib.sha256((h.ROOT/'Fortran'/n/'DRIVER.exe').read_bytes()).hexdigest() for n in ['vlidort_standard','vlidort_qs']}
(OUT/'manifest.json').write_text(json.dumps(dict(hashes=hashes,points=3,wavelength=477,height_m=550,sza=60,raa=90),indent=2))
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--angles',type=float,nargs='+',default=[89.99,89.99999,90.,90.00001,90.0001,90.001,90.01])
parser.add_argument('--append',action='store_true')
args=parser.parse_args()
rows=[]
if args.append:
 with (OUT/'angle_sweep.csv').open() as stream:rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(stream)]
with h.isolated(hashes) as work:
 p=h.scene(550,90,90,False)
 for angle in args.angles:
  p.viewing_geometry.vza=angle
  results={}
  for mode in ['qs','standard_native_fo','aircraft_combined']:
   folder=c.generate(p,work,mode,3)
   if angle==90 and mode=='qs':
    schema={}
    for f in (folder/'INPUT').rglob('*'):
     if f.is_file():schema[str(f.relative_to(folder/'INPUT'))]=f.read_text()[:3500]
    (OUT/'input_schema.json').write_text(json.dumps(schema,indent=2))
   log=OUT/f'v{angle:.8f}_{mode}.log'
   result=c.run(folder,log);results[mode]=result
   np.savez(OUT/f'v{angle:.8f}_{mode}.npz',grid=p.vertical_grid,**result)
  # The QS full-run log reports its native first order separately.
  text=(OUT/f'v{angle:.8f}_qs.log').read_text()
  fo=float(re.search(r'Radiance first order:\s*([\d.Ee+\-]+)',text)[1])
  row=dict(angle=angle,fo=fo)
  for mode,result in results.items():
   row[mode+'_I']=float(result['radiance'][0])
   row[mode+'_MS']=float(result['radiance'][0])-fo
   row[mode+'_I_error_pct']=100*(result['radiance'][0]/results['qs']['radiance'][0]-1)
   for key in ['absorption','haze']:
    row[mode+'_'+key+'_error_pct']=100*float(np.max(np.abs(result[key]-results['qs'][key]))/np.max(np.abs(results['qs'][key])))
  rows.append(row);h.write_csv(OUT/'angle_sweep.csv',rows);print(json.dumps(row),flush=True)
  if angle==90:
   folder=c.generate(p,work,'aircraft_combined',3,ss=True)
   single=c.run(folder,OUT/'v90_native_fo.log')
   np.testing.assert_allclose(single['radiance'][0],fo,rtol=1e-10,atol=1e-14)
print('ANGLE SWEEP COMPLETE',flush=True)
