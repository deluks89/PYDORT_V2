"""Check AC/QS at and near the surface horizon, without perturbing the requested geometry."""
import hashlib,json
import numpy as np
from compare_limb_altitudes import make_scene
import validate_limb_hybrid as h
import validate_combined_ac as c

OUT=h.ROOT/'Examples/testing/limb_altitude_comparison/surface_continuity'
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 p=make_scene()
 hashes={n:hashlib.sha256((h.ROOT/'Fortran'/n/'DRIVER.exe').read_bytes()).hexdigest() for n in ['vlidort_qs','vlidort_standard']}
 # Extra heights are intentional continuity diagnostics; the altitude comparison
 # itself still uses only existing layer boundaries on the unmodified default grid.
 cases=[(0,90),(0,89.99999),(0,90.00001),(0.001,90),(0.01,90),(0.1,90),(1,90),
        (0.001,89.99999),(0.001,90.00001),(1,89.999),(1,90.001),(0,180),(1,0)]
 (OUT/'manifest.json').write_text(json.dumps(dict(hashes=hashes,cases=cases,grid=p.vertical_grid.tolist()),indent=2))
 rows=[]
 with h.isolated(hashes) as work:
  for height,vza in cases:
   p.viewing_geometry.observer_height_agl=height;p.viewing_geometry.vza=vza
   for mode in ['qs','aircraft_combined']:
    tag=f'h{height:.10g}_v{vza:.10g}_{mode}'
    folder=c.generate(p,work,mode,3)
    result=c.run(folder,OUT/(tag+'.log'))
    p.results_rt=h.rt.ResultsRT(p);p.results_rt._read_raw_results((folder/'OUTPUT').as_posix()+'/')
    row=dict(height_m=height,vza=vza,mode=mode,radiance=float(result['radiance'][0]),
      o4_scd=float(p.results_rt.get_scds('o4').iloc[0]),cpu_seconds=result['cpu_seconds'])
    rows.append(row);h.write_csv(OUT/'comparison.csv',rows)
    np.savez(OUT/(tag+'.npz'),**result)
    print(json.dumps(row),flush=True)
 check(rows)

def check(rows):
 for r in rows:
  if r['mode']=='qs':continue
  q=next(q for q in rows if q['mode']=='qs' and q['height_m']==r['height_m'] and q['vza']==r['vza'])
  for key in ['radiance','o4_scd']:
   error=abs(r[key]/q[key]-1)
   assert error<1e-4,(r['height_m'],r['vza'],key,error)
 for height,angles in [(0,[89.99999]),(0.001,[89.99999,90.00001]),(1,[89.999,90.001])]:
  for angle in angles:
   for key in ['radiance','o4_scd']:
    def value(mode,h,v):
     return next(r[key] for r in rows if r['mode']==mode and r['height_m']==h and r['vza']==v)
    delta_ac=value('aircraft_combined',height,angle)-value('aircraft_combined',height,90)
    delta_qs=value('qs',height,angle)-value('qs',height,90)
    assert abs(delta_ac-delta_qs)/abs(value('qs',height,90))<1e-5,(height,angle,key)
 for height in [0.001,0.01,0.1,1]:
  for key in ['radiance','o4_scd']:
   def value(mode,h):
    return next(r[key] for r in rows if r['mode']==mode and r['height_m']==h and r['vza']==90)
   assert abs((value('aircraft_combined',height)-value('aircraft_combined',0))-(value('qs',height)-value('qs',0)))/value('qs',0)<1e-5,(height,key)
 print('PASS: angular and altitude increments track QS; no artificial horizon jump.',flush=True)
 print('PASS: all surface/horizon cases agree with QS within 0.01%.',flush=True)
if __name__=='__main__':
 import argparse,csv
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--check-existing',action='store_true')
 if parser.parse_args().check_existing:
  with (OUT/'comparison.csv').open() as f: rows=list(csv.DictReader(f))
  for r in rows:
   for key in r:
    if key!='mode': r[key]=float(r[key])
  check(rows)
 else: main()
