"""Compare shared-work AC against its previous binary64 executable and QS."""
from pathlib import Path
import argparse, hashlib, json, shutil, os
import numpy as np
import validate_limb_hybrid as h
import validate_combined_ac as c

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--reference',type=Path,default=h.ROOT/'tmp/ac_shared_reference/DRIVER.exe')
 p.add_argument('--output',type=Path,default=h.ROOT/'Examples/testing/limb_hybrid/shared_ms')
 p.add_argument('--case',type=float,nargs=3,action='append')
 p.add_argument('--points',type=int,nargs='+',default=[3,5,9])
 p.add_argument('--wavelengths',type=int,default=8)
 p.add_argument('--streams',type=int,default=4)
 p.add_argument('--fine-points',type=int,default=4)
 p.add_argument('--jacobian-tolerance',type=float,default=1e-5,
  help='Peak-normalized change allowed against the previous binary (eight-stream binary64 check: 5e-5).')
 p.add_argument('--vector',action='store_true');p.add_argument('--delta-m',action='store_true')
 p.add_argument('--property-jacobian',action='store_true')
 p.add_argument('--property-fd',action='store_true')
 p.add_argument('--profile',action='store_true')
 p.add_argument('--allow-model-change',action='store_true',help='Report old/new differences without requiring numerical equivalence.')
 p.add_argument('--forward',action='store_true');p.add_argument('--fd',action='store_true')
 args=p.parse_args()
 if args.profile:os.environ['PYDORT_AC_SHARED_PROFILE']='1'
 args.output.mkdir(parents=True,exist_ok=True)
 hashes={n:hashlib.sha256((h.ROOT/'Fortran'/n/'DRIVER.exe').read_bytes()).hexdigest() for n in ['vlidort_standard','vlidort_qs']}
 (args.output/'manifest.json').write_text(json.dumps(dict(arguments=vars(args)|{'output':str(args.output),'reference':str(args.reference)},
  hashes=hashes,reference_hash=hashlib.sha256(args.reference.read_bytes()).hexdigest()),indent=2))
 rows=[]
 for height,vza,raa in args.case or [(15000,91,90),(550,90,90)]:
  case=f'h{height:.10g}_v{vza:.10g}_a{raa:.10g}'
  with h.isolated(hashes) as work:
   scene=h.scene(height,vza,raa,args.vector)
   scene.settings_rt.n_halfspace_streams=args.streams
   scene.settings_rt.n_full_sphericity_fine_layers=args.fine_points
   if args.property_jacobian:scene.settings_rt.jacobian_list=['absorption','haze','haze_asym0']
   scene.settings_rt.do_deltam_phase_function_scaling=args.delta_m
   scene.settings_rt.do_radiance_jacobians=not args.forward
   if args.wavelengths>1:
    scene.wavelengths=np.linspace(470,490,args.wavelengths)
    scene.surface_albedo.change_wavelengths(scene.wavelengths);scene.aerosols.run_models(scene.wavelengths)
   qs=c.run(c.generate(scene,work,'qs',3),args.output/(case+'_qs.log'))
   np.savez(args.output/(case+'_qs.npz'),grid=scene.vertical_grid,**qs)
   for points in args.points:
    results={}
    for name,exe in [('previous',args.reference),('shared',h.ROOT/'Fortran/vlidort_standard/DRIVER.exe')]:
     folder=c.generate(scene,work,'aircraft_combined',points);shutil.copy2(exe,folder/'DRIVER.exe')
     result=c.run(folder,args.output/f'{case}_p{points}_{name}.log');results[name]=result
     np.savez(args.output/f'{case}_p{points}_{name}.npz',grid=scene.vertical_grid,**result)
    row=dict(case=case,points=points,qs_cpu=qs['cpu_seconds'],previous_cpu=results['previous']['cpu_seconds'],
             shared_cpu=results['shared']['cpu_seconds'],previous_ms=results['previous']['ms_seconds'],shared_ms=results['shared']['ms_seconds'])
    for key in ['radiance']+scene.settings_rt.jacobian_list:
     if key not in qs:continue
     row[key+'_change']=float(np.max(np.abs(results['shared'][key]-results['previous'][key]))/max(np.max(np.abs(results['previous'][key])),1e-30))
     for name in results:
      row[name+'_'+key+'_qs_error']=float(np.max(np.abs(results[name][key]-qs[key]))/max(np.max(np.abs(qs[key])),1e-30))
    rows.append(row);h.write_csv(args.output/'comparison.csv',rows);print(json.dumps(row),flush=True)
    for key in ['radiance']+scene.settings_rt.jacobian_list:
     if key+'_change' in row and not args.allow_model_change:
      tolerance=1e-8 if key=='radiance' else args.jacobian_tolerance
      assert row[key+'_change']<tolerance, (case,points,key,row[key+'_change'])
    if args.property_fd:
     assert args.property_jacobian and args.wavelengths==1
     c.property_difference(scene,work,points,results['shared'],args.output/f'{case}_p{points}_property_fd')
    if args.fd:
     assert args.wavelengths==1
     h.finite_difference(folder,results['shared'],args.output/f'{case}_p{points}_fd')
 print('Shared AC comparison complete.',flush=True)
if __name__=='__main__':main()
