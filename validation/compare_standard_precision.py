"""Isolated binary64/128 Standard speed and derivative comparison against QS.

Production sources and executables are never replaced. --build copies sources
and removes only the real-kind promotion flag in the isolated makefile.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import validate_limb_hybrid as h
import validate_combined_ac as c

ROOT = h.ROOT
BUILD = ROOT/'tmp/standard_binary64/Fortran'
OUTPUT = ROOT/'Examples/testing/standard_precision'


def build():
    from build_standard_binary64 import main as build_isolated
    build_isolated()


def peak_error(a,b):
    return float(np.max(np.abs(a-b))/max(np.max(np.abs(b)),1e-30))


def main():
    global OUTPUT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUTPUT)
    parser.add_argument('--case',action='append',choices=['limb_scalar','toa_near_nadir_scalar','toa_off60_vector'])
    parser.add_argument('--standard128-exe',type=Path,default=ROOT/'tmp/standard_binary128/DRIVER.exe',
                        help='Explicit quadruple-precision reference executable')
    parser.add_argument('--build',action='store_true')
    parser.add_argument('--repeats',type=int,default=2)
    parser.add_argument('--wavelengths',type=int,default=16)
    parser.add_argument('--skip-fd',action='store_true')
    args=parser.parse_args()
    OUTPUT=args.output
    if args.build: build()
    OUTPUT.mkdir(parents=True,exist_ok=True)
    executables={'standard64':BUILD/'vlidort_standard/DRIVER.exe',
                 'standard128':args.standard128_exe,
                 'qs':ROOT/'Fortran/vlidort_qs/DRIVER.exe'}
    hashes={name:hashlib.sha256(p.read_bytes()).hexdigest() for name,p in executables.items()}
    production={name:hashlib.sha256((ROOT/'Fortran'/name/'DRIVER.exe').read_bytes()).hexdigest()
                for name in ('vlidort_standard','vlidort_qs')}
    (OUTPUT/'manifest.json').write_text(json.dumps(dict(hashes=hashes,arguments=vars(args)|{'output':str(OUTPUT),'standard128_exe':str(args.standard128_exe)},
        precision='Standard O2 binary64 vs O2 binary128; QS production O3 binary64',
        scene='15 layers, 4 halfspace streams, SZA60, RAA90, albedo0.2, AOD0.1, scale height2km, SSA0.95, g0.7',
        cases='limb: 15km/VZA91 internal native FO baseline; TOA:70km/VZA179 scalar and VZA120 vector, public Standard, delta-M enabled'),indent=2))
    rows=[]; errors=[]
    for case,height,vza,vector,delta,mode in [
        ('limb_scalar',15000,91,False,False,'standard_native_fo'),
        ('toa_near_nadir_scalar',70000,179,False,True,'standard'),
        ('toa_off60_vector',70000,120,True,True,'standard')]:
      if args.case and case not in args.case: continue
      with h.isolated(production) as work:
        p=h.scene(height,vza,90,vector)
        p.settings_rt.do_deltam_phase_function_scaling=delta
        p.wavelengths=np.linspace(470,490,args.wavelengths)
        p.surface_albedo.change_wavelengths(p.wavelengths)
        p.aerosols.run_models(p.wavelengths)
        for linearized in (False,True):
          p.settings_rt.do_radiance_jacobians=linearized
          results={}
          for backend in ('standard128','standard64','qs'):
            folder=c.generate(p,work,'qs' if backend=='qs' else mode,3)
            shutil.copy2(executables[backend],folder/'DRIVER.exe')
            for stale in (folder/'OUTPUT').glob('*.out'): stale.unlink()
            for repeat in range(args.repeats):
              result=c.run(folder,OUTPUT/f'{case}_{linearized}_{backend}_{repeat}.log')
              results[backend]=result
              row=dict(case=case,linearized=linearized,backend=backend,repeat=repeat,
                       wavelengths=args.wavelengths,**{k:result[k] for k in ('cpu_seconds','seconds','fo_seconds','ms_seconds')})
              rows.append(row); h.write_csv(OUTPUT/'timing.csv',rows)
              print(json.dumps(row),flush=True)
            np.savez(OUTPUT/f'{case}_{linearized}_{backend}.npz',grid=p.vertical_grid,**result)
          for key in ('radiance','absorption','haze'):
            if key not in results['qs']: continue
            row=dict(case=case,linearized=linearized,quantity=key,
                     binary64_vs_128_peak_relative=peak_error(results['standard64'][key],results['standard128'][key]),
                     binary64_vs_qs_peak_relative=peak_error(results['standard64'][key],results['qs'][key]),
                     binary128_vs_qs_peak_relative=peak_error(results['standard128'][key],results['qs'][key]))
            errors.append(row); h.write_csv(OUTPUT/'accuracy.csv',errors)
            print(json.dumps(row),flush=True)
    if not args.skip_fd:
      fd_rows=[]
      with h.isolated(production) as work:
        p=h.scene(70000,179,90,False)
        p.settings_rt.do_deltam_phase_function_scaling=True
        for backend in ('standard128','standard64','qs'):
          folder=c.generate(p,work,'qs' if backend=='qs' else 'standard',3)
          shutil.copy2(executables[backend],folder/'DRIVER.exe')
          base=c.run(folder,OUTPUT/f'fd_{backend}_analytic.log')
          config=folder/'INPUT/configuration.inp'
          lines=config.read_text().splitlines(); lines[14]='False ! independent forward FD'
          config.write_text('\n'.join(lines)+'\n')
          path=folder/'INPUT/aerosol_extinction_profiles.inp'
          header='\n'.join(path.read_text().splitlines()[:3])
          values=np.loadtxt(path,skiprows=3,ndmin=2)
          for layer in (3,len(values)-1):
            for step in (1e-6,1e-7):
              def perturb(mult):
                with tempfile.TemporaryDirectory(prefix='fd_',dir=work) as temp:
                  isolated=Path(temp)
                  shutil.copytree(folder,isolated,dirs_exist_ok=True)
                  changed=values.copy(); changed[layer,0]+=mult*step
                  np.savetxt(isolated/'INPUT/aerosol_extinction_profiles.inp',changed,header=header,comments='')
                  return c.run(isolated,OUTPUT/f'fd_{backend}_{layer}_{step}_{mult}.log')['radiance'][0]
              # FD jobs use independent folders; timing measurements above remain sequential.
              with ThreadPoolExecutor(max_workers=3) as pool:
                intensities=list(pool.map(perturb,(1,2,3)))
              fd=(-5*intensities[0]+8*intensities[1]-3*intensities[2])/(2*step)
              analytic=base['haze'][layer]
              row=dict(backend=backend,layer_top_down=layer,step=step,analytic=float(analytic),
                       finite_difference=float(fd),relative_error=float((analytic-fd)/max(abs(fd),1e-12)))
              fd_rows.append(row); h.write_csv(OUTPUT/'finite_difference.csv',fd_rows)
              print(json.dumps(row),flush=True)
    summaries=[]
    for case in dict.fromkeys(r['case'] for r in rows):
      for linearized in (False,True):
        times={b:float(np.median([r['cpu_seconds'] for r in rows
               if r['case']==case and r['linearized']==linearized and r['backend']==b]))
               for b in executables}
        summaries.append(dict(case=case,linearized=linearized,**times,
            qs_over_standard64=times['qs']/times['standard64'],
            standard128_over_64=times['standard128']/times['standard64']))
    (OUTPUT/'summary.json').write_text(json.dumps(summaries,indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4))
    labels=list(dict.fromkeys(r['case'] for r in rows))
    for ax,linearized in zip(axes,(False,True)):
        for i,backend in enumerate(executables):
            ax.bar(np.arange(len(labels))+(i-1)*.25,
                   [next(r[backend] for r in summaries if r['case']==case and r['linearized']==linearized)
                    for case in labels],width=.25,label=backend)
        ax.set_xticks(np.arange(len(labels)),labels,rotation=15)
        ax.set_ylabel('Median process CPU time (seconds)')
        ax.set_title('Radiance + Jacobians' if linearized else 'Radiance only')
        ax.legend();ax.grid(axis='y',alpha=.25)
    fig.tight_layout();fig.savefig(OUTPUT/'timing.png',dpi=160);plt.close(fig)
    assert all(hashlib.sha256(path.read_bytes()).hexdigest()==hashes[key] for key,path in executables.items())
    print('Precision comparison complete; production binaries unchanged.',flush=True)

if __name__=='__main__': main()
