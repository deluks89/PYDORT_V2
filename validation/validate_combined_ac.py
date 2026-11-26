"""Validate the single-executable AC backend against QS and a spherical-FO baseline.

The baseline uses exactly the same native spherical FO, with ordinary local
Standard MS; it does not silently substitute QS for unsupported limb rays.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import subprocess
import time
import numpy as np
import validate_limb_hybrid as h
from profile_limb_hybrid import process_cpu


def run(folder, log):
    start = time.perf_counter()
    process = subprocess.Popen([str(folder/'DRIVER.exe')], cwd=folder,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    text, _ = process.communicate(timeout=900)
    elapsed = time.perf_counter()-start
    cpu = process_cpu(process)
    log.write_text(text)
    if process.returncode:
        raise RuntimeError(f'{log}: {text[-3000:]}')
    result = {'radiance':np.loadtxt(folder/'OUTPUT/radiances.out').reshape(-1)}
    for path in (folder/'OUTPUT').glob('jacs_radiance_*.out'):
        result[path.stem.removeprefix('jacs_radiance_')]=np.loadtxt(path).reshape(-1)
    for key, value in result.items():
        assert np.all(np.isfinite(value)), (log, key)
    result.update(seconds=elapsed, cpu_seconds=cpu)
    timings = re.findall(r'AC native seconds FO/MS:\s*([\d.Ee+\-]+)\s+([\d.Ee+\-]+)', text)
    result['fo_seconds'] = sum(float(t[0]) for t in timings)
    result['ms_seconds'] = sum(float(t[1]) for t in timings)
    return result


def generate(p, work, mode, points, ss=False):
    h.generate(p, work, mode, points, ss)
    if mode != 'qs':
        p._stage_standard_backend_input()
    return work/('vlidort_qs' if mode == 'qs' else 'vlidort_standard')


def property_difference(p,work,points,base,output):
    """Perturb HG directly, independently of the exported optical Jacobians."""
    output.mkdir(parents=True,exist_ok=True)
    value=p.aerosols.haze.asymmetry_parameter.value
    p.settings_rt.do_radiance_jacobians=False
    values=[]; step=1e-5
    try:
        for sign in [-1,1]:
            p.aerosols.haze.asymmetry_parameter.value=value+sign*step
            p.aerosols.run_models(p.wavelengths)
            folder=generate(p,work,'aircraft_combined',points)
            values.append(run(folder,output/('property_fd_'+str(sign)+'.log'))['radiance'][0])
    finally:
        p.aerosols.haze.asymmetry_parameter.value=value
        p.settings_rt.do_radiance_jacobians=True
        p.aerosols.run_models(p.wavelengths)
    fd=(values[1]-values[0])/(2*step)
    analytic=np.sum(base['haze_asym0'])
    passed=bool(np.isclose(analytic,fd,rtol=2e-3,atol=1e-7))
    (output/'property_fd.json').write_text(json.dumps(dict(analytic=analytic,fd=fd,passed=passed),indent=2))
    assert passed,(analytic,fd)


def plot(output, rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cases = list(dict.fromkeys(r['case'] for r in rows))
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for label in ['baseline'] + [f'ac{n}' for n in sorted(set(r['points'] for r in rows if r['model']=='ac'))]:
        selected = [r for r in rows if r['label']==label]
        for ax, key in zip(axes, ['radiance_error_pct','haze_peak_error_pct','cpu_seconds']):
            display='Spherical FO + local MS' if label=='baseline' else 'AC, '+label[2:]+' points/leg'
            ax.plot([r['case'] for r in selected], [r[key] for r in selected], 'o-', label=display)
            ax.tick_params(axis='x', labelrotation=65)
            ax.grid(alpha=.25)
    for ax, title in zip(axes, ['Radiance error vs QS (%)','Aerosol Jacobian peak error (%)','Process CPU time (s)']):
        ax.set_title(title)
        ax.legend(handlelength=4)
    axes[2].plot(cases, [next(r['qs_cpu_seconds'] for r in rows if r['case']==c) for c in cases], 'ks--', label='QS')
    axes[2].legend(handlelength=4)
    axes[0].set_yscale('symlog',linthresh=.1)
    axes[1].set_yscale('symlog',linthresh=.1)
    fig.tight_layout(); fig.savefig(output/'comparison.png', dpi=160); plt.close(fig)
    if 'qs_haze' not in np.load(output/(rows[0]['case']+'_'+rows[0]['label']+'.npz')):
        return
    fig, axes = plt.subplots(len(cases), 2, figsize=(11, 3*len(cases)), squeeze=False)
    for i, case in enumerate(cases):
        chosen=[r for r in rows if r['case']==case]
        data=np.load(output/(case+'_'+chosen[0]['label']+'.npz'))
        z=(data['grid'][1:]+data['grid'][:-1])[::-1]/2000
        for ax, key in zip(axes[i], ['absorption','haze']):
            ax.plot(data['qs_'+key][:len(z)],z,'k-',lw=2,label='QS')
            for row in chosen:
                d=np.load(output/(case+'_'+row['label']+'.npz'))
                ax.plot(d['model_'+key][:len(z)],z,'--',label=row['label'])
            ax.set(title=case,xlabel='dI/d'+key,ylabel='Altitude (km)')
            ax.grid(alpha=.25); ax.legend(handlelength=4)
    fig.tight_layout(); fig.savefig(output/'jacobian_profiles.png',dpi=150); plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick',action='store_true')
    parser.add_argument('--vector',action='store_true')
    parser.add_argument('--forward',action='store_true')
    parser.add_argument('--property-jacobian',action='store_true')
    parser.add_argument('--property-fd',action='store_true')
    parser.add_argument('--streams',type=int,default=4)
    parser.add_argument('--delta-m',action='store_true')
    parser.add_argument('--points',type=int,nargs='+',default=[3,5])
    parser.add_argument('--case',type=float,nargs=3,action='append')
    parser.add_argument('--wavelengths',type=int,default=1)
    parser.add_argument('--fd',action='store_true')
    parser.add_argument('--public-api',action='store_true')
    parser.add_argument('--output',type=Path,default=h.ROOT/'Examples/testing/limb_hybrid/combined')
    args=parser.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    hashes={name:hashlib.sha256((h.ROOT/'Fortran'/name/'DRIVER.exe').read_bytes()).hexdigest()
            for name in ['vlidort_qs','vlidort_standard']}
    (args.output/'manifest.json').write_text(json.dumps(dict(arguments=vars(args)|{'output':str(args.output)},
        executable_sha256=hashes,baseline='Native spherical FO + ordinary local Standard MS; horizontal MS uses 89.99999 deg'),indent=2))
    cases=args.case or ([(15000,91,90),(550,90,90)] if args.quick else
        [(15000,89,90),(15000,90,90),(15000,91,90),(15000,93,90),(550,90,90),(550,90.5,90),
         (15000,91,0),(15000,91,180)])
    rows=[]
    for height,vza,raa in cases:
        case=f'h{height:.10g}_v{vza:.10g}_a{raa:.10g}'
        with h.isolated(hashes) as work:
            p=h.scene(height,vza,raa,args.vector)
            p.settings_rt.do_deltam_phase_function_scaling=args.delta_m
            p.settings_rt.n_halfspace_streams=args.streams
            p.settings_rt.do_radiance_jacobians=not args.forward
            if args.property_jacobian:
                p.settings_rt.jacobian_list=['absorption','haze','haze_asym0']
            if args.wavelengths>1:
                p.wavelengths=np.linspace(470,490,args.wavelengths)
                p.surface_albedo.change_wavelengths(p.wavelengths)
                p.aerosols.run_models(p.wavelengths)
            qs=run(generate(p,work,'qs',5),args.output/(case+'_qs.log'))
            for model,points in [('baseline',0)]+[('ac',n) for n in args.points]:
                label='baseline' if model=='baseline' else f'ac{points}'
                mode='standard_native_fo' if model=='baseline' else 'aircraft_combined'
                folder=generate(p,work,mode,max(points,2))
                result=run(folder,args.output/(case+'_'+label+'.log'))
                row=dict(case=case,model=model,label=label,points=points,
                         radiance_error_pct=100*(result['radiance'][0]/qs['radiance'][0]-1),
                         **{key:result[key] for key in ['seconds','cpu_seconds','fo_seconds','ms_seconds']},
                         qs_seconds=qs['seconds'],qs_cpu_seconds=qs['cpu_seconds'])
                for key in p.settings_rt.jacobian_list:
                    row[key+'_peak_error_pct']=(100*np.max(np.abs(result[key]-qs[key]))/np.max(np.abs(qs[key]))
                                                if key in qs else float('nan'))
                rows.append(row)
                np.savez(args.output/(case+'_'+label+'.npz'),grid=p.vertical_grid,
                         **{'qs_'+k:v for k,v in qs.items()},**{'model_'+k:v for k,v in result.items()})
                h.write_csv(args.output/'comparison.csv',rows)
                print(json.dumps(row),flush=True)
            if args.public_api:
                from unittest.mock import patch
                from contextlib import redirect_stdout
                import io
                p.settings_rt.sphericity_correction=0
                p.settings_rt.n_aircraft_source_points=args.points[-1]
                calls=[]
                original=h.rt.run_subprocess
                def counted(executable):
                    calls.append(executable)
                    return original(executable)
                with patch.object(h.rt,'run_subprocess',counted), redirect_stdout(io.StringIO()):
                    p.run_rt(run_aerosol_models=False)
                assert len(calls)==1 and 'vlidort_standard' in calls[0],calls
                public=np.loadtxt(work/'vlidort_standard/OUTPUT/radiances.out').reshape(-1)
                np.testing.assert_allclose(public,result['radiance'],rtol=1e-12)
                (args.output/(case+'_public_api.json')).write_text(json.dumps(dict(driver_calls=calls)))
            if args.property_fd:
                assert args.property_jacobian and args.wavelengths==1
                property_difference(p,work,args.points[-1],result,args.output/(case+'_property_fd'))
            if args.fd:
                assert args.wavelengths==1
                folder=generate(p,work,'aircraft_combined',args.points[-1])
                h.finite_difference(folder,result,args.output/(case+'_fd'))
    plot(args.output,rows)
    print('Combined AC validation complete.',flush=True)

if __name__=='__main__':
    main()
