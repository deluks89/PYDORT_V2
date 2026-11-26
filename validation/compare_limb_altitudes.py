"""Compare horizontal limb radiance and O4 SCD at 20 default-grid boundaries."""
from pathlib import Path
from contextlib import redirect_stdout
import hashlib, io, json, time
import numpy as np
import validate_limb_hybrid as h
import validate_combined_ac as c

OUT=h.ROOT/'Examples/testing/limb_altitude_comparison'
MODES=[('qs','QS'),('standard','VLIDORT (FO SC only)'),('aircraft_combined','AC (3 points)')]

def plot(rows):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 styles={'qs':('black','o','-'),'standard':('tab:blue','s','--'),'aircraft_combined':('tab:orange','^','-.')}
 fig,axs=plt.subplots(2,3,figsize=(15,9),constrained_layout=True)
 reference=[r for r in rows if r['mode']=='qs']
 for mode,label in MODES:
  rr=[r for r in rows if r['mode']==mode]
  z=np.array([r['height_m']/1000 for r in rr])
  color,marker,ls=styles[mode]
  for ax,key,scale in [(axs[0,0],'radiance',1),(axs[0,1],'o4_scd',1e43),(axs[0,2],'cpu_seconds',1),(axs[1,2],'wall_seconds',1)]:
   ax.plot([r[key]/scale for r in rr],z,color=color,marker=marker,ls=ls,label=label,ms=4)
  if mode!='qs':
   for ax,key in [(axs[1,0],'radiance'),(axs[1,1],'o4_scd')]:
    error=[100*(r[key]/next(q[key] for q in reference if q['height_m']==r['height_m'])-1) for r in rr]
    ax.plot(error,z,color=color,marker=marker,ls=ls,label=label,ms=4)
 for ax,title in zip(axs.flat,['Normalized radiance (sr⁻¹)','O₄ SCD (10⁴³ molecule² cm⁻⁵)','Native process CPU time (s)','Radiance difference from QS (%)','O₄ SCD difference from QS (%)','Input generation + driver wall time (s)']):
  ax.set(xlabel=title,ylabel='Aircraft altitude (km)')
  ax.grid(alpha=.25);ax.legend(handlelength=4,fontsize=8)
 fig.suptitle('Horizontal limb: VZA 90°, SZA 60°, RAA 90°, 477 nm\nDefault 82-layer grid; aircraft exactly on layer boundaries',fontsize=14)
 fig.savefig(OUT/'comparison.png',dpi=180)
 fig.savefig(OUT/'comparison.pdf')
 plt.close(fig)

def make_scene():
 p=h.rt.Pydort()
 p.aerosols.remove_all();p.aerosols.haze=h.rt.OpticalAerosol()
 grid=p.vertical_grid.copy()
 a=p.aerosols.haze
 a.extinction_multiplier.reference_wavelength=477.
 a.extinction_multiplier.set_constant(1.)
 a.profile.values[:]=np.exp(-a.profile.get_layer_middles()/2000)
 a.profile.values[:]*=.1/np.sum(a.profile.values*np.diff(grid)/1000)
 a.asymmetry_parameter.set_constant(.7);a.single_scattering_absorption.set_constant(.05)
 p.wavelengths=[477.]
 p.surface_albedo.set_constant(.2)
 p.viewing_geometry.sza=60;p.viewing_geometry.saa=0;p.viewing_geometry.vaa=90;p.viewing_geometry.vza=90
 p.settings_rt.do_radiance_jacobians=True;p.settings_rt.jacobian_list=['absorption']
 p.settings_rt.n_aircraft_source_points=3
 p.aerosols.run_models(p.wavelengths)
 return p

def main():
 OUT.mkdir(parents=True,exist_ok=True)
 start=time.perf_counter()
 p=make_scene()
 grid=p.vertical_grid.copy()
 a=p.aerosols.haze
 eligible=grid[(grid>=0)&(grid<=15000)]
 heights=np.array([eligible[np.argmin(abs(eligible-t))] for t in np.linspace(0,15000,20)])
 assert len(np.unique(heights))==20 and heights[0]==0 and heights[-1]==15000
 assert np.all(np.isin(heights,grid))
 hashes={n:hashlib.sha256((h.ROOT/'Fortran'/n/'DRIVER.exe').read_bytes()).hexdigest() for n in ['vlidort_qs','vlidort_standard']}
 manifest=dict(heights_m=heights.tolist(),vertical_grid_m=grid.tolist(),layers=len(grid)-1,
  wavelength_nm=477,sza=60,raa=90,vza=90,aerosol_AOD=.1,aerosol_scale_height_m=2000,
  aerosol_SSA=.95,aerosol_g=.7,aerosol_reference_wavelength_nm=float(a.profile_wavelength),
  albedo=.2,streams=p.settings_rt.n_halfspace_streams,delta_m=p.settings_rt.do_deltam_phase_function_scaling,
  fine_points=p.settings_rt.n_full_sphericity_fine_layers,ac_points=3,vector=False,
  gases=p.trace_gases.get_name_list(),hashes=hashes,initialization_seconds=time.perf_counter()-start,
  timing='CPU=driver process user+kernel time; native_wall includes launch delay; wall adds input generation; initialization separate',
  scd='Absorption-Jacobian box AMFs weighted by O4 partial vertical columns; not a spectral DOAS fit')
 (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2))
 print(json.dumps(manifest),flush=True)
 rows=[]
 with h.isolated(hashes) as work:
  for i,height in enumerate(heights):
   p.viewing_geometry.observer_height_agl=float(height)
   # Rotate order to avoid systematically favoring the first/last model.
   modes=MODES[i%3:]+MODES[:i%3]
   for mode,label in modes:
    t=time.perf_counter()
    folder=c.generate(p,work,mode,3)
    generation=time.perf_counter()-t
    tag=f'h{height:g}_{mode}'
    result=c.run(folder,OUT/(tag+'.log'))
    p.results_rt=h.rt.ResultsRT(p)
    p.results_rt._read_raw_results((folder/'OUTPUT').as_posix()+'/')
    scd=float(p.results_rt.get_scds('o4').iloc[0])
    radiance=float(result['radiance'][0])
    manual=float(np.dot(p.get_o4_profile()*np.diff(grid)*100, -result['absorption'][::-1]/radiance))
    np.testing.assert_allclose(scd,manual,rtol=1e-12)
    assert np.isfinite(scd) and scd>0 and radiance>0
    row=dict(height_m=float(height),mode=mode,radiance=radiance,o4_scd=scd,
      cpu_seconds=result['cpu_seconds'],native_wall_seconds=result['seconds'],
      input_seconds=generation,wall_seconds=generation+result['seconds'],
      fo_seconds=result['fo_seconds'],ms_seconds=result['ms_seconds'])
    rows.append(row);h.write_csv(OUT/'comparison.csv',rows)
    np.savez(OUT/(tag+'.npz'),grid=grid,o4_profile=p.get_o4_profile(),o4_scd=scd,**result)
    print(json.dumps(row),flush=True)
 assert len(rows)==60
 plot(rows)
 print('PASS: 60 simulations, 20 exact default-grid boundaries; SCDs verified against public API.',flush=True)

if __name__=='__main__':main()
