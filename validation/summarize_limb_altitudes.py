"""Summarize the saved limb altitude comparison, including the corrected surface endpoint."""
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from compare_limb_altitudes import OUT,MODES,plot

with (OUT/'comparison.csv').open() as f: rows=list(csv.DictReader(f))
for r in rows:
 for k in r:
  if k!='mode':r[k]=float(r[k])
assert len(rows)==60
plot(rows)
summary={}
fig,axs=plt.subplots(1,2,figsize=(10,6),constrained_layout=True)
for mode,label in MODES:
 rr=sorted([r for r in rows if r['mode']==mode],key=lambda r:r['height_m'])
 above=[r for r in rr if r['height_m']>0]
 stats=dict(mean_cpu_s=float(np.mean([r['cpu_seconds'] for r in rr])),
  median_cpu_s=float(np.median([r['cpu_seconds'] for r in rr])),
  mean_wall_s=float(np.mean([r['wall_seconds'] for r in rr])),
  total_cpu_s=sum(r['cpu_seconds'] for r in rr),total_wall_s=sum(r['wall_seconds'] for r in rr))
 for key in ['radiance','o4_scd']:
  err=[]
  for r in rr:
   q=next(q for q in rows if q['mode']=='qs' and q['height_m']==r['height_m'])
   err.append(100*(r[key]/q[key]-1))
  stats[key+'_surface_error_pct']=err[0]
  stats[key+'_max_abs_error_pct']=max(abs(e) for e in err)
  stats[key+'_max_abs_error_above_surface_pct']=max(abs(e) for e in err[1:])
  if mode!='qs':
   ax=axs[0 if key=='radiance' else 1]
   ax.plot(err[1:],[r['height_m']/1000 for r in above],
    's--' if mode=='standard' else '^-.',label=label,ms=4)
 summary[mode]=stats
for ax,title in zip(axs,['Radiance difference from QS (%)','O₄ SCD difference from QS (%)']):
 ax.set(xlabel=title,ylabel='Aircraft altitude (km)')
 ax.axvline(0,color='gray',lw=.7);ax.grid(alpha=.25);ax.legend(handlelength=4)
fig.suptitle('Elevated aircraft: detail of model differences\nSurface endpoint retained in comparison.png; excluded only from this zoom')
fig.savefig(OUT/'differences_above_surface.png',dpi=180)
plt.close(fig)
(OUT/'summary.json').write_text(json.dumps(summary,indent=2))
m=json.loads((OUT/'manifest.json').read_text())
text=['# Limb altitude comparison','',
'VZA 90°, SZA 60°, relative azimuth 90°, 477 nm; original 82-layer default grid.',
'Twenty aircraft heights lie exactly on existing layer boundaries; no grid or altitude offsets.',
'Scalar, four streams, delta-M enabled, four fine points, AC three source points.',
'Aerosol AOD 0.1 at 477 nm, exponential 2 km scale height, SSA 0.95, HG g=0.7; albedo 0.2.',
'Default pressure, temperature and trace-gas profiles, with gas absorption enabled.',
'O4 SCDs use the public get_scds method and were independently checked against the same box-AMF sum.',
'They are Jacobian-derived SCDs, not fitted differential SCDs. Units: molecule^2 cm^-5.',
'Normalized radiance units: sr^-1.','',
'## Corrected surface endpoint','',
'AC now includes surface-tangent rays in the atmospheric limb path, preserving MS at 0 m / VZA 90.',
'The table includes all 20 heights, including the corrected surface endpoint.',
'Near-surface and near-horizon checks against QS are saved in surface_continuity/.',
'All 13 pairs passed; explicit angular and altitude increments also track QS.',
'At ground level, looking below the horizon sees the surface rather than sky: both solvers retain this physical jump.','',
'| Model | Max radiance error (%) | Max O4 SCD error (%) | Mean process CPU (s) | Mean input + native wall (s) |',
'|---|---:|---:|---:|---:|']
for mode,label in MODES:
 s=summary[mode]
 text.append(f"| {label} | {s['radiance_max_abs_error_pct']:.6g} | {s['o4_scd_max_abs_error_pct']:.6g} | {s['mean_cpu_s']:.4f} | {s['mean_wall_s']:.4f} |")
text+=['','Process CPU includes native input/output; native wall time also includes environment launch delay.',
'Input generation is timed separately in comparison.csv and added to plotted end-to-end wall time.',
f"One-time scene initialization: {m['initialization_seconds']:.3f} s (excludes interpreter startup/imports).",
'Each model/altitude was measured once; execution order rotated by altitude. No concurrency was used.',
'Standard is the actual standard backend (sphericity_correction=1), not QS fallback or the research native-FO baseline.',
'Executable hashes, settings and exact grids are recorded in manifest.json.',
'Logs and NPZ arrays retain each native result; comparison.csv contains all 60 measurements.']
(OUT/'README.md').write_text('\n'.join(text)+'\n',encoding='utf8')
print(json.dumps(summary,indent=2))

# Surface-neighborhood diagnostics: separate sky continuity from the physical
# surface/sky jump at a ground-level observer.
folder=OUT/'surface_continuity'
with (folder/'comparison.csv').open() as f: cr=list(csv.DictReader(f))
for r in cr:
 for k in r:
  if k!='mode':r[k]=float(r[k])
errors=[]
for r in cr:
 if r['mode']!='aircraft_combined':continue
 q=next(q for q in cr if q['mode']=='qs' and q['height_m']==r['height_m'] and q['vza']==r['vza'])
 errors.append(dict(height_m=r['height_m'],vza=r['vza'],
  radiance_error_pct=100*(r['radiance']/q['radiance']-1),
  o4_scd_error_pct=100*(r['o4_scd']/q['o4_scd']-1)))
(folder/'summary.json').write_text(json.dumps(dict(geometry_pairs=len(errors),max_radiance_error_pct=max(abs(r['radiance_error_pct']) for r in errors),max_o4_scd_error_pct=max(abs(r['o4_scd_error_pct']) for r in errors)),indent=2))
with (folder/'errors.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(errors[0]));w.writeheader();w.writerows(errors)
fig,axs=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
for ax,key,title in [(axs[0,0],'radiance_error_pct','Horizontal radiance'),(axs[0,1],'o4_scd_error_pct','Horizontal O₄ SCD')]:
 rr=sorted([r for r in errors if r['vza']==90],key=lambda r:r['height_m'])
 ax.plot([r['height_m'] for r in rr],[r[key] for r in rr],'o-')
 ax.set_xscale('symlog',linthresh=.001)
 ax.set(xlabel='Height above surface (m)',ylabel='AC − QS (%)',title=title)
for ax,height in [(axs[1,0],.001),(axs[1,1],0)]:
 for mode,label in [('qs','QS'),('aircraft_combined','AC')]:
  rr=sorted([r for r in cr if r['mode']==mode and r['height_m']==height and abs(r['vza']-90)<.0001],key=lambda r:r['vza'])
  ax.plot([1e5*(r['vza']-90) for r in rr],[r['radiance'] for r in rr],
   'o-' if mode=='qs' else 'x--',label=label)
 ax.set(xlabel='VZA − 90° (10⁻⁵ degrees)',ylabel='Normalized radiance (sr⁻¹)',
  title='1 mm above ground: continuous sky view' if height else 'At ground: physical sky-to-surface jump')
 ax.legend(handlelength=4)
for ax in axs.flat:ax.grid(alpha=.25)
fig.savefig(folder/'continuity.png',dpi=180)
plt.close(fig)
