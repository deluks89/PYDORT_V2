"""Analyze off-horizontal limb comparisons and empirical 1%/5% applicability."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from compare_limb_angles import OUT,ANGLES
from compare_limb_altitudes import MODES as BASE_MODES
MODES=[(mode, "AC (3 points/leg)" if mode=="aircraft_combined" else label) for mode,label in BASE_MODES]
import validate_limb_hybrid as h

FIELDS=['radiance','o4_scd','absorption','haze']
LABELS=['Radiance','O₄ SCD','Absorption Jacobian','Aerosol extinction Jacobian']
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--metrics-only',action='store_true')
 args=parser.parse_args()
 manifest=json.loads((OUT/'manifest.json').read_text());heights=np.array(manifest['heights_m'],float)
 rows=[json.loads(f.read_text()) for f in OUT.glob('v*_h*_*.json')]
 lookup={(r['vza'],r['height_m'],r['mode']):r for r in rows}
 arrays={}
 for key,r in lookup.items():
  if r['status']=='ok':
   with np.load(OUT/f"v{r['vza']:g}_h{r['height_m']:g}_{r['mode']}.npz") as a:
    arrays[key]={p:a[p][::-1].copy() for p in ['absorption','haze']}
 errors=[];summary=[]
 for angle in ANGLES:
  for height in heights:
   q=lookup.get((angle,height,'qs'))
   for mode,label in MODES:
    r=lookup.get((angle,height,mode))
    if r is None:continue
    entry=dict(vza=angle,height_m=height,mode=mode,status=r['status'],path_type=r['path_type'])
    if r['status']=='ok' and q and q['status']=='ok':
     for field in ['radiance','o4_scd']:
      entry[field+'_error_pct']=100*(r[field]/q[field]-1)
     for field in ['absorption','haze']:
      ref=arrays[angle,height,'qs'][field];v=arrays[angle,height,mode][field];j=int(np.argmax(abs(v-ref)))
      entry[field+'_error_pct']=100*np.max(abs(v-ref))/np.max(abs(ref))
      entry[field+'_worst_layer_bottom_m']=manifest['vertical_grid_m'][j]
      entry[field+'_worst_layer_top_m']=manifest['vertical_grid_m'][j+1]
     entry['all_fields_max_error_pct']=max(abs(entry[f+'_error_pct']) for f in FIELDS)
     entry['observables_max_error_pct']=max(abs(entry[f+'_error_pct']) for f in ['radiance','o4_scd'])
     entry['within_1pct']=entry['all_fields_max_error_pct']<=1
     entry['within_5pct']=entry['all_fields_max_error_pct']<=5
    errors.append(entry)
 for angle in ANGLES:
  for mode,label in MODES:
   rr=[r for r in errors if r['vza']==angle and r['mode']==mode]
   valid=[r for r in rr if 'all_fields_max_error_pct' in r]
   times=[lookup[angle,r['height_m'],mode] for r in valid]
   s=dict(vza=angle,mode=mode,completed=len(rr),valid=len(valid),
     unsupported=sum(r['status']=='unsupported' for r in rr),failed=sum(r['status']=='failed' for r in rr),
     heights_within_1pct=[r['height_m'] for r in valid if r['within_1pct']],
     heights_within_5pct=[r['height_m'] for r in valid if r['within_5pct']],
     observable_heights_within_1pct=[r['height_m'] for r in valid if r['observables_max_error_pct']<=1],
     observable_heights_within_5pct=[r['height_m'] for r in valid if r['observables_max_error_pct']<=5])
   if valid:
    for f in FIELDS:
     worst=max(valid,key=lambda r:abs(r[f+'_error_pct']))
     s[f+'_max_abs_error_pct']=abs(worst[f+'_error_pct']);s[f+'_worst_height_m']=worst['height_m']
     s[f+'_heights_within_1pct']=[r['height_m'] for r in valid if abs(r[f+'_error_pct'])<=1]
     s[f+'_heights_within_5pct']=[r['height_m'] for r in valid if abs(r[f+'_error_pct'])<=5]
    s['mean_cpu_seconds']=float(np.mean([r['cpu_seconds'] for r in times]))
    s['mean_concurrent_wall_seconds']=float(np.mean([r['wall_seconds'] for r in times]))
    s['matched_qs_cpu_ratio']=sum(lookup[angle,r['height_m'],'qs']['cpu_seconds'] for r in valid)/sum(r['cpu_seconds'] for r in times)
   summary.append(s)
 def write_csv(path,records):
  fields=list(dict.fromkeys(k for r in records for k in r))
  with path.open('w',newline='',encoding='utf8') as f:
   w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
   w.writerows({k:json.dumps(v) if isinstance(v,list) else v for k,v in r.items()} for r in records)
 write_csv(OUT/'errors.csv',errors);write_csv(OUT/'measurements.csv',rows);write_csv(OUT/'angle_summary.csv',summary)
 (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
 print(json.dumps(summary,indent=2),flush=True)
 if args.metrics_only:return
 assert len(rows)==len(ANGLES)*len(heights)*3,'Incomplete sweep: use --metrics-only while running.'
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.colors import LogNorm,SymLogNorm
 from matplotlib.backends.backend_pdf import PdfPages
 from matplotlib.ticker import SymmetricalLogLocator
 errlookup={(r['vza'],r['height_m'],r['mode']):r for r in errors}
 x=np.arange(len(ANGLES));xe=np.arange(len(ANGLES)+1)-.5;z=heights/1000
 ze=np.r_[z[0]-(z[1]-z[0])/2,(z[1:]+z[:-1])/2,z[-1]+(z[-1]-z[-2])/2]
 fig,axs=plt.subplots(2,4,figsize=(17,10),constrained_layout=True)
 cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#c9c9c9')
 for row,mode in enumerate(['standard','aircraft_combined']):
  for col,field in enumerate(FIELDS):
   matrix=np.array([[abs(errlookup.get((a,hh,mode),{}).get(field+'_error_pct',np.nan)) for a in ANGLES] for hh in heights])
   vmax=max(5,max(abs(r.get(field+'_error_pct',0)) for r in errors if r['mode']!='qs'))
   im=axs[row,col].pcolormesh(xe,ze,np.maximum(matrix,1e-5),cmap=cmap,norm=LogNorm(vmin=1e-4,vmax=vmax))
   axs[row,col].set_xticks(x,ANGLES);axs[row,col].set(xlabel='VZA (degrees)',ylabel='Aircraft altitude (km)',ylim=(0,15),
     title=('VLIDORT FO SC only' if row==0 else 'AC (3 points/leg)')+'\n'+LABELS[col])
   cb=fig.colorbar(im,ax=axs[row,col],shrink=.85);cb.set_label('Absolute difference from QS (%)')
   cb.set_ticks([v for v in [.0001,.001,.01,.1,1,5,10,100] if v<=vmax])
 fig.suptitle('Off-horizontal comparison: 20 boundary-aligned heights, 477 nm\nGray = unsupported/failed; Jacobian differences normalized by each QS profile peak',fontsize=14)
 fig.savefig(OUT/'comparison.png',dpi=170,bbox_inches='tight');fig.savefig(OUT/'comparison.pdf',bbox_inches='tight');plt.close(fig)
 # Empirical availability/accuracy classes: 0 unavailable, 1 >5%, 2 <=5%, 3 <=1%.
 from matplotlib.colors import ListedColormap,BoundaryNorm
 fig,axs=plt.subplots(1,2,figsize=(11,7),constrained_layout=True)
 for ax,mode in zip(axs,['standard','aircraft_combined']):
  mat=np.array([[0 if 'all_fields_max_error_pct' not in errlookup.get((a,hh,mode),{}) else
    (3 if errlookup[a,hh,mode]['within_1pct'] else (2 if errlookup[a,hh,mode]['within_5pct'] else 1))
    for a in ANGLES] for hh in heights])
  im=ax.pcolormesh(xe,ze,mat,cmap=ListedColormap(['#bdbdbd','#d95f5f','#f2cf63','#60ae78']),norm=BoundaryNorm([-.5,.5,1.5,2.5,3.5],4))
  ax.set_xticks(x,ANGLES);ax.set(xlabel='VZA (degrees)',ylabel='Aircraft altitude (km)',ylim=(0,15),title='VLIDORT FO SC only' if mode=='standard' else 'AC (3 points/leg)')
 cb=fig.colorbar(im,ax=list(axs),ticks=[0,1,2,3],shrink=.7);cb.ax.set_yticklabels(['Unavailable','Any field >5%','All fields ≤5%','All fields ≤1%'])
 fig.suptitle('Empirical applicability for this atmosphere and these sampled geometries\nAll four criteria: radiance, O₄ SCD, absorption and aerosol Jacobians')
 fig.savefig(OUT/'applicability.png',dpi=180,bbox_inches='tight');plt.close(fig)
 fig,axs=plt.subplots(1,2,figsize=(11,7),constrained_layout=True)
 for ax,mode in zip(axs,['standard','aircraft_combined']):
  mat=np.array([[0 if 'observables_max_error_pct' not in errlookup.get((a,hh,mode),{}) else
   (3 if errlookup[a,hh,mode]['observables_max_error_pct']<=1 else (2 if errlookup[a,hh,mode]['observables_max_error_pct']<=5 else 1))
   for a in ANGLES] for hh in heights])
  im=ax.pcolormesh(xe,ze,mat,cmap=ListedColormap(['#bdbdbd','#d95f5f','#f2cf63','#60ae78']),norm=BoundaryNorm([-.5,.5,1.5,2.5,3.5],4))
  ax.set_xticks(x,ANGLES);ax.set(xlabel='VZA (degrees)',ylabel='Aircraft altitude (km)',ylim=(0,15),title='VLIDORT FO SC only' if mode=='standard' else 'AC (3 points/leg)')
 cb=fig.colorbar(im,ax=list(axs),ticks=[0,1,2,3],shrink=.7);cb.ax.set_yticklabels(['Unavailable','I or O₄ >5%','I and O₄ ≤5%','I and O₄ ≤1%'])
 fig.suptitle('Applicability for radiance and O₄ SCD only — Jacobian profile accuracy not required')
 fig.savefig(OUT/'observable_applicability.png',dpi=180,bbox_inches='tight');plt.close(fig)
 styles={'qs':('black','-'),'standard':('tab:blue','--'),'aircraft_combined':('tab:orange','-.')}
 grid=np.array(manifest['vertical_grid_m'])/1000
 units={'absorption':'∂(I/E☉)/∂τ_abs (sr⁻¹)','haze':'∂(I/E☉)/∂β_aer (km sr⁻¹)'}
 for angle in ANGLES:
  folder=OUT/f'vza_{angle}';folder.mkdir(exist_ok=True)
  fig,axs=plt.subplots(2,4,figsize=(16,10),constrained_layout=True)
  for mode,label in MODES:
   color,ls=styles[mode];rr=[lookup[angle,hh,mode] for hh in heights]
   for ax,field,title in zip(axs[0],['radiance','o4_scd','cpu_seconds','wall_seconds'],
     ['Normalized radiance (sr⁻¹)','O₄ SCD (10⁴³ molecule² cm⁻⁵)','Native process CPU (s)','Concurrent wall + export (s)']):
    scale=1e43 if field=='o4_scd' else 1
    ax.plot([r.get(field,np.nan)/scale for r in rr],z,color=color,ls=ls,marker='o',ms=3,label=label)
    ax.set_xlabel(title)
   if mode!='qs':
    for ax,field in zip(axs[1],FIELDS):
     ax.plot([errlookup[angle,hh,mode].get(field+'_error_pct',np.nan) for hh in heights],z,color=color,ls=ls,marker='o',ms=3,label=label)
     ax.set_xlabel(LABELS[FIELDS.index(field)]+' difference (%)')
     ax.set_xscale('symlog',linthresh=.01)
  for ax in axs.flat:ax.set_ylabel('Aircraft altitude (km)');ax.grid(alpha=.2);ax.legend(handlelength=4,fontsize=8)
  fig.suptitle(f'VZA {angle}° — radiance, O₄ SCD, Jacobians and timings\nMissing Standard values denote unsupported turning rays; no QS fallback counted')
  fig.savefig(folder/'comparison.png',dpi=170,bbox_inches='tight');plt.close(fig)
  profile_rows=[]
  for field in ['absorption','haze']:
   data={mode:np.array([arrays.get((angle,hh,mode),{}).get(field,np.full(82,np.nan)) for hh in heights]) for mode,label in MODES}
   q=data['qs'];peaks=np.nanmax(abs(q),axis=1)
   delta={mode:100*(data[mode]-q)/peaks[:,None] for mode in ['standard','aircraft_combined']}
   rawmax=max(np.nanmax(abs(d)) for d in data.values() if np.any(np.isfinite(d)))
   common=max(np.nanmax(abs(d)) for d in delta.values() if np.any(np.isfinite(d)))
   acmax=max(np.nanmax(abs(delta['aircraft_combined'])),1e-12)
   for ymax,suffix in [(20,''),(70,'_full_atmosphere')]:
    fig,axs=plt.subplots(2,3,figsize=(15,9),constrained_layout=True)
    cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('#c9c9c9')
    for col,(mode,label) in enumerate(MODES):
     im=axs[0,col].pcolormesh(ze,grid,data[mode].T,cmap=cm,norm=SymLogNorm(rawmax*1e-4,vmin=-rawmax,vmax=rawmax))
     axs[0,col].set_title(label)
    fig.colorbar(im,ax=list(axs[0]),label=units[field]+' (symlog)',shrink=.8)
    for col,mode in enumerate(['standard','aircraft_combined']):
     im=axs[1,col].pcolormesh(ze,grid,delta[mode].T,cmap=cm,vmin=-common,vmax=common)
     axs[1,col].set_title(('VLIDORT' if col==0 else 'AC')+' − QS: common scale')
    fig.colorbar(im,ax=list(axs[1,:2]),label='Difference / QS peak (%)',shrink=.8)
    im=axs[1,2].pcolormesh(ze,grid,delta['aircraft_combined'].T,cmap=cm,vmin=-acmax,vmax=acmax)
    axs[1,2].set_title('AC − QS: magnified');fig.colorbar(im,ax=axs[1,2],label='Difference / QS peak (%)',shrink=.8)
    for ax in axs.flat:
     ax.set(xlim=(0,15),ylim=(0,ymax),xlabel='Aircraft altitude (km)',ylabel='Perturbed layer altitude (km)')
     ax.plot(z,z,'k:',lw=.7)
     if angle>90:
      tangent=(6371+z)*np.sin(np.deg2rad(angle))-6371
      ax.plot(z,np.where(tangent>=0,tangent,np.nan),'k--',lw=.8)
    fig.suptitle(f'VZA {angle}°: '+('Absorption optical-depth' if field=='absorption' else 'Aerosol extinction')+' Jacobians\nDotted: observer altitude; dashed: above-surface tangent altitude')
    fig.savefig(folder/(field+'_heatmaps'+suffix+'.png'),dpi=150,bbox_inches='tight');plt.close(fig)
   for mode,label in MODES:
    for i,hh in enumerate(heights):
     for j in range(82):
      profile_rows.append(dict(vza=angle,height_m=hh,mode=mode,jacobian=field,layer_bottom_m=1000*grid[j],
        layer_top_m=1000*grid[j+1],value=data[mode][i,j],qs_value=q[i,j]))
  write_csv(folder/'profiles.csv',profile_rows)
  with PdfPages(folder/'all_profiles.pdf') as pdf:
   for i,hh in enumerate(heights):
    fig,axs=plt.subplots(1,2,figsize=(11,8),constrained_layout=True)
    for ax,field in zip(axs,['absorption','haze']):
     for mode,label in MODES:
      if (angle,hh,mode) not in arrays:continue
      color,ls=styles[mode];ax.stairs(arrays[angle,hh,mode][field],grid,orientation='horizontal',baseline=None,color=color,ls=ls,label=label)
     peak=np.max(abs(arrays[angle,hh,'qs'][field]))
     ax.set_xscale('symlog',linthresh=peak*1e-4)
     lo,hi=ax.get_xlim();loc=SymmetricalLogLocator(base=10,linthresh=peak*1e-4);loc.set_params(numticks=5)
     ax.set_xticks([v for v in loc.tick_values(lo,hi) if lo<=v<=hi and (v==0 or abs(v)>=peak*1e-3)])
     ax.set(xlabel=units[field],ylabel='Perturbed layer altitude (km)',ylim=(0,70))
     ax.axhline(hh/1000,color='gray',ls=':')
     tangent=(6371+hh/1000)*np.sin(np.deg2rad(angle))-6371
     if angle>90 and tangent>=0:ax.axhline(tangent,color='gray',ls='--')
     ax.grid(alpha=.2);ax.legend(handlelength=4,fontsize=8)
    fig.suptitle(f'VZA {angle}°, aircraft {hh:g} m — full atmosphere'+('\nStandard unsupported: ray turns above surface' if lookup[angle,hh,'standard']['status']=='unsupported' else ''))
    pdf.savefig(fig,bbox_inches='tight');plt.close(fig)
 # Human-readable applicability report; lists refer only to tested heights.
 text=['# Off-horizontal limb comparison','',
  'Same 82-layer default atmosphere and 20 boundary-aligned aircraft heights as the horizontal comparison.',
  '477 nm; SZA 60°, RAA 90°; AOD 0.1, 2 km aerosol scale height, SSA 0.95, g=0.7, albedo 0.2.',
  'Scalar, four streams, delta-M enabled, four fine points, AC three angular source points per ray leg.',
  'Turning rays use two legs (up to six angular samples total); degenerate legs can use fewer. Curved-path quadrature integrates these solutions without additional MS solves.',
  'No production solver changes were made for the sweep.','',
  '## Reading the errors and applicability map','',
  'Radiance/O4 errors are relative to QS. Jacobian errors use the maximum absolute layer difference divided by the QS profile peak at that aircraft height.',
  'The applicability.png classes require all four quantities to pass; observable_applicability.png requires radiance and O4 only. Both describe this atmosphere, wavelength and sampled geometries only.',
  'Unsupported Standard cells are not substituted with QS. Its public interface would fall back to QS for downward turning rays.',
  'impact_altitude_m is the minimum radius of the infinite viewing line minus Earth radius; it is an actual tangent height only for turning rays.',
  'A failing Standard case is distinct from an unsupported geometry; measurements.csv records both statuses.',
  'CPU times are native process user+kernel time. Solvers for each geometry ran concurrently in separate folders; wall times are not serial benchmarks.',
  'No two QS cases were run concurrently. Matched QS CPU ratios use only heights where the corresponding model and QS both succeeded.','',
  '| VZA | Model | Supported successful cases | Max I (%) | Max O4 (%) | Max absorption J (%) | Max aerosol J (%) | All four ≤1% | All four ≤5% |',
  '|---|---|---:|---:|---:|---:|---:|---:|---:|']
 for s in summary:
  if s['mode']=='qs':continue
  text.append('| '+ ' | '.join([str(s['vza']),'Standard' if s['mode']=='standard' else 'AC',str(s['valid'])+'/20']+
    [f"{s.get(f+'_max_abs_error_pct',float('nan')):.5g}" for f in FIELDS]+[str(len(s['heights_within_1pct'])),str(len(s['heights_within_5pct']))])+' |')
 text+=['','Exact passing heights are listed in angle_summary.csv and summary.json. Full layer profiles and figures are in vza_*/.',
  'comparison.png shows error magnitudes; applicability.png shows the joint threshold result. Per-angle comparison.png includes signed radiance/O4 differences and timings.',
  'A gray cell means no supported successful result; it must not be interpreted as small error.']
 (OUT/'README.md').write_text('\n'.join(text)+'\n',encoding='utf8')
 print('ANGLE COMPARISON PLOTS COMPLETE',flush=True)
if __name__=='__main__':main()
