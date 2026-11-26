"""Plot and summarize the saved full-profile limb Jacobian comparison."""
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm,TwoSlopeNorm
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import SymmetricalLogLocator
from compare_limb_jacobians import OUT
from compare_limb_altitudes import MODES
import validate_limb_hybrid as h

def main():
 m=json.loads((OUT/'manifest.json').read_text())
 heights=np.array(m['heights_m'],float);grid=np.array(m['vertical_grid_m'],float)/1000
 z=(grid[:-1]+grid[1:])/2
 data={}
 for mode,label in MODES:
  for key in ['absorption','haze']:
   data[mode,key]=np.array([np.load(OUT/f'h{a:g}_{mode}.npz')[key][::-1] for a in heights])
   assert data[mode,key].shape==(20,82)
 metrics=[];profiles=[]
 for key in ['absorption','haze']:
  q=data['qs',key];peak=np.max(abs(q),axis=1)
  for mode,label in MODES:
   values=data[mode,key]
   for i,height in enumerate(heights):
    delta=values[i]-q[i];worst=int(np.argmax(abs(delta)))
    mask=abs(q[i])>=.01*peak[i]
    metrics.append(dict(jacobian=key,mode=mode,height_m=height,
     qs_peak=float(peak[i]),peak_error_pct=float(100*np.max(abs(delta))/peak[i]),
     worst_layer_bottom_m=grid[worst]*1000,worst_layer_top_m=grid[worst+1]*1000,
     rms_error_pct=float(100*np.sqrt(np.mean(delta**2))/peak[i]),
     local_error_where_qs_above_1pct_peak_pct=float(100*np.max(abs(delta[mask]/q[i,mask])))))
    for j in range(len(z)):
     profiles.append(dict(jacobian=key,mode=mode,height_m=height,
      layer_bottom_m=grid[j]*1000,layer_top_m=grid[j+1]*1000,
      value=values[i,j],qs_value=q[i,j],difference=delta[j],peak_normalized_difference_pct=100*delta[j]/peak[i]))
 h.write_csv(OUT/'errors.csv',metrics);h.write_csv(OUT/'all_profiles.csv',profiles)
 with (OUT/'timings.csv').open() as f: timings=list(csv.DictReader(f))
 for r in timings:
  for k in r:
   if k!='mode':r[k]=float(r[k])
 styles={'qs':('black','-'),'standard':('tab:blue','--'),'aircraft_combined':('tab:orange','-.')}
 units={'absorption':'∂(I/E☉)/∂τ_abs (sr⁻¹)','haze':'∂(I/E☉)/∂β_aer (km sr⁻¹)'}
 titles={'absorption':'Absorption optical-depth Jacobian','haze':'Aerosol extinction Jacobian'}
 x=heights/1000;xe=np.r_[x[0]-(x[1]-x[0])/2,(x[1:]+x[:-1])/2,x[-1]+(x[-1]-x[-2])/2]
 for key in ['absorption','haze']:
  q=data['qs',key];peak=np.max(abs(q),axis=1)
  rawmax=max(np.max(abs(data[mode,key])) for mode,label in MODES)
  norm=SymLogNorm(linthresh=rawmax*1e-4,vmin=-rawmax,vmax=rawmax)
  errors={mode:100*(data[mode,key]-q)/peak[:,None] for mode in ['standard','aircraft_combined']}
  common=max(np.max(abs(v)) for v in errors.values())
  acmax=max(np.max(abs(errors['aircraft_combined'])),1e-12)
  for ymax,suffix in [(20,''),(70,'_full_atmosphere')]:
   fig,axs=plt.subplots(2,3,figsize=(15,9),constrained_layout=True)
   for col,(mode,label) in enumerate(MODES):
    im=axs[0,col].pcolormesh(xe,grid,data[mode,key].T,cmap='RdBu_r',norm=norm,shading='flat')
    axs[0,col].set_title(label)
   fig.colorbar(im,ax=list(axs[0,:]),label=units[key]+'; symmetric logarithmic scale',shrink=.85)
   for col,mode in enumerate(['standard','aircraft_combined']):
    im=axs[1,col].pcolormesh(xe,grid,errors[mode].T,cmap='RdBu_r',vmin=-common,vmax=common)
    axs[1,col].set_title(('VLIDORT' if col==0 else 'AC')+' − QS: common scale')
   fig.colorbar(im,ax=list(axs[1,:2]),label='Difference / QS profile peak (%)',shrink=.85)
   im=axs[1,2].pcolormesh(xe,grid,errors['aircraft_combined'].T,cmap='RdBu_r',vmin=-acmax,vmax=acmax)
   axs[1,2].set_title('AC − QS: magnified scale')
   fig.colorbar(im,ax=axs[1,2],label='Difference / QS profile peak (%)',shrink=.85)
   for ax in axs.flat:
    ax.set(xlim=(0,15),ylim=(0,ymax),xlabel='Aircraft altitude (km)',ylabel='Perturbed layer altitude (km)')
    ax.plot(x,x,'k:',lw=.8)
   fig.suptitle(titles[key]+' — 477 nm, VZA 90°, default 82-layer grid\nDotted line: perturbed layer altitude = aircraft altitude',fontsize=14)
   fig.savefig(OUT/(key+'_heatmaps'+suffix+'.png'),dpi=160,bbox_inches='tight')
   plt.close(fig)
 def profile_panel(ax,i,key):
  for mode,label in MODES:
   color,ls=styles[mode]
   ax.stairs(data[mode,key][i],grid,orientation='horizontal',baseline=None,color=color,ls=ls,lw=1.6,label=label)
  peak=np.max(abs(data['qs',key][i]))
  ax.set_xscale('symlog',linthresh=peak*1e-4)
  locator=SymmetricalLogLocator(base=10,linthresh=peak*1e-4)
  locator.set_params(numticks=5)
  lo,hi=ax.get_xlim()
  ticks=[t for t in locator.tick_values(lo,hi) if lo<=t<=hi and (t==0 or abs(t)>=peak*1e-3)]
  ax.set_xticks(ticks)
  ax.tick_params(axis='x',labelsize=9)
  ax.axhline(heights[i]/1000,color='gray',ls=':',lw=1)
  ax.axvline(0,color='gray',lw=.6)
  ax.set(xlabel=units[key],ylabel='Perturbed layer altitude (km)')
  ax.grid(alpha=.2)
  ax.legend(handlelength=4,fontsize=8)
 selected=[0,6,12,19]
 fig,axs=plt.subplots(2,4,figsize=(16,10),constrained_layout=True)
 for col,i in enumerate(selected):
  for row,key in enumerate(['absorption','haze']):
   profile_panel(axs[row,col],i,key);axs[row,col].set_ylim(0,20)
   axs[row,col].set_title(f'{titles[key]}\nAircraft {heights[i]/1000:g} km')
 fig.suptitle('Selected Jacobian profiles; symmetric-logarithmic x axes; dotted horizontal line = observer')
 fig.savefig(OUT/'selected_profiles.png',dpi=180,bbox_inches='tight');plt.close(fig)
 with PdfPages(OUT/'all_profiles.pdf') as pdf:
  for i in range(len(heights)):
   fig,axs=plt.subplots(1,2,figsize=(11,8),constrained_layout=True)
   for ax,key in zip(axs,['absorption','haze']):
    profile_panel(ax,i,key);ax.set_ylim(0,70);ax.set_title(titles[key])
   fig.suptitle(f'Aircraft at {heights[i]:g} m; VZA 90°, 477 nm — full atmosphere')
   pdf.savefig(fig,bbox_inches='tight');plt.close(fig)
 fig,axs=plt.subplots(1,3,figsize=(14,6),constrained_layout=True)
 for mode,label in MODES:
  color,ls=styles[mode]
  if mode!='qs':
   for ax,key in zip(axs[:2],['absorption','haze']):
    rr=[r for r in metrics if r['mode']==mode and r['jacobian']==key]
    ax.plot([r['peak_error_pct'] for r in rr],x,color=color,ls=ls,marker='o',ms=3,label=label)
  rr=sorted([r for r in timings if r['mode']==mode],key=lambda r:r['height_m'])
  axs[2].plot([r['cpu_seconds'] for r in rr],x,color=color,ls=ls,marker='o',ms=3,label=label)
 for ax,title in zip(axs,['Absorption: max |ΔJ| / QS peak (%)','Aerosol: max |ΔJ| / QS peak (%)','Process CPU time (s)']):
  ax.set(xlabel=title,ylabel='Aircraft altitude (km)');ax.grid(alpha=.25);ax.legend(handlelength=4,fontsize=8)
 for ax in axs[:2]:ax.set_xscale('log')
 fig.suptitle('Full-profile Jacobian accuracy and computation time')
 fig.savefig(OUT/'comparison.png',dpi=180,bbox_inches='tight');fig.savefig(OUT/'comparison.pdf');plt.close(fig)
 fig,axs=plt.subplots(1,2,figsize=(10,6),constrained_layout=True)
 for mode,label in MODES:
  rr=sorted([r for r in timings if r['mode']==mode],key=lambda r:r['height_m'])
  color,ls=styles[mode]
  for ax,key in zip(axs,['cpu_seconds','wall_seconds']):
   ax.plot([r[key] for r in rr],x,color=color,ls=ls,marker='o',ms=3,label=label)
 for ax,title in zip(axs,['Native process CPU time (s)','Input generation + native wall time (s)']):
  ax.set(xlabel=title,ylabel='Aircraft altitude (km)');ax.grid(alpha=.25);ax.legend(handlelength=4,fontsize=8)
 fig.suptitle('Radiance with absorption and aerosol extinction Jacobians')
 fig.savefig(OUT/'timings.png',dpi=180,bbox_inches='tight');plt.close(fig)
 summary={}
 for mode,label in MODES:
  rr=[r for r in timings if r['mode']==mode]
  summary[mode]=dict(mean_cpu_seconds=float(np.mean([r['cpu_seconds'] for r in rr])),
   mean_wall_seconds=float(np.mean([r['wall_seconds'] for r in rr])))
  for key in ['absorption','haze']:
   r=max([r for r in metrics if r['mode']==mode and r['jacobian']==key],key=lambda r:r['peak_error_pct'])
   summary[mode][key]=r
 (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
 text=['# Limb Jacobian comparison','',
 'Same 20 aircraft heights on original default-grid boundaries as the radiance/O4 comparison.',
 '477 nm, VZA 90, SZA 60, RAA 90, scalar, four streams, delta-M enabled, four fine points; AC three angular source points.',
 'AOD 0.1 at 477 nm, 2 km exponential aerosol profile, SSA 0.95, g=0.7, surface albedo 0.2; default trace gases retained.',
 'Absorption: derivative of normalized radiance with respect to layer absorption optical depth.',
 'Aerosol: derivative with respect to haze extinction coefficient (km^-1), holding its SSA and phase function fixed.',
 'These are derivatives of radiance, not derivatives of O4 SCD or finite perturbation radiance differences.','',
 'Difference normalization is per aircraft height and parameter: 100*(J_model-J_QS)/max_layer(abs(J_QS)).',
 'No division by a vanishing local QS Jacobian is used in the main plots.',
 'errors.csv additionally reports local relative differences only where |J_QS| is at least 1% of its profile peak.',
 'All 82 layers are retained in NPZ files and all_profiles.csv; native NPZ arrays are top-down, CSV and plots bottom-up.',
 'Heatmaps use actual layer boundaries. Raw-value panels share symmetric-log color scales across models.',
 'Difference panels for Standard and AC share a linear color scale; the third panel magnifies AC separately.',
 'Heatmap zooms cover 0–20 km; *_full_atmosphere.png and all_profiles.pdf include the full 0–70 km atmosphere.',
 'The PDF contains one full-profile page per aircraft height. Selected profiles use symmetric-logarithmic x axes.','',
 '| Model | Worst absorption peak error (%) | Worst aerosol peak error (%) | Mean CPU (s) | Mean input + driver wall (s) |',
 '|---|---:|---:|---:|---:|']
 for mode,label in MODES:
  s=summary[mode];text.append(f"| {label} | {s['absorption']['peak_error_pct']:.6g} | {s['haze']['peak_error_pct']:.6g} | {s['mean_cpu_seconds']:.4f} | {s['mean_wall_seconds']:.4f} |")
 text+=['','Timing includes both requested Jacobian families. Process CPU excludes launch delay; wall time includes it.',
 'One run per model/height, rotated execution order. Initialization is recorded separately in manifest.json.',
 'All 60 runs were finite. Radiance and absorption Jacobians reproduced the prior absorption-only run within numerical tolerance.',
 'Production solver code was not changed for this comparison.']
 (OUT/'README.md').write_text('\n'.join(text)+'\n',encoding='utf8')
 print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
