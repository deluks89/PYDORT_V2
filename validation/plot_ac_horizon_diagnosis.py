"""Plot the completed AC horizon diagnosis."""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1]
out=root/'Examples/testing/limb_hybrid/horizon_accuracy_audit'
def table(name):
 with (out/name).open() as f:return [{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]
angles=table('angle_sweep.csv');refine=table('vertical_refinement.csv')
fig,ax=plt.subplots(2,2,figsize=(12,8))
styles=[('aircraft_combined','AC: three points','o-','#bd3039'),('standard_native_fo','Standard MS + spherical FO','s--','#2868a6')]
for mode,label,style,color in styles:
 for j,key in enumerate(['I_error_pct','haze_error_pct']):
  ax[0,j].plot([r['angle']-90 for r in angles],[r[mode+'_'+key] for r in angles],style,color=color,label=label,lw=1.8)
  ax[1,j].plot([r['layer_spacing_m'] for r in refine],[abs(r[mode+'_'+key]) for r in refine],style,color=color,label=label,lw=1.8)
for a in ax[0]:
 a.set_xscale('symlog',linthresh=1e-5)
 a.set_xlabel('VZA minus 90° [degrees]')
 a.axvline(0,color='0.5',ls=':',lw=1)
for a in ax[1]:
 a.set_xscale('log');a.set_yscale('log');a.invert_xaxis()
 a.set_xlabel('Lowest-layer subdivision spacing [m]')
for a in ax[:,0]:a.set_ylabel('Radiance difference from QS [%]')
for a in ax[:,1]:a.set_ylabel('Aerosol Jacobian error [% of profile peak]')
ax[1,0].set_ylabel('Absolute radiance difference from QS [%]')
ax[0,0].set_title('Horizon transition: coarse vertical grid')
ax[0,1].set_title('Jacobian discontinuity at the second-leg transition')
ax[1,0].set_title('VZA = 90°: same optical profiles, finer layers')
ax[1,1].set_title('VZA = 90°: derivatives mapped to original layers')
for a in ax.flat:
 a.grid(alpha=.25,which='both');a.legend(handlelength=4,fontsize=8)
fig.suptitle('AC accuracy diagnosis: observer 550 m, 477 nm, SZA 60°, RAA 90°\nThree angular source points retained throughout',fontsize=13)
fig.tight_layout(rect=(0,0,1,.94))
fig.savefig(out/'horizon_accuracy.png',dpi=170)
print(out/'horizon_accuracy.png')
