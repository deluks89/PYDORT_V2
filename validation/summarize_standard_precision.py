"""Summarize saved precision benchmarks without rerunning the models."""
from pathlib import Path
import csv
import json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]/'Examples/testing/standard_precision'

def main():
    timings=[];accuracy=[]
    for directory in (ROOT,ROOT/'satellite',ROOT/'batch128'):
        if not (directory/'timing.csv').exists():continue
        rows=list(csv.DictReader((directory/'timing.csv').open()))
        for case,linearized in dict.fromkeys((r['case'],r['linearized']) for r in rows):
            selected=[r for r in rows if r['case']==case and r['linearized']==linearized]
            times={b:float(np.median([float(r['cpu_seconds']) for r in selected if r['backend']==b]))
                   for b in ('standard64','standard128','qs') if any(r['backend']==b for r in selected)}
            if len(times)<3:continue
            timings.append(dict(case=case,linearized=linearized=='True',wavelengths=int(selected[0]['wavelengths']),
                **times,speedup_vs_qs=times['qs']/times['standard64'],precision_speedup=times['standard128']/times['standard64']))
            if linearized!='True':continue
            a=np.load(directory/f'{case}_True_standard64.npz')
            b=np.load(directory/f'{case}_True_standard128.npz')
            grid=a['grid'];alt=((grid[1:]+grid[:-1])/2000)[::-1]
            for key in ('absorption','haze'):
                ref=b[key].reshape(-1,len(alt));value=a[key].reshape(ref.shape)
                diff=np.abs(value-ref);peak=np.max(np.abs(ref))
                relative=diff/np.maximum(np.abs(ref),peak*1e-12)
                index=np.unravel_index(np.argmax(diff),diff.shape)
                significant=np.abs(ref)>.01*peak
                accuracy.append(dict(case=case,wavelengths=int(selected[0]['wavelengths']),quantity=key,
                    peak_relative=float(diff.max()/peak),max_local_relative=float(relative.max()),
                    max_relative_above_one_percent_peak=float(relative[significant].max()),
                    largest_difference_altitude_km=float(alt[index[1]]),largest_difference_layer=int(index[1]),
                    largest_difference_wavelength_index=int(index[0]),reference=float(ref[index]),binary64=float(value[index])))
    report=dict(timings=timings,local_precision_errors=accuracy)
    (ROOT/'combined_summary.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
