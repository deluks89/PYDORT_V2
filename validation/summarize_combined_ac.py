"""Verify saved numerical equivalence checks and summarize combined-AC timings."""
from pathlib import Path
import csv
import json
import numpy as np

root=Path(__file__).resolve().parents[1]/'Examples/testing/limb_hybrid'
rows=list(csv.DictReader((root/'combined_final/comparison.csv').open()))
errors=[]
for row in rows:
    if row['model']!='ac': continue
    filename=row['case']+'_ac3.npz'
    current=np.load(root/'combined_final'/filename)
    original=np.load(root/'combined_suite'/filename)
    for key in ['radiance','absorption','haze']:
        a=current['model_'+key];b=original['model_'+key]
        error=float(np.max(np.abs(a-b))/max(np.max(np.abs(b)),1e-20))
        assert error<1e-8,(filename,key,error)
        errors.append(error)
    calls=json.loads((root/'combined_final'/(row['case']+'_public_api.json')).read_text())['driver_calls']
    assert len(calls)==1,calls
forward=np.load(root/'combined_forward8/h15000_v91_a90_ac3.npz')
linear=np.load(root/'combined_batch8/h15000_v91_a90_ac3.npz')
np.testing.assert_allclose(forward['model_radiance'],linear['model_radiance'],rtol=1e-12,atol=1e-14)
np.testing.assert_allclose(forward['qs_radiance'],linear['qs_radiance'],rtol=1e-12,atol=1e-14)
# The earlier MS-only prototype is optional historical evidence.
prototype=[]
for row in rows:
    if row['model']!='ac': continue
    path=root/'corrected_suite'/(row['case']+'_p5.npz')
    if not path.exists(): continue
    old=np.load(path); combined=np.load(root/'combined_suite'/(row['case']+'_ac5.npz'))
    for key in ['radiance','absorption','haze']:
        expected=old['fo_'+key]+old['ms_'+key]
        error=float(np.max(np.abs(combined['model_'+key]-expected))/max(np.max(np.abs(expected)),1e-20))
        assert error<1e-7,(row['case'],key,error)
        prototype.append(error)
ac=[r for r in rows if r['model']=='ac']
result=dict(geometry_cases=len(ac),single_driver_call_checks=len(ac),
            fo_workspace_optimization_max_relative_change=max(errors),
            forward_matches_linearized=True,
            old_two_executable_prototype_max_relative_change=max(prototype) if prototype else None,
            single_wavelength_cpu_savings_percent=[100*(1-float(r['cpu_seconds'])/float(r['qs_cpu_seconds'])) for r in ac])
(root/'combined_final/verification.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
