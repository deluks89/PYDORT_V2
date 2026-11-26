"""Independent quadrature checks for curved AC source weights and derivatives."""
from pathlib import Path
import os, subprocess, json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'tmp/ac_chord_integrals';OUT.mkdir(parents=True,exist_ok=True)
driver="""program check
use aircraft_chord_quadrature_m
use MSSeg_PARS_m,only:fpk
implicit none
type(chord_rule)::q
real(fpk)::tau,top,bot,b,te,dl,kr,ki,nodes(3)
complex(fpk)::wd,wu,kd,ku,td,tu,forward_down,forward_up
integer::stat,which
do
 read(*,*,iostat=stat)tau,top,bot,b,te,dl,kr,ki,nodes,which
 if(stat/=0)exit
 call make_chord_rule(tau,top,bot,6371.0_fpk,b,te,dl,nodes,which,60.0_fpk,q)
 call chord_weights(cmplx(kr,ki,fpk),tau,q,wd,wu,kd,ku,td,tu)
 write(*,'(12ES26.16E3)')real(wd),aimag(wd),real(wu),aimag(wu),real(kd),aimag(kd), &
 real(ku),aimag(ku),real(td),aimag(td),real(tu),aimag(tu)
 call chord_weights(cmplx(kr,ki,fpk),tau,q,forward_down,forward_up,kd,ku,td,tu,.false.)
 if(forward_down/=wd.or.forward_up/=wu)error stop 'forward/linearized integral mismatch'
enddo
end program
"""
(OUT/'check.f90').write_text(driver)
env=dict(os.environ);env['TEMP']=env['TMP']=str(OUT)
subprocess.run(['C:/mingw64/bin/gfortran.exe','-O2','-static','-I'+str(ROOT/'Fortran/vlidort_standard/mod'),'-J'+str(OUT),
 str(ROOT/'Fortran/vlidort_standard/aircraft_chord_quadrature.f90'),str(OUT/'check.f90'),'-o',str(OUT/'check.exe')],env=env,cwd=OUT,check=True)
cases=[]
r=6371.
for geom in ['tangent_down','tangent_up','vertical','partial']:
 top,bot=2.,0.
 if geom.startswith('tangent'):
  b=r+.55;tmax=np.sqrt((r+top-b)*(r+top+b))
  if geom.endswith('down'):te,dl=0.,tmax
  else:te,dl=-tmax,tmax
 elif geom=='vertical':b=0.;te=r;dl=2.
 else:
  b=r-1.;te=np.sqrt((r+.55-b)*(r+.55+b));dl=np.sqrt((r+2-b)*(r+2+b))-te
 endpoints=[abs(te)/np.hypot(b,te),abs(te+dl)/np.hypot(b,te+dl)]
 nodes=np.linspace(min(endpoints),max(endpoints),3)
 if geom=='vertical':nodes=np.array([.98,.99,1.])
 for tau in [1e-8,.01,1.,10.]:
  for k in [.3,1.2+.2j,20.]:
   for which in [1,2,3]:
    cases.append((tau,top,bot,b,te,dl,complex(k),nodes,which,geom))
data=''.join(' '.join(map(str,[tau,top,bot,b,te,dl,k.real,k.imag,*nodes,which]))+'\n' for tau,top,bot,b,te,dl,k,nodes,which,_ in cases)
result=subprocess.run([str(OUT/'check.exe')],input=data,text=True,capture_output=True,check=True,cwd=OUT)
values=np.loadtxt(result.stdout.splitlines());gx,gw=np.polynomial.legendre.leggauss(768)
maxerror=0.
for row,(tau,top,bot,b,te,dl,k,nodes,which,geom) in zip(values,cases):
 d=.5*dl*(gx+1);weights=.5*dl*gw
 rr=np.hypot(b,te+d);mu=np.abs(te+d)/rr
 f=((r+top)-rr)/(top-bot);a=d/(top-bot)
 basis=np.ones_like(mu)
 for j in range(3):
  if j!=which-1:basis*= (mu-nodes[j])/(nodes[which-1]-nodes[j])
 w=weights/(top-bot)*basis*np.exp(-tau*a)
 ed=np.exp(-k*tau*f);eu=np.exp(-k*tau*(1-f))
 ref=np.array([np.sum(tau*w*ed),np.sum(tau*w*eu),np.sum(-tau*tau*w*f*ed),
  np.sum(-tau*tau*w*(1-f)*eu),np.sum(w*ed*(1-tau*a-k*tau*f)),
  np.sum(w*eu*(1-tau*a-k*tau*(1-f)))])
 got=row[::2]+1j*row[1::2]
 maxerror=max(maxerror,float(np.max(np.abs(got-ref))))
 np.testing.assert_allclose(got,ref,rtol=2e-7,atol=2e-9,err_msg=str((geom,tau,k,which)))
(OUT/'result.json').write_text(json.dumps(dict(cases=len(cases),maximum_absolute_difference=maxerror,passed=True),indent=2))
print('PASS',len(cases),'curved integral / derivative cases; max absolute error',maxerror,flush=True)
