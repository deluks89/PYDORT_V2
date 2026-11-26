module aircraft_chord_quadrature_m
! Integrate an existing diffuse source solution along a spherical chord.
! No RT solve is performed here. Geometry and angular interpolation are independent
! of the differentiated optical properties; derivatives are of the integrand.
use MSSeg_PARS_m, only: fpk
implicit none
private
public :: chord_rule, make_chord_rule, chord_weights
type chord_rule
 real(fpk),allocatable :: fraction(:),remaining(:),weight(:)
end type
contains
real(fpk) function angular_weight(mu,nodes,which) result(w)
real(fpk),intent(in)::mu,nodes(:)
integer,intent(in)::which
integer::j,k,first,last,count
real(fpk)::x
count=size(nodes);w=0
if(count==1)then
 w=1;return
endif
x=max(nodes(1),min(nodes(count),mu))
j=1
do while(j<count-1)
 if(x<=nodes(j+1))exit
 j=j+1
enddo
first=max(1,min(count-2,j-1));last=min(count,first+2)
if(count==2)then
 first=1;last=2
endif
if(which<first.or.which>last)return
w=1
do k=first,last
 if(k/=which)w=w*(x-nodes(k))/(nodes(which)-nodes(k))
enddo
end function

subroutine make_chord_rule(tau,top,bottom,radius,impact,exit_t,length,nodes,which,kbound,rule)
real(fpk),intent(in)::tau,top,bottom,radius,impact,exit_t,length,nodes(:),kbound
integer,intent(in)::which
type(chord_rule),intent(out)::rule
real(fpk),parameter::gx(8)=[ &
 -0.96028985649753623168_fpk,-0.79666647741362673959_fpk, &
 -0.52553240991632898582_fpk,-0.18343464249564980494_fpk, &
  0.18343464249564980494_fpk, 0.52553240991632898582_fpk, &
  0.79666647741362673959_fpk, 0.96028985649753623168_fpk]
real(fpk),parameter::gw(8)=[ &
 0.10122853629037625915_fpk,0.22238103445337447054_fpk, &
 0.31370664587788728734_fpk,0.36268378337836198297_fpk, &
 0.36268378337836198297_fpk,0.31370664587788728734_fpk, &
 0.22238103445337447054_fpk,0.10122853629037625915_fpk]
real(fpk)::h,rt,tail,cuts(size(nodes)+2),d,t,r,mu,sgn,tmp,a,b,r0,r1,span,mid,half
integer::nc,i,j,k,nseg,npoints,index,pass
h=top-bottom;rt=radius+top
! Neglected tail is below exp(-36); the same integrand derivatives are retained.
tail=min(length,36.0_fpk*h/tau)
cuts(1:2)=[0.0_fpk,tail];nc=2
sgn=1
if(exit_t+0.5_fpk*length<0)sgn=-1
! Split at stencil boundaries so angular interpolation is polynomial on each interval.
if(size(nodes)>3)then
 do i=2,size(nodes)-1
  if(nodes(i)>=1)cycle
  d=sgn*impact*nodes(i)/sqrt(1-nodes(i)**2)-exit_t
  if(d>0.and.d<tail)then
   nc=nc+1;cuts(nc)=d
  endif
 enddo
endif
do i=2,nc
 tmp=cuts(i);j=i-1
 do while(j>=1)
  if(cuts(j)<=tmp)exit
  cuts(j+1)=cuts(j);j=j-1
 enddo
 cuts(j+1)=tmp
enddo
npoints=0
do pass=1,2
 index=0
 do i=1,nc-1
  a=cuts(i);b=cuts(i+1)
  r0=sqrt(impact**2+(exit_t+a)**2);r1=sqrt(impact**2+(exit_t+b)**2)
  span=abs((b-a)*(2*exit_t+a+b)/(r0+r1))/h
  nseg=max(1,ceiling((tau*(b-a)/h+kbound*tau*span)/4.0_fpk))
  if(pass==1)then
   npoints=npoints+8*nseg;cycle
  endif
  do j=1,nseg
   half=(b-a)/(2*nseg);mid=a+(2*j-1)*half
   do k=1,8
    index=index+1;d=mid+half*gx(k);t=exit_t+d
    r=sqrt(impact**2+t*t);mu=abs(t)/r
    rule%fraction(index)=max(0.0_fpk,min(1.0_fpk,((rt-impact)*(rt+impact)-t*t)/((rt+r)*h)))
    rule%remaining(index)=d/h
    rule%weight(index)=half*gw(k)/h*angular_weight(mu,nodes,which)*exp(-tau*d/h)
   enddo
  enddo
 enddo
 if(pass==1)allocate(rule%fraction(npoints),rule%remaining(npoints),rule%weight(npoints))
enddo
end subroutine

subroutine chord_weights(k,tau,rule,wd,wu,kd,ku,td,tu,derivatives)
complex(fpk),intent(in)::k
real(fpk),intent(in)::tau
type(chord_rule),intent(in)::rule
complex(fpk),intent(out)::wd,wu,kd,ku,td,tu
logical,intent(in),optional::derivatives
logical::need_derivatives
complex(fpk)::ed,eu,bd,bu
real(fpk)::f,a,w,x
integer::i
need_derivatives=.true.
if(present(derivatives))need_derivatives=derivatives
wd=0;wu=0;kd=0;ku=0;td=0;tu=0
do i=1,size(rule%weight)
 f=rule%fraction(i);a=rule%remaining(i);w=rule%weight(i);x=tau*f
 if(aimag(k)==0)then
  ed=exp(-real(k,fpk)*x);eu=exp(-real(k,fpk)*(tau-x))
 else
  ed=exp(-k*x);eu=exp(-k*(tau-x))
 endif
 bd=tau*w*ed;bu=tau*w*eu
 wd=wd+bd;wu=wu+bu
 if(.not.need_derivatives)cycle
 kd=kd-x*bd;ku=ku-(tau-x)*bu
 td=td+w*ed*(1-tau*a-k*x)
 tu=tu+w*eu*(1-tau*a-k*(tau-x))
enddo
end subroutine
end module aircraft_chord_quadrature_m
