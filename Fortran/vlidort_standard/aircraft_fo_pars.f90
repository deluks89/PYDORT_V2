! Build the shared geometry/FO kernels with the Standard driver's dimensions.
! This executable never links the QS MS solver. Keeping a separate FOMS module
! here avoids allocating QS-sized optical/Jacobian workspaces on every FO call.
module FOMS_pars_m
   use VLIDORT_PARS_m
   implicit none
   integer, parameter :: ffp=fpk, maxmessages=10
   integer, parameter :: maxdisos=MAXSTREAMS, maxlays=MAXLAYERS
   integer, parameter :: maxmoms=MAXMOMENTS_INPUT, maxJacobs=MAX_ATMOSWFS
   integer, parameter :: maxfine=MAXFINELAYERS, maxsegments=2*MAXLAYERS, maxgeoms=1
   character(len=*), parameter :: DWFL1_2='(A,I3,2(A,L1))'
   real(fpk), parameter :: small=1.0e-12_fpk, dtr=pie/180.0_fpk
end module FOMS_pars_m
