
! ###############################################################
! #                                                             #
! #                     VLIDORT_QS_1p4                          #
! #                      ("VQS_1p4")                            #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                   (Quasi - Spherical)                       #
! #                    -       -                                #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Authors :     Robert. J. D. Spurr (1)                      #
! #                Matt Christi                                 #
! #                                                             #
! #  Address (1) : RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                                                             #
! #  Tel:          (617) 492 1183                               #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  This Version :   VQS_1p4                                   #
! #  Release Date :   30 June 2023                              #
! #                                                             #
! #  Previous VQS Versions under Standard GPL 3.0:              #
! #  --------------------------------------------               #
! #                                                             #
! #      1.0   F90, limited release May    2012                 #
! #      1.1   F90, limited release June   2015                 #
! #      1.2   F90, limited release March  2016                 #
! #      1.3   F90, limited release March  2020                 #
! #                                                             #
! #  Features Summary of Recent VQS Versions:                   #
! #  ---------------------------------------                    #
! #                                                             #
! #      FO + spher-enhanced MS for Stokes vectors     (1.0)    #
! #      Profile Jacobians added                       (1.1)    #
! #      MS code upgraded to VLIDORT 2.7 standards     (1.2)    #
! #      FO/MS codes upgraded to VLIDORT 2.8 standards (1.3)    #
! #                                                             #
! #  Features Summary of This VQS Version                       #
! #  ------------------------------------                       #
! #                                                             #
! #  1.4, released 30 June 2023.                                #
! #    ==> Formal VQS package release                           #
! #    ==> Formal VQS regular and linearized prof. master subs  #
! #    ==> New Fortran Type Structures for I/O                  #
! #    ==> MS computation performance enhancements              #
! #    ==> Source code made thread-safe for OpenMP              #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 1.4 of the VQS_1p4 software library.            #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VQS Copyright (c) 2012-2023.                               #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VQS_1p4 ( Version 1.4 )                    #
! #                                                                 #
! # VQS_1p4 is free software: you can redistribute it               #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VQS_1p4 is distributed in the hope that it will be              #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VQS_1p4        #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ###############################################################
! #                                                             #
! # GEOMETRY WORKHORSE ROUTINES for FOMS2020                    #
! #                                                             #
! #            FindSunAngle_D                                   #
! #            FindSunPaths_D                                   #
! #            FindSunPaths_T                                   #
! #            Geometry_Exceptions                              #
! #            GAULEG_NG                                        #
! #                                                             #
! ###############################################################

!  March 2020 update

!  Based on "Titan-Probe" FOMS code developed by R. Spurr in May 2012
!  This is a recapitulation and upgrade (VLIDORT Version 2.7)
!  All references to "Probe" have been replaced by "Attrex", now FOMS2020

module FOMS_Geometry_Routines_m

use FOMS_pars_m

!  Everything Public

public

contains


subroutine FindSunAngle_D (Radius,SolarDirection,CumAngle,theta,stheta,ctheta,DirSun)

   use FOMS_pars_m, only : ffp, zero, one

!  Find the solar anlge along the LOS path, for given radius and cumulative angle from BOA
!    SolarDirection is defined at BOA, with azimuth relative to the LOS direction.

   implicit none

!  Inputs

   real(ffp) , Intent(in)    :: Radius,SolarDirection(3),CumAngle

!  Outputs

   real(ffp) , Intent(out)   :: theta,stheta,ctheta
   logical   , Intent(InOut) :: DirSun

!  Local

   real(ffp) :: px(3),b

   px(1) = Radius * sin(CumAngle)
   px(2) = zero
   px(3) = Radius * cos(CumAngle)
   b = DOT_PRODUCT(px,SolarDirection)
   ctheta = -b/Radius
   DirSun = (ctheta.ge.zero)
   stheta = sqrt(one-ctheta*ctheta)
   theta  = acos(ctheta)

!  Done

   return
end subroutine FindSunAngle_D

subroutine FindSunPaths_D (Maxlayers,Radstart,Radii,thstart,sthstart,N,sunpaths)

   use FOMS_pars_m, only : ffp, zero, one

!  Sunpaths for the Direct-sun illumination
!  Starting point is Radstart on the LOS path, with solar angle thstart, in layer N

   implicit none

!  inputs

   INTEGER   , Intent(In)   :: maxlayers, N
   real(ffp) , Intent(In)   :: Radstart,Radii(0:maxlayers)
   real(ffp) , Intent(In)   :: thstart,sthstart

!  Output

   real(ffp), Intent(InOut) :: Sunpaths(maxlayers)

!  Local

   integer    :: n1, k
   real(ffp)  :: sth0, th0, sth1, th1, ks1

!  Layer boundary upper

   N1 = N - 1

!  special case

   if ( sthstart .eq. zero ) then
      sunpaths(n) =  radii(N1) - Radstart
      do k = n1, 1, -1
         sunpaths(k) = radii(k-1)-radii(k) 
      enddo
      return
   endif

!  First layer

   sth0 = sthstart
   th0  = thstart
   sth1 = sth0*Radstart/radii(N1)
   th1  = dasin(sth1)
   ks1  = th0-th1
   sunpaths(n) = dsin(ks1)*Radstart/sth1

!  Other layers to TOA

   sth0 = sth1
   th0  = th1
   do k = n1, 1, -1
      sth1 = sth0*radii(k)/radii(k-1)
      th1  = dasin(sth1)
      ks1  = th0-th1
      sunpaths(k) = dsin(ks1)*radii(k)/sth1 
      sth0 = sth1
      th0  = th1
   enddo

!  Done

   return
end subroutine FindSunPaths_D

subroutine FindSunPaths_T (Maxlayers,nlayers,Radstart,Radii,thstart,sthstart,N,sunpaths,NT,fail,message)

   use FOMS_pars_m, only : ffp, zero, one, pie

!  Sunpaths for the Tangent-height illumination
!  Starting point is Radstart on the LOS path, with solar angle thstart, in layer N

   implicit none

!  inputs

   INTEGER   , Intent(In)   :: maxlayers, nlayers, N
   real(ffp) , Intent(In)   :: Radstart,Radii(0:maxlayers)
   real(ffp) , Intent(In)   :: thstart,sthstart

!  Output

   INTEGER   , Intent(InOut)   :: NT
   real(ffp) , Intent(InOut)   :: Sunpaths(maxlayers)

!  Exception handling

   logical      , intent(inout) :: fail
   character*(*), intent(inout) :: message

!  Local

   logical    :: trawl
   integer    :: n1, k
   real(ffp)  :: sth0, th0, sth1, th1, ks1, tanr

!  Layer boundary upper

   N1 = N - 1

!  tangent height

   NT = N
   tanr = sthstart * Radstart

!  Exception handling: tangent height is below ground (No illumination case)

   if ( tanr .le. radii(nlayers) ) then
      fail = .true. ; message = 'Solar beam cut off by Ground' ; return
   endif

!  Find which layer NT contains tangent height

   k = n1 ; trawl = .true.
   do while (k.ge.n1.and.trawl)
      trawl = (radii(k).gt.tanr) ; k = k + 1
   enddo
   nt = k-1 !; write(*,*)n,nt

!  Distances for layers N and below to NT

   if ( nt.gt.n ) then
      th0  = pie - thstart ; sth0 = sthstart
      sth1 = sth0*Radstart/radii(n)
      th1  = dasin(sth1) ; ks1  = th0-th1
      sunpaths(n) = 2.0d0 * dsin(ks1)*Radstart/sth1
      sth0 = sth1 ; th0 = th1
      do k = n+1,nt-1
        sth1 = sth0*radii(k-1)/radii(k)
        th1  = dasin(sth1) ; ks1  = th0-th1
        sunpaths(k) = 2.0d0 * dsin(ks1)*radii(k)/sth0
        sth0 = sth1 ; th0 = th1
      enddo
      sth1 = 1.0d0 ; ks1 = 0.5*pie - th0
      sunpaths(nt) = 2.0d0 * dsin(ks1)*radii(nt-1)
   else if ( nt.eq.n ) then
      sunpaths(n) = - 2.0d0 * Radstart * dcos(thstart)
   endif

!  Rest of layer n up to the upper boundary

   th0 = pie - thstart ; sth0 = sthstart
   sth1 = sth0*Radstart/radii(N1)
   th1  = dasin(sth1) ; ks1  = th0-th1
   sunpaths(n) = sunpaths(n) + dsin(ks1)*Radstart/sth1
   sth0 = sth1 ; th0 = th1

!  Trawl up from layers above n, to TOA

   do k = n1, 1, -1
      sth1 = sth0*radii(k)/radii(k-1)
      th1  = dasin(sth1)
      ks1  = th0-th1
      sunpaths(k) = dsin(ks1)*radii(k)/sth1 
      sth0 = sth1
      th0  = th1
   enddo

!  Done

   return
end subroutine FindSunPaths_T

SUBROUTINE GAULEG_NG(X1,X2,X,W,N,NMAX)

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input/Output

      INTEGER  , intent(in)  :: N,NMAX
      REAL(ffp), intent(in)  :: X1, X2
      REAL(ffp), intent(out) :: X(NMAX),W(NMAX)

      INTEGER     :: I, M, J
      REAL(ffp)   :: XM,XL,P1,P2,P3,PP,Z,Z1
      REAL(ffp), PARAMETER :: EPS = 3.0D-14

      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)

      DO I=1,M
            Z=DCOS(3.141592654D0*(I-.25D0)/(N+.5D0))
            Z1 = 0.0d0
            DO WHILE (DABS(Z-Z1).GT.EPS)
                  P1=1.D0
                  P2=0.D0
                  DO J=1,N
                        P3=P2
                        P2=P1
                        P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
                  ENDDO
                  PP=N*(Z*P1-P2)/(Z*Z-1.D0)
                  Z1=Z
                  Z=Z1-P1/PP
            ENDDO
            X(I)=XM-XL*Z
            X(N+1-I)=XM+XL*Z
            W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
            W(N+1-I)=W(I)
      ENDDO
      RETURN
END SUBROUTINE GAULEG_NG

subroutine Geometry_Exceptions ( fail, phrase, trace, status )
   logical, intent(inout) :: fail
   character*(*)          :: phrase
   character*(*)          :: trace
   integer, intent(inout) :: status
   if ( fail ) then
      status = 2
      trace  = adjustl(trim(phrase))
      return
   else
      status = 1
   endif 
end subroutine Geometry_Exceptions

!  End module

end module FOMS_Geometry_Routines_m
