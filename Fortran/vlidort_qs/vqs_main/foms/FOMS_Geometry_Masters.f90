
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
! # GEOMETRY MASTER ROUTINES for the FOMS                       #
! #                                                             #
! #            FOMS_Geometry_Master (calling:--)                #
! #              FOMS_sphergeom_Initial                         #
! #              FOMS_sphergeom_Quads                           #
! #              FOMS_sphergeom_Solar                           #
! #                                                             #
! #            Write_FOMS_Geometry                              #
! #                                                             #
! ###############################################################

!  Based on "Titan-Probe" FOMS code developed by R. Spurr in May 2012
!  This is a recapitulation and upgrade (VLIDORT Version 2.7)
!  All references to "Probe" have been replaced by "Attrex", now FOMS

!  March 2020 Update

module FOMS_Geometry_Masters_m

use FOMS_pars_m
use FOMS_Geometry_Routines_m

!  note Private/public distinction

public  FOMS_Geometry_Master, Write_FOMS_Geometry
private FOMS_sphergeom_Initial, FOMS_sphergeom_Quads, FOMS_sphergeom_Solar

contains

subroutine FOMS_Geometry_Master &
       ( nlayers, heights, eradius, height_Camera,            & ! input 
         alpha_Camera, theta_Camera, Relaz_Camera,            & ! Input
         nfine_input, doCrit, Acrit, extinction,              & ! Input
         Regime_idx, Exact_condition, Raycon, nsegments,      & ! Output
         seglayer_map, segupdown_map, layermask_map, radii,   & ! Output
         Lospaths, Rsegs, alpha, sina, cosa, cota, nfinedivs, & ! Output
         wfine, rsegsfine, alphafine, csqfine, cotfine,       & ! Output
         sunpaths, ntraverse, sunpathsfine, ntraversefine,    & ! Output
         cosscat, theta, Relaz, thetafine, relazfine,         & ! Output
         fail, status, message, trace )                         ! Output

!  Input:
!  ======

!    * Given a triplet of angles at the Instrument ("Camera"), and given height of
!      the Instrument, the earth radius and a stratification of optically
!      uniform layers with known height levels and layer extinction values.
!      Angles in degress, heights etc in [km].
!    * Also given some control on the fine-level sublayering, and a critical
!      value for the extinction (Acrit).

!  Output:
!  ======

!  Description required, but see comments below

   use FOMS_pars_m

   implicit none

!  Input arguments
!  ===============

!  Layer control

   integer  , intent(in)  :: nlayers

!  Radius + heights, Camera height

   real(ffp), intent(in)  :: eradius, heights (0:maxlays)
   real(ffp), intent(in)  :: Height_Camera

!  input angles at the Camera Instrument (Degrees)

   real(ffp), intent(in)  :: alpha_Camera, theta_Camera, Relaz_Camera

!  template number of fine divisions

   integer  , intent(in)  :: nfine_input

!  flag for using Critical Attenuation + threshold

   logical  , intent(in)  :: doCrit
   real(ffp), intent(in)  :: ACrit

!  Extinctions

   real(ffp), intent(in)  :: Extinction(maxlays)

!  Outputs
!  -------

!  Regime index, number of segments, Exactness condition, Ray constant
!    Regime_idx = 1, Upwelling from the ground
!    Regime_idx = 2, Upwelling through the limb from TOA (tangent height)
!    Regime_idx = 3, Downwelling from TOA
!    Exact condition = Camera is positioned at one of level boundaries

   integer  , intent(out)  :: nsegments, Regime_idx
   real(ffp), intent(out)  :: Raycon
   logical  , intent(out)  :: Exact_condition

!  Segment to layer map, segment up/down switch, Layer-mask (1=up,2=dn)
!    Bookkeeping on Line-of-sight path

   integer  , intent(out)  :: seglayer_map  (maxsegments)
   logical  , intent(out)  :: segupdown_map (maxsegments)
   logical  , intent(out)  :: layermask_map (maxlays,2)

!  Radii, Lospaths

   real(ffp), intent(out)  :: radii    (0:maxlays)
   real(ffp), intent(out)  :: Lospaths (maxsegments)
   real(ffp), intent(out)  :: Rsegs    (0:maxsegments)

!  LOS angles and trig ratios thereof
!    Convention, range always 0 to 90 degrees 

   real(ffp), intent(out)  :: alpha    (0:maxsegments)
   real(ffp), intent(out)  :: cosa     (0:maxsegments)
   real(ffp), intent(out)  :: sina     (0:maxsegments)
   real(ffp), intent(out)  :: cota     (0:maxsegments)

!  Quadrature values of radii and Los-angles.
!   Segment fine-layering

   integer  , intent(out)  :: nfinedivs (maxsegments)
   real(ffp), intent(out)  :: wfine     (maxsegments,maxfine)
   real(ffp), intent(out)  :: csqfine   (maxsegments,maxfine)
   real(ffp), intent(out)  :: cotfine   (maxsegments,maxfine)
   real(ffp), intent(out)  :: alphafine (maxsegments,maxfine)
   real(ffp), intent(out)  :: rsegsfine (maxsegments,maxfine)

!  Solar paths to segment boundaries

   integer  , intent(out)  :: ntraverse  (0:maxsegments)
   real(ffp), intent(out)  :: sunpaths   (0:maxsegments,maxlays)

!  Solar paths for Fine levels

   integer  , intent(out)  :: ntraversefine(maxsegments,maxfine)
   real(ffp), intent(out)  :: sunpathsfine (maxsegments,maxlays,maxfine)

!  scattering angle cosine, segment-boundary and fine-grid angles

   real(ffp), intent(out)  :: cosscat

!  Solar and azimuth angles at segment boundaries and at segment
!  fine-levels (Radians)

   real(ffp), intent(out)  :: theta  (0:maxsegments)
   real(ffp), intent(out)  :: relaz  (0:maxsegments)
   real(ffp), intent(out)  :: thetafine (maxsegments,maxfine)
   real(ffp), intent(out)  :: Relazfine (maxsegments,maxfine)

!  Exception handling

   logical      , intent(out) :: fail
   integer      , intent(out) :: status
   character*(*), intent(out) :: message
   character*(*), intent(out) :: trace

!  Local Critical values

   logical    :: Crit_map    (maxsegments)
   real(ffp)  :: alpha_crits (0:maxsegments)
   real(ffp)  :: Rad_crits   (0:maxsegments)

!  Initial calculations for the LOS path

   call FOMS_sphergeom_Initial &
       ( nlayers, heights, eradius, height_Camera, alpha_Camera, & ! Input
         Regime_idx, Exact_condition, nsegments, seglayer_map,   & ! Output
         segupdown_map, layermask_map, radii, Raycon, Lospaths,  & ! Output
         Rsegs, alpha, sina, cosa, cota )                          ! Output

!  Set up the fine-level quadrature along the LOS path segments.

   call FOMS_sphergeom_Quads &
       ( nsegments, nfine_input, doCrit, Acrit,             & ! input
         seglayer_map, segupdown_map, Raycon,               & ! Input
         Lospaths, rsegs, alpha, sina, cosa, extinction,    & ! Input
         Crit_map, Alpha_Crits, Rad_crits, nfinedivs,       & ! Output
         wfine, rsegsfine, alphafine, csqfine, cotfine )      ! output

!  Obtain Solar paths and angles, for segment boundaries and fine-levels.

   call FOMS_sphergeom_Solar &
       ( nlayers, nsegments, regime_idx, theta_Camera, Relaz_Camera, & ! input
         seglayer_map, segupdown_map, Exact_condition,               & ! Input
         radii, rsegs, alpha, nfinedivs, rsegsfine, alphafine,       & ! Input
         sunpaths, ntraverse, sunpathsfine, ntraversefine,           & ! Output
         cosscat, theta, Relaz, thetafine, relazfine,                & ! Output
         fail, status, message, trace )                                ! Output

!  Finish

   return
end subroutine FOMS_Geometry_Master

!

subroutine Write_FOMS_Geometry &
       ( Filename, Regime_idx, nsegments, eradius, cosscat, & 
         seglayer_map, segupdown_map, Lospaths, Rsegs,      & 
         alpha, theta, Relaz, nfinedivs, rsegsfine,         & 
         sunpaths, ntraverse, sunpathsfine, ntraversefine )  

!  This is just a write-routine for the geometry.
!    All arguments are input.

   use FOMS_pars_m, only : maxsegments,maxlays,maxfine,dtr

   implicit none

!  Input arguments
!  ===============

!  Earth radius

   real(ffp), intent(in)   :: eradius

!  number of segments, regime index

   integer  , intent(in)  :: nsegments, Regime_idx

!  Segment to layer map, segment up/down switch

   integer  , intent(in)  :: seglayer_map  (maxsegments)
   logical  , intent(in)  :: segupdown_map (maxsegments)

!  Segment radii, Lospaths

   real(ffp), intent(in)  :: Lospaths (maxsegments)
   real(ffp), intent(in)  :: Rsegs    (0:maxsegments)

!  LOS angles 

   real(ffp), intent(in)  :: alpha    (0:maxsegments)

!  Quadrature values for segment radii

   integer  , intent(in)  :: nfinedivs (maxsegments)
   real(ffp), intent(in)  :: rsegsfine (maxsegments,maxfine)

!  Solar paths to segment boundaries

   integer  , intent(in)  :: ntraverse  (0:maxsegments)
   real(ffp), intent(in)  :: sunpaths   (0:maxsegments,maxlays)

!  Solar paths for Fine levels

   integer  , intent(in)  :: ntraversefine(maxsegments,maxfine)
   real(ffp), intent(in)  :: sunpathsfine (maxsegments,maxlays,maxfine)

!  scattering angle and associated segment-boundary angles

   real(ffp), intent(in)  :: cosscat
   real(ffp), intent(in)  :: theta  (0:maxsegments)
   real(ffp), intent(in)  :: relaz  (0:maxsegments)

!  Filename for output

   character*(*), intent(in) :: Filename

!  Local

   integer :: s, s1, j, n
   integer, parameter :: lu = 67

   open(lu,file=Trim(Filename),status='replace')

   write(lu,*)'Surface Critical Angle (deg) = ',asin(eradius/rsegs(nsegments))/dtr
   write(lu,*)'Path  Scattering Angle (deg) = ',acos(cosscat)/dtr

!  write-up for LOS-path angles and distances at Segment boundaries

   s = 0; write(lu,22)' Regime Index = ',Regime_idx,'         ; TOA angles = ', alpha(s)/dtr,theta(s)/dtr, Relaz(s)/dtr
   write(lu,'(/a/)')'Segment#  Seg/Layer Seg/Up/dn     Lospath [km]    LosAngle       SunAngle       RelAZM'

   do s = 1, nsegments
      write(lu,23)s, seglayer_map(s), segupdown_map(s), lospaths(s), &
               alpha(s)/dtr,theta(s)/dtr, Relaz(s)/dtr
   enddo
22 format(/a,i2,a,T46,3f15.7)
23 format(2(1x,i4,5x),3x,L2,5x,4f15.7)

!   write(lu,'(/a/)')'Segment#  FineLevel NTraverse  SegHts[km]   Sunpaths [km] ------>'
   write(lu,'(/a/)')' Sg Fn nt  SegHt   Sunpaths [km] ------>'

!  write-up for Solar paths to segement boundaries and fine-levels
!    Rob Fix 9/7/15. Format 24/25 changed to get more SF, formerly F7.2/F5.2, now all F8.3

   do s = 0, nsegments
      write(lu,24)s,'  -',ntraverse(s),rsegs(s)-eradius,(sunpaths(s,n),n=1,ntraverse(s))
      s1 = s + 1
      if (s.ne.nsegments) then
        do j = 1, nfinedivs(s1)
           write(lu,25)s1,j,ntraversefine(s1,j),rsegsfine(s1,j)-eradius,(sunpathsfine(s1,n,j),n=1,ntraversefine(s1,j))
        enddo
      endif
   enddo
24 format(i3,A3,i3,f8.3,2x,200f8.3)
25 format(3i3,f8.3,2x,200f8.3)

   close(lu)

!  Finish

   return
end subroutine Write_FOMS_Geometry

!
 
subroutine FOMS_sphergeom_Initial                                       &
   ( nlayers, heights, eradius, height_Camera, alpha_Camera,                & ! Input
     Regime_idx, Exact_condition, nsegments, seglayer_map, segupdown_map,   & ! output
     layermask_map, radii, Raycon, Lospaths, Rsegs, alpha, sina, cosa, cota ) ! Output

   use FOMS_pars_m, only : ffp, maxlays, maxsegments, pie, dtr, zero, one

   implicit none

!  Completely stand-alone geometry routine for Initial LOS (Line of Sight) path setup

!    starting inputs are - Camera value of VZA (alpha_Camera), in degrees
!                        - height grid, earth radius, Camera height in [km]

!  Convention. Alpha_Camera > 90 is downwelling

!  inputs
!  ------

!  Layer control and altitudes

   integer  , intent(in)   :: nlayers
   real(ffp), intent(in)   :: eradius, heights (0:maxlays)

!  Camera height and Camera angle

   real(ffp), intent(in)   :: height_Camera
   real(ffp), intent(in)   :: alpha_Camera

!  Outputs
!  -------

!  Regime index, number of segments

   integer, intent(out) :: nsegments, Regime_idx
   logical, intent(out) :: Exact_condition

!  Segment to layer map, segment up/down switch, Layer-mask (1=up, 2_dn)

   integer  , intent(out)  :: seglayer_map  (maxsegments)
   logical  , intent(out)  :: segupdown_map (maxsegments)
   logical  , intent(out)  :: layermask_map (maxlays,2)

!  Radii, Ray constant, Lospaths [km]

   real(ffp), intent(out)  :: radii    (0:maxlays)
   real(ffp), intent(out)  :: Raycon
   real(ffp), intent(out)  :: Lospaths (maxsegments)

!  LOS angles and trig ratios thereof
!    Convention, range always 0 to 90 degrees 

   real(ffp), intent(out)  :: Rsegs    (0:maxsegments)
   real(ffp), intent(out)  :: alpha    (0:maxsegments)
   real(ffp), intent(out)  :: cosa     (0:maxsegments)
   real(ffp), intent(out)  :: sina     (0:maxsegments)
   real(ffp), intent(out)  :: cota     (0:maxsegments)

!  Local
!  -----

   logical         :: sud1(maxsegments)
   integer         :: n, nCamera, s, s1, nt, ston1(maxsegments)
   real(ffp)       :: Camera_radius, alpha_R, sa, cum, losp1(maxsegments)
   real(ffp)       :: salpha1, calpha1, alpha_Crit, ct, ct1, xx, yy, ht, hsegs(0:maxsegments)

!  debug output

   integer, parameter  :: du = 39
   !logical, parameter :: do_debug_output = .true.
   logical, parameter :: do_debug_output = .false.

!  Zero output

      nsegments = 0 ; Regime_idx = 0 ; seglayer_map = 0 
      segupdown_map = .false. ; layermask_map = .false.

      Raycon = zero ; radii = zero   ; Lospaths = zero  ; rsegs = zero
      alpha = zero  ; cota = zero    ; cosa = zero      ; sina = zero

!  Radii

      do n = 0, nlayers
        radii(n) = eradius + heights(n)
      enddo
      Camera_radius = eradius + height_Camera
      Exact_condition = .false.

!  Critical angle for Nadir-view upwelling to Camera

      alpha_Crit = dasin(eradius/Camera_radius)/dtr

!  Regime index
!     1 = Start at ground, nadir-view upwelling
!     2 = Start at TOA, limb view through tangent height
!     3 = Start at TOA, nadir-view downwelling

      alpha_R    = alpha_Camera * dtr
      if ( alpha_Camera .ge. 90.0_ffp ) then
         Regime_idx = 3
         alpha_R    = pie - alpha_R
      else if ( alpha_Camera.lt.alpha_Crit ) then
         Regime_idx = 1
      else
         Regime_idx = 2
      endif

!  initial values of Cosine/Sine

      if ( alpha_Camera .eq. 90.0_ffp ) then
         salpha1 = 1.0_ffp
         calpha1 = 0.0_ffp
      else
         salpha1 = dsin(alpha_R)
         calpha1 = dcos(alpha_R)  
      endif

!  Ray constant

      Raycon = salpha1 * Camera_radius

!  Regime 1 LOS geometries; Upwelling from BOA
!  -------------------------------------------

      if ( Regime_idx .eq. 1 ) then

!  Find Camera layer, Set exactness Condition

         nCamera = 0
         do while (radii(nCamera).ge.Camera_radius)
           nCamera = nCamera + 1
         enddo
         if (radii(nCamera-1).eq.Camera_radius ) Exact_condition = .true.
         layermask_map(nCamera:nlayers,1) = .true. ! All up-layers, BOA to NCamera

!  Loop from BOA to probe, assigning LOS angles using Sine rule
 
         s = 0
         do n = nlayers, nCamera, -1
            sa = Raycon  / radii (n)
            sina(s)  = sa         ; cosa(s)  = dsqrt(one-sa*sa)
            alpha(s) = dasin(sa)  ; cota(s)  = cosa(s) / sina(s) ; hsegs(s) = heights(n)
            if (s.gt.0) Lospaths(s) = radii(n)*cosa(s) - radii(n+1)*cosa(s-1)
            s = s + 1 ; seglayer_map(s) = n ; segupdown_map(s) = .true.
         enddo
         sina(s)  = salpha1 ;  cosa(s)  = calpha1
         alpha(s) = alpha_R ; cota(s)  = cosa(s) / sina(s)
         nsegments = s ; hsegs(s) = height_Camera
         Lospaths(s) = Camera_radius*calpha1 - radii(nCamera)*cosa(s-1)

!  debug output

         if ( do_debug_output ) then
            open(du,file='Reg1_debug.dat',status='replace')
            write(du,'(i4,4x,2f12.7)')0,alpha(0)/dtr,heights(nlayers)
            sa = Camera_radius*calpha1 - radii(nlayers)*cosa(0)
            do s = 1, nsegments
               write(du,'(2i4,3f12.7)') s,seglayer_map(s),alpha(s)/dtr,hsegs(s),Lospaths(s)
            enddo
            write(du,*)'Check total path',sa,sum(Lospaths(1:nsegments))
            close(du)
         endif

      endif

!  Regime 2 LOS geometries: Limb path through tangent height
!  ---------------------------------------------------------

!  Start Segment counting from Camera, then invert when done

      if ( Regime_idx .eq. 2 ) then

!  Find Camera layer, Set exactness Condition

         nCamera = 0
         do while (radii(nCamera).ge.Camera_radius)
           nCamera = nCamera + 1
         enddo
         if (radii(nCamera-1).eq.Camera_radius ) Exact_condition = .true.

!  Find tangent layer

         nt = 0
         do while (radii(nt).gt.Raycon)
           nt = nt + 1
         enddo
         ht = Raycon - eradius
         layermask_map(1:nt,2)      = .true.  ! All dn-layers from TOA to NT
         layermask_map(nCamera:nt,1) = .true.  ! All up-layers from NT to NCamera

!  initialize at probe

         ct1 = Camera_radius*calpha1
         s = 0 ; sina(s) = salpha1 ; cosa(s) = calpha1 ; hsegs(0) = height_Camera

!  Camera layer is tangent layer. 2 segments either side of tangent height (1 up, 1 down)

         if ( nt.eq.nCamera ) then
            s = s + 1 ; sina(s) = one ; cosa(s) = zero
            losp1(s) = ct1 ; hsegs(s) = ht
            ston1(s) = nCamera ; sud1(s) = .true.
            s = s + 1
            sina(s) = salpha1
            cosa(s) = calpha1
            losp1(s) = ct1 ; hsegs(s) = heights(nCamera-1)
            ston1(s) = nCamera ; sud1(s) = .false.
         endif

!  Camera layer is above tangent layer

         if ( nt.gt.nCamera ) then

!  ..( Camera layer, Upwelling )

            s = s + 1
            sina(s) = Raycon / radii(nCamera)
            cosa(s) = sqrt(one-sina(s)*sina(s))
            ct = radii(nCamera)*cosa(s)
            losp1(s) = ct1 - ct ; hsegs(s) = heights(nCamera) ; ct1 = ct
            ston1(s) = nCamera ; sud1(s) = .true.

!  ..( Camera layer to tangent layer, intervening levels, Upwelling )

            do n = nCamera + 1, nt - 1
               s = s + 1
               sina(s) = Raycon / radii(n)
               cosa(s) = sqrt(one-sina(s)*sina(s))
               ct = radii(n)*cosa(s)
               losp1(s) = ct1 - ct ; hsegs(s) = heights(n) ; ct1 = ct
               ston1(s) = n ; sud1(s) = .true.
            enddo

!  ..( tangent layer, 2 segments, 1 up 1 down)

            s = s + 1 ;  sina(s) = one ; cosa(s) = zero
            losp1(s) = ct1
            ston1(s) = nt ; sud1(s) = .true. ; hsegs(s) = ht

            s = s + 1
            sina(s) = Raycon / radii(nt-1)
            cosa(s) = sqrt(one-sina(s)*sina(s))
            losp1(s) = ct1
            ston1(s) = nt ; sud1(s) = .false. ; hsegs(s) = heights(nt-1)

         endif

!  Trawl from upper level of Tangent layer to TOA, Downwelling

         do n = nt - 1, 1, -1
            s = s + 1
            sina(s) = Raycon / radii(n-1)
            cosa(s) = sqrt(one-sina(s)*sina(s))
            ct = radii(n-1)*cosa(s)
            losp1(s) = ct - ct1 
            ston1(s) = n ; sud1(s) = .false. ; hsegs(s) = heights(n-1) ; ct1 = ct
         enddo

!  Inversion section
!  -----------------

!  Invert Cos and sin ratios, also segment boundary heights 

         nsegments = s
         do s = 0, nsegments/2
            s1 = nsegments - s
            xx = sina(s)  ; yy = sina(s1)  ; sina(s)  = yy ; sina(s1)  = xx
            xx = cosa(s)  ; yy = cosa(s1)  ; cosa(s)  = yy ; cosa(s1)  = xx
            xx = hsegs(s) ; yy = hsegs(s1) ; hsegs(s) = yy ; hsegs(s1) = xx
         enddo

!  invert maps and lospaths; set alpha values, cotangents

         do s = 1, nsegments
            s1 = nsegments + 1 - s
            segupdown_map(s) = sud1 (s1)
            seglayer_map(s)  = ston1(s1)
            lospaths(s)      = losp1(s1)
            alpha(s-1)       = dasin(sina(s-1))
            cota(s-1)        = cosa(s-1) / sina(s-1)
         enddo
         alpha(nsegments) = alpha_R
         cota(nsegments)  = calpha1 / salpha1

!  Debug output

         if ( do_debug_output ) then
            open(du,file='Reg2_debug.dat',status='replace')
            write(du,'(i4,4x,2f12.7)')0,alpha(0)/dtr,hsegs(0)
            cum = pie - alpha(0) - alpha_R
            sa = sqrt(Camera_radius*Camera_radius+radii(0)*radii(0)-&
                      2.0d0*Camera_radius*radii(0)*cos(cum))
            do s = 1, nsegments-1
               write(du,'(2i4,3f12.7)') s,seglayer_map(s),alpha(s)/dtr,hsegs(s),Lospaths(s)
            enddo
            s = nsegments
            write(du,'(2i4,3f12.7)') s,seglayer_map(s),alpha(s)/dtr,height_Camera,Lospaths(s)
            write(du,*)'Check total path',sa,sum(Lospaths(1:nsegments))
            close(du)
         endif

      endif

!  Regime 3 LOS geometries; Downwelling from TOA
!  ---------------------------------------------

      if ( Regime_idx .eq. 3 ) then

!   Find probe layer, Set exactness Condition

         nCamera = 0
         do while (radii(nCamera).gt.Camera_radius)
           nCamera = nCamera + 1
         enddo
         if (radii(nCamera).eq.Camera_radius ) Exact_condition = .true.
         layermask_map(1:nCamera,2) = .true.  ! All dn-layers from TOA to NCamera

!  Trawl down from BOA to AirCraft

         s = 0
         do n = 1, nCamera
            sa = Raycon  / radii (n-1)
            sina(s)  = sa         ; cosa(s)  = dsqrt(one-sa*sa)
            alpha(s) = dasin(sa)  ; cota(s)  = cosa(s) / sina(s) ; hsegs(s) = heights(n-1)
            if(s.gt.0)Lospaths(s) = radii(n-2)*cosa(s-1) - radii(n-1)*cosa(s)
            s = s + 1 ; seglayer_map(s) = n ; segupdown_map(s) = .false.
         enddo
         sina(s)  = salpha1 ; cosa(s)  = calpha1
         alpha(s) = alpha_R ; cota(s)  = cosa(s) / sina(s)
         nsegments = s ; hsegs(s) = height_Camera
         Lospaths(s) = radii(nCamera-1)*cosa(s-1) - Camera_radius*calpha1

!  debug output

         if ( do_debug_output ) then
            open(du,file='Reg3_debug.dat',status='replace')
            write(du,'(i4,4x,2f12.7)')0,alpha(0)/dtr,heights(0)
            sa = - Camera_radius*calpha1 + radii(0)*cosa(0)
            do s = 1, nsegments
               write(du,'(2i4,3f12.7)') s,seglayer_map(s),alpha(s)/dtr,hsegs(s),Lospaths(s)
            enddo
            write(du,*)'Check total path',sa,sum(Lospaths(1:nsegments))
            close(du)
         endif

      endif

!  Segment Radii (all Regimes)

      do s = 0, nsegments
         Rsegs(s) = eradius + hsegs(s)
      enddo

!  debug stop

      if ( do_debug_output ) Stop 'Sphergeom_initial, geometry master - DEBUG STOP!!'

!  Done

   return
end subroutine FOMS_sphergeom_Initial

!

subroutine FOMS_sphergeom_Quads                           &
       ( nsegments, nfine_input, doCrit, Acrit,             & ! input
         seglayer_map, segupdown_map, Raycon,               & ! Input
         Lospaths, rsegs, alpha, sina, cosa, extinction,    & ! Input
         Crit_map, Alpha_Crits, Rad_crits, nfinedivs,       & ! Output
         wfine, rsegsfine, alphafine, csqfine, cotfine )      ! output

   use FOMS_pars_m, only : ffp, maxlays, maxfine, maxsegments, dtr, zero, one

   implicit none

!  Completely stand-alone geometry routine for LOS path Quadrature

!    starting inputs are - path Segment information, layer extinctions
!                        - Criticality flag and threshold (ACrit, DoCrit)

!    Convention. Alpha_Camera > 90 is downwelling

!  inputs
!  ------

!  number of segments

   integer  , intent(in)   :: nsegments

!  template number of fine divisions

   integer  , intent(in)   :: nfine_input

!  flag for using Critical Attenuation + threshold

   logical  , intent(in)   :: doCrit
   real(ffp), intent(in)   :: ACrit

!  Segment to layer map, segment up/down switch

   integer  , intent(in)   :: seglayer_map  (maxsegments)
   logical  , intent(in)   :: segupdown_map (maxsegments)

!  Extinction, Alphas, Lospaths, Ray constant

   real(ffp), intent(in)   :: raycon
   real(ffp), intent(in)   :: extinction (maxlays)

   real(ffp), intent(in)   :: Lospaths (maxsegments)
   real(ffp), intent(in)   :: Rsegs    (0:maxsegments)
   real(ffp), intent(in)   :: alpha    (0:maxsegments)
   real(ffp), intent(in)   :: sina     (0:maxsegments)
   real(ffp), intent(in)   :: cosa     (0:maxsegments)

!  Outputs
!  -------

!  Critical values

   logical  , intent(out)  :: Crit_map (maxsegments)
   real(ffp), intent(out)  :: alpha_crits    (0:maxsegments)
   real(ffp), intent(out)  :: Rad_crits      (0:maxsegments)

!  Quadrature values

   integer  , intent(out)  :: nfinedivs (maxsegments)
   real(ffp), intent(out)  :: wfine     (maxsegments,maxfine)
   real(ffp), intent(out)  :: csqfine   (maxsegments,maxfine)
   real(ffp), intent(out)  :: cotfine   (maxsegments,maxfine)
   real(ffp), intent(out)  :: alphafine (maxsegments,maxfine)
   real(ffp), intent(out)  :: rsegsfine (maxsegments,maxfine)

!  Local
!  -----

   integer   :: s, n, j
   real(ffp) :: tanac, vsign, dcrit, bcrit, bval, b10, atten, od
   real(ffp) :: tfine(maxfine), afine(maxfine), a1, a2, csfine

!  debug output

   integer, parameter  :: du = 39
   !logical, parameter :: do_debug_output = .true.
   logical, parameter :: do_debug_output = .false.

!  Initialize to no criticality

   Crit_map = .false.
   alpha_Crits(0:nsegments) = alpha(0:nsegments)
   Rad_Crits  (0:nsegments) = Rsegs(0:nsegments)

!  Initialize quadrature values

   nfinedivs = nfine_input ; wfine = zero ; rsegsfine = zero
   csqfine = zero ; cotfine = zero ; alphafine = zero

!  If criticality not flagged, assign nfinedivs and quadratures, then Skip Criticality
!  Rob Fix 9/2/15, bval should not be bigger than the allowed dimensioning (MAXFINE)

   if ( .not. doCrit ) then
      do s = 1, nsegments
         n = seglayer_map(s)
         od = lospaths(s) * extinction(n)
!         Atten = exp ( - od ) ; bval = -log10(atten) 
         Atten = exp ( - od ) ; bval = min(-log10(atten),real(maxfine,ffp))
         nfinedivs(s) = max(nint(bval),nfine_input)
         vsign = one ; if ( segupdown_map(s) ) vsign = -one
         a1 = alpha(s-1) ; a2 = alpha(s)
         call gauleg_ng (a1,a2,tfine,afine,nfinedivs(s),maxfine)
         do j = 1, nfinedivs(s)
            csfine = one / sin(tfine(j))
            rsegsfine(s,j) = raycon * csfine
            alphafine(s,j) = tfine(j)
            wfine(s,j)     = vsign * afine(j)
            cotfine(s,j)   = cos(tfine(j)) * csfine
            csqfine(s,j)   = csfine * csfine
         enddo
      enddo
      go to 67
   endif

!  Criticality examination
!   Rob Fix 9/2/15. Vsign was wrong way round !!

   bcrit = -log(Acrit) ; b10 = -log10(acrit)
   do s = 1, nsegments
      n = seglayer_map(s)
      od = lospaths(s) * extinction(n)
      if ( od .gt. bcrit ) then
         Crit_map(s) = .true.
         dcrit = bcrit / extinction(n)
         vsign = one ; if ( segupdown_map(s) ) vsign = -one   ! Changed sign here 9/2/15
         tanac = vsign * rsegs(s) * sina(s) / ( dcrit +vsign*rsegs(s)*cosa(s) )
         alpha_Crits(s-1) = atan(tanac)
         Rad_Crits(s-1)   = rsegs(s) * sina(s) / sin(alpha_Crits(s-1))
         nfinedivs(s) = max(nint(b10),nfine_input)    
         a1 = alpha_Crits(s-1) ; a2 = alpha(s)
!write(91,*)s,n,a1*180.0d0/acos(-1.0d0),a2*180.0d0/acos(-1.0d0),a2-a1
      else
         Atten = exp ( - od ) ; bval = -log10(atten)
         nfinedivs(s) = max(nint(bval),nfine_input)
         vsign = one ; if ( segupdown_map(s) ) vsign = -one
         a1 = alpha(s-1) ; a2 = alpha(s)
!write(91,*)s,n,a1*180.0d0/acos(-1.0d0),a2*180.0d0/acos(-1.0d0),a2-a1
      endif
!      write(*,*)s,n,od,Rad_Crits(s-1),nfinedivs(s),a1*180.0d0/acos(-1.0d0),a2*180.0d0/acos(-1.0d0)
      call gauleg_ng (a1,a2,tfine,afine,nfinedivs(s),maxfine)
      do j = 1, nfinedivs(s)
         csfine = one / sin(tfine(j))
         rsegsfine(s,j) = raycon * csfine
         alphafine(s,j) = tfine(j)
         wfine(s,j)     = vsign * afine(j)
         cotfine(s,j)   = cos(tfine(j)) * csfine
         csqfine(s,j)   = csfine * csfine
!write(91,*)s,j,afine(j),wfine(s,j),tfine(j)*180.0d0/acos(-1.0d0)
      enddo
   enddo

!  Continuation point

67 continue

!  debug output

   if ( do_debug_output ) then
      Open(du,file='Sphergeom_quads.dat',status='replace')
      do s = 1, nsegments
         write(du,'(2i4,l2,4f12.5)')s,nfinedivs(s),segupdown_map(s),alpha_Crits(s-1)/dtr,(alphafine(s,j)/dtr,j=1,3)
      enddo
      close(du)
   endif

!  debug stop

   if ( do_debug_output ) Stop 'Sphergeom_Quads, geometry master - DEBUG STOP!!'

!  Finish

   return
end subroutine FOMS_sphergeom_Quads

subroutine FOMS_sphergeom_Solar &
       ( nlayers, nsegments, regime_idx, theta_Camera, Relaz_Camera,& ! input
         seglayer_map, segupdown_map, Exact_condition,              & ! Input
         radii, rsegs, alpha, nfinedivs, rsegsfine, alphafine,      & ! Input
         sunpaths, ntraverse, sunpathsfine, ntraversefine,          & ! Output
         cosscat, theta, Relaz, thetafine, relazfine,               & ! Output
         fail, status, message, trace )                               ! Output

   use FOMS_pars_m, only : ffp, maxlays, maxfine, maxsegments, dtr, zero, one, Pie

   implicit none

!  Completely stand-alone geometry routine for LOS path Quadrature

!    starting inputs are - Segment information, layer extinctions
!                        - LOSpath segment radii/angles at boundaries and fine levels

!  NB. There is no criticality assumed for SOLAR paths here.
!      [ Note 8 June 2015 ]

!    Convention. Alpha_Camera > 90 is downwelling

!  inputs
!  ------

!  number of segments, Regime index, exactness condition

   integer  , intent(in)  :: nlayers, nsegments, regime_idx
   logical  , intent(in)  :: Exact_condition
   integer  , intent(in)  :: nfinedivs(maxsegments)

!  sun angles at Camera

   real(ffp), intent(in)  :: theta_Camera
   real(ffp), intent(in)  :: Relaz_Camera

!  Segment to layer map, segment up/down switch

   integer  , intent(in)  :: seglayer_map  (maxsegments)
   logical  , intent(in)  :: segupdown_map (maxsegments)

!  LOSpath segment angles, radii

   real(ffp), intent(in)  :: radii    (0:maxlays)
   real(ffp), intent(in)  :: alpha    (0:maxsegments)
   real(ffp), intent(in)  :: rsegs    (0:maxsegments)

!  Quadrature values (angles/radii)

   real(ffp), intent(in)  :: alphafine (maxsegments,maxfine)
   real(ffp), intent(in)  :: rsegsfine (maxsegments,maxfine)

!  outputs
!  -------

!  Solar paths to segment boundaries

   integer  , intent(out)  :: ntraverse  (0:maxsegments)
   real(ffp), intent(out)  :: sunpaths   (0:maxsegments,maxlays)

!  Solar paths for Fine levels

   integer  , intent(out)  :: ntraversefine(maxsegments,maxfine)
   real(ffp), intent(out)  :: sunpathsfine (maxsegments,maxlays,maxfine)

!  scattering angle cosine, segment-boundary and fine-grid solar/azimuth angles

   real(ffp), intent(out)  :: cosscat
   real(ffp), intent(out)  :: theta  (0:maxsegments)
   real(ffp), intent(out)  :: relaz  (0:maxsegments)
   real(ffp), intent(out)  :: thetafine (maxsegments,maxfine)
   real(ffp), intent(out)  :: Relazfine (maxsegments,maxfine)

!  Exception handling

   logical      , intent(out) :: fail
   integer      , intent(out) :: status
   character*(*), intent(out) :: message
   character*(*), intent(out) :: trace

!  local

   integer       :: s, n, j, nt
   real(ffp)     :: SolarDirection(3), term1, term2, vsign, radius
   real(ffp)     :: salpha_Camera,  calpha_Camera, Relaz_Camera_R,  sRelaz_Camera
   real(ffp)     :: theta_Camera_R, stheta_Camera, ctheta_Camera, cRelaz_Camera
   real(ffp)     :: ctheta, stheta, cRelaz, CumAngle, sthf, cthf
   logical       :: DirSun, dsunf

!  debug output

   integer, parameter :: du = 39
   !logical, parameter :: do_debug_output = .true.
   logical, parameter :: do_debug_output = .false.

!  Initialise output

   ntraverse = 0      ; ntraversefine = 0
   sunpaths  = zero   ; sunpathsfine  = zero

   Relaz     = zero   ; theta     = zero ; cosscat = zero
   Relazfine = zero   ; thetafine = zero

   fail = .false. ; message = ' ' ; trace = ' ' ; status = 0

!  Nominal traverse paths for Full illumination

   if (Regime_idx.eq.1) then
      ntraverse(0) = nlayers
      do s = 1, nsegments
         n = seglayer_map(s)
         ntraverse(s) = n - 1
         do j = 1, nfinedivs(s)
            ntraversefine(s,j) = n
         enddo
      enddo
      if ( .not. Exact_condition ) ntraverse(nsegments) = seglayer_map(nsegments)
   else if (Regime_idx.eq.3) then
      ntraverse(0) = 0
      do s = 1, nsegments
         n = seglayer_map(s)
         ntraverse(s) = n
         do j = 1, nfinedivs(s)
            ntraversefine(s,j) = n
         enddo
      enddo
      if ( .not. Exact_condition ) ntraverse(nsegments) = seglayer_map(nsegments)
   else if (Regime_idx.eq.2) then
      ntraverse(0) = 0
      do s = 1, nsegments
         if ( .not.segupdown_map(s) )then
            n = seglayer_map(s)
            ntraverse(s) = n
            do j = 1, nfinedivs(s)
               ntraversefine(s,j) = n
            enddo
         else
            n = seglayer_map(s)
            ntraverse(s) = n - 1
            do j = 1, nfinedivs(s)
               ntraversefine(s,j) = n
            enddo
         endif
      enddo
      if ( .not. Exact_condition ) ntraverse(nsegments) = seglayer_map(nsegments)
   endif

!   do s = 0, nsegments
!      write(*,*)s,segupdown_map(s),Exact_condition,ntraverse(s),rsegs(s)-6371.0d0
!   enddo

!  Angles at Camera

   salpha_Camera  = dsin(alpha(nsegments))
   calpha_Camera  = dcos(alpha(nsegments))

   theta_Camera_R = theta_Camera * DTR
   stheta_Camera  = dsin(theta_Camera_R)
   ctheta_Camera  = dcos(theta_Camera_R)

   Relaz_Camera_R = Relaz_Camera * dtr
   cRelaz_Camera  = dcos(Relaz_Camera_R)
   sRelaz_Camera  = dsin(Relaz_Camera_R)

!  define Unit solar vector at Camera

   SolarDirection(1) = - stheta_Camera * cRelaz_Camera
   SolarDirection(2) = - stheta_Camera * sRelaz_Camera
   SolarDirection(3) = - ctheta_Camera

!  Cosine of scattering angle at Camera
!    Note sign change for Regime_Idx = 3 (vsign = -1)

   term1 = salpha_Camera * stheta_Camera * cRelaz_Camera
   term2 = calpha_Camera * ctheta_Camera
   vsign = one; if ( Regime_idx .eq. 3 ) vsign = -one
   cosscat = - vsign * term2 + term1 

!  regime 3; start loop over all segments. DOWNWELLING from TOA

   if ( Regime_idx .eq. 3 )  then
      Cumangle =  alpha(nsegments) - alpha(0) ; Radius = rsegs(0)
      call FindSunAngle_D(Radius,SolarDirection,CumAngle,theta(0),stheta,ctheta,Dirsun)
      if ( .not. Dirsun ) then
         call FindSunPaths_T (maxlays,nlayers,Radius,Radii,theta(0),stheta,1, &
                              sunpaths(0,:),NT,fail,message)
         call Geometry_Exceptions(fail,'Solar, Regime 3 at TOA',trace,status); if ( fail ) return
      endif
      do s = 1, nsegments
         n = seglayer_map(s)
         Cumangle =  alpha(nsegments) - alpha(s) ; Radius = rsegs(s)
         call FindSunAngle_D(Radius,SolarDirection,CumAngle,theta(s),stheta,ctheta,Dirsun)
         if ( Dirsun ) then
            call FindSunPaths_D(maxlays,Radius,Radii,theta(s),stheta,n,sunpaths(s,:))
         else
            call FindSunPaths_T (maxlays,nlayers,Radius,Radii,theta(s),stheta,n, &
                                 sunpaths(s,:),NT,fail,message)
            call Geometry_Exceptions(fail,'Solar, Regime 3 down from TOA',trace,status); if ( fail ) return
         endif
         do j = 1, nfinedivs(s)
            CumAngle = alpha(nsegments) - alphafine(s,j) ; Radius = rsegsfine(s,j)
            call FindSunAngle_D(Radius,SolarDirection,CumAngle,thetafine(s,j),sthf,cthf,Dsunf)
            if ( Dsunf ) then
               call FindSunPaths_D(maxlays,Radius,Radii,thetafine(s,j),sthf,n,sunpathsfine(s,:,j))
            else
               call FindSunPaths_T(maxlays,nlayers,Radius,Radii,thetafine(s,j),sthf,n, &
                                   sunpathsfine(s,:,j),NT,fail,message)
               call Geometry_Exceptions(fail,'Solar, Regime 3 fine layers',trace,status); if ( fail ) return
            endif
         enddo
      enddo
   endif

!  regime 1; start loop over all segments. UPWELLING from BOA

   if ( Regime_idx .eq. 1 )  then
      n = nlayers
      Cumangle = alpha(0) - alpha(nsegments); radius = Rsegs(0)
      call FindSunAngle_D(Radius,SolarDirection,CumAngle,theta(0),stheta,ctheta,Dirsun)
      if ( Dirsun ) then
         call FindSunPaths_D(maxlays,Radius,Radii,theta(0),stheta,ntraverse(0),sunpaths(0,:))
      else
         call FindSunPaths_T (maxlays,nlayers,Radius,Radii,theta(0),stheta,n, &
                              sunpaths(0,:),ntraverse(0),fail,message)
         call Geometry_Exceptions(fail,'Solar, Regime 1 at BOA',trace,status); if ( fail ) return
      endif
      do s = 1, nsegments
         n = seglayer_map(s)
         Cumangle = alpha(s) - alpha(nsegments); radius = rsegs(s)
         call FindSunAngle_D(Radius,SolarDirection,CumAngle,theta(s),stheta,ctheta,Dirsun)
         if ( Dirsun ) then
            call FindSunPaths_D(maxlays,Radius,Radii,theta(s),stheta,ntraverse(s),sunpaths(s,:))
         else
            call FindSunPaths_T (maxlays,nlayers,Radius,Radii,theta(s),stheta,n,&
                                 sunpaths(s,:),ntraverse(s),fail,message)
            call Geometry_Exceptions(fail,'Solar, Regime 1 up from BOA',trace,status); if ( fail ) return
         endif
         do j = 1, nfinedivs(s)
            CumAngle = alphafine(s,j) - alpha(nsegments) ; Radius = rsegsfine(s,j)
            call FindSunAngle_D(Radius,SolarDirection,CumAngle,thetafine(s,j),sthf,cthf,Dsunf)
            if ( Dsunf ) then
               call FindSunPaths_D(maxlays,Radius,Radii,thetafine(s,j),sthf,ntraversefine(s,j),sunpathsfine(s,:,j))
            else
               call FindSunPaths_T(maxlays,nlayers,Radius,Radii,thetafine(s,j),sthf,n, &
                                   sunpathsfine(s,:,j),ntraversefine(s,j),fail,message)
               call Geometry_Exceptions(fail,'Solar, Regime 1 fine layers',trace,status); if ( fail ) return
            endif
         enddo
      enddo
   endif

!  regime 2. Upwelling to Camera, Limb path through tangent height

   if ( Regime_idx .eq. 2 )  then
      Cumangle =  Pie - alpha(nsegments) - alpha(0) ; radius = Rsegs(0)
      call FindSunAngle_D(Radius,SolarDirection,CumAngle,theta(0),stheta,ctheta,Dirsun)
      if ( .not. Dirsun ) then
         call FindSunPaths_T (maxlays,nlayers,Radius,Radii,theta(0),stheta,1, &
                              sunpaths(0,:),NT,fail,message)
         call Geometry_Exceptions(fail,'Solar, Regime 2 at TOA',trace,status); if ( fail ) return
      endif
      do s = 1, nsegments
         n = seglayer_map(s) ; radius = rsegs(s)
         if (      segupdown_map(s) ) Cumangle =  alpha(s) - alpha(nsegments) 
         if ( .not.segupdown_map(s) ) Cumangle =  Pie - alpha(nsegments) - alpha(s)
         call FindSunAngle_D(Radius,SolarDirection,CumAngle,theta(s),stheta,ctheta,Dirsun)

         if ( Dirsun ) then
            call FindSunPaths_D(maxlays,Radius,Radii,theta(s),stheta,n,sunpaths(s,:))
         else
            call FindSunPaths_T(maxlays,nlayers,Radius,Radii,theta(s),stheta,n, &
                                sunpaths(s,:),NT,fail,message)
            call Geometry_Exceptions(fail,'Solar, Regime 2 down from TOA',trace,status); if ( fail ) return
         endif
         do j = 1, nfinedivs(s)
            if ( segupdown_map(s)      ) CumAngle = alphafine(s,j) - alpha(nsegments)
            if ( .not.segupdown_map(s) ) Cumangle =  Pie - alpha(nsegments) - alphafine(s,j)
            Radius = rsegsfine(s,j)
            call FindSunAngle_D(Radius,SolarDirection,CumAngle,thetafine(s,j),sthf,cthf,Dsunf)
            if ( Dsunf ) then
               call FindSunPaths_D(maxlays,Radius,Radii,thetafine(s,j),sthf,n,sunpathsfine(s,:,j))
            else
               call FindSunPaths_T(maxlays,nlayers,Radius,Radii,thetafine(s,j),sthf,n, &
                                   sunpathsfine(s,:,j),NT,fail,message)
               call Geometry_Exceptions(fail,'Solar, Regime 2 fine layers',trace,status); if ( fail ) return
            endif
         enddo
      enddo
   endif

!   Azimuth angles. [Simply set using the constant scattering angle]

   do s = 1, nsegments
      if (      segupdown_map(s) ) vsign =  one 
      if ( .not.segupdown_map(s) ) vsign = -one
      crelaz = one
      if ( sin(theta(s-1)).ne.zero ) then
         cRelaz = (cosscat+vsign*cos(theta(s-1))*cos(alpha(s-1))) / &
                              sin(theta(s-1))/sin(alpha(s-1))
      endif
      if ( crelaz.gt.one )  crelaz = one
      if ( crelaz.lt.-one ) crelaz = -one
      if ( relaz_Camera .eq. zero    .and. crelaz .gt. zero ) crelaz = one
      if ( relaz_Camera .eq. 180.0d0 .and. crelaz .lt. zero ) crelaz = -one
      Relaz(s-1) = acos(cRelaz)
      do j = 1, nfinedivs(s)
         cRelaz = (cosscat+vsign*cos(thetafine(s,j))*cos(alphafine(s,j))) / &
                                 sin(thetafine(s,j))/sin(alphafine(s,j))
         if ( crelaz.gt.one )  crelaz = one
         if ( crelaz.lt.-one ) crelaz = -one
         if ( relaz_Camera .eq. zero    .and. crelaz .gt. zero ) crelaz = one
         if ( relaz_Camera .eq. 180.0d0 .and. crelaz .lt. zero ) crelaz = -one
         Relazfine(s,j) = acos(cRelaz)
      enddo
   enddo
   Relaz(nsegments) = relaz_Camera_R
   Theta(nsegments) = theta_Camera_R

!  debug output

   if ( do_debug_output ) then
      open(du,file='SolarBeams.dat',status='replace')
      write(du,*)'Regime, Cos/deg scatangle', Regime_idx, cosscat, acos(cosscat)/dtr
      s = 0 ; write(du,*)s,Rsegs(s),alpha(s)/dtr,theta(s)/dtr, Relaz(s)/dtr
      do s = 1, nsegments
         do j = 1, nfinedivs(s)
            write(du,*)s,j,Rsegsfine(s,j),alphafine(s,j)/dtr,thetafine(s,j)/dtr, Relazfine(s,j)/dtr
         enddo
         write(du,*)s,Rsegs(s),alpha(s)/dtr,theta(s)/dtr, Relaz(s)/dtr
      enddo
      close(du)
   endif

!  Warning message

   if ( status .eq. 1 ) then
      message = 'Non-direct sun encountered along path'
   endif

!  debug stop

   if ( do_debug_output ) Stop 'Solar beams, geometry master - DEBUG STOP!!'

!  Normal return

   return

!  Done

   return
end subroutine FOMS_sphergeom_Solar

!  End module

end module FOMS_Geometry_Masters_m
