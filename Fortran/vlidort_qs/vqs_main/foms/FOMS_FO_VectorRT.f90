
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
! # FO ( single-scatter) Vector RT for FOMS                     #
! #                                                             #
! #        FOMS_FO_VectorRT                                     #
! #                                                             #
! ###############################################################

!  Based on "Titan-Probe" FOMS code developed by R. Spurr in May 2012
!  This is a recapitulation and upgrade (VLIDORT Version 2.7)
!  All references to "Probe" have been replaced by "Attrex", now FOMS

!  March 2020 Update
!  8/19/21. Update to Mark2 (first update for 2.5 years)

module FOMS_FO_VectorRT_m

   use FOMS_pars_m

!  For a given wavelength, this routine will calculate FO Stokes Vector at Planet Probe
!     -- For the Atmospheric Single-scatter and Surface Direct-Beam solar sources,

!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.
!       This will perform Enhanced-PS (LOS-path sphericity) calculations

!  7/21/21. Mk2. Addition of the solar beam attenuation at the camera position
!  12/1/21. Mk2. Added subroutine "FOMS_FO_Write_Input", and control for it (do_debug_FO_input)

!  4/14/23. Mk3. Add "OutputDirect-flux" arrays as Output variables (1 extra line)
!                  - Replace endsuntrans output ==> now local variable
!                  - Add theta inputs (necessary for flux calculations) and Flux vector
!                  - Add 2 Control flags for flux (same as in the MS routine) 

!  4/24/24. Add "OutputDirect-fluxes_SB" (SEGMENT BOUNDARY) arrays as Output variables (1 extra line)
!      -- This output was formerly written to file (fort.91) as part of a flux debug

public  :: FOMS_FO_VectorRT
private :: FOMS_FO_Write_Input

contains

subroutine FOMS_FO_VectorRT &
   ( do_debug_FO_input, do_debug_flux, do_Calculate_flux,              & ! Inputs (4/14/23 flux added)
     do_sunlight, do_deltam_scaling, Regime_idx, nsegments, nlayers,   & ! Inputs
     nstokes, nmoments, seglayer_map, segupdown_map, layermask_map,    & ! Inputs 
     nfinedivs, ntraverse, radii, Raycon, cota, theta, sunpaths,       & ! Inputs  (4/14/23 theta added)
     Mu0, Genspher_up, Genspher_dn, Rotations_up, Rotations_dn,        & ! Inputs
     ntraversefine, wfine, csqfine, cotfine, rsegsfine, sunpathsfine,  & ! Inputs
     reflec, extinction, omega, truncfac, greekmat, flux, Fluxvec,     & ! Inputs (Optical)
     stokes, lostrans, sfi_factors, sfi_opdeps,                        & ! Outputs
     OutputDirect_Actinic_flux, OutputDirect_Regular_flux,             & ! Flux Outputs. Added, 4/14/23
     OutputDirect_Actinic_fluxes_SB, OutputDirect_Regular_fluxes_SB )    ! Flux Outputs. Added, 4/24/24

!  Stand alone routine for SS field with Solar sources alone

!  7/21/21. Mk2. Addition of the solar beam attenuation at the camera position
!    - This is necessary for calculating the direct downelling fluxes at the camera
!    - Here, it is calculated anyway, so simple to output

!  4/14/23. Mk3. Add "OutputDirect-flux" arrays as Output variables (1 extra line)
!                  - Replace endsuntrans output ==> now local variable
!                  - Add theta inputs (necessary for flux calculations) and Flux vector
!                  - Add 2 Control flags for flux (same as in the MS routine) 

!  4/24/24. Add "OutputDirect-fluxes_SB" (SEGMENT BOUNDARY) arrays as Output variables (1 extra line)
!      -- This output was formerly written to file (fort.91) as part of a flux debug

   use FOMS_pars_m

   implicit none

!  Inputs
!  ======

!  12/1/21 (Mick mod). Added input debug flag to input arguments

   logical  , intent(in)  :: do_debug_FO_input

!  4/14/23. Add control flags for Flux calculations

   logical  , intent(in)  :: do_debug_flux
   logical  , intent(in)  :: do_Calculate_flux

!  Other flags

   logical  , intent(in)  :: do_sunlight
   logical  , intent(in)  :: do_deltam_scaling

!  number of segments, number of layers, regime index

   integer  , intent(in)  :: nsegments, nlayers, Regime_idx

!  Number of Stokes components, moments

   INTEGER  , Intent(in)  :: nstokes
   INTEGER  , Intent(in)  :: nmoments

!  Segment to layer map, segment up/down switch, layer mask (1=up, 2=dn)

   integer  , intent(in)  :: seglayer_map  (maxsegments)
   logical  , intent(in)  :: segupdown_map (maxsegments)
   logical  , intent(in)  :: layermask_map (maxlays,2)

!  Number of fine divisions

   integer  , intent(in)  :: nfinedivs (maxsegments)

!  Solar path geometry to segment boundaries
!    -- 4/14/23. theta values added

   integer  , intent(in)  :: ntraverse  (0:maxsegments)
   real(ffp), intent(in)  :: theta      (0:maxsegments)
   real(ffp), intent(in)  :: sunpaths   (0:maxsegments,maxlays)

!  Ray constant, radii, Cotangents

   real(ffp), Intent(in)  :: Raycon, cota(0:maxsegments)
   real(ffp), intent(in)  :: radii  (0:maxlays)

!  Mu0 = cos(theta_boa), required for surface term (Regime 1 only)

   real(ffp), Intent(in)  :: Mu0

!  spehrical functions and rotation angles

   real(ffp), intent(in) :: GENSPHER_UP(0:maxmoms,4)
   real(ffp), intent(in) :: GENSPHER_DN(0:maxmoms,4)
   real(ffp), intent(in) :: ROTATIONS_UP(4), ROTATIONS_DN(4)

!  Quadrature values

   real(ffp), intent(in)  :: wfine     (maxsegments,maxfine)
   real(ffp), intent(in)  :: csqfine   (maxsegments,maxfine)
   real(ffp), intent(in)  :: cotfine   (maxsegments,maxfine)
   real(ffp), intent(in)  :: rsegsfine (maxsegments,maxfine)

!  Fine level Solar path geometry

   integer  , intent(in)  :: ntraversefine(maxsegments,maxfine)
   real(ffp), intent(in)  :: sunpathsfine (maxsegments,maxlays,maxfine)

!  Atmosphere optical inputs

   real(ffp), Intent(in) :: EXTINCTION  ( maxlays )
   real(ffp), Intent(in) :: OMEGA       ( maxlays )
   real(ffp), Intent(in) :: TRUNCFAC    ( maxlays )
   real(ffp), Intent(in) :: GREEKMAT    ( maxlays,0:maxmoms, 4, 4 )

!  Solar Flux and Surface reflectivity (Could be the albedo: Regime 1 only)
!    -- 4/14/23. Added Fluxvec input

   real(ffp), Intent(in) :: REFLEC(4,4), FLUX, FLUXVEC(4)

!  outputs
!  -------

!  7/21/21. Mk2. Add endsuntrans, the direct-sun transmittance at the camera
!  4/14/23. Mk2. Here, removed endsuntrans, L_endsuntrans. Set these later on as local variables
!   real(ffp), Intent(out) :: endsuntrans

!  Stokes output and SFI factors for integration

   real(ffp), Intent(out) :: stokes(4)
   real(ffp), Intent(out) :: lostrans    ( maxsegments )
   real(ffp), Intent(out) :: sfi_factors ( maxsegments, maxfine )
   real(ffp), Intent(out) :: sfi_opdeps  ( maxsegments, maxfine )

!  4/14/23. Mk2. Add Output direct flux arrays as Output. Set endsuntrans as local variable
!    -- 5/1/24. Must be Intent(InOut), otherwise you lose it

   real(ffp), Intent(inout) :: OutputDirect_Actinic_flux(4)
   real(ffp), Intent(inout) :: OutputDirect_Regular_flux(4)

!  4/24/24.  Add direct fluxes_SB arrays as Output. (Formerly written to dump file fort.91)
!    -- 5/1/24. Must be Intent(InOut), otherwise you lose it

   real(ffp), Intent(inout) :: OutputDirect_Actinic_fluxes_SB(4,0:maxsegments)
   real(ffp), Intent(inout) :: OutputDirect_Regular_fluxes_SB(4,0:maxsegments)

!  LOCAL
!  -----

!  Attenuations

   real(ffp)  :: attenuations     (0:maxsegments)
   real(ffp)  :: attenuationsfine (maxsegments,maxfine)

!  Scattering

   real(ffp)  :: tms          (maxlays)
   real(ffp)  :: exactscat_up (maxlays,4,4)
   real(ffp)  :: exactscat_dn (maxlays,4,4)

!  Source function integration results

   real(ffp)  :: sources       ( maxsegments, 4 )
   real(ffp)  :: multiplier    ( maxsegments )

!  Help. 4/14/23. endsuntrans variables now local. Also need "muend"

   integer    :: n, ns, nt, j, s, js, s0
   real(ffp)  :: argum(maxfine), tran(maxfine), endsuntrans, muend
   real(ffp)  :: sumd, sum1, cons, help, tran_1, kn, ke, m4, term1(4)
   real(ffp)  :: cot_1, cot_2, suntau(0:maxsegments)
   real(ffp)  :: cumsource_ss(4), cumsource_db(4), fmat(6)

!  4/14/23. Additional local flux variables (debug only)

   real(ffp)  :: Actinic_flux(2)
   real(ffp)  :: Regular_flux(2)

!  cutoff parameter

   real(ffp), parameter  :: cutoff = 88.0_ffp

!  Check inputs if desired
!    12/1/21 (Mick mod). Added this section

   if ( do_debug_FO_input ) then
     call FOMS_FO_Write_Input &
      ( do_sunlight, do_deltam_scaling, Regime_idx, nsegments, nlayers,  & ! Inputs
        nstokes, nmoments, seglayer_map, segupdown_map, layermask_map,   & ! Inputs 
        nfinedivs, ntraverse, radii, Raycon, cota, sunpaths,             & ! Inputs 
        Mu0, Genspher_up, Genspher_dn, Rotations_up, Rotations_dn,       & ! Inputs
        ntraversefine, wfine, csqfine, cotfine, rsegsfine, sunpathsfine, & ! Inputs
        reflec, extinction, omega, truncfac, greekmat, flux)               ! Inputs
   endif

!  Zero the output

   Stokes      = zero ; lostrans   = zero
   sfi_factors = zero ; sfi_opdeps = zero

!  7/21/21, 8/18/21. Mk2. Initialize endsuntrans (the direct-sun transmittance at camera)
!  4/14/23. replace endsuntrans with Output direct flux, and zero them here

   OutputDirect_Regular_flux = zero
   OutputDirect_Actinic_flux = zero

!  4/24/24.  Initialize Output direct fluxes_SB arrays

   OutputDirect_Regular_fluxes_SB = zero
   OutputDirect_Actinic_fluxes_SB = zero

!  Zero the local sources and scattering matrices

   sources      = zero ; fmat         = zero
   exactscat_up = zero ; exactscat_dn = zero

!  shorthand

   ns = nstokes

!  TMS factors

   do n = 1, nlayers
      if ( do_deltam_scaling ) then
         help = one - truncfac(n) * omega(n)
         tms(n) = omega(n) / help
      else
         tms(n) = omega(n)
      endif
   enddo

!  Scattering functions
!  --------------------

!  Scalar only

   if ( nstokes .eq. 1 ) then
      do n = 1, nlayers
         if ( layermask_map(n,1) ) then
           sum1 = DOT_PRODUCT(GenSpher_up(0:nmoments,1), Greekmat(n,0:nmoments,1,1))
           exactscat_up(n,1,1) = sum1 * tms(n)
         endif
         if ( layermask_map(n,2) ) then
           sum1 = DOT_PRODUCT(GenSpher_dn(0:nmoments,1), Greekmat(n,0:nmoments,1,1))
           exactscat_dn(n,1,1) = sum1 * tms(n)
         endif
      enddo
   endif

!  Vector with Sunlight

   if ( nstokes .gt. 1 .and. do_sunlight ) then
      do n = 1, nlayers
         if ( layermask_map(n,1) ) then
            fmat(1) = DOT_PRODUCT(GenSpher_up(0:nmoments,1), Greekmat(n,0:nmoments,1,1))
            fmat(2) = DOT_PRODUCT(GenSpher_up(0:nmoments,2), Greekmat(n,0:nmoments,1,2))
            exactscat_up(n,1,1) = + fmat(1)
            exactscat_up(n,2,1) = - fmat(2) * Rotations_up(3)
            exactscat_up(n,3,1) = + fmat(2) * Rotations_up(4)
            exactscat_up(n,1:ns,1) = tms(n) * exactscat_up(n,1:ns,1) 
         endif
         if ( layermask_map(n,2) ) then
            fmat(1) = DOT_PRODUCT(GenSpher_dn(0:nmoments,1), Greekmat(n,0:nmoments,1,1))
            fmat(2) = DOT_PRODUCT(GenSpher_dn(0:nmoments,2), Greekmat(n,0:nmoments,1,2))
            exactscat_dn(n,1,1) = + fmat(1)
            exactscat_dn(n,2,1) = - fmat(2) * Rotations_dn(3)
            exactscat_dn(n,3,1) = + fmat(2) * Rotations_dn(4)
            exactscat_dn(n,1:ns,1) = tms(n) * exactscat_dn(n,1:ns,1) 
         endif
!write(56,*)n,exactscat_up(n,1:ns,1), exactscat_dn(n,1:ns,1) 
      enddo
   endif

!  Vector General case, USE FULL 4X4 MATRIX
!    CODE INTRODUCED BUT NOT TESTED, 05 OCTOBER 2010. PLACEHOLDER HERE....
!    ** See FO_VLIDORT/FirstOrder_Vector_RTCalcs for detailed Code....

   if ( nstokes .gt. 1 .and. .not. do_sunlight ) then
!    P L A C E H O L D E R
   endif

!  Attenuations
!  ============

!  Initialize

   Attenuations = zero ; Attenuationsfine = zero ; Suntau = zero 

!  Attenuations to End points (including TOA).

   do s = 0, nsegments
      nt = ntraverse(s) ; sumd = DOT_PRODUCT(extinction(1:nt),sunpaths(s,1:nt))
      suntau(s) = sumd  ; If (sumd .lt. cutoff ) Attenuations(s) = exp( - sumd )
   enddo

!  7/21/21, 8/19/21. Mk2. Endsuntrans is just the end-point Attenuation
!  4/14/23. Derive Output flux only at final segment boundary point (s = nsegments)
!   -- Need the attenuations (endsuntrans), muend and Flux/Fluxvec
!   -- 5/1/24. Not restricted to Regime 3

!   if ( do_Calculate_flux .and. Regime_idx .eq. 3 ) then
   if ( do_Calculate_flux ) then
     endsuntrans = Attenuations(nsegments)
     muend      = cos ( theta(nsegments) )
     OutputDirect_Actinic_flux(1:nstokes) = Flux * endsuntrans * Fluxvec(1:nstokes) 
     OutputDirect_Regular_flux(1:nstokes) = Pi4 * OutputDirect_Actinic_flux(1:nstokes) * muend
   endif

!  Old code, pre 4/14/23
!   endsuntrans = Attenuations(nsegments)

!  fine-layer attenuations

   do s = 1, nsegments
      n = seglayer_map(s)
      do j = 1, nfinedivs(s)
         nt = ntraversefine(s,j) ; sumd = DOT_PRODUCT(extinction(1:nt),sunpathsfine(s,1:nt,j))
         if (sumd .lt. cutoff ) Attenuationsfine(s,j) = exp( - sumd )
      enddo
   enddo

!  4/14/23. Compute flux at top and bottom, debug output (same style as MS routine) to fort.91
!     Use direct flux to the segment boundary points.
!  -- 5/1/24. Not restricted to Regime 3

!   if ( do_Calculate_flux .and. do_debug_flux .and. Regime_idx.eq.3 ) then
   if ( do_Calculate_flux .and. do_debug_flux ) then

!  Top of first segment

     endsuntrans = Attenuations(0)
     muend       = cos (theta(0))
     Actinic_flux(1) = Flux * endsuntrans * Fluxvec(1) 
     Regular_flux(1) = Pi4 * Actinic_flux(1) * muend

!  4/24/24. Copy these fluxes to the output "Fluxes_SB" arrays

     OutputDirect_Actinic_fluxes_SB(1,0) = Actinic_flux(1)
     OutputDirect_Regular_fluxes_SB(1,0) = Regular_flux(1)

!  Main segment loop

     do s = 1, nsegments

!  Bottom of segment

       endsuntrans = Attenuations(s)
       muend       = cos (theta(s))
       Actinic_flux(2) = Flux * endsuntrans * Fluxvec(1) 
       Regular_flux(2) = Pi4 * Actinic_flux(2) * muend

!  4/24/24. Copy these fluxes to the output "Fluxes_SB" arrays
!   --replaces the debug write-up to fort.91

       OutputDirect_Actinic_fluxes_SB(1,s) = Actinic_flux(2)
       OutputDirect_Regular_fluxes_SB(1,s) = Regular_flux(2)

!  write-up top/bottom (Rob debug)

!       if ( s.eq.1) write(91,'(/A/)')' segm. layer    ActFlux @TOP   ActFlux @BOT    RegFlux @TOP   RegFlux @BOT'
!       write(91,455)s,s,Actinic_flux(1),Actinic_flux(2),Regular_flux(1),Regular_flux(2)
!455   format(i4,3x,i4,2(1x,1p2e15.5))

!  re-set the segment-top values

       Actinic_flux(1) = Actinic_flux(2) 
       Regular_flux(1) = Regular_flux(2)

!  End segment loop

     enddo

!  Check
!     do s = 0, nsegments
!       write(991,455)s,s,OutputDirect_Actinic_fluxes_SB(1,s),OutputDirect_Regular_fluxes_SB(1,s)
!     enddo

!  End Flux output and debug clause

   endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Layer integrated multipliers
!  ============================

!  Regime 1. Upwelling from BOA

   if ( regime_idx .eq. 1 ) then
      do s = 1, nsegments
         n = seglayer_map(s) ; js = nfinedivs(s)
         cot_2 = cota(s) ; cot_1 = cota(s-1)
         kn = extinction(n) ;  ke = raycon * kn 

!  LOS transmittance

         cons = raycon * ( cot_2 - cot_1 ) ; tran_1 = kn * cons
         if ( tran_1 .lt. cutoff ) lostrans(s) = exp ( - tran_1 )

!  Multiplier

         do j = 1, nfinedivs(s)
            argum(j) = Raycon * ( cot_2 - cotfine(s,j) )  !  Up
            tran(j)  = exp ( - kn * argum(j) )
            sfi_opdeps (s,j) = kn * ( radii(n-1) - rsegsfine(s,j) )
            sfi_factors(s,j) = ke * csqfine(s,j) * tran(j) * wfine(s,j)
         enddo
         multiplier(s) = DOT_PRODUCT(attenuationsfine(s,1:js),sfi_factors(s,1:js))

!  End

      enddo
   endif

!  Regime 2. Limb path

   if ( regime_idx .eq. 2 ) then

!  tangent layer
!      nt = maxval(seglayer_map(1:nsegments))

!  Start segment loop

      do s = 1, nsegments

!  Start

         n = seglayer_map(s) ; js = nfinedivs(s)
         kn = extinction(n)  ; ke = raycon * kn
         cot_2 = cota(s) ; cot_1 = cota(s-1)

!  Downwelling vs. Upwelling segment

         if ( .not.segUpdown_map(s) ) then
            cons = raycon * ( cot_1 - cot_2 ) 
            do j = 1, nfinedivs(s)
               argum(j) = Raycon * ( cotfine(s,j) - cot_2 )
            enddo
         else
            cons = raycon * ( cot_2 - cot_1 ) 
            do j = 1, nfinedivs(s)
               argum(j) = Raycon * ( cot_2 - cotfine(s,j) )
            enddo
         endif

!  LOS transmittance

         tran_1 = kn * cons
         if ( tran_1 .lt. cutoff ) lostrans(s) = exp ( - tran_1 )

!  Multiplier

         do j = 1, nfinedivs(s)
            tran(j)  = exp ( - kn * argum(j) )
            sfi_opdeps (s,j) = kn * ( radii(n-1) - rsegsfine(s,j) )
            sfi_factors(s,j) = ke * csqfine(s,j) * tran(j) * wfine(s,j)
         enddo
         multiplier(s) = DOT_PRODUCT(attenuationsfine(s,1:js),sfi_factors(s,1:js))

!  End

      enddo
   endif

!  Regime 3. Downwelling from TOA

   if ( regime_idx .eq. 3 ) then
      do s = 1, nsegments

!  Start

         n = seglayer_map(s) ; js = nfinedivs(s)
         kn = extinction(n)  ;  ke = raycon * kn
         cot_2 = cota(s) ; cot_1 = cota(s-1)

!  LOS transmittance

         cons = raycon * ( cot_1 - cot_2 ) ; tran_1 = kn * cons
         if ( tran_1 .lt. cutoff ) lostrans(s) = exp ( - tran_1 )

!  Multiplier

         do j = 1, nfinedivs(s)
            argum(j) = Raycon * ( cotfine(s,j) - cot_2 )
            tran(j)  = exp ( - kn * argum(j) )
            sfi_opdeps (s,j) = kn * ( radii(n-1) - rsegsfine(s,j) )
            sfi_factors(s,j) = ke * csqfine(s,j) * tran(j) * wfine(s,j)
         enddo
         multiplier(s) = DOT_PRODUCT(attenuationsfine(s,1:js),sfi_factors(s,1:js))

!   End

      enddo
   endif

!  Layer sources
!  =============

!   PLACEHOLDER - non-sunlight case

!  Regime 1. Upwelling from BOA

   if ( regime_idx .eq. 1 .and. do_sunlight ) then
      do s = 1, nsegments
         n = seglayer_map(s)
         sources(s,1:ns) = exactscat_up(n,1:ns,1) * multiplier(s)
      enddo
   endif

!  Regime 2. Limb path

   if ( regime_idx .eq. 2 .and. do_sunlight ) then
      do s = 1, nsegments
         n = seglayer_map(s)
         if ( .not.segUpdown_map(s) ) then
            sources(s,1:ns) = exactscat_dn(n,1:ns,1) * multiplier(s)
         else
            sources(s,1:ns) = exactscat_up(n,1:ns,1) * multiplier(s)
         endif
      enddo
   endif

!  Regime 3. Downwelling from TOA

   if ( regime_idx .eq. 3 .and. do_sunlight ) then
      do s = 1, nsegments
         n = seglayer_map(s)
         sources(s,1:ns) = exactscat_dn(n,1:ns,1) * multiplier(s)
      enddo
   endif

!  Source function integration
!  ===========================

!  RT recursion

   CUMSOURCE_SS = zero
   do s = 1, nsegments
      CUMSOURCE_SS(1:ns) = LOSTRANS(S) * CUMSOURCE_SS(1:ns) + SOURCES(S,1:ns)
   enddo
   STOKES(1:ns) = FLUX * CUMSOURCE_SS(1:ns)

!  OLD CODE
!   dub = 668 ; if (regime_idx.eq.3 ) dub = 669
!   if ( Regime_idx .eq. 1 ) CUMSOURCE_DB(1:ns) = four * Mu0 * REFLEC(1:ns,1) * attenuations(0)
!   write(dub,'(a,f10.4)')'Recursion Output for viewing angle = ',atan(one/cota(nsegments))*180.0d0/acos(-1.0d0)
!   CUMSOURCE_SS = zero
!   do s = 1, nsegments
!      CUMSOURCE_SS(1:ns) = LOSTRANS(S) * CUMSOURCE_SS(1:ns) + SOURCES(S,1:ns)
!if (s.gt.20)write(*,*)'FO',s,lostrans(s),sources(s,1)
!     if ( Regime_idx .eq. 1 ) then
!        CUMSOURCE_DB(1:ns) = LOSTRANS(S) * CUMSOURCE_DB(1:ns)
!        write(dub,667)Regime_idx,s,seglayer_map(s),segUpdown_map(s),cota(s),atan(one/cota(s))*180.0d0/acos(-1.0d0),&
!                lostrans(s),multiplier(s),sources(s,1),flux*cumsource_ss(1),flux*cumsource_db(1)
!     else
!        write(dub,667)Regime_idx,s,seglayer_map(s),segUpdown_map(s),cota(s),atan(one/cota(s))*180.0d0/acos(-1.0d0),&
!                lostrans(s),multiplier(s),sources(s,1),flux*cumsource_ss(1)        
!     endif
!   enddo
!   STOKES(1:ns) = FLUX * CUMSOURCE_SS(1:ns)

!667 format(i2,2i4,1x,l2,1x,2f12.7,1p5e16.8)

!  Add surface term for Regime 1. Valid for BRDF as well.
!  Assumes natural light only (unpolarized solar beam)

   CUMSOURCE_DB = zero ; M4 = zero ; s0 =  0
   if ( Regime_idx .eq. 1 ) then
      M4 = four * MU0 ; CUMSOURCE_DB(1:ns) = M4 * REFLEC(1:ns,1) * attenuations(s0)
      do s = 1, nsegments
         CUMSOURCE_DB(1:ns) = LOSTRANS(S) * CUMSOURCE_DB(1:ns)
      enddo
      STOKES(1:ns) = STOKES(1:ns) + FLUX * CUMSOURCE_DB(1:ns)
   endif

!  Finish

   return
end subroutine FOMS_FO_VectorRT

! PRIVATE ROUTINE

subroutine FOMS_FO_Write_Input &
   ( do_sunlight, do_deltam_scaling, Regime_idx, nsegments, nlayers,  & ! Inputs
     nstokes, nmoments, seglayer_map, segupdown_map, layermask_map,   & ! Inputs 
     nfinedivs, ntraverse, radii, Raycon, cota, sunpaths,             & ! Inputs 
     Mu0, Genspher_up, Genspher_dn, Rotations_up, Rotations_dn,       & ! Inputs
     ntraversefine, wfine, csqfine, cotfine, rsegsfine, sunpathsfine, & ! Inputs
     reflec, extinction, omega, truncfac, greekmat, flux)               ! Inputs (Optical)

!  Check inputs if desired
!mick mod 12/1/2021  - added this subroutine (but not tested yet)

   use FOMS_pars_m

   implicit none

!  Inputs
!  ======

!  Standard
!  --------

!  Flags

   logical  , intent(in)  :: do_sunlight
   logical  , intent(in)  :: do_deltam_scaling

!  Number of segments, number of layers, regime index

   integer  , intent(in)  :: nsegments, nlayers, Regime_idx

!  Number of Stokes components, moments

   integer  , intent(in)  :: nstokes
   integer  , intent(in)  :: nmoments

!  Segment to layer map, segment up/down switch, layer mask (1=up, 2=dn)

   integer  , intent(in)  :: seglayer_map  (maxsegments)
   logical  , intent(in)  :: segupdown_map (maxsegments)
   logical  , intent(in)  :: layermask_map (maxlays,2)

!  Number of fine divisions

   integer  , intent(in)  :: nfinedivs (maxsegments)

!  Solar path geometry to segment boundaries

   integer  , intent(in)  :: ntraverse  (0:maxsegments)
   real(ffp), intent(in)  :: sunpaths   (0:maxsegments,maxlays)

!  Ray constant, radii, Cotangents

   real(ffp), intent(in)  :: Raycon, cota(0:maxsegments)
   real(ffp), intent(in)  :: radii  (0:maxlays)

!  Mu0 = cos(theta_boa), required for surface term (Regime 1 only)

   real(ffp), intent(in)  :: Mu0

!  spehrical functions and rotation angles

   real(ffp), intent(in)  :: genspher_up(0:maxmoms,4)
   real(ffp), intent(in)  :: genspher_dn(0:maxmoms,4)
   real(ffp), intent(in)  :: rotations_up(4), rotations_dn(4)

!  Quadrature values

   real(ffp), intent(in)  :: wfine     (maxsegments,maxfine)
   real(ffp), intent(in)  :: csqfine   (maxsegments,maxfine)
   real(ffp), intent(in)  :: cotfine   (maxsegments,maxfine)
   real(ffp), intent(in)  :: rsegsfine (maxsegments,maxfine)

!  Fine level Solar path geometry

   integer  , intent(in)  :: ntraversefine(maxsegments,maxfine)
   real(ffp), intent(in)  :: sunpathsfine (maxsegments,maxlays,maxfine)

!  Atmosphere optical inputs

   real(ffp), intent(in)  :: extinction  ( maxlays )
   real(ffp), intent(in)  :: omega       ( maxlays )
   real(ffp), intent(in)  :: truncfac    ( maxlays )
   real(ffp), intent(in)  :: greekmat    ( maxlays,0:maxmoms, 4, 4 )

!  Solar Flux and Surface reflectivity (Could be the albedo: Regime 1 only)

   real(ffp), intent(in)  :: flux, reflec(4,4)

!  Local variables
!  ---------------

   integer :: outunit
   integer :: lay,mom,nf,ngreek,o1,o2,s,ss
!   integer :: geo,i,j,k,sza,ulev,ura,uva,uog, ngeoms

   integer :: gmask(8)
   data gmask / 1, 2, 5, 6, 11, 12, 15, 16 /

   integer :: gmask2(8,2)

!  Open output file

   OUTUNIT = 101
   OPEN (OUTUNIT,file = 'Debug\FOMS_FO_Write_Input.dbg',status = 'replace')

!  Only output non-zero entries of GreekMat

   if ( nstokes.eq.1 ) then
     ngreek = 1
   elseif ( nstokes.eq.3 ) then
     ngreek = 5
   elseif ( nstokes.eq.4 ) then
     ngreek = 8
   endif

!  Define non-zero GreekMat components (a 2-D 4x4 array):
!                O1               O2
   gmask2(1,1) = 1; gmask2(1,2) = 1
   if (nstokes.gt.1) then
     gmask2(2,1) = 1; gmask2(2,2) = 2
     gmask2(3,1) = 2; gmask2(3,2) = 1
     gmask2(4,1) = 2; gmask2(4,2) = 2
     gmask2(5,1) = 3; gmask2(5,2) = 3
     if (nstokes.gt.3) then
       gmask2(6,1) = 3; gmask2(6,2) = 4
       gmask2(7,1) = 4; gmask2(7,2) = 3
       gmask2(8,1) = 4; gmask2(8,2) = 4
     endif
   endif

!  Write all VQS inputs to file

   WRITE(OUTUNIT,*)
   WRITE(OUTUNIT,'(A)') '---------------'
   WRITE(OUTUNIT,'(A)') 'Standard Inputs'
   WRITE(OUTUNIT,'(A)') '---------------'

!  Control

   WRITE(OUTUNIT,*)
   WRITE(OUTUNIT,DWFL) 'do_sunlight       = ',do_sunlight
   WRITE(OUTUNIT,DWFL) 'do_deltam_scaling = ',do_deltam_scaling

!  Number of segments, number of layers, regime index

   WRITE(OUTUNIT,*)
   WRITE(OUTUNIT,DWFI) 'nsegments  = ',nsegments
   WRITE(OUTUNIT,DWFI) 'nlayers    = ',nlayers
   WRITE(OUTUNIT,DWFI) 'Regime_idx = ',Regime_idx

!  Number of Stokes components, moments

   WRITE(OUTUNIT,*)
   WRITE(OUTUNIT,DWFI) 'nstokes    = ',nstokes
   WRITE(OUTUNIT,DWFI) 'nmoments   = ',nmoments

!  Segment to layer map, segment up/down switch, layer mask (1=up, 2=dn)

   WRITE(OUTUNIT,*)
   do s=1,nsegments
     WRITE(OUTUNIT,DWFI1) 's = ',s,' seglayer_map(s)  = ',seglayer_map(s)
   enddo

   WRITE(OUTUNIT,*)
   do s=1,nsegments
     WRITE(OUTUNIT,DWFL1) 's = ',s,' segupdown_map(s) = ',segupdown_map(s)
   enddo

   WRITE(OUTUNIT,*)
   do lay=1,nlayers
     WRITE(OUTUNIT,DWFL1_2) 'lay = ',lay,&
                           ' layermask_map(lay,1) = ',layermask_map(lay,1),&
                           ' layermask_map(lay,2) = ',layermask_map(lay,2)
   enddo

!  Number of fine divisions

   WRITE(OUTUNIT,*)
   do s=1,nsegments
     WRITE(OUTUNIT,DWFI1) 's = ',s,' nfinedivs(s) = ',nfinedivs(s)
   enddo

!  Solar path geometry to segment boundaries

   WRITE(OUTUNIT,*)
   do s=0,nsegments
     WRITE(OUTUNIT,DWFI1) 's = ',s,' ntraverse(s) = ',ntraverse(s)
   enddo

   WRITE(OUTUNIT,*)
   do lay=1,nlayers
     do s=0,nsegments
       WRITE(OUTUNIT,DWFR2) 'lay = ',lay,' s = ',s,' sunpaths(s,lay) = ',sunpaths(s,lay)
     enddo
   enddo

!  Ray constant, radii, Cotangents

   WRITE(OUTUNIT,*)
   WRITE(OUTUNIT,DWFR) 'Raycon = ',Raycon

   WRITE(OUTUNIT,*)
   do s=0,nsegments
     WRITE(OUTUNIT,DWFR1) 's = ',s,' cota(s) = ',cota(s)
   enddo

   WRITE(OUTUNIT,*)
   do lay=0,nlayers
     WRITE(OUTUNIT,DWFR1) 'lay = ',lay,' radii(lay) = ',radii(lay)
   enddo

!  Mu0 = cos(theta_boa), required for surface term (Regime 1 only)

   WRITE(OUTUNIT,*)
   WRITE(OUTUNIT,DWFR) 'Mu0 = ',Mu0

!  Spherical functions and rotation angles

   if (Regime_idx.eq.1 .or. Regime_idx.eq.2) then
     WRITE(OUTUNIT,*)
     do o1=1,nstokes
       do mom=0,nmoments
         WRITE(OUTUNIT,DWFR2) 'mom = ',mom,' o1 = ',o1,' genspher_up(mom,o1) = ',genspher_up(mom,o1)
       enddo
     enddo

     WRITE(OUTUNIT,*)
     do o1=1,nstokes
       WRITE(OUTUNIT,DWFR1) 'o1 = ',o1,' rotations_up(o1) = ',rotations_up(o1)
     enddo
   endif

   if (Regime_idx.eq.2 .or. Regime_idx.eq.3) then
     WRITE(OUTUNIT,*)
     do o1=1,nstokes
       do mom=0,nmoments
         WRITE(OUTUNIT,DWFR2) 'mom = ',mom,' o1 = ',o1,' genspher_dn(mom,o1) = ',genspher_dn(mom,o1)
       enddo
     enddo

     WRITE(OUTUNIT,*)
     do o1=1,nstokes
       WRITE(OUTUNIT,DWFR1) 'o1 = ',o1,' rotations_dn(o1) = ',rotations_dn(o1)
     enddo
   endif

!  Quadrature values

   WRITE(OUTUNIT,*)
   do s=1,nsegments
     do nf=1,nfinedivs(s)
       WRITE(OUTUNIT,DWFR2) 's = ',s,' nf = ',nf,' wfine(s,nf)    = ',wfine(s,nf)
     enddo
   enddo

   WRITE(OUTUNIT,*)
   do s=1,nsegments
     do nf=1,nfinedivs(s)
       WRITE(OUTUNIT,DWFR2) 's = ',s,' nf = ',nf,' csqfine(s,nf)  = ',csqfine(s,nf)
     enddo
   enddo

   WRITE(OUTUNIT,*)
   do s=1,nsegments
     do nf=1,nfinedivs(s)
       WRITE(OUTUNIT,DWFR2) 's = ',s,' nf = ',nf,' cotfine(s,nf)  = ',cotfine(s,nf)
     enddo
   enddo

   WRITE(OUTUNIT,*)
   do s=1,nsegments
     do nf=1,nfinedivs(s)
       WRITE(OUTUNIT,DWFR2) 's = ',s,' nf = ',nf,' rsegsfine(s,nf) = ',rsegsfine(s,nf)
     enddo
   enddo

!  Fine level Solar path geometry
!mick fix 6/6/2023 - swapped 1st & 2nd dims of "sunpathsfine"

   WRITE(OUTUNIT,*)
   do s=1,nsegments
     do nf=1,nfinedivs(s)
       WRITE(OUTUNIT,DWFI2) 's = ',s,' nf = ',nf,' ntraversefine(s,nf) = ',ntraversefine(s,nf)
     enddo
   enddo

   WRITE(OUTUNIT,*)
   do s=1,nsegments
     do lay=1,nlayers
       do nf=1,nfinedivs(s)
         WRITE(OUTUNIT,DWFR3) 's = ',s,' lay = ',lay,' nf = ',nf,&
                              ' sunpathsfine(s,lay,nf) = ',sunpathsfine(s,lay,nf)
       enddo
     enddo
   enddo

!  Atmosphere optical inputs

   WRITE(OUTUNIT,*)
   do lay=1,nlayers
     WRITE(OUTUNIT,DWFR1) 'lay = ',lay,' extinction(lay) = ',extinction(lay)
   enddo

   WRITE(OUTUNIT,*)
   do lay=1,nlayers
     WRITE(OUTUNIT,DWFR1) 'lay = ',lay,' omega(lay)      = ',omega(lay)
   enddo

   WRITE(OUTUNIT,*)
   do lay=1,nlayers
     WRITE(OUTUNIT,DWFR1) 'lay = ',lay,' truncfac(lay)   = ',truncfac(lay)
   enddo

   do ss=1,ngreek
     o1 = gmask2(ss,1); o2 = gmask2(ss,2); WRITE(OUTUNIT,*)
     do lay=1,nlayers
       do mom=0,nmoments
         WRITE(OUTUNIT,DWFR5)  'ss = ',ss,' o1 = ',o1,' o2 = ',o2,' lay = ',lay,' mom = ',mom,&
           ' greekmat(lay,mom,o1,o2) = ',greekmat(lay,mom,o1,o2)
       enddo
     enddo
   enddo

!  Solar Flux and Surface reflectivity (Could be the albedo: Regime 1 only)

   WRITE(OUTUNIT,*)
   WRITE(OUTUNIT,DWFR) 'Flux = ',Flux

   WRITE(OUTUNIT,*)
   do o1=1,4
     do o2=1,4
       WRITE(OUTUNIT,DWFR2) 'o1 = ',o1,' o2 = ',o2,&
                           ' reflec(o1,o2) = ',reflec(o1,o2)
     enddo
   enddo

!  Close output file

   CLOSE (OUTUNIT)

end subroutine FOMS_FO_Write_Input

! End module

end module FOMS_FO_VectorRT_m
