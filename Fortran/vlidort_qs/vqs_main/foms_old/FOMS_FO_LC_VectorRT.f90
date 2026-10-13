
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
! #       Delivers Column Jacobians in addition                 #
! #                                                             #
! #        FOMS_FO_LC_VectorRT                                  #
! #                                                             #
! ###############################################################

!  Based on "Titan-Probe" FOMS code developed by R. Spurr in May 2012
!  This is a recapitulation and upgrade (VLIDORT Version 2.7)
!  All references to "Probe" have been replaced by "FOMS"

!  March 2020 Update (regular code)
!  8/19/21. Update to Mark2 (first update for 2.5 years)
!  Mk2a. 9/21/21. This module is new (Column Jacobians)

module FOMS_FO_LC_VectorRT_m

   use FOMS_pars_m

!  For a given wavelength, this routine will calculate FO Stokes Vector at Planet Probe
!     -- For the Atmospheric Single-scatter and Surface Direct-Beam solar sources,
!     -- This subroutine delivers column Jacobians in addition

!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.
!       This will perform Enhanced-PS (LOS-path sphericity) calculations

!  7/21/21. Mk2. Addition of the solar beam attenuation at the camera position
!  8/19/21. Mk2, Add also the Jacobians of this attenuation
!  9/21/21. Mk2a. New module for column Jacobians
!  12/1/21. Mk2. Added subroutine "FOMS_FO_LC_Write_Input", and control for it (do_debug_FO_input)

!  5/30/22. Mk2. Add "OutputDirect-flux" arrays as Output variables (2 extra lines)
!                  - Replace endsuntrans & LC_Endsuntrans output ==> now local variables
!                  - Add theta inputs (necessary for flux calculations) and Flux vector
!                  - Add 2 Control flags for flux (same as in the MS routine) 

!  4/24/24. Add "OutputDirect-fluxes_SB" (SEGMENT BOUNDARY) arrays as Output variables (1 extra line)
!      -- This output was formerly written to file (fort.91) as part of a flux debug
!      -- Do same with LC Jacobians of these quantities (formerly written to flux debug file fort.92)

public  :: FOMS_FO_LC_VectorRT
private :: FOMS_FO_LC_Write_Input

contains

subroutine FOMS_FO_LC_VectorRT &
   ( do_debug_FO_input, do_debug_flux, do_Calculate_flux,               & ! Inputs (5/30/22 flux added)
     do_sunlight, do_deltam_scaling, Regime_idx, nsegments, nlayers,    & ! Inputs
     nstokes, nmoments, seglayer_map, segupdown_map, layermask_map,     & ! Inputs 
     nfinedivs, ntraverse, radii, Raycon, cota, theta, sunpaths,        & ! Inputs  (5/30/22 theta added)
     Mu0, Genspher_up, Genspher_dn, Rotations_up, Rotations_dn,         & ! Inputs
     ntraversefine, wfine, csqfine, cotfine, rsegsfine, sunpathsfine,   & ! Inputs
     do_Column_Wfs, Lvaryflags, Lvarynums, Lvarymoms,                   & ! Lin Control
     reflec, extinction, omega, truncfac, greekmat, flux, Fluxvec,      & ! Inputs (Optical)
     L_extinction, L_omega, L_truncfac, L_greekmat,                     & ! Inputs (Optical, linearized)
     stokes, lostrans, sfi_factors, sfi_opdeps,                         & ! Outputs
     LC_stokes, L_lostrans, L_sfi_factors, L_sfi_opdeps,                & ! Outputs (Linearized)
     OutputDirect_Actinic_flux, OutputDirect_Regular_flux,                  & ! Flux    Outputs. Added, 5/30/22
     LC_OutputDirect_Actinic_flux, LC_OutputDirect_Regular_flux,            & ! LC Flux Outputs  Added, 5/30/22
     OutputDirect_Actinic_fluxes_SB, OutputDirect_Regular_fluxes_SB,        & ! SB Flux    Outputs. Added, 4/24/24
     LC_OutputDirect_Actinic_fluxes_SB, LC_OutputDirect_Regular_fluxes_SB )   ! SB LC Flux Outputs. Added, 4/24/24

!  Stand alone routine for SS field with Solar sources alone

!  7/21/21. Mk2. Addition of the solar beam attenuation at the camera position
!    - This is necessary for calculating the direct downelling fluxes at the camera
!    - Here, it is calculated anyway, so simple to output. Include also Jacobians of this.

!  9/21/21. Mk2a. New module for column Jacobians

!  5/30/22. Mk2. Add "OutputDirect-flux" arrays as Output variables (2 extra lines)
!                  - Replace endsuntrans & LC_Endsuntrans output ==> now local variables
!                  - Add theta inputs (necessary for flux calculations) and Flux vector
!                  - Add 2 Control flags for flux (same as in the MS routine) 

!  4/24/24. Add "OutputDirect-fluxes_SB" (SEGMENT BOUNDARY) arrays as Output variables (1 extra line)
!      -- This output was formerly written to file (fort.91) as part of a flux debug
!      -- Do same with LC Jacobians of these quantities (formerly written to flux debug file fort.92)

   use FOMS_pars_m

   implicit none         

!  Inputs
!  ======

!  12/1/21 (Mick mod). Added input debug flag to input arguments

   logical  , intent(in)  :: do_debug_FO_input

!  5/30/22. Add control flags for Flux calculations

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
!    -- 5/30/22. theta values added

   integer  , intent(in)  :: ntraverse  (0:maxsegments)
   real(ffp), intent(in)  :: theta      (0:maxsegments)
   real(ffp), intent(in)  :: sunpaths   (0:maxsegments,maxlays)

!  Ray constant, radii, Cotangents

   real(ffp), Intent(in)  :: Raycon, cota(0:maxsegments)
   real(ffp), intent(in)  :: radii  (0:maxlays)

!  Mu0 = cos(theta_boa), required for surface term (Regime 1 only)

   real(ffp), Intent(in)  :: Mu0

!  spehrical functions and rotation angles

   REAL(ffp), intent(in) :: GENSPHER_UP(0:maxmoms,4)
   REAL(ffp), intent(in) :: GENSPHER_DN(0:maxmoms,4)
   REAL(ffp), intent(in) :: ROTATIONS_UP(4), ROTATIONS_DN(4)

!  Quadrature values

   real(ffp), intent(in)  :: wfine     (maxsegments,maxfine)
   real(ffp), intent(in)  :: csqfine   (maxsegments,maxfine)
   real(ffp), intent(in)  :: cotfine   (maxsegments,maxfine)
   real(ffp), intent(in)  :: rsegsfine (maxsegments,maxfine)

!  Fine level Solar path geometry

   integer  , intent(in)  :: ntraversefine(maxsegments,maxfine)
   real(ffp), intent(in)  :: sunpathsfine (maxsegments,maxlays,maxfine)

!  Atmosphere optical inputs

   REAL(ffp), Intent(in) :: EXTINCTION  ( maxlays )
   REAL(ffp), Intent(in) :: OMEGA       ( maxlays )
   REAL(ffp), Intent(in) :: TRUNCFAC    ( maxlays )
   REAL(ffp), Intent(in) :: GREEKMAT    ( maxlays,0:maxmoms, 4, 4 )

!  Solar Flux and Surface reflectivity (Could be the albedo: Regime 1 only)
!    -- 5/30/22. Added Fluxvec input

   REAL(ffp), Intent(in) :: REFLEC(4,4), FLUX, FLUXVEC(4)

!  Jacobian control

   LOGICAL, Intent(in) ::  Do_Column_WFs
   LOGICAL, Intent(in) ::  Lvaryflags(maxlays)
   INTEGER, Intent(in) ::  Lvarynums (maxlays)
   LOGICAL, Intent(in) ::  Lvarymoms (maxlays,maxJacobs)

!  Linearized optical inputs

   REAL(ffp), Intent(in) :: L_EXTINCTION  ( maxlays, maxJacobs )
   REAL(ffp), Intent(in) :: L_OMEGA       ( maxlays, maxJacobs )
   REAL(ffp), Intent(in) :: L_TRUNCFAC    ( maxlays, maxJacobs )
   REAL(ffp), Intent(in) :: L_GREEKMAT    ( maxlays,0:maxmoms, 4, 4, maxJacobs )

!  outputs
!  -------

!  7/21/21. Mk2. Add endsuntrans, the direct-sun transmittance at the camera
!  8/19/21. Mk2. Add L_endsuntrans
!  5/30/22. Mk2. Here, removed endsuntrans, L_endsuntrans. Set these later on as local variables
!   real(ffp), Intent(out) :: endsuntrans
!   real(ffp), Intent(out) :: LC_endsuntrans ( maxJacobs )

!  Stokes output and SFI factors for integration

   real(ffp), Intent(out) :: stokes(4)
   real(ffp), Intent(out) :: lostrans    ( maxsegments )
   real(ffp), Intent(out) :: sfi_factors ( maxsegments, maxfine )
   real(ffp), Intent(out) :: sfi_opdeps  ( maxsegments, maxfine )

!  Linearized Output

   real(ffp), Intent(out) :: LC_Stokes      ( 4, maxJacobs )
   real(ffp), Intent(out) :: L_lostrans     ( maxsegments, maxJacobs )
   real(ffp), Intent(out) :: L_sfi_factors  ( maxsegments, maxfine, maxJacobs )
   real(ffp), Intent(out) :: L_sfi_opdeps   ( maxsegments, maxfine, maxJacobs )

!  5/30/22. Mk2. Add Output direct flux arrays as Output. Set endsuntrans & LC_Endsuntrans as local variables

   real(ffp), Intent(out) :: OutputDirect_Actinic_flux(4)
   real(ffp), Intent(out) :: OutputDirect_Regular_flux(4)
   real(ffp), Intent(out) :: LC_OutputDirect_Actinic_flux ( maxJacobs, 4 )
   real(ffp), Intent(out) :: LC_OutputDirect_Regular_flux ( maxJacobs, 4 )

!  4/24/24.  Add direct fluxes_SB arrays as Output. (Formerly written to dump file fort.91)
!      -- Do same with LC Jacobians of these quantities (formerly written to flux debug file fort.92)

   real(ffp), Intent(out) :: OutputDirect_Actinic_fluxes_SB(4,0:maxsegments)
   real(ffp), Intent(out) :: OutputDirect_Regular_fluxes_SB(4,0:maxsegments)
   real(ffp), Intent(out) :: LC_OutputDirect_Actinic_fluxes_SB(4,0:maxsegments,maxJacobs)
   real(ffp), Intent(out) :: LC_OutputDirect_Regular_fluxes_SB(4,0:maxsegments,maxJacobs)

!  LOCAL
!  -----

!  Attenuations

   real(ffp)  :: attenuations     (0:maxsegments)
   real(ffp)  :: attenuationsfine (maxsegments,maxfine)
   real(ffp)  :: L_attenuations     (0:maxsegments, maxJacobs)
   real(ffp)  :: L_attenuationsfine (maxsegments,maxfine, maxJacobs)

!  Scattering

   real(ffp)  :: tms          (maxlays)
   real(ffp)  :: exactscat_up (maxlays,4,4)
   real(ffp)  :: exactscat_dn (maxlays,4,4)

   real(ffp)  :: L_tms          (maxlays,maxJacobs)
   real(ffp)  :: L_exactscat_up (maxlays,4,4,maxJacobs)
   real(ffp)  :: L_exactscat_dn (maxlays,4,4,maxJacobs)

!  Source function integration results

   real(ffp)  :: sources       ( maxsegments, 4 )
   real(ffp)  :: multiplier    ( maxsegments )
   real(ffp)  :: L_multiplier    ( maxsegments, maxJacobs )
   real(ffp)  :: L_sources       ( maxsegments, 4, maxJacobs )

!  Help. 5/30/22. endsuntrans variables now local. Also need "muend"

   logical    :: Qvary(maxlays)
   integer    :: n, ns, nt, j, js, k, s, s0, q, Qnums(maxlays), n_column_wfs
   real(ffp)  :: argum(maxfine), tran(maxfine), endsuntrans
   real(ffp)  :: L_help, L_sum, L_tran, L_suntau, L_sumd, LC_endsuntrans, muend
   real(ffp)  :: sumd, sum1, cons, help, tran_1, kn, ke, m4, term1(4), L_kn, L_kn_reg
   real(ffp)  :: cot_1, cot_2, suntau(0:maxsegments), fmat(6), L_fmat(6), L_exact(4)
   real(ffp)  :: cumsource_ss(0:maxsegments,4), cumsource_db(0:maxsegments,4)
   real(ffp)  :: L_cumsource_ss ( 4, maxJacobs ), L_cumsource_db ( 4, maxJacobs )

!  5/30/22. Additional local flux variables (debug only)

   real(ffp)  :: Actinic_flux(2), LC_Actinic_flux ( maxJacobs, 2 )
   real(ffp)  :: Regular_flux(2), LC_Regular_flux ( maxJacobs, 2 )

!  cutoff parameter

   real(ffp), parameter  :: cutoff = 88.0_ffp

!  Check inputs if desired
!    12/1/21 (Mick mod). Added this section

   if ( do_debug_FO_input ) then
     call FOMS_FO_LC_Write_Input &
      ( do_sunlight, do_deltam_scaling, Regime_idx, nsegments, nlayers,  & ! Inputs
        nstokes, nmoments, seglayer_map, segupdown_map, layermask_map,   & ! Inputs 
        nfinedivs, ntraverse, radii, Raycon, cota, sunpaths,             & ! Inputs 
        Mu0, Genspher_up, Genspher_dn, Rotations_up, Rotations_dn,       & ! Inputs
        ntraversefine, wfine, csqfine, cotfine, rsegsfine, sunpathsfine, & ! Inputs
        do_Column_Wfs, Lvaryflags, Lvarynums, Lvarymoms,                 & ! Lin Control
        reflec, extinction, omega, truncfac, greekmat, flux,             & ! Inputs
        L_extinction, L_omega, L_truncfac, L_greekmat)                     ! Inputs (Optical, linearized)
   endif

!  Zero the output

   Stokes      = zero ; lostrans   = zero 
   sfi_factors = zero ; sfi_opdeps = zero

   LC_Stokes     = zero ; L_lostrans   = zero
   L_sfi_factors = zero ; L_sfi_opdeps = zero

!  7/21/21, 8/18/21. Mk2. Initialize endsuntrans and LC_Endsuntrans (the direct-sun transmittance at camera)
!  5/30/22. replace endsuntrans & LC_Endsuntrans with Output direct flux, and zero them here

   OutputDirect_Regular_flux = zero ; LC_OutputDirect_Regular_flux = zero
   OutputDirect_Actinic_flux = zero ; LC_OutputDirect_Actinic_flux = zero

!  Zero the local sources and scattering matrices

   sources      = zero ; fmat         = zero
   exactscat_up = zero ; exactscat_dn = zero

   L_sources      = zero
   L_exactscat_up = zero ; L_exactscat_dn = zero

!  4/24/24.  Initialize Output direct fluxes_SB arrays

   OutputDirect_Regular_fluxes_SB = zero
   OutputDirect_Actinic_fluxes_SB = zero
   LC_OutputDirect_Regular_fluxes_SB = zero
   LC_OutputDirect_Actinic_fluxes_SB = zero

!  shorthand

   ns = nstokes

!  Local Column Jacobian control

   Qvary = .false. ; QNums = 0
   if ( do_column_wfs ) then
      Qvary(1:nlayers) = Lvaryflags(1:nlayers)
      QNums(1:nlayers) = Lvarynums (1:nlayers)
   endif
   n_column_wfs = maxval(QNums(1:nlayers))

!  TMS factors and linearization

   do n = 1, nlayers
      if ( do_deltam_scaling ) then
         help = one - truncfac(n) * omega(n)
         tms(n) = omega(n) / help
         if ( Qvary(n) ) then
            do q = 1, Qnums(n)
               L_help = - L_truncfac(n,q)*omega(n) - truncfac(n) * L_omega(n,q)
               L_tms(n,q) = ( L_omega(n,q) - tms(n)*L_help ) / help
            enddo
         endif
      else
         tms(n) = omega(n)
         if ( Qvary(n) ) then
            do q = 1, Qnums(n)
               L_tms(n,q) = L_omega(n,q)
            enddo
         endif
      endif
   enddo

!  Scattering functions and linearization
!  --------------------------------------

!  Scalar only

   if ( nstokes .eq. 1 ) then
      do n = 1, nlayers
         if ( layermask_map(n,1) ) then
            sum1 = DOT_PRODUCT(GenSpher_up(0:nmoments,1), Greekmat(n,0:nmoments,1,1))
            exactscat_up(n,1,1) = sum1 * tms(n)
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  L_exactscat_up(n,1,1,q) = sum1 * L_tms(n,q)
                  if ( Lvarymoms(n,q) ) then
                     L_sum = DOT_PRODUCT(GenSpher_up(0:nmoments,1), L_Greekmat(n,0:nmoments,1,1,q))
                     L_exactscat_up(n,1,1,q) = L_exactscat_up(n,1,1,q) + L_sum * tms(n)
                  endif
               enddo
            endif
         endif
         if ( layermask_map(n,2) ) then
            sum1 = DOT_PRODUCT(GenSpher_dn(0:nmoments,1), Greekmat(n,0:nmoments,1,1))
            exactscat_dn(n,1,1) = sum1 * tms(n)
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  L_exactscat_dn(n,1,1,q) = sum1 * L_tms(n,q)
                  if ( Lvarymoms(n,q) ) then
                     L_sum = DOT_PRODUCT(GenSpher_dn(0:nmoments,1), L_Greekmat(n,0:nmoments,1,1,q))
                     L_exactscat_dn(n,1,1,q) = L_exactscat_dn(n,1,1,q) + L_sum * tms(n)
                  endif
               enddo
            endif
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
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  L_exact(1:ns) = exactscat_up(n,1:ns,1) * L_tms(n,q)
                  if ( Lvarymoms(n,q) ) then
                     L_fmat(1) = DOT_PRODUCT(GenSpher_up(0:nmoments,1), L_Greekmat(n,0:nmoments,1,1,q))
                     L_fmat(2) = DOT_PRODUCT(GenSpher_up(0:nmoments,2), L_Greekmat(n,0:nmoments,1,2,q))
                     L_exactscat_up(n,1,1,q) = + L_fmat(1)
                     L_exactscat_up(n,2,1,q) = - L_fmat(2) * Rotations_up(3)
                     L_exactscat_up(n,3,1,q) = + L_fmat(2) * Rotations_up(4)
                     L_exactscat_up(n,1:ns,1,q) = tms(n) * L_exactscat_up(n,1:ns,1,q) + L_exact(1:ns)
                  else
                     L_exactscat_up(n,1:ns,1,q) = L_exact(1:ns)
                  endif
               enddo
            endif
            exactscat_up(n,1:ns,1) = tms(n) * exactscat_up(n,1:ns,1) 
!            write(*,*)n,exactscat_up(n,3,1),Rotations_up(4),Rotations_up(3)
         endif
         if ( layermask_map(n,2) ) then
            fmat(1) = DOT_PRODUCT(GenSpher_dn(0:nmoments,1), Greekmat(n,0:nmoments,1,1))
            fmat(2) = DOT_PRODUCT(GenSpher_dn(0:nmoments,2), Greekmat(n,0:nmoments,1,2))
            exactscat_dn(n,1,1) = + fmat(1)
            exactscat_dn(n,2,1) = - fmat(2) * Rotations_dn(3)
            exactscat_dn(n,3,1) = + fmat(2) * Rotations_dn(4)
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  L_exact(1:ns) = exactscat_dn(n,1:ns,1) * L_tms(n,q)
                  if ( Lvarymoms(n,q) ) then
                     L_fmat(1) = DOT_PRODUCT(GenSpher_dn(0:nmoments,1), L_Greekmat(n,0:nmoments,1,1,q))
                     L_fmat(2) = DOT_PRODUCT(GenSpher_dn(0:nmoments,2), L_Greekmat(n,0:nmoments,1,2,q))
                     L_exactscat_dn(n,1,1,q) = + L_fmat(1)
                     L_exactscat_dn(n,2,1,q) = - L_fmat(2) * Rotations_dn(3)
                     L_exactscat_dn(n,3,1,q) = + L_fmat(2) * Rotations_dn(4)
                     L_exactscat_dn(n,1:ns,1,q) = tms(n) * L_exactscat_dn(n,1:ns,1,q) + L_exact(1:ns)
                  else
                     L_exactscat_dn(n,1:ns,1,q) = L_exact(1:ns)
                  endif
               enddo
            endif
            exactscat_dn(n,1:ns,1) = tms(n) * exactscat_dn(n,1:ns,1) 
 !           write(*,*)n,exactscat_dn(n,3,1),Rotations_dn(4),Rotations_dn(3)
         endif
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

   Attenuations = ZERO     ; Attenuationsfine = zero ;  Suntau = zero 
   L_Attenuations = zero   ; L_Attenuationsfine  = zero

!  Attenuations to End points (including TOA).

   do s = 0, nsegments
      nt = ntraverse(s) ; sumd = DOT_PRODUCT(extinction(1:nt),sunpaths(s,1:nt))
      suntau(s) = sumd  ; If (sumd .lt. cutoff ) Attenuations(s) = exp( - sumd )
      if ( do_column_wfs ) then
         do q = 1, n_column_wfs
            L_suntau = zero
            do k = 1, nlayers
               if (  k.le.ntraverse(s) ) then
                  L_suntau = L_suntau + L_extinction(k,q) * sunpaths(s,k)
               endif
               L_Attenuations(s,q) = - Attenuations(s) * L_suntau
            enddo
         enddo
      endif
   enddo

!  7/21/21, 8/19/21. Mk2. Endsuntrans is just the end-point Attenuation
!   -- Add the derivatives here.

!  5/30/22. Derive Output flux only at final segment boundary point (s = nsegments)
!   -- Need the attenuations (endsuntrans), muend and Flux/Fluxvec

   if ( do_Calculate_flux .and. Regime_idx .eq. 3 ) then
     endsuntrans = Attenuations(nsegments)
     muend      = cos ( theta(nsegments) )
     OutputDirect_Actinic_flux(1:nstokes) = Flux * endsuntrans * Fluxvec(1:nstokes) 
     OutputDirect_Regular_flux(1:nstokes) = Pi4 * OutputDirect_Actinic_flux(1:nstokes) * muend
     if ( do_Column_wfs ) then
       do q = 1, n_column_wfs
           LC_endsuntrans = L_Attenuations(nsegments,q)
           LC_OutputDirect_Actinic_flux(q,1:nstokes) = Flux * LC_endsuntrans * Fluxvec(1:nstokes)
           LC_OutputDirect_Regular_flux(q,1:nstokes) = Pi4 * LC_OutputDirect_Actinic_flux(q,1:nstokes) * muend
       enddo
     endif
   endif

!  Old code, pre 5/30/22
!   endsuntrans = Attenuations(nsegments)
!   if ( do_column_wfs ) then
!      do q = 1, n_column_wfs
!          LC_endsuntrans(q) = L_Attenuations(nsegments,q)
!      enddo
!   endif

!  fine-layer attenuations

   do s = 1, nsegments
      n = seglayer_map(s)
      do j = 1, nfinedivs(s)
         nt = ntraversefine(s,j) ; sumd = DOT_PRODUCT(extinction(1:nt),sunpathsfine(s,1:nt,j))
         if (sumd .lt. cutoff ) Attenuationsfine(s,j) = exp( - sumd )
         if ( do_column_wfs ) then
            do q = 1, n_column_wfs
               L_sumd = zero
               do k = 1, nlayers
                  L_sumd = L_sumd + L_extinction(k,q) * sunpathsfine(s,k,j)
               enddo
               L_Attenuationsfine(s,j,q) = - Attenuationsfine(s,j) * L_sumd 
            enddo
         endif
      enddo
   enddo

!  5/30/22. Compute flux at top and bottom, debug output (same style as MS routine) to fort.91 and fort.92
!     Use direct flux to the segment boundary points.

   if ( do_Calculate_flux .and. do_debug_flux .and. Regime_idx.eq.3 ) then

!  Top of first segment

     endsuntrans = Attenuations(0)
     muend       = cos (theta(0))
     Actinic_flux(1) = Flux * endsuntrans * Fluxvec(1) 
     Regular_flux(1) = Pi4 * Actinic_flux(1) * muend
     if ( do_column_wfs ) then
       do q = 1, n_column_wfs
         LC_endsuntrans = L_Attenuations(0,q)
         LC_Actinic_flux(q,1) = Flux * LC_endsuntrans * Fluxvec(1)
         LC_Regular_flux(q,1) = Pi4 * LC_Actinic_flux(q,1) * muend
       enddo
     endif

!  4/24/24. Copy these fluxes to the output "Fluxes_SB" arrays

     OutputDirect_Actinic_fluxes_SB(1,0) = Actinic_flux(1)
     OutputDirect_Regular_fluxes_SB(1,0) = Regular_flux(1)
     if ( do_column_wfs ) then
       do q = 1, n_column_wfs
         LC_OutputDirect_Actinic_fluxes_SB(1,0,q) = LC_Actinic_flux(q,1) 
         LC_OutputDirect_Regular_fluxes_SB(1,0,q) = LC_Regular_flux(q,1) 
       enddo
     endif

!  Main segment loop

     do s = 1, nsegments

!  Bottom of segment

       endsuntrans = Attenuations(s)
       muend       = cos (theta(s))
       Actinic_flux(2) = Flux * endsuntrans * Fluxvec(1) 
       Regular_flux(2) = Pi4 * Actinic_flux(2) * muend
       if ( do_column_wfs ) then
         do q = 1, n_column_wfs
           LC_endsuntrans = L_Attenuations(s,q)
           LC_Actinic_flux(q,2) = Flux * LC_endsuntrans * Fluxvec(1)
           LC_Regular_flux(q,2) = Pi4 * LC_Actinic_flux(q,2) * muend
         enddo
       endif

!  4/24/24. Copy these fluxes to the output "Fluxes_SB" arrays
!   -- replaces the debug write-up to fort.91 and fort.92

       OutputDirect_Actinic_fluxes_SB(1,s) = Actinic_flux(2)
       OutputDirect_Regular_fluxes_SB(1,s) = Regular_flux(2)
       if ( do_column_wfs ) then
         do q = 1, n_column_wfs
           LC_OutputDirect_Actinic_fluxes_SB(1,s,q) = LC_Actinic_flux(q,2) 
           LC_OutputDirect_Regular_fluxes_SB(1,s,q) = LC_Regular_flux(q,2) 
         enddo
       endif

!  write-up top/bottom.  Rob debug

       !if ( s.eq.1) write(91,'(/A/)')' segm. layer    ActFlux @TOP   ActFlux @BOT    RegFlux @TOP   RegFlux @BOT'
       !write(91,455)s,s,Actinic_flux(1),Actinic_flux(2),Regular_flux(1),Regular_flux(2)
       !if ( do_Column_Wfs ) then
       !   write(92,'(/A/)')' segm. layer    ActFlux @TOP   ActFlux @BOT    RegFlux @TOP   RegFlux @BOT'
       !   write(92,455)s,s,Actinic_flux(1),Actinic_flux(2),Regular_flux(1),Regular_flux(2) 
       !   write(92,'(/A/)')' Column Jacobians, ActFlux @TOP ==>'
       !   write(92,456)k,LC_Actinic_flux(1:n_column_wfs,1)
       !   write(92,'(/A/)')' Column Jacobians, ActFlux @BOT ==>'
       !   write(92,456)k,LC_Actinic_flux(1:n_column_wfs,2)
       !   write(92,'(/A/)')' Column Jacobians, RegFlux @TOP ==>'
       !   write(92,456)k,LC_Regular_flux(1:n_column_wfs,1)
       !   write(92,'(/A/)')' Column Jacobians, RegFlux @BOT ==>'
       !   write(92,456)k,LC_Regular_flux(1:n_column_wfs,2)
       !endif

!  format for this Flux debug
!455   format(i4,3x,i4,2(1x,1p2e15.5))
!456   format(T10,1p12e15.5)

!  re-set the segment-top values

       Actinic_flux(1) = Actinic_flux(2) 
       Regular_flux(1) = Regular_flux(2)
       if ( do_Column_Wfs ) then
         LC_Actinic_flux(1:n_column_wfs,1) = LC_Actinic_flux(1:n_column_wfs,2) 
         LC_Regular_flux(1:n_column_wfs,1) = LC_Regular_flux(1:n_column_wfs,2)
       endif

!  End segment loop

     enddo

!  End debug clause

   endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Layer integrated Multipliers
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
         if ( do_column_wfs ) then
            do q = 1, n_column_wfs
               L_lostrans(s,q) = - L_extinction(n,q) * cons * lostrans(s)
            enddo
         endif

!  Multiplier

         do j = 1, nfinedivs(s)
            argum(j) = Raycon * ( cot_2 - cotfine(s,j) )  !  Up
            tran(j)  = exp ( - kn * argum(j) )
            sfi_opdeps (s,j) = kn * ( radii(n-1) - rsegsfine(s,j) )
            sfi_factors(s,j) = ke * csqfine(s,j) * tran(j) * wfine(s,j)
         enddo
         multiplier(s) = DOT_PRODUCT(attenuationsfine(s,1:js),sfi_factors(s,1:js))

         if ( do_column_wfs ) then
            do q = 1, n_column_wfs
               L_multiplier(s,q) = DOT_PRODUCT(L_attenuationsfine(s,1:js,q),sfi_factors(s,1:js)) 
               L_kn = L_extinction(n,q) ; L_kn_reg = L_kn/kn
               do j = 1, nfinedivs(s)
                  L_tran = - argum(j) * L_kn
                  L_sfi_opdeps (s,j,q) = sfi_opdeps (s,j) * L_kn_reg
                  L_sfi_factors(s,j,q) = sfi_factors(s,j) * ( L_kn_reg - argum(j) * L_kn )
               enddo
               L_multiplier(s,q)  = L_multiplier(s,q) + &
                                  DOT_PRODUCT(attenuationsfine(s,1:js),L_sfi_factors(s,1:js,q))
            enddo
         endif

!  end

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
         if ( do_column_wfs ) then
            do q = 1, n_column_wfs
               L_lostrans(s,q) = - L_extinction(n,q) * cons * lostrans(s)
            enddo
         endif

!write(*,*)S,lostrans(s),raycon, cons
!  Multiplier

         do j = 1, nfinedivs(s)
            tran(j)  = exp ( - kn * argum(j) )
            sfi_opdeps (s,j) = kn * ( radii(n-1) - rsegsfine(s,j) )
            sfi_factors(s,j) = ke * csqfine(s,j) * tran(j) * wfine(s,j)
         enddo
         multiplier(s) = DOT_PRODUCT(attenuationsfine(s,1:js),sfi_factors(s,1:js))
         if ( do_column_wfs ) then
            do q = 1, n_column_wfs
               L_multiplier(s,q) = DOT_PRODUCT(L_attenuationsfine(s,1:js,q),sfi_factors(s,1:js))
               L_kn = L_extinction(n,q) ; L_kn_reg = L_kn/kn
               do j = 1, nfinedivs(s)
                  L_tran = - argum(j) * L_kn
                  L_sfi_opdeps (s,j,q) = sfi_opdeps (s,j) * L_kn_reg
                  L_sfi_factors(s,j,q) = sfi_factors(s,j) * ( L_kn_reg - argum(j) * L_kn )
               enddo
               L_multiplier(s,q)  = L_multiplier(s,q) + &
                                  DOT_PRODUCT(attenuationsfine(s,1:js),L_sfi_factors(s,1:js,q))
            enddo
         endif

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
         if ( do_column_wfs ) then
            do q = 1, n_column_wfs
               L_lostrans(s,q) = - L_extinction(n,q) * cons * lostrans(s)
            enddo
         endif

!  Multiplier

         do j = 1, nfinedivs(s)
            argum(j) = Raycon * ( cotfine(s,j) - cot_2 )
            tran(j)  = exp ( - kn * argum(j) )
            sfi_opdeps (s,j) = kn * ( radii(n-1) - rsegsfine(s,j) )
            sfi_factors(s,j) = ke * csqfine(s,j) * tran(j) * wfine(s,j)
         enddo
         multiplier(s) = DOT_PRODUCT(attenuationsfine(s,1:js),sfi_factors(s,1:js))
         if ( do_column_wfs ) then
            do q = 1, n_column_wfs
               L_multiplier(s,q) = DOT_PRODUCT(L_attenuationsfine(s,1:js,q),sfi_factors(s,1:js)) 
               L_kn = L_extinction(n,q) ; L_kn_reg = L_kn/kn
               do j = 1, nfinedivs(s)
                  L_tran = - argum(j) * L_kn
                  L_sfi_opdeps (s,j,q) = sfi_opdeps (s,j) * L_kn_reg
                  L_sfi_factors(s,j,q) = sfi_factors(s,j) * ( L_kn_reg - argum(j) * L_kn )
               enddo
               L_multiplier(s,q)  = L_multiplier(s,q) + &
                                  DOT_PRODUCT(attenuationsfine(s,1:js),L_sfi_factors(s,1:js,q))
            enddo
         endif

!  End

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
         if ( do_column_wfs ) then
            do q = 1, n_column_wfs
               L_sources(s,1:ns,q) = exactscat_up(n,1:ns,1) * L_multiplier(s,q)
               L_sources(s,1:ns,q) = L_sources(s,1:ns,q) + &
                             L_exactscat_up(n,1:ns,1,q) * multiplier(s)
            enddo
         endif
      enddo
   endif

!  Regime 2. Limb path

   if ( regime_idx .eq. 2 .and. do_sunlight ) then
      do s = 1, nsegments
         n = seglayer_map(s)
         if ( .not.segUpdown_map(s) ) then
            sources(s,1:ns) = exactscat_dn(n,1:ns,1) * multiplier(s)
            if ( do_column_wfs ) then
               do q = 1, n_column_wfs
                  L_sources(s,1:ns,q) = exactscat_dn(n,1:ns,1) * L_multiplier(s,q)
                  L_sources(s,1:ns,q) = L_sources(s,1:ns,q) + &
                               L_exactscat_dn(n,1:ns,1,q) * multiplier(s)
               enddo
            endif
         else
            sources(s,1:ns) = exactscat_up(n,1:ns,1) * multiplier(s)
            if ( do_column_wfs ) then
               do q = 1, n_column_wfs
                  L_sources(s,1:ns,q) = exactscat_up(n,1:ns,1) * L_multiplier(s,q)
                  L_sources(s,1:ns,q) = L_sources(s,1:ns,q) + &
                               L_exactscat_up(n,1:ns,1,q) * multiplier(s)
               enddo
            endif
         endif
      enddo
   endif

!  Regime 3. Downwelling from TOA

   if ( regime_idx .eq. 3 .and. do_sunlight ) then
      do s = 1, nsegments
         n = seglayer_map(s)
         sources(s,1:ns) = exactscat_dn(n,1:ns,1) * multiplier(s)
         if ( do_column_wfs ) then
            do q = 1, n_column_wfs
               L_sources(s,1:ns,q) = exactscat_dn(n,1:ns,1) * L_multiplier(s,q)
               L_sources(s,1:ns,q) = L_sources(s,1:ns,q) + &
                             L_exactscat_dn(n,1:ns,1,q) * multiplier(s)
            enddo
         endif
      enddo
   endif

!  Source function integration
!  ===========================

!  RT recursion

   CUMSOURCE_SS(0,:) = zero
   do s = 1, nsegments
      CUMSOURCE_SS(s,1:ns) = LOSTRANS(S) * CUMSOURCE_SS(s-1,1:ns) + SOURCES(S,1:ns)
   enddo
   STOKES(1:ns) = FLUX * CUMSOURCE_SS(nsegments,1:ns)
   
!  Add surface term for Regime 1. Valid for BRDF as well.
!  Assumes natural light only (unpolarized solar beam)

   CUMSOURCE_DB(0,:) = zero ; M4 = zero ; s0 =  0
   if ( Regime_idx .eq. 1 ) then
      M4 = four * MU0 ; CUMSOURCE_DB(0,1:ns) = M4 * REFLEC(1:ns,1) * attenuations(s0)
      do s = 1, nsegments
         CUMSOURCE_DB(s,1:ns) = LOSTRANS(S) * CUMSOURCE_DB(s-1,1:ns)
      enddo
      STOKES(1:ns) = STOKES(1:ns) + FLUX * CUMSOURCE_DB(nsegments,1:ns)
   endif

!  Add surface term for Regime 1. Valid for BRDF as well.
!  Assumes natural light only (unpolarized solar beam)

!  Column Wfs (atmospheric term). All Regimes

   if ( do_column_wfs ) then
      L_CUMSOURCE_SS = zero
      do s = 1, nsegments
         n = seglayer_map(s)
         do q = 1, n_column_wfs
            L_cumsource_ss(1:ns,q) = L_SOURCES(S,1:ns,Q) +   LOSTRANS(S)   * L_CUMSOURCE_SS(1:ns,Q) &
                                                         + L_LOSTRANS(S,Q) *   CUMSOURCE_SS(S-1,1:ns)
         enddo
      enddo
      do q = 1, n_column_wfs
         LC_STOKES(1:ns,Q) = FLUX * L_CUMSOURCE_SS(1:ns,Q)
      enddo
   endif

!  Add Column Wfs (Surface direct beam term, Regime 1 only). Valid for BRDF as well.
!  Assumes natural light only (unpolarized solar beam)

   L_CUMSOURCE_DB = zero ; s0 = 0
   if ( Regime_idx .eq. 1 .and. do_column_wfs ) then
      Term1(1:ns) = M4 * reflec(1:ns,1)
      do q = 1, n_column_wfs
         L_CUMSOURCE_DB(1:ns,q) = Term1(1:ns) * L_attenuations(s0,q)
      enddo
      do s = 1, nsegments
         n = seglayer_map(s)
         do q = 1, n_column_wfs
               L_cumsource_db(1:ns,q) = LOSTRANS(S) * L_CUMSOURCE_db(1:ns,Q) &
                                  + L_LOSTRANS(S,Q) *   CUMSOURCE_db(S-1,1:ns)
         enddo
      enddo
      do q = 1, n_column_wfs
         LC_STOKES(1:ns,Q) = LC_STOKES(1:ns,Q) + FLUX * L_CUMSOURCE_db(1:ns,Q)
      enddo
   endif

!  Finish

   return
end subroutine FOMS_FO_LC_VectorRT

! PRIVATE ROUTINE

subroutine FOMS_FO_LC_Write_Input &
   ( do_sunlight, do_deltam_scaling, Regime_idx, nsegments, nlayers,  & ! Inputs
     nstokes, nmoments, seglayer_map, segupdown_map, layermask_map,   & ! Inputs 
     nfinedivs, ntraverse, radii, Raycon, cota, sunpaths,             & ! Inputs 
     Mu0, Genspher_up, Genspher_dn, Rotations_up, Rotations_dn,       & ! Inputs
     ntraversefine, wfine, csqfine, cotfine, rsegsfine, sunpathsfine, & ! Inputs
     do_Column_Wfs, Lvaryflags, Lvarynums, Lvarymoms,                 & ! Lin Control
     reflec, extinction, omega, truncfac, greekmat, flux,             & ! Inputs (Optical)
     L_extinction, L_omega, L_truncfac, L_greekmat)                     ! Inputs (Optical, linearized)

!  Check inputs if desired
!mick mod 12/1/2021  - added this subroutine

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

!  Linearized
!  ----------

!  Jacobian control

   logical, intent(in)    :: do_Column_Wfs
   logical, intent(in)    :: Lvaryflags(maxlays)
   integer, intent(in)    :: Lvarynums (maxlays)
   logical, intent(in)    :: Lvarymoms (maxlays,maxJacobs)

!  Linearized optical inputs

   real(ffp), intent(in)  :: L_extinction  ( maxlays, maxjacobs )
   real(ffp), intent(in)  :: L_omega       ( maxlays, maxjacobs )
   real(ffp), intent(in)  :: L_truncfac    ( maxlays, maxjacobs )
   real(ffp), intent(in)  :: L_greekmat    ( maxlays,0:maxmoms, 4, 4, maxJacobs )

!  Local variables
!  ---------------

   integer :: outunit, n_wfs
   integer :: lay,mom,nf,ngreek,o1,o2,s,ss,wf
!   integer :: geo,i,j,k,ulev,ura,uva,uog,sza,ngeoms


   integer :: gmask(8)
   data gmask / 1, 2, 5, 6, 11, 12, 15, 16 /

   integer :: gmask2(8,2)

!  Open output file

   OUTUNIT = 101
   OPEN (OUTUNIT,file = 'Debug/FOMS_FO_LC_Write_Input.dbg',status = 'replace')

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

   n_wfs = 0
   if (Do_Column_WFs) n_wfs = maxval(Lvarynums)

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

   WRITE(OUTUNIT,*)
   do o1=1,4
     do mom=0,nmoments
       WRITE(OUTUNIT,DWFR2) 'mom = ',mom,' o1 = ',o1,' genspher_up(mom,o1) = ',genspher_up(mom,o1)
     enddo
   enddo

   WRITE(OUTUNIT,*)
   do o1=1,4
     do mom=0,nmoments
       WRITE(OUTUNIT,DWFR2) 'mom = ',mom,' o1 = ',o1,' genspher_dn(mom,o1) = ',genspher_dn(mom,o1)
     enddo
   enddo

   WRITE(OUTUNIT,*)
   do o1=1,4
     WRITE(OUTUNIT,DWFR1) 'o1 = ',o1,' rotations_up(o1) = ',rotations_up(o1)
   enddo

   WRITE(OUTUNIT,*)
   do o1=1,4
     WRITE(OUTUNIT,DWFR1) 'o1 = ',o1,' rotations_dn(o1) = ',rotations_dn(o1)
   enddo

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

   WRITE(OUTUNIT,*)
   do s=1,nsegments
     do nf=1,nfinedivs(s)
       WRITE(OUTUNIT,DWFI2) 's = ',s,' nf = ',nf,' ntraversefine(s,nf) = ',ntraversefine(s,nf)
     enddo
   enddo

   WRITE(OUTUNIT,*)
   do lay=1,nlayers
     do s=1,nsegments
       do nf=1,nfinedivs(s)
         WRITE(OUTUNIT,DWFR3) 'lay = ',lay,' s = ',s,' nf = ',nf,&
                             ' sunpathsfine(lay,s,nf) = ',sunpathsfine(lay,s,nf)
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

!  Linearized Inputs
!  -----------------

   WRITE(OUTUNIT,*)
   WRITE(OUTUNIT,'(A)') '-----------------'
   WRITE(OUTUNIT,'(A)') 'Linearized Inputs'
   WRITE(OUTUNIT,'(A)') '-----------------'

   WRITE(OUTUNIT,*)
   WRITE(OUTUNIT,DWFL) 'do_Column_Wfs = ',do_Column_Wfs

   if (do_Column_Wfs) then

!  Jacobian control

     WRITE(OUTUNIT,*)
     do lay=1,nlayers
       WRITE(OUTUNIT,DWFL1) 'lay = ',lay,' Lvaryflags(lay) = ',Lvaryflags(lay)
     enddo

     WRITE(OUTUNIT,*)
     do lay=1,nlayers
       WRITE(OUTUNIT,DWFI1) 'lay = ',lay,' Lvarynums(lay) = ',Lvarynums(lay)
     enddo

     WRITE(OUTUNIT,*)
     do wf=1,n_wfs
       do lay=1,nlayers
         WRITE(OUTUNIT,DWFL2) 'wf = ',wf,' lay = ',lay,&
                             ' Lvarymoms(lay,wf) = ',Lvarymoms(lay,wf)
       enddo
     enddo

!  Optical inputs

     WRITE(OUTUNIT,*)
     do wf=1,n_wfs
       do lay=1,nlayers
         WRITE(OUTUNIT,DWFR2) 'wf = ',wf,' lay = ',lay,&
                             ' L_extinction(lay,wf) = ',L_extinction(lay,wf)
       enddo
     enddo

     WRITE(OUTUNIT,*)
     do wf=1,n_wfs
       do lay=1,nlayers
         WRITE(OUTUNIT,DWFR2) 'wf = ',wf,' lay = ',lay,&
                             ' L_omega(lay,wf)      = ',L_omega(lay,wf)
       enddo
     enddo

     WRITE(OUTUNIT,*)
     do wf=1,n_wfs
       do lay=1,nlayers
         WRITE(OUTUNIT,DWFR2) 'wf = ',wf,' lay = ',lay,&
                             ' L_truncfac(lay,wf)   = ',L_truncfac(lay,wf)
       enddo
     enddo

     do wf=1,n_wfs
       do ss=1,ngreek
         o1 = gmask2(ss,1); o2 = gmask2(ss,2); WRITE(OUTUNIT,*)
         do lay=1,nlayers
           do mom=0,nmoments
             WRITE(OUTUNIT,DWFR6) 'wf = ',wf,' ss = ',ss,' o1 = ',o1,' o2 = ',o2,' lay = ',lay,' mom = ',mom,&
               ' L_greekmat(lay,mom,o1,o2,wf) = ',L_greekmat(lay,mom,o1,o2,wf)
           enddo
         enddo
       enddo
     enddo
   endif

!  Close output file

   CLOSE (OUTUNIT)

end subroutine FOMS_FO_LC_Write_Input

! End module

end module FOMS_FO_LC_VectorRT_m
