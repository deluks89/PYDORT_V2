
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

!  4/24/24. Added FO (Direct) and MS (Diffuse) Flux outputs at all segment boundaries
!    -- Extra lines of output in FOMS_FO_VectorRT and FOMS_MS_VectorRT_Mk2/3
!    -- Extra outputs are all written to type structure arrays (new in vqs_output_def)

      MODULE vqs_masters_m

      PUBLIC  :: VQS_MASTER
      !PRIVATE ::

      CONTAINS

      SUBROUTINE VQS_MASTER ( do_debug_input, &
        VQS_FixIn, &
        VQS_ModIn, &
        VQS_Out )

!  Parameter files

      USE MSSeg_PARS_m
      USE FOMS_pars_m

!  I/O Type structures

      USE VQS_Inputs_def_m
      USE VQS_Outputs_def_m

!  FOMS modules, Basic set

      USE FOMS_Geometry_Masters_m
      USE FOMS_VectorSpherfuncs_m
      USE FOMS_FO_VectorRT_m
      USE FOMS_MS_VectorRT_Mk2_Alt_m
      USE FOMS_MS_VectorRT_Mk3_Exp_m

      USE VQS_INPUTS_m
      USE VQS_WRITEMODULES_m

      IMPLICIT NONE

!  VQS I/O
!  =======

      LOGICAL, INTENT (IN) :: DO_DEBUG_INPUT

!  VQS input structures

      TYPE(VQS_Fixed_Inputs), INTENT (IN)       :: VQS_FixIn
      TYPE(VQS_Modified_Inputs), INTENT (INOUT) :: VQS_ModIn

!  VQS output structure

      TYPE(VQS_Outputs), INTENT (OUT)           :: VQS_Out

!  Inputs (Local Versions)
!  ======

!  General
!  -------

!  Dimensions

      integer     :: Method_Index

      integer     :: nlayers
      integer     :: nstokes
      integer     :: nstreams

!  FO Geometry
!  -----------

!  Control

      logical     :: do_debug_FO_geometry

!  Geophysical

      real(fpk)   :: FO_pradius, FO_heights ( 0:maxlays )

!  Camera

      real(fpk)   :: height_Camera
      real(fpk)   :: theta_Camera
      real(fpk)   :: relaz_Camera
      real(fpk)   :: alpha_Camera

!  Number of fine divisions

      integer     :: nfine_input

!  Flag for using Critical Attenuation + threshold

      logical     :: doCrit
      real(fpk)   :: ACrit

!  Filename for Geometry output

      character(len=100):: Geo_Filename

!  Both FO & MS codes
!  ------------------

!  Control

      logical     :: do_deltam_scaling
      logical     :: do_calculate_fluxes
      logical     :: do_debug_fluxes

!  Atmosphere optical inputs

      real(fpk)   :: FO_omega ( maxlays )
      real(fpk)   :: FO_greekmat ( maxlays, 0:maxmoms, maxstokes, maxstokes )

!  Surface reflectivity (Regime 1 only)

      real(fpk)   :: FO_reflec ( maxstokes, maxstokes )

!  FO code only
!  ------------

!  Control

      logical     :: do_sunlight
      logical     :: do_debug_FO_input

!  Solar Flux

      real(fpk)   :: solarflux, fluxvec ( maxstokes )

!  Atmosphere optical inputs

      integer     :: FO_nmoments
      real(fpk)   :: FO_extinction ( maxlays )
      real(fpk)   :: FO_truncfac   ( maxlays )
      real(fpk)   :: extinction_m  ( maxlays )

!  MS code only
!  ------------

!  Control

      logical     :: do_MSRT, RT_PROGRESS
      logical     :: do_debug_MS_input

      logical     :: do_plane_parallel
      logical     :: do_Rayleigh_only
      logical     :: do_lambertian_surface

!  Control (not currently active)

      logical     :: do_sleave
      logical     :: do_sl_isotropic

!  Solar Flux

      real(fpk)   :: flux_factor

!  Atmosphere optical inputs

      real(fpk)   :: FO_deltaus ( maxlays )

!  VQS Diagnostic (RT Solutions Use Only)

      integer     :: atmos_index

!  Outputs (Local Versions)
!  =======

!  Exception handling
!  ------------------

!  Exception handling for input checking.
!     Message Length should be at least 120 Characters

      integer      :: status_inputcheck
      integer      :: ncheckmessages

      character (len=120), dimension(0:max_messages)  :: checkmessages
      character (len=120), dimension(0:max_messages)  :: actions

!  Exception handling for model calculation.
!mick note 6/6/2023 - this arrangement is a little different than VLIDORT

!      integer      :: status_calculation
      integer      :: nmessages
      character (len=120), dimension(0:max_messages)  :: messages

      logical :: verbose = .false.
      !logical :: verbose = .true.

!  Internal
!  ========

!  Optical properties
!  ------------------

      real(fpk)  :: dnm4, tms

!  FO Geometry Output
!  ------------------

!  Regime index, number of segments

      integer    :: nsegments, Regime_idx
      logical    :: Exact_condition

!  Segment to layer map, segment up/down switch, Layer-mask (1=up,2=dn)

      integer    :: seglayer_map  ( maxsegments )
      logical    :: segupdown_map ( maxsegments )
      logical    :: layermask_map ( maxlays, 2 )

!  Radii, Ray constant, Lospaths

      real(fpk)  :: radii    ( 0:maxlays )
      real(fpk)  :: Raycon
      real(fpk)  :: Lospaths ( maxsegments )
      real(fpk)  :: Rsegs    ( 0:maxsegments )

!  LOS angles and trig ratios thereof
!    Convention, range always 0 to 90 degrees 

      real(fpk)  :: alpha    ( 0:maxsegments )
      real(fpk)  :: cosa     ( 0:maxsegments )
      real(fpk)  :: sina     ( 0:maxsegments )
      real(fpk)  :: cota     ( 0:maxsegments )

!  Critical values

!      logical    :: Crit_map    ( maxsegments )
!      real(fpk)  :: alpha_crits ( 0:maxsegments )
!      real(fpk)  :: Rad_crits   ( 0:maxsegments )

!  Quadrature values

      integer    :: nfinedivs ( maxsegments)
      real(fpk)  :: wfine     ( maxsegments, maxfine )
      real(fpk)  :: csqfine   ( maxsegments, maxfine )
      real(fpk)  :: cotfine   ( maxsegments, maxfine )
      real(fpk)  :: alphafine ( maxsegments, maxfine )
      real(fpk)  :: rsegsfine ( maxsegments, maxfine )

!  Solar Paths to segment boundaries

      integer    :: ntraverse ( 0:maxsegments )
      real(fpk)  :: sunpaths  ( 0:maxsegments, maxlays )

!  Solar Paths to fine levels

      integer    :: ntraversefine ( maxsegments, maxfine )
      real(fpk)  :: sunpathsfine  ( maxsegments, maxlays, maxfine )

!  Scattering angle cosine, segment-boundary and fine-grid angles

      real(fpk)  :: cosscat
      real(fpk)  :: theta     ( 0:maxsegments )
      real(fpk)  :: relaz     ( 0:maxsegments )
      real(fpk)  :: thetafine ( maxsegments, maxfine )
      real(fpk)  :: relazfine ( maxsegments, maxfine )

!  Exception handling

      logical            :: fail_geo
      integer            :: status_geo
      character(len=100) :: message_geo
      character(len=100) :: trace_geo

!  Spherical Functions
!  -------------------

!  Set starter flag before first call

      logical    :: starter

!  Solar zenith angle cosine

      real(fpk)  :: Mu0

!  Vsign

      real(fpk), parameter :: vsign_up =   one
      real(fpk), parameter :: vsign_dn = - one

!  Help

      real(fpk)  :: GSHELP(7,maxmoms)

!  Number of moments and spherical functions

      real(fpk)  :: Genspher_up ( 0:maxmoms, maxstokes )
      real(fpk)  :: Genspher_dn ( 0:maxmoms, maxstokes )
      real(fpk)  :: Rotations_up ( maxstokes )
      real(fpk)  :: Rotations_dn ( maxstokes )

!  FO code output
!  --------------

!  Stokes vector and integration factors

      real(fpk)  :: lostrans    ( maxsegments )
      real(fpk)  :: sfi_factors ( maxsegments, maxfine )
      real(fpk)  :: sfi_opdeps  ( maxsegments, maxfine )

!  MS code output
!  --------------

!  Exception handling
!mick fix 6/6/2023 - changed extent of messages_MS & actions_MS from 25 to 0:MAX_MESSAGES

      logical            :: fail_MS
      integer            :: status_MS
      character(len=100) :: messages_MS ( 0:max_messages )
      character(len=100) :: actions_ms  ( 0:max_messages )
      integer            :: nmessages_MS

!  Misc
!  ----

!  Help variables

      integer :: n, l, o1, o2, mo1, mo2, nstreams_2

!  Timing

      real    :: e1, e2

! -----------------------------------------
! SFC - January, 30th 2018
! Parameter to avoid output formatting error
      real(fpk), parameter :: smallest = 1.0e-50_fpk
! -----------------------------------------

!  Start program
!  =============

!  Initialize some variables
!  -------------------------

!  Input checks

      VQS_Out%Status%STATUS_INPUTCHECK = VQS_SUCCESS

      NCHECKMESSAGES = 0
      CHECKMESSAGES  = ' '
      ACTIONS        = ' '

!  Model calculation

      VQS_Out%Status%STATUS_CALCULATION = VQS_SUCCESS

      NMESSAGES = 0
      MESSAGES  = ' '

!  Check input dimensions
!  ----------------------

      CALL VQS_CHECK_INPUT_DIMS &
      ( VQS_FixIn, VQS_ModIn, &
        STATUS_INPUTCHECK, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )

!  Exception handling

      IF ( STATUS_INPUTCHECK .EQ. VQS_SERIOUS ) THEN
        VQS_Out%Status%STATUS_INPUTCHECK = STATUS_INPUTCHECK
        VQS_Out%Status%NCHECKMESSAGES    = NCHECKMESSAGES
        VQS_Out%Status%CHECKMESSAGES     = CHECKMESSAGES
        VQS_Out%Status%ACTIONS           = ACTIONS
        RETURN
      ENDIF

!  Check regular inputs
!  --------------------

      CALL VQS_CHECK_INPUT &
      ( VQS_FixIn, VQS_ModIn, &
        STATUS_INPUTCHECK, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )

!  Exception handling
!  mick note 6/6/2023 - since this is the last check sub, if all is well, 
!                       define remaining three TS versions of the input check
!                       vars here in an ELSE section

      IF ( STATUS_INPUTCHECK .EQ. VQS_SERIOUS ) THEN
        VQS_Out%Status%STATUS_INPUTCHECK = STATUS_INPUTCHECK
        VQS_Out%Status%NCHECKMESSAGES    = NCHECKMESSAGES
        VQS_Out%Status%CHECKMESSAGES     = CHECKMESSAGES
        VQS_Out%Status%ACTIONS           = ACTIONS
        RETURN
      ELSE
        VQS_Out%Status%NCHECKMESSAGES    = NCHECKMESSAGES
        VQS_Out%Status%CHECKMESSAGES     = CHECKMESSAGES
        VQS_Out%Status%ACTIONS           = ACTIONS
      ENDIF

!  ====================================
!  BEGIN COPY INPUTS TO LOCAL VARIABLES
!  ====================================

!  -----------------------
!  Standard Inputs - Fixed
!  -----------------------

!  VQS_Fixed_Boolean
!  -----------------

      do_MSRT               = VQS_FixIn%Bool%DO_MSRT
      RT_PROGRESS           = VQS_FixIn%Bool%DO_RT_PROGRESS
      do_calculate_fluxes   = VQS_FixIn%Bool%DO_CALCULATE_FLUXES

      do_sunlight           = VQS_FixIn%Bool%DO_SUNLIGHT
      do_plane_parallel     = VQS_FixIn%Bool%DO_PLANE_PARALLEL
      do_lambertian_surface = VQS_FixIn%Bool%DO_LAMBERTIAN_SURFACE

      do_debug_FO_geometry  = VQS_FixIn%Bool%DO_DEBUG_FO_GEOMETRY
      do_debug_FO_input     = VQS_FixIn%Bool%DO_DEBUG_FO_INPUT
      do_debug_MS_input     = VQS_FixIn%Bool%DO_DEBUG_MS_INPUT
      do_debug_fluxes       = VQS_FixIn%Bool%DO_DEBUG_FLUXES

      !(not currently active)
      !do_sleave            = VQS_FixIn%Bool%DO_SLEAVE
      !do_sl_isotropic      = VQS_FixIn%Bool%DO_SL_ISOTROPIC

!  tempo
!       do_calculate_fluxes = .true.
!       do_debug_fluxes     = .true.

!  VQS_Fixed_Control
!  -----------------

      Method_Index = VQS_FixIn%Cont%METHOD_INDEX

      nlayers      = VQS_FixIn%Cont%NLAYERS
      nstokes      = VQS_FixIn%Cont%NSTOKES
      nstreams     = VQS_FixIn%Cont%NSTREAMS
      nfine_input  = VQS_FixIn%Cont%NFINE_INPUT

      doCrit       = VQS_FixIn%Cont%DOCRIT
      ACrit        = VQS_FixIn%Cont%ACRIT

!  VQS_Fixed_Sunrays
!  -----------------

      flux_factor = VQS_FixIn%Sunrays%FLUX_FACTOR

!  VQS_Fixed_UserValues
!  --------------------

      height_Camera = VQS_FixIn%UserVal%HEIGHT_CAMERA

      theta_Camera  = VQS_FixIn%UserVal%THETA_CAMERA
      relaz_Camera  = VQS_FixIn%UserVal%RELAZ_CAMERA
      alpha_Camera  = VQS_FixIn%UserVal%ALPHA_CAMERA

      Geo_Filename  = VQS_FixIn%UserVal%GEO_FILENAME

!  VQS_Fixed_Geophys
!  -----------------

      FO_pradius    = VQS_FixIn%Geophys%PLANET_RADIUS
      FO_heights(0:nlayers) = VQS_FixIn%Geophys%HEIGHT_GRID(0:nlayers)

!  VQS_Fixed_Optical
!  -----------------

      atmos_index   = VQS_FixIn%Optical%ATMOS_INDEX

!  --------------------------
!  Standard Inputs - Modified
!  --------------------------

!  VQS_Modified_Boolean
!  --------------------

      do_Rayleigh_only  = VQS_ModIn%MBool%DO_RAYLEIGH_ONLY
      do_deltam_scaling = VQS_ModIn%MBool%DO_DELTAM_SCALING

!  VQS_Modified_Control
!  --------------------
!  (These type structure inputs used directly)

!  VQS_Modified_Optical
!  --------------------
!  (These type structure inputs used directly)


!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  VQS input debug

      IF ( do_debug_input ) THEN
        CALL VQS_DEBUG_INPUT_MASTER()
      END IF

!  Initialize outputs
!  ------------------

!  Main outputs (radiances and fluxes)

      VQS_Out%Main%StokesTot(1:nstokes)   = ZERO
      VQS_Out%Main%StokesFO (1:nstokes)   = ZERO
      VQS_Out%Main%StokesMS (1:nstokes)   = ZERO

      VQS_Out%Main%ActFluxDiff(1:nstokes) = ZERO
      VQS_Out%Main%RegFluxDiff(1:nstokes) = ZERO
      VQS_Out%Main%ActFluxDir (1:nstokes) = ZERO
      VQS_Out%Main%RegFluxDir (1:nstokes) = ZERO

!  Prepare additional solar inputs
!  -------------------------------

      fluxvec = ZERO; fluxvec(1) = ONE
      solarflux = flux_factor / (4.0_fpk * pie) ! FO

!  Prepare optical properties
!  --------------------------

!  Zero FO inputs

      FO_Extinction = ZERO ; FO_deltaus  = ZERO ; FO_omega  = ZERO
      FO_Greekmat   = ZERO ; FO_truncfac = ZERO ; FO_Reflec = ZERO

!  Set FO optical properties and miscellaneous inputs

      do n = 1, nlayers
        FO_extinction(n) = VQS_FixIn%Optical%DELTAU_VERT_INPUT(n) / &
                          (VQS_FixIn%Geophys%HEIGHT_GRID(n-1) - VQS_FixIn%Geophys%HEIGHT_GRID(n))
        FO_deltaus(n)    = VQS_FixIn%Optical%DELTAU_VERT_INPUT(n)
        FO_omega(n)      = VQS_ModIn%MOptical%OMEGA_TOTAL_INPUT(n)
        do L = 0, VQS_ModIn%MCont%NGREEK_MOMENTS_INPUT
          do o1 = 1, 4 ; mo1 = 4*(o1-1) ; do o2 = 1, 4 ; mo2 = mo1 + o2
            FO_greekmat(N,L,o1,o2) = VQS_FixIn%Optical%GREEKMAT_TOTAL_INPUT(L,N,mo2)
          enddo ; enddo
        enddo
      enddo
      FO_nmoments    = VQS_ModIn%MCont%NGREEK_MOMENTS_INPUT
      FO_reflec(1,1) = VQS_FixIn%Optical%LAMBERTIAN_ALBEDO

!  Set truncation factor
!    - Only set this if multiple scattering is being done

      if ( do_MSRT .and. .not.do_Rayleigh_only ) then
        nstreams_2 = 2 * nstreams ; dnm4 = dble(2*nstreams_2+1)
        do n = 1, nlayers
          FO_truncfac(n) = FO_greekmat(n,nstreams_2,1,1)/ dnm4
          tms = one-FO_omega(n)*FO_truncfac(n)
          extinction_m(n) = FO_extinction(n) * tms
        enddo
      else
        do n = 1, nlayers
          extinction_m(n) = FO_extinction(n)
        enddo
      endif

!  Prepare geometry
!  ----------------

!  FO Geometry call

      call cpu_time(e1)
      call FOMS_Geometry_Master &
       ( nlayers, FO_heights, FO_pradius, height_Camera,      & ! Input
         alpha_Camera, theta_Camera, relaz_Camera,            & ! Input
         nfine_input, doCrit, Acrit, FO_extinction,           & ! Input
         Regime_idx, Exact_condition, Raycon, nsegments,      & ! Output
         seglayer_map, segupdown_map, layermask_map, radii,   & ! Output
         Lospaths, Rsegs, alpha, sina, cosa, cota, nfinedivs, & ! Output
         wfine, rsegsfine, alphafine, csqfine, cotfine,       & ! Output
         sunpaths, ntraverse, sunpathsfine, ntraversefine,    & ! Output
         cosscat, theta, relaz, thetafine, relazfine,         & ! Output
         fail_geo, status_geo, message_geo, trace_geo )         ! Output
      call cpu_time(e2)
      VQS_ModIn%MTiming%GEO_TIME = e2 - e1

!  Exception handling on this FO geometry

      if ( fail_geo ) then
        VQS_Out%Status%status_calculation = VQS_SERIOUS
        VQS_Out%Status%messages(nmessages+1) = Trim(message_geo)
        VQS_Out%Status%messages(nmessages+2) = Trim(trace_geo)
        VQS_Out%Status%nmessages = nmessages+2
        return
      else
        if ( status_geo .eq. 1 ) then
          VQS_Out%Status%status_calculation = VQS_SERIOUS
          VQS_Out%Status%messages(nmessages+1) = ' SZA > 90; cannot proceed with MS calculation, Status = 1'
          VQS_Out%Status%NMessages = nmessages+1
          return
        endif
      endif

!      write(*,'(a,3f10.5,a,i2)')&
!       ' Geometry Angles = ', theta_Camera, alpha_Camera, relaz_Camera,&
!       '     Regime index (1=Nadir-up, 2=Limbpath-up, 3=ToaDown) = ',Regime_idx

!  Write geometry (debug only)

      if ( do_debug_FO_geometry ) then
        call Write_FOMS_Geometry &
          ( Geo_Filename, Regime_idx, nsegments, FO_pradius, cosscat, & ! Inputs
            seglayer_map, segupdown_map, Lospaths, Rsegs,             & ! Inputs
            alpha, theta, relaz, nfinedivs, rsegsfine,                & ! Inputs
            sunpaths, ntraverse, sunpathsfine, ntraversefine )          ! Inputs
      endif

!  Prepare Spherical functions
!  ---------------------------

!  Set starter flag before first call
!    -- 12/10/20. Think there was a fault here in the use of the starter flag
!    -- 12/10/20. Think there was a fault with the VSIGN_DN value in regime 2, should be negative.

      call cpu_time(e1)
        starter = .true.
      if ( Regime_idx .eq. 1 ) then
        call FOMS_VectorSpherfuncs &
           ( STARTER, FO_NMOMENTS, NSTOKES, DO_SUNLIGHT, VSIGN_UP, & ! Inputs
             THETA_Camera, ALPHA_Camera, RELAZ_Camera, COSSCAT,    & ! Inputs
             ROTATIONS_UP, GSHELP, GENSPHER_UP )                     ! Outputs
      else if ( Regime_idx .eq. 3 ) then
        call FOMS_VectorSpherfuncs &
           ( STARTER, FO_NMOMENTS, NSTOKES, DO_SUNLIGHT, VSIGN_DN, & ! Inputs
             THETA_Camera, ALPHA_Camera, RELAZ_Camera, COSSCAT,    & ! Inputs
             ROTATIONS_DN, GSHELP, GENSPHER_DN )                     ! Outputs
      else
        call FOMS_VectorSpherfuncs &
           ( STARTER, FO_NMOMENTS, NSTOKES, DO_SUNLIGHT, VSIGN_UP, & ! Inputs
             THETA_Camera, ALPHA_Camera, RELAZ_Camera, COSSCAT,    & ! Inputs
             ROTATIONS_UP, GSHELP, GENSPHER_UP )                     ! Outputs
        starter = .true.
        call FOMS_VectorSpherfuncs &
           ( STARTER, FO_NMOMENTS, NSTOKES, DO_SUNLIGHT, -VSIGN_DN, & ! Inputs
             THETA_Camera, ALPHA_Camera, RELAZ_Camera, COSSCAT,     & ! Inputs
             ROTATIONS_DN, GSHELP, GENSPHER_DN )                      ! Outputs
      endif

      if (Verbose) write(*,'(/1x,a)') ' - Doing FO RT'

!  Do FO radiative transfer
!  ========================

!  7/21/21. Mk2. Addition of solar beam attenuation at camera position ("endsuntrans")
!    - Necessary for calculating direct downwelling fluxes at camera
!    - Always calculated, so simple to output automaticaaly

!mick mod 12/1/2021 - added "do_debug_FO_input" to inputs

!  4/14/23. Mk3. Flux Calculations now inside the FO routine (just as in MS routine below)
!    - Add "Direct-flux" arrays as Output variables (1 extra lines)
!    - remove endsuntrans output
!    - Add theta inputs (necessary for flux calculations) and Flux vector
!    - Add 2 Control flags for Fluxes (same as in the MS routine) 

!mick fix 6/6/2023 - replaced "FO_extinction" with "extinction_m" in FO call below

!  4/24/24. Add "OutputDirect-fluxes_AllSB" (SEGMENT BOUNDARY) arrays as Output variables (1 extra line)
!      -- This output was formerly written to file (fort.91) as part of a flux debug
!      -- Use the type structure arrays directly (newly defined in vqs_outputs_def)

      if ( Regime_idx == 1 ) Mu0 = cos(theta(0)) !  Set Mu0; only required for Regime 1
      call FOMS_FO_VectorRT ( &
         do_debug_FO_input, do_debug_fluxes, do_calculate_fluxes,          & ! Inputs (5/30/22 flags added)
         do_sunlight, do_deltam_scaling, Regime_idx, nsegments, nlayers,   & ! Inputs
         nstokes, FO_nmoments, seglayer_map, segupdown_map, layermask_map, & ! Inputs 
         nfinedivs, ntraverse, radii, Raycon, cota, theta, sunpaths,       & ! Inputs (5/30/22 theta added)
         Mu0, Genspher_up, Genspher_dn, Rotations_up, Rotations_dn,        & ! Inputs
         ntraversefine, wfine, csqfine, cotfine, rsegsfine, sunpathsfine,  & ! Inputs
         FO_reflec, extinction_m, FO_omega, FO_truncfac, FO_greekmat,      & ! Inputs Optical
         solarflux, fluxvec,                                               & ! Inputs Solar
         VQS_Out%Main%StokesFO, lostrans, sfi_factors, sfi_opdeps,         & ! Outputs
         VQS_Out%Main%ActFluxDir, VQS_Out%Main%RegFluxDir,                 & ! Flux Outputs. Added, 4/14/23
         VQS_Out%Main%ActFluxDir_AllSB, VQS_Out%Main%RegFluxDir_AllSB )      ! Flux Outputs. Added, 4/24/24
      call cpu_time(e2)
      VQS_ModIn%MTiming%SS_TIME = e2 - e1

!  10/15/23. Introduce progress control

      if ( RT_PROGRESS ) then
         write(*,'(A,I1)')' -- Done  VQS     Regular FO calculation, Regime ',Regime_idx
      endif

!  Small number leveling

      if ( abs(VQS_Out%Main%StokesFO(1)) .lt. smallest ) VQS_Out%Main%StokesFO(1) = zero

!  Define initial value of Stokes vector

      VQS_Out%Main%StokesTot = VQS_Out%Main%StokesFO

!  Do MS radiative transfer
!  ========================

!  Only if flagged
!    Rob Fix 5/30/18. Add npert/eps to argument list (debug only)

      if ( do_MSRT ) then

!  Progress

        if (Verbose) write(*,'(/1x,a)') ' - Doing MS RT'

!  Default surface flags

        do_sleave       = .false.
        do_sl_isotropic = .false.

!  04/05/23. Method Index (set by hand)
!     1 = Fast (Mk3 implementation), 2 = Slow (Mk2 implementation)

!mick mod 12/1/2021 - added "do_debug_MS_input" to inputs

!  4/24/24. Add "OutputDiffuse-fluxes_AllSB" (SEGMENT BOUNDARY) arrays as Output variables (1 extra line)
!      -- This output was formerly written to file (fort.81) as part of a flux debug
!      -- Use the type structure arrays directly (newly defined in vqs_outputs_def)

!  8/26/24. AllSB output line in the calling statement, was in the wrong order (AAARGH !!!!!!!!!)
!    -- NOW        VQS_Out%Main%ActFluxDiff_AllSB, VQS_Out%Main%RegFluxDiff_AllSB,    &
!    -- WAS        VQS_Out%Main%RegFluxDiff_AllSB, VQS_Out%Main%ActFluxDiff_AllSB,    &

        call cpu_time(e1)

        if ( Method_Index.eq. 1 ) then

          call FOMS_MS_VectorRT_Mk3_Exp ( do_debug_MS_input, RT_Progress,       & ! Inputs Monitoring
             do_plane_parallel, do_Rayleigh_only, do_deltam_scaling,            & ! Inputs RT control flags
             do_lambertian_surface, do_sleave, do_sl_isotropic,                 & ! Inputs (surface flags)
             do_calculate_fluxes, do_debug_fluxes, Regime_idx,                  & ! Inputs (flux-flags/Regime)
             nstreams, nlayers, nstokes, FO_nmoments, FO_pradius, FO_heights,   & ! Inputs (integers/heights)
             nsegments, seglayer_map, segupdown_map, alpha, theta, relaz,       & ! Inputs (segments/geometry)
             nfinedivs, alphafine, thetafine, relazfine,                        & ! Inputs (fine-division control)
             lostrans, sfi_factors, sfi_opdeps,                                 & ! Inputs Optical
             FO_reflec, FO_deltaus, FO_omega, FO_greekmat, flux_factor,         & ! Inputs Optical
             VQS_Out%Main%StokesMS,          VQS_Out%Main%StokesTot,            & ! Output Stokes
             VQS_Out%Main%RegFluxDiff,       VQS_Out%Main%ActFluxDiff,          & ! Output Fluxes
             VQS_Out%Main%ActFluxDiff_AllSB, VQS_Out%Main%RegFluxDiff_AllSB,    & ! Output Fluxes_SB, Added, 4/24/24
             fail_MS, status_MS, nmessages_MS, messages_MS, actions_MS )          ! Exception handling

        else if ( Method_Index .eq. 2 ) then

          call FOMS_MS_VectorRT_Mk2_Alt ( do_debug_MS_input, RT_Progress,       & ! Inputs Monitoring
             do_plane_parallel, do_Rayleigh_only, do_deltam_scaling,            & ! Inputs RT control flags
             do_lambertian_surface, do_sleave, do_sl_isotropic,                 & ! Inputs (surface flags)
             do_calculate_fluxes, do_debug_fluxes, Regime_idx,                  & ! Inputs (flux-flags/Regime)
             nstreams, nlayers, nstokes, FO_nmoments, FO_pradius, FO_heights,   & ! Inputs (integers/heights)
             nsegments, seglayer_map, segupdown_map, alpha, theta, Relaz,       & ! Inputs (segments/geometry)
             nfinedivs, alphafine, thetafine, relazfine,                        & ! Inputs (fine-division control)
             lostrans, sfi_factors, sfi_opdeps,                                 & ! Inputs Optical
             FO_reflec, FO_deltaus, FO_omega, FO_greekmat, flux_factor,         & ! Inputs Optical
             VQS_Out%Main%StokesMS,          VQS_Out%Main%StokesTot,            & ! Output Stokes
             VQS_Out%Main%RegFluxDiff,       VQS_Out%Main%ActFluxDiff,          & ! Output Fluxes
             VQS_Out%Main%ActFluxDiff_AllSB, VQS_Out%Main%RegFluxDiff_AllSB,    & ! Output Fluxes_SB, Added, 4/24/24
             fail_MS, status_MS, nmessages_MS, messages_MS, actions_MS )          ! Exception handling

        endif

        call cpu_time(e2)
        VQS_ModIn%MTiming%MS_TIME = e2 - e1 

!  Exception handling. Add Actions

        if ( fail_MS ) then
          VQS_Out%Status%Status_Calculation = VQS_SERIOUS
          do n = 1, nmessages_MS
            VQS_Out%Status%messages(2*n-1+nmessages) = Trim(messages_MS(n))
            VQS_Out%Status%messages(2*n  +nmessages) = Trim(actions_MS(n))
          enddo
          nmessages = nmessages + 2*nmessages_MS
          if ( Method_Index.eq.1 ) VQS_Out%Status%messages(nmessages+1) = '@@@@@ Running with the Mk3 Fast MS Code @@@@@'
          if ( Method_Index.eq.2 ) VQS_Out%Status%messages(nmessages+1) = '@@@@@ Running with the Mk2 Slow MS Code @@@@@'
          nmessages = nmessages + 1
          VQS_Out%Status%NMessages = nmessages
          return
        endif

      else

        fail_MS = .false.

        VQS_Out%Main%StokesMS    = ZERO
        VQS_Out%Main%RegFluxDiff = ZERO
        VQS_Out%Main%ActFluxDiff = ZERO

        VQS_ModIn%MTiming%MS_TIME = 0.0

!  End MSRT calculation

      endif

!  ==========================================================
!  BEGIN COPY LOCAL VARIABLES TO OUTPUTS (pure OUT variables)
!  ==========================================================

!  Diagnostic information

      VQS_Out%Main%Regime_idx = Regime_idx
      VQS_Out%Main%NSegments  = NSegments

!  Exception handling for model calculation.
!mick note 6/6/2023 - this arrangement is a little different than VLIDORT

      VQS_Out%Status%NMESSAGES = NMESSAGES
      VQS_Out%Status%MESSAGES  = MESSAGES

!  ===================================
!  END COPY LOCAL VARIABLES TO OUTPUTS
!  ===================================

!  Finish

      RETURN

      CONTAINS

      SUBROUTINE VQS_DEBUG_INPUT_MASTER()

        CALL VQS_WRITE_STD_INPUT ( VQS_FixIn, VQS_ModIn )

      END SUBROUTINE VQS_DEBUG_INPUT_MASTER

      END SUBROUTINE VQS_MASTER

!  End Module

      END MODULE vqs_masters_m
