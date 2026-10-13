
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p3                         #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
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
! #  This Version :   VLIDORT_2p8p3                             #
! #  Release Date :   31 March 2021                             #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released        August 2014                 #
! #      2.8   F90, released        May    2017                 #
! #      2.8.1 F90, released        August 2019                 # 
! #      2.8.2 F90, limited release May    2020                 # 
! #                                                             #
! #  Features Summary of Recent VLIDORT Versions:               #
! #  -------------------------------------------                #
! #                                                             #
! #      NEW: TOTAL COLUMN JACOBIANS         (2.4)              #
! #      NEW: BPDF Land-surface KERNELS      (2.4R)             #
! #      NEW: Thermal Emission Treatment     (2.4RT)            #
! #      Consolidated BRDF treatment         (2.4RTC)           #
! #      f77/f90 Release                     (2.5)              #
! #      External SS / New I/O Structures    (2.6)              #
! #                                                             #
! #      SURFACE-LEAVING / BRDF-SCALING      (2.7)              #
! #      TAYLOR Series / OMP THREADSAFE      (2.7)              #
! #      New Water-Leaving Treatment         (2.8)              #
! #      LBBF & BRDF-Telescoping, enabled    (2.8)              #
! #      Several Performance Enhancements    (2.8)              #
! #      Water-leaving coupled code          (2.8.1)            #
! #      Planetary problem, media properties (2.8.1)            #
! #      Doublet geometry post-processing    (2.8.2)            #
! #      Reduction zeroing, dynamic memory   (2.8.2)            #
! #                                                             #
! #  Features Summary of This VLIDORT Version                   #
! #  ----------------------------------------                   #
! #                                                             #
! #   2.8.3, released 31 March 2021.                            #
! #     ==> Green's function RT solutions (Nstokes = 1 or 3)    #
! #     ==> Sphericity Corrections using MS source terms        #
! #     ==> BRDF upgrades, including new snow reflectance       #
! #     ==> SLEAVE Upgrades, extended water-leaving treatment   #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.3 of the VLIDORT_2p8 software library.      #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2024.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VLIDORT_2p8p3 ( Version 2.8.3 )            #
! #                                                                 #
! # VLIDORT_2p8p3 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p3 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p3  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #       AC_LPS_SPHERICITY_CORRECTION_MASTER (master, public)  #
! #          MS2pt_AC_LPS_SPHERCORR_V1                          #
! #                                                             #
! ###############################################################

MODULE VLIDORT_AC_LPS_SPHERCORR_ROUTINES_m

!  HERE IS THE HISTORY
!  ===================

!  2-point correction history --
!    V1: 11/20/19. Originally Coded 20-29 November 2019 for VLIDORT Version 2.8.1
!    V2: 12/18/19. Extension to include BOA downwelling situation as well as TOA upwelling.
!    V3: 03/01/20. Renamed, along with addition of 3-point and multi-point Corrections

!  1/31/21. Developed for Official VLIDORT Version 2.8.3 Example
!     R. Spurr. RT Solutions Inc.

!  11/28/22. Upgrade to Version 2.8.3. 
!    -- Gather together sphericity correction code all in one module
!    -- Use single master routine to call required LPS correction (new)

!  1/31/25. Introduce sphericity corrections for airplane viewing. 2-point only.
!    -- 2-point, 3-point and Multi-point formalism still in place (however)
!    -- Aircraft must be at one of the levels in the atmosphere
!    -- Geometry at the aircraft zenith not the same as LOS geometry at BOA or TOA
!         ** necessary to find BOA/TOA geometries.

!  HERE IS THE VLIDORT IMPLEMENTATION
!  ==================================

!  This is the LPS linearization of the sphericity correction

!  Module files for VLIDORT.

      USE VLIDORT_PARS_m
      USE VLIDORT_IO_DEFS_m
      USE VLIDORT_LIN_IO_DEFS_m
      USE VLIDORT_LPS_MASTERS_m
 
!  11/28/22. Upgrade to Version 2.8.3.
!    -- Geometry conversion subroutines are public (needed elsewhere)

!  11/20/23. Version 2.8.3. Sphericity Correction variables are now in type structures
!    -- Remove MSType argument from Master I/O (becomes a local variable)
!    -- 3-point SC now uses weighting factor --> call it Version 2

!  1/31/25. 2-point AC sphericity correction - need 2 geometry subroutines

      USE VLIDORT_AC_SPHERCORR_ROUTINES_m, Only : AC_to_BOA_geometry_conversion, AC_to_BOATOA_geometry_conversion

!  11/28/22. Upgrade to Version 2.8.3.
!    -- NEW LPS_master routine is public

public  ::  AC_LPS_SPHERICITY_CORRECTION_MASTER

!  All MSType routines are private

private :: MS2pt_AC_LPS_SPHERCORR_V1

contains

SUBROUTINE AC_LPS_SPHERICITY_CORRECTION_MASTER  ( do_debug_input, & 
        VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,          & ! VLIDORT Standard   Inputs
        VLIDORT_LinFixIn, VLIDORT_LinModIn, VLIDORT_LinSup, & ! VLIDORT Linearized Inputs
        STANDARD_FO_STOKES     , STANDARD_MS_STOKES     , STANDARD_STOKES     , & ! Standard-VLIDORT     (STOKES)
        CORRECTED_FO_STOKES    , CORRECTED_MS_STOKES    , CORRECTED_STOKES    , & ! Sphericity-corrected (STOKES)
        STANDARD_FO_PROFILEWF  , STANDARD_MS_PROFILEWF  , STANDARD_PROFILEWF  , & ! Standard-VLIDORT     (ProfileWFs)
        CORRECTED_FO_PROFILEWF , CORRECTED_MS_PROFILEWF , CORRECTED_PROFILEWF , & ! Sphericity-corrected (ProfileWFs)
        STANDARD_FO_SURFACEWF  , STANDARD_MS_SURFACEWF  , STANDARD_SURFACEWF  , & ! Standard-VLIDORT     (SurfaceWFs)
        CORRECTED_FO_SURFACEWF , CORRECTED_MS_SURFACEWF , CORRECTED_SURFACEWF , & ! Sphericity-corrected (SurfaceWFs)
        VLIDORT_Out, VLIDORT_LinOut, Twilight_Flag, Fail, Main_Message, Local_Message, Local_Action )  ! Exceptions and errors

!  Implicit none

      IMPLICIT NONE

!  Debug input flag

      LOGICAL, intent(in) :: do_debug_input

!  11/20/23. MSType is now a local variables, copied from type-structure SC input
!    -- Sphericity type ('2' = 2-point, '3' = 3-point, 'M' = Multi-point)
!      CHARACTER*1, intent(in) :: MSTYPE

!  VLIDORT structures (including supplements i/o)

      TYPE(VLIDORT_Fixed_Inputs)   , Intent(in)    :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), Intent(inout) :: VLIDORT_ModIn
      TYPE(VLIDORT_Sup_InOut)      , Intent(inout) :: VLIDORT_Sup

!  VLIDORT Linearized Input structures (including supplements i/o)

      TYPE(VLIDORT_Fixed_LinInputs)   , Intent(in)    :: VLIDORT_LinFixIn
      TYPE(VLIDORT_Modified_LinInputs), Intent(inout) :: VLIDORT_LinModIn
      TYPE(VLIDORT_LinSup_InOut)      , Intent(inout) :: VLIDORT_LinSup

!  VLIDORT Output Structures

      TYPE(VLIDORT_Outputs)        , Intent(inout) :: VLIDORT_Out
      TYPE(VLIDORT_LinOutputs)     , Intent(inout) :: VLIDORT_LinOut

!  Standard output with explanations

      DOUBLE PRECISION, INTENT(OUT) :: STANDARD_FO_STOKES (4)  ! 1. Standard  First-order (SS/DB) STOKES for AC geometry
      DOUBLE PRECISION, INTENT(OUT) :: STANDARD_MS_STOKES (4)  ! 2. Standard  Multiple-scattering STOKES for AC geometry
      DOUBLE PRECISION, INTENT(OUT) :: STANDARD_STOKES    (4)  ! 3. Standard  Total STOKES: 3 = 1 + 2
      DOUBLE PRECISION, INTENT(OUT) :: CORRECTED_FO_STOKES(4)  ! 4. Corrected First-order (SS/DB) STOKES for BOA geometry
      DOUBLE PRECISION, INTENT(OUT) :: CORRECTED_MS_STOKES(4)  ! 5. Corrected Multiple-scattering STOKES for 2-pt interpolated geometries
      DOUBLE PRECISION, INTENT(OUT) :: CORRECTED_STOKES   (4)  ! 6. Corrected Total STOKES: 6 = 4 + 5

!  Column Jacobian output with explanations

      DOUBLE PRECISION, INTENT(OUT) :: STANDARD_FO_PROFILEWF  ( MAX_ATMOSWFS, MAXLAYERS, 4 ) ! 1. Standard  First-order (SS/DB) LP Jac for AC geometry
      DOUBLE PRECISION, INTENT(OUT) :: STANDARD_MS_PROFILEWF  ( MAX_ATMOSWFS, MAXLAYERS, 4 ) ! 2. Standard  Multiple-scattering LP Jac for AC geometry
      DOUBLE PRECISION, INTENT(OUT) :: STANDARD_PROFILEWF     ( MAX_ATMOSWFS, MAXLAYERS, 4 ) ! 3. Standard  Total LP Jac: 3 = 1 + 2
      DOUBLE PRECISION, INTENT(OUT) :: CORRECTED_FO_PROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, 4 ) ! 4. Corrected First-order (SS/DB) LP Jac for BOA geometry
      DOUBLE PRECISION, INTENT(OUT) :: CORRECTED_MS_PROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, 4 ) ! 5. Corrected Multiple-scattering LP Jac for 2-pt interpolated geometries
      DOUBLE PRECISION, INTENT(OUT) :: CORRECTED_PROFILEWF    ( MAX_ATMOSWFS, MAXLAYERS, 4 ) ! 6. Corrected Total LP Jac: 6 = 4 + 5

!  Surface Jacobian output with explanations

      DOUBLE PRECISION, INTENT(OUT) :: STANDARD_FO_SURFACEWF  ( MAX_SURFACEWFS, 4 ) ! 1. Standard  First-order (SS/DB) LS Jac for AC geometry
      DOUBLE PRECISION, INTENT(OUT) :: STANDARD_MS_SURFACEWF  ( MAX_SURFACEWFS, 4 ) ! 2. Standard  Multiple-scattering LS Jac for AC geometry
      DOUBLE PRECISION, INTENT(OUT) :: STANDARD_SURFACEWF     ( MAX_SURFACEWFS, 4 ) ! 3. Standard  Total LS Jac: 3 = 1 + 2
      DOUBLE PRECISION, INTENT(OUT) :: CORRECTED_FO_SURFACEWF ( MAX_SURFACEWFS, 4 ) ! 4. Corrected First-order (SS/DB) LS Jac for BOA geometry
      DOUBLE PRECISION, INTENT(OUT) :: CORRECTED_MS_SURFACEWF ( MAX_SURFACEWFS, 4 ) ! 5. Corrected Multiple-scattering LS Jac for 2-pt interpolated geometries
      DOUBLE PRECISION, INTENT(OUT) :: CORRECTED_SURFACEWF    ( MAX_SURFACEWFS, 4 ) ! 6. Corrected Total LS Jac: 6 = 4 + 5

!  Twilight flag and exception handling

      LOGICAL      , intent(OUT) :: Twilight_Flag, FAIL
      CHARACTER*(*), intent(OUT) :: Main_Message, Local_Message, Local_Action

!  11/20/23. Version 2.8.3. MSType is now a local variable, copied from type-structure SC input
!    -- Sphericity type ('2' = 2-point, '3' = 3-point, 'M' = Multi-point)

      CHARACTER*1      :: MSTYPE

!  11/20/23. Version 2.8.3. 3-point Sphericity Correction mid-geometry weighting factor
!    -- New Variable, first introduced for 2.8.5

      DOUBLE PRECISION :: SPHERCORR_3PTFACTOR

!  Start Code
!  ==========

!  11/20/23. MSType is now a local variable, copied from type-structure SC input
!    -- Sphericity type ('2' = 2-point, '3' = 3-point, 'M' = Multi-point)

      MSType = VLIDORT_FixIn%Cont%TS_SPHERCORR_TYPE

!  11/20/23. Version 2.8.3. 3-point Sphericity Correction mid-geometry weighting factor
!    -- New Variable, first introduced for 2.8.5

      SPHERCORR_3PTFACTOR = ZERO
      if ( MSType .eq. '3' ) SPHERCORR_3PTFACTOR = VLIDORT_FixIn%Cont%TS_SPHERCORR_3PTFACTOR

!  Call to dedicated LPS_sphericity correction routine (either 2-point, 3-point or multi-point)
!   -- 11/20/23. Version 2.8.3. 3-point correction has the 3-pt weighting factor argument (Version 2)
!   -- 1/31/25. 2-point only for the AC scenarios.

      if ( MSType .eq. '2' ) then
        CALL MS2pt_AC_LPS_SPHERCORR_V1 ( do_debug_input, &
            VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup, VLIDORT_LinFixIn, VLIDORT_LinModIn, VLIDORT_LinSup, & ! VLIDORT Inputs
            STANDARD_FO_STOKES     , STANDARD_MS_STOKES     , STANDARD_STOKES     , & ! Standard-VLIDORT      (Stokes)
            CORRECTED_FO_STOKES    , CORRECTED_MS_STOKES    , CORRECTED_STOKES    , & ! Sphericity-corrected  (Stokes)
            STANDARD_FO_PROFILEWF  , STANDARD_MS_PROFILEWF  , STANDARD_PROFILEWF  , & ! Standard-VLIDORT      (ProfileWFs)
            CORRECTED_FO_PROFILEWF , CORRECTED_MS_PROFILEWF , CORRECTED_PROFILEWF , & ! Sphericity-corrected  (ProfileWFs)
            STANDARD_FO_SURFACEWF  , STANDARD_MS_SURFACEWF  , STANDARD_SURFACEWF  , & ! Standard-VLIDORT      (SurfaceWFs)
            CORRECTED_FO_SURFACEWF , CORRECTED_MS_SURFACEWF , CORRECTED_SURFACEWF , & ! Sphericity-corrected  (SurfaceWFs)
            VLIDORT_Out, VLIDORT_LinOut, Twilight_Flag, Fail, Local_Message, Local_Action )                  ! Exceptions and errors
      else if ( MSType .eq. '3' ) then
!  Placeholder
      else if ( MSType .eq. 'M' ) then
!  Placeholder
      endif

!  Exception handling

      IF ( Fail ) THEN
        Main_message = 'Failure from Call to MS'//MStype//'pt_LPS_SPHERICITY_CORRECTION'
      endif

!  return

      return
END SUBROUTINE AC_LPS_SPHERICITY_CORRECTION_MASTER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    P R I V A T E    R O U T I N E S
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine MS2pt_AC_LPS_SPHERCORR_V1 ( do_debug_input, & 
            VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup, VLIDORT_LinFixIn, VLIDORT_LinModIn, VLIDORT_LinSup, & ! VLIDORT Inputs
            STANDARD_FO_STOKES     , STANDARD_MS_STOKES     , STANDARD_STOKES     , & ! Standard-VLIDORT      (Stokes)
            CORRECTED_FO_STOKES    , CORRECTED_MS_STOKES    , CORRECTED_STOKES    , & ! Sphericity-corrected  (Stokes)
            STANDARD_FO_PROFILEWF  , STANDARD_MS_PROFILEWF  , STANDARD_PROFILEWF  , & ! Standard-VLIDORT      (ProfileWFs)
            CORRECTED_FO_PROFILEWF , CORRECTED_MS_PROFILEWF , CORRECTED_PROFILEWF , & ! Sphericity-corrected  (ProfileWFs)
            STANDARD_FO_SURFACEWF  , STANDARD_MS_SURFACEWF  , STANDARD_SURFACEWF  , & ! Standard-VLIDORT      (SurfaceWFs)
            CORRECTED_FO_SURFACEWF , CORRECTED_MS_SURFACEWF , CORRECTED_SURFACEWF , & ! Sphericity-corrected  (SurfaceWFs)
            VLIDORT_Out, VLIDORT_LinOut, Twilight_Flag, Fail, Local_Message, Local_Action ) ! Exceptions and errors

!  Implicit none

      IMPLICIT NONE

!  VLIDORT Input structures (including supplements i/o)

      TYPE(VLIDORT_Fixed_Inputs)   , Intent(in)    :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), Intent(inout) :: VLIDORT_ModIn
      TYPE(VLIDORT_Sup_InOut)      , Intent(inout) :: VLIDORT_Sup

!  VLIDORT Linearized Input structures (including supplements i/o)

      TYPE(VLIDORT_Fixed_LinInputs)   , Intent(in)    :: VLIDORT_LinFixIn
      TYPE(VLIDORT_Modified_LinInputs), Intent(inout) :: VLIDORT_LinModIn
      TYPE(VLIDORT_LinSup_InOut)      , Intent(inout) :: VLIDORT_LinSup

!  VLIDORT Output Structures

      TYPE(VLIDORT_Outputs)        , Intent(inout) :: VLIDORT_Out
      TYPE(VLIDORT_LinOutputs)     , Intent(inout) :: VLIDORT_LinOut

!  debug input flag

      LOGICAL, intent(in) :: do_debug_input

!  Output
!  ------
!mick fix 12/20/2022 - changed intent of main outputs to "inout"

!  Standard output with explanations

      DOUBLE PRECISION, INTENT(OUT) :: STANDARD_FO_STOKES (4)  ! 1. Standard  First-order (SS/DB) STOKES for AC geometry
      DOUBLE PRECISION, INTENT(OUT) :: STANDARD_MS_STOKES (4)  ! 2. Standard  Multiple-scattering STOKES for AC geometry
      DOUBLE PRECISION, INTENT(OUT) :: STANDARD_STOKES    (4)  ! 3. Standard  Total STOKES: 3 = 1 + 2
      DOUBLE PRECISION, INTENT(OUT) :: CORRECTED_FO_STOKES(4)  ! 4. Corrected First-order (SS/DB) STOKES for BOA geometry
      DOUBLE PRECISION, INTENT(OUT) :: CORRECTED_MS_STOKES(4)  ! 5. Corrected Multiple-scattering STOKES for 2-pt interpolated geometries
      DOUBLE PRECISION, INTENT(OUT) :: CORRECTED_STOKES   (4)  ! 6. Corrected Total STOKES: 6 = 4 + 5

!  Column Jacobian output with explanations

      DOUBLE PRECISION, INTENT(OUT) :: STANDARD_FO_PROFILEWF  ( MAX_ATMOSWFS, MAXLAYERS, 4 ) ! 1. Standard  First-order (SS/DB) LP Jac for AC geometry
      DOUBLE PRECISION, INTENT(OUT) :: STANDARD_MS_PROFILEWF  ( MAX_ATMOSWFS, MAXLAYERS, 4 ) ! 2. Standard  Multiple-scattering LP Jac for AC geometry
      DOUBLE PRECISION, INTENT(OUT) :: STANDARD_PROFILEWF     ( MAX_ATMOSWFS, MAXLAYERS, 4 ) ! 3. Standard  Total LP Jac: 3 = 1 + 2
      DOUBLE PRECISION, INTENT(OUT) :: CORRECTED_FO_PROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, 4 ) ! 4. Corrected First-order (SS/DB) LP Jac for BOA geometry
      DOUBLE PRECISION, INTENT(OUT) :: CORRECTED_MS_PROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, 4 ) ! 5. Corrected Multiple-scattering LP Jac for 2-pt interpolated geometries
      DOUBLE PRECISION, INTENT(OUT) :: CORRECTED_PROFILEWF    ( MAX_ATMOSWFS, MAXLAYERS, 4 ) ! 6. Corrected Total LP Jac: 6 = 4 + 5

!  Surface Jacobian output with explanations

      DOUBLE PRECISION, INTENT(OUT) :: STANDARD_FO_SURFACEWF  ( MAX_SURFACEWFS, 4 ) ! 1. Standard  First-order (SS/DB) LS Jac for AC geometry
      DOUBLE PRECISION, INTENT(OUT) :: STANDARD_MS_SURFACEWF  ( MAX_SURFACEWFS, 4 ) ! 2. Standard  Multiple-scattering LS Jac for AC geometry
      DOUBLE PRECISION, INTENT(OUT) :: STANDARD_SURFACEWF     ( MAX_SURFACEWFS, 4 ) ! 3. Standard  Total LS Jac: 3 = 1 + 2
      DOUBLE PRECISION, INTENT(OUT) :: CORRECTED_FO_SURFACEWF ( MAX_SURFACEWFS, 4 ) ! 4. Corrected First-order (SS/DB) LS Jac for BOA geometry
      DOUBLE PRECISION, INTENT(OUT) :: CORRECTED_MS_SURFACEWF ( MAX_SURFACEWFS, 4 ) ! 5. Corrected Multiple-scattering LS Jac for 2-pt interpolated geometries
      DOUBLE PRECISION, INTENT(OUT) :: CORRECTED_SURFACEWF    ( MAX_SURFACEWFS, 4 ) ! 6. Corrected Total LS Jac: 6 = 4 + 5

!  Twilight flag and exception handling

      LOGICAL      , intent(INOUT) :: Twilight_Flag, FAIL
      CHARACTER*(*), intent(INOUT) :: Local_Message, Local_Action

!  Local
!  -----

!  Variables for doing the sphericity calculation

      INTEGER   :: N, G, G1, G2, G3, Q, NC, K
      REAL(fpk) :: MU, MU0, MU1, D01, F0, F1, TRANS, SMSST(4), AMSST(MAXLAYERS,4), COSSCAT
      REAL(fpk) :: CUMSOURCE(0:MAXLAYERS,4), LP_CUMSOURCE(4), LS_CUMSOURCE(MAX_SURFACEWFS,4)
      REAL(fpk) :: LP_TRANS, LP_SMSST(MAX_ATMOSWFS,MAXLAYERS,4), LP_AMSST(MAX_ATMOSWFS,MAXLAYERS,MAXLAYERS,4)
      REAL(fpk) :: LS_SMSST(MAX_SURFACEWFS,4), LS_AMSST(MAX_SURFACEWFS,MAXLAYERS,4)

!  Proxies

      LOGICAL   :: DO_DNWELLING, DO_UPWELLING
      LOGICAL   :: DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION, LAYER_VARY_FLAG ( MAXLAYERS )
      INTEGER   :: NLAYERS, NS, NSTOKES, NSWF, ACLEVEL, NGEOMS_LOCAL, LAYER_VARY_NUMBER ( MAXLAYERS )
      REAL(fpk) :: EARTH_RADIUS, BOA_HEIGHT, AC_HEIGHT, TOA_HEIGHT
      REAL(fpk) :: OBSGEOMS_BOA(3), OBSGEOMS_AC(3), OBSGEOMS_TOA(3)

!  Miscellaneous variables

      LOGICAL :: Verbose

!  Start Code
!  ==========

!  Geometry indices
!    G1 = BOAGEOM ; G2 = ACGEOM ; G3 = TOAGEOM ( downwelling only)

      G1 = 1 ; G2 = 2 ; G3 = 3

!  Initialize

      CORRECTED_FO_STOKES = zero ; CORRECTED_FO_PROFILEWF = zero ; CORRECTED_FO_SURFACEWF = zero
      CORRECTED_MS_STOKES = zero ; CORRECTED_MS_PROFILEWF = zero ; CORRECTED_MS_SURFACEWF = zero
      CORRECTED_STOKES    = zero ; CORRECTED_PROFILEWF    = zero ; CORRECTED_SURFACEWF    = zero
      STANDARD_FO_STOKES  = zero ; STANDARD_FO_PROFILEWF  = zero ; STANDARD_FO_SURFACEWF  = zero
      STANDARD_MS_STOKES  = zero ; STANDARD_MS_PROFILEWF  = zero ; STANDARD_MS_SURFACEWF  = zero
      STANDARD_STOKES     = zero ; STANDARD_PROFILEWF     = zero ; STANDARD_SURFACEWF     = zero

!  Set exceptions

      FAIL          = .false.
      Twilight_Flag = .false.
      Local_Message = ' '
      Local_Action  = ' '

!  Set Proxies (saved values), NLAYERS, Directions.....

      ns           = VLIDORT_FixIn%Cont%TS_nstokes
      nlayers      = VLIDORT_FixIn%Cont%TS_nlayers
      DO_DNWELLING = VLIDORT_FixIn%Bool%TS_DO_DNWELLING
      DO_UPWELLING = VLIDORT_FixIn%Bool%TS_DO_UPWELLING
      EARTH_RADIUS = VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS   

!  Linearization proxies

      DO_PROFILE_LINEARIZATION = VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION
      DO_SURFACE_LINEARIZATION = VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION
      LAYER_VARY_FLAG          = VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG
      LAYER_VARY_NUMBER        = VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER
      NSWF                     = VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS

!  Set debug output flag

      Verbose   = .false.
!      Verbose   = .true.

!  Set the DO_MSSTS flag for the 2-point sphericity correction
!  There are very strict conditions for this specialist option, as follows :==>
!     1. Either DO_UPWELLING or DO_DNWELLING must be set, Not Both !!!!
!     2. DO_FULLRAD_MODE must be set
!     3. DO_OBSERVATION_GEOMETRY must be set, with N_USER_OBSGEOMS = 2
!     4. DO_FOCORR and DO_FOCORR_OUTGOING must both be set
!     5a. Upwelling  : N_USER_LEVELS = 1, and USER_LEVELS(1) = 0.0            [ TOA output only ]
!     5b. Downwelling: N_USER_LEVELS = 1, and USER_LEVELS(1) = Real(nlayers)  [ BOA output only ]
!     5c. Up or Down  : N_USER_LEVELS = 1, and USER_LEVELS(1) = Real(ACLEVEL)  [  AC output (up or down) ]

!  These checks have been implemented inside VLIDORT.
!      VLIDORT_FixIn%Bool%TS_DO_MSSTS = .true.

!  Geometry Conversion
!  ===================

!  AC Level. Note that input geometry is given for the AirCraft location

      ACLEVEL = INT(VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1))

!  inputs to conversion routines (AirCraft geometry and height, height grid, earth radius)
!    -- Original VLIDORT User_obsgeom input is for AC, will get over written

      BOA_HEIGHT       = VLIDORT_FixIn%Chapman%TS_height_grid(NLAYERS)
      TOA_HEIGHT       = VLIDORT_FixIn%Chapman%TS_height_grid(0)
      AC_HEIGHT        = VLIDORT_FixIn%Chapman%TS_height_grid(ACLEVEL)
      OBSGEOMS_AC(1:3) = VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1:3)

!  Conversion. Different for Upwelling and downwelling.

      IF ( DO_UPWELLING ) THEN
        Call AC_to_BOA_geometry_conversion ( &
          EARTH_RADIUS, BOA_HEIGHT, AC_HEIGHT,                       & ! input AC/BOA heights 
          OBSGEOMS_AC(1), OBSGEOMS_AC(2), OBSGEOMS_AC(3),            & ! input AC Geometry
          cosscat, OBSGEOMS_BOA(1), OBSGEOMS_BOA(2), OBSGEOMS_BOA(3) ) ! output BOA geometry, scattering angle
      ELSE IF ( DO_DNWELLING ) THEN
        Call AC_to_BOATOA_geometry_conversion ( &
          EARTH_RADIUS, BOA_HEIGHT, AC_HEIGHT, TOA_HEIGHT,            & ! Input AC/BOA/TOA heights
          OBSGEOMS_AC(1), OBSGEOMS_AC(2), OBSGEOMS_AC(3),             & ! Input AC Geometry
          cosscat, OBSGEOMS_BOA(1), OBSGEOMS_BOA(2), OBSGEOMS_BOA(3), & ! Output BOA geometry, scattering angle
          OBSGEOMS_TOA(1), OBSGEOMS_TOA(2), OBSGEOMS_TOA(3) )           ! Output TOA geometry
      ENDIF

!  Debug

      if ( Verbose ) then
         write(*,*)ACOS(COSSCAT)/DEG_TO_RAD
         write(*,*)OBSGEOMS_BOA(1), OBSGEOMS_BOA(2), OBSGEOMS_BOA(3)
         write(*,*)OBSGEOMS_AC(1) , OBSGEOMS_AC(2) , OBSGEOMS_AC(3)
         if ( do_dnwelling ) write(*,*)OBSGEOMS_TOA(1), OBSGEOMS_TOA(2), OBSGEOMS_TOA(3)
      endif

!  Reset the VLIDORT second geometrical input (TOA)

      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(2,1:3) = OBSGEOMS_TOA(1:3)
      VLIDORT_ModIn%MSunrays%TS_SZANGLES(2)                = OBSGEOMS_TOA(1)
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(2)     = OBSGEOMS_TOA(2)
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(2)            = OBSGEOMS_TOA(3)

!  Reset the VLIDORT Obervational geometry inputs (BOA and AC). Downwelling, add TOA
!  re-set geometry numbers. 2 or 3 Geometries

      VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS = G2
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(G1,1:3) = OBSGEOMS_BOA(1:3)
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(G2,1:3) = OBSGEOMS_AC(1:3)
      IF ( DO_DNWELLING ) THEN
        VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS = G3
        VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(G3,1:3) = OBSGEOMS_TOA(1:3)
      ENDIF
      NGEOMS_LOCAL = VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS

!  reset bookkeeping angles (safety)

      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES      = NGEOMS_LOCAL
      VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES = NGEOMS_LOCAL
      VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS  = NGEOMS_LOCAL
      DO G = 1, NGEOMS_LOCAL
        VLIDORT_ModIn%MSunrays%TS_SZANGLES(G)            = VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(G,1)
        VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(G) = VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(G,2)
        VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(G)        = VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(G,3)
      ENDDO

!  DOWNWELLING ONLY. 
!  Twilight condition at TOA (SZA > 90); set dummy output flag and skip calculation. 

      if ( DO_DNWELLING .and. OBSGEOMS_TOA(G1).gt.90.0 ) then
         Twilight_flag = .true. ; return
      endif

!  Call to VLIDORT

      CALL VLIDORT_LPS_master ( do_debug_input, &
          VLIDORT_FixIn,    & ! INPUTS
          VLIDORT_ModIn,    & ! INPUTS (possibly modified)
          VLIDORT_Sup,      & ! INPUTS/OUTPUTS
          VLIDORT_Out,      & ! OUTPUTS
          VLIDORT_LinFixIn, & ! INPUTS
          VLIDORT_LinModIn, & ! INPUTS (possibly modified)
          VLIDORT_LinSup,   & ! INPUTS/OUTPUTS
          VLIDORT_LinOut )    ! OUTPUTS

!  Exception handling (simply done here)

      if ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK  .eq. VLIDORT_SERIOUS .or. &
           VLIDORT_Out%Status%TS_STATUS_CALCULATION .eq. VLIDORT_SERIOUS ) then
         Local_message = 'Some errors arising from VLIDORT CALL in MS2pt_AC_LPS_SPHERCORR_V1'
         Local_Action  = 'In Driver, use call to VLIDORT_WRITE_STATUS to examine errors'
         FAIL = .true. ; return
      Endif

!  Store Radiance results for CORRECTED FO VLIDORT calculations using BOA Geometry
!    ==> OUTPUT LEVEL IS ALWAYS AT THE AIRCRAFT

      if ( DO_UPWELLING ) then
         CORRECTED_FO_STOKES(1:ns) = VLIDORT_Sup%SS%TS_STOKES_SS(1,G1,1:ns,UPIDX) &
                                   + VLIDORT_Sup%SS%TS_STOKES_DB(1,1:ns,G1)
         IF ( DO_PROFILE_LINEARIZATION ) THEN
            DO K = 1, NLAYERS ; IF ( LAYER_VARY_FLAG(K) ) THEN ; DO Q = 1, LAYER_VARY_NUMBER(K)
               CORRECTED_FO_PROFILEWF(Q,K,1:NS) = VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_SS(Q,K,1,G1,1:ns,UPIDX) &
                                                + VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_DB(Q,K,1,1:ns,G1) 
            ENDDO ; ENDIF ; ENDDO
         ENDIF
         IF ( DO_SURFACE_LINEARIZATION ) THEN
            DO Q = 1, NSWF
               CORRECTED_FO_SURFACEWF(Q,1:NS) = VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(Q,1,1:ns,G1) 
            ENDDO
         ENDIF
      else if ( DO_DNWELLING ) then
         CORRECTED_FO_STOKES(1:ns)  = VLIDORT_Sup%SS%TS_STOKES_SS(1,G1,1:ns,DNIDX)
         IF ( DO_PROFILE_LINEARIZATION ) THEN
            DO K = 1, NLAYERS ; IF ( LAYER_VARY_FLAG(K) ) THEN ; DO Q = 1, LAYER_VARY_NUMBER(K)
               CORRECTED_FO_PROFILEWF(Q,K,1:NS) = VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_SS(Q,K,1,G1,1:ns,DNIDX)
            ENDDO ; ENDIF ; ENDDO
         ENDIF
      endif

!  Store Radiance results for STANDARD FO+MS VLIDORT calculation using AC geometry.
!    ==> OUTPUT LEVEL IS ALWAYS AT THE AIRCRAFT

      if ( DO_UPWELLING ) then
         STANDARD_FO_STOKES(1:ns)  = VLIDORT_Sup%SS%TS_STOKES_SS(1,G2,1:ns,UPIDX) & 
                                      + VLIDORT_Sup%SS%TS_STOKES_DB(1,G2,1:ns) 
         STANDARD_MS_STOKES(1:ns)  = VLIDORT_Out%Main%TS_STOKES(1,G2,1:ns,UPIDX)  - STANDARD_FO_STOKES(1:ns)
         STANDARD_STOKES   (1:ns)  = VLIDORT_Out%Main%TS_STOKES(1,G2,1:ns,UPIDX)
         IF ( DO_PROFILE_LINEARIZATION ) THEN
            DO K = 1, NLAYERS ; IF ( LAYER_VARY_FLAG(K) ) THEN ; DO Q = 1, LAYER_VARY_NUMBER(K)
               STANDARD_FO_PROFILEWF(Q,K,1:ns) = VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_SS(Q,K,1,G2,1:ns,UPIDX) &
                                               + VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_DB(Q,K,1,G2,1:ns) 
               STANDARD_MS_PROFILEWF(Q,K,1:ns) = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,K,1,G2,1:ns,UPIDX) &
                                               - STANDARD_FO_PROFILEWF(Q,K,1:ns)
               STANDARD_PROFILEWF(Q,K,1:ns)    = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,K,1,G2,1:ns,UPIDX) 
            ENDDO ; ENDIF ; ENDDO
         ENDIF
         IF ( DO_SURFACE_LINEARIZATION ) THEN
            DO Q = 1, NSWF
               STANDARD_FO_SURFACEWF(Q,1:ns) = VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(Q,1,G2,1:ns) 
               STANDARD_MS_SURFACEWF(Q,1:ns) = VLIDORT_LinOut%Surf%TS_SURFACEWF(Q,1,G2,1:ns,UPIDX) - STANDARD_FO_SURFACEWF(Q,1:ns)
               STANDARD_SURFACEWF(Q,1:ns)    = VLIDORT_LinOut%Surf%TS_SURFACEWF(Q,1,G2,1:ns,UPIDX) 
            ENDDO
         ENDIF
      else if ( DO_DNWELLING ) then
         STANDARD_FO_STOKES(1:ns)  = VLIDORT_Sup%SS%TS_STOKES_SS(1,G2,1:ns,DNIDX)
         STANDARD_MS_STOKES(1:ns)  = VLIDORT_Out%Main%TS_STOKES(1,G2,1:ns,DNIDX)  - STANDARD_FO_STOKES(1:ns)
         STANDARD_STOKES   (1:ns)  = VLIDORT_Out%Main%TS_STOKES(1,G2,1:ns,DNIDX)
         IF ( DO_PROFILE_LINEARIZATION ) THEN
            DO K = 1, NLAYERS ; IF ( LAYER_VARY_FLAG(K) ) THEN ; DO Q = 1, LAYER_VARY_NUMBER(K)
               STANDARD_FO_PROFILEWF(Q,K,1:ns) = VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_SS(Q,K,1,G2,1:ns,dnIDX)
               STANDARD_MS_PROFILEWF(Q,K,1:ns) = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,K,1,G2,1:ns,DNIDX) &
                                               - STANDARD_FO_PROFILEWF(Q,K,1:ns)
               STANDARD_PROFILEWF(Q,K,1:ns)    = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,K,1,G2,1:ns,DNIDX) 
            ENDDO ; ENDIF ; ENDDO
         ENDIF
         IF ( DO_SURFACE_LINEARIZATION ) THEN
            DO Q = 1, NSWF
               STANDARD_MS_SURFACEWF(Q,1:ns) = VLIDORT_LinOut%Surf%TS_SURFACEWF(Q,1,G2,1:ns,DNIDX)
               STANDARD_SURFACEWF(Q,1:ns)    = VLIDORT_LinOut%Surf%TS_SURFACEWF(Q,1,G2,1:ns,DNIDX) 
            ENDDO
         ENDIF
      endif

!  perform 2-POINT SPHERICITY CORRECTION, Radiative Transfer
!  =========================================================

!  initialize

      AMSST = zero ; LP_AMSST = zero ; LS_AMSST = zero
      SMSST = zero ; LP_SMSST = zero ; LS_SMSST = zero

!  set surface MS source term, using BOA geometry, as LOS originates at the surface (UPWELLING ONLY)

      if ( DO_UPWELLING ) then
         SMSST(1:NS)  = VLIDORT_Out%Main%TS_SURF_MSSTS(1,1:NS)
         IF ( DO_PROFILE_LINEARIZATION ) THEN
            DO K = 1, NLAYERS ; IF ( LAYER_VARY_FLAG(K) ) THEN ; DO Q = 1, LAYER_VARY_NUMBER(K)
               LP_SMSST(Q,K,1:NS)  = VLIDORT_LinOut%Prof%TS_LP_SURF_MSSTS(Q,K,1,1:NS)
            ENDDO ; ENDIF ; ENDDO
         ENDIF
         IF ( DO_SURFACE_LINEARIZATION ) THEN
            DO Q = 1, NSWF
               LS_SMSST(Q,1:NS)  = VLIDORT_LinOut%Surf%TS_LS_SURF_MSSTS(Q,1,1:NS)
            ENDDO
         ENDIF
      endif

!  Alternative: get the Layer MS source terms by end-point interpolation with Cos(VZA)
!     -- Always works, however may not be accurate around azimuth-flip conditions
!     -- Los angle comes from the FO calculation PATHGEOMS, only for BOA geometry G1
!     -- Upwelling  , use geometries G1 and G2 (BOA and AC) for interpolating
!     -- Downwelling, use geometries G2 and G3 (AC and TOA) for interpolating

      if ( DO_UPWELLING  ) then
        Mu0 = COS(OBSGEOMS_BOA(2)*DEG_TO_RAD)
        Mu1 = COS(OBSGEOMS_AC(2) *DEG_TO_RAD)
        D01 = ONE / ( Mu1 - Mu0 )
        do N = ACLEVEL + 1, NLAYERS
          Mu = 0.5_fpk * ( COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N)) + COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N-1)) )
          F0 = ( Mu1 - Mu ) * D01 ; F1 = ONE - F0
          AMSST(N,1:ns) = F0 * VLIDORT_Out%Main%TS_LAYER_MSSTS(G1,1:ns,N) + F1 * VLIDORT_Out%Main%TS_LAYER_MSSTS(G2,1:ns,N)
          IF ( DO_PROFILE_LINEARIZATION ) THEN
            DO K = 1, NLAYERS ; IF ( LAYER_VARY_FLAG(K) ) THEN ; DO Q = 1, LAYER_VARY_NUMBER(K)
              LP_AMSST(Q,K,N,1:ns) = F0 * VLIDORT_LinOut%Prof%TS_LP_LAYER_MSSTS(Q,K,G1,1:ns,N) &
                                   + F1 * VLIDORT_LinOut%Prof%TS_LP_LAYER_MSSTS(Q,K,G2,1:ns,N)  
            ENDDO ; ENDIF ; ENDDO
          ENDIF
          IF ( DO_SURFACE_LINEARIZATION ) THEN
            DO Q = 1, NSWF
              LS_AMSST(Q,N,1:ns) = F0 * VLIDORT_LinOut%Surf%TS_LS_LAYER_MSSTS(Q,G1,1:ns,N) &
                                 + F1 * VLIDORT_LinOut%Surf%TS_LS_LAYER_MSSTS(Q,G2,1:ns,N)  
            ENDDO
          ENDIF
        ENDDO
      ENDIF

      if ( DO_DNWELLING  ) then
        Mu0 = COS(OBSGEOMS_AC(2) *DEG_TO_RAD)
        Mu1 = COS(OBSGEOMS_TOA(2)*DEG_TO_RAD)
        D01 = ONE / ( Mu1 - Mu0 )
        do N = 1, ACLEVEL
          Mu = 0.5_fpk * ( COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N)) + COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N-1)) )
          F0 = ( Mu1 - Mu ) * D01 ; F1 = ONE - F0
          AMSST(N,1:ns) = F0 * VLIDORT_Out%Main%TS_LAYER_MSSTS(G2,1:ns,N) + F1 * VLIDORT_Out%Main%TS_LAYER_MSSTS(G3,1:ns,N)
          IF ( DO_PROFILE_LINEARIZATION ) THEN
            DO K = 1, NLAYERS ; IF ( LAYER_VARY_FLAG(K) ) THEN ; DO Q = 1, LAYER_VARY_NUMBER(K)
              LP_AMSST(Q,K,N,1:ns) = F0 * VLIDORT_LinOut%Prof%TS_LP_LAYER_MSSTS(Q,K,G2,1:ns,N) &
                                   + F1 * VLIDORT_LinOut%Prof%TS_LP_LAYER_MSSTS(Q,K,G3,1:ns,N)  
            ENDDO ; ENDIF ; ENDDO
          ENDIF
          IF ( DO_SURFACE_LINEARIZATION ) THEN
            DO Q = 1, NSWF
              LS_AMSST(Q,N,1:ns) = F0 * VLIDORT_LinOut%Surf%TS_LS_LAYER_MSSTS(Q,G2,1:ns,N) &
                                 + F1 * VLIDORT_LinOut%Surf%TS_LS_LAYER_MSSTS(Q,G3,1:ns,N)  
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  Radiative Transfer recursion for spherically-corrected MS fields
!  ================================================================

!  Perform the radiative transfer recursion for MS_SPHERCORR fields, the spherically-corrected MS field
!   ==> Start with the surface surface terms SMSST, LP_SMSST, LS_SMSST, Upwelling scenario

      if ( DO_UPWELLING ) then
         CUMSOURCE(0,1:NS) = SMSST(1:NS) ; NC = 0
         DO N = NLAYERS, ACLEVEL+1, -1
            NC = NC + 1 ; TRANS = VLIDORT_Out%Main%TS_LOSTRANS(G1,N)
            CUMSOURCE(NC,1:NS) = TRANS * CUMSOURCE(NC-1,1:NS) + AMSST(N,1:ns)
         ENDDO
         CORRECTED_MS_STOKES(1:NS) = CUMSOURCE(NC,1:NS) 
         IF ( DO_PROFILE_LINEARIZATION ) THEN
           DO K = 1, NLAYERS ; IF ( LAYER_VARY_FLAG(K) ) THEN ; DO Q = 1, LAYER_VARY_NUMBER(K)
             LP_CUMSOURCE(1:NS) = LP_SMSST(Q,K,1:NS) ; NC = 0
             DO N = NLAYERS, ACLEVEL+1, -1
               NC = NC + 1 ; TRANS  = VLIDORT_Out%Main%TS_LOSTRANS(G1,N)
               LP_TRANS = zero ; if ( K.eq.N) LP_TRANS = VLIDORT_LinOut%Prof%TS_LP_LOSTRANS(Q,G1,N)
               LP_CUMSOURCE(1:NS) = TRANS * LP_CUMSOURCE(1:NS) + LP_TRANS * CUMSOURCE(NC-1,1:NS) + LP_AMSST(Q,K,N,1:NS)
             ENDDO
             CORRECTED_MS_PROFILEWF(Q,K,1:NS) = LP_CUMSOURCE(1:NS)
           ENDDO ; ENDIF ; ENDDO
         ENDIF
         IF ( DO_SURFACE_LINEARIZATION ) THEN
           DO Q = 1, NSWF
             LS_CUMSOURCE(Q,1:NS) = LS_SMSST(Q,1:NS)
             DO N = NLAYERS, ACLEVEL+1, -1
               NC = NC + 1 ; TRANS  = VLIDORT_Out%Main%TS_LOSTRANS(G1,N)
               LS_CUMSOURCE(Q,1:NS) = TRANS * LS_CUMSOURCE(Q,1:NS) + LS_AMSST(Q,N,1:NS)
             ENDDO
             CORRECTED_MS_SURFACEWF(Q,1:NS) = LS_CUMSOURCE(Q,1:NS)
           ENDDO
         ENDIF
      endif

!  Downwelling case

      if ( DO_DNWELLING ) THEN
         CUMSOURCE(0,1:NS) = ZERO
         DO N = 1, ACLEVEL
            TRANS = VLIDORT_Out%Main%TS_LOSTRANS(1,N)
            CUMSOURCE(N,1:NS) = TRANS * CUMSOURCE(N-1,1:NS) + AMSST(N,1:NS) 
         ENDDO
!  The recursion ends at the aircraft, not at NLAYERS (which is uncomputed).
         CORRECTED_MS_STOKES(1:NS) = CUMSOURCE(ACLEVEL,1:NS)
         IF ( DO_PROFILE_LINEARIZATION ) THEN
           DO K = 1, NLAYERS ; IF ( LAYER_VARY_FLAG(K) ) THEN ; DO Q = 1, LAYER_VARY_NUMBER(K)
             LP_CUMSOURCE(1:NS) = ZERO
             DO N = 1, ACLEVEL
               TRANS    = VLIDORT_Out%Main%TS_LOSTRANS(G1,N)
               LP_TRANS = zero ; if ( K.eq.N) LP_TRANS = VLIDORT_LinOut%Prof%TS_LP_LOSTRANS(Q,G1,N)
               LP_CUMSOURCE(1:NS) = TRANS * LP_CUMSOURCE(1:NS) + LP_TRANS * CUMSOURCE(N-1,1:NS) + LP_AMSST(Q,K,N,1:NS)
             ENDDO
             CORRECTED_MS_PROFILEWF(Q,K,1:NS) = LP_CUMSOURCE(1:NS)
           ENDDO ; ENDIF ; ENDDO
         ENDIF
         IF ( DO_SURFACE_LINEARIZATION ) THEN
           DO Q = 1, NSWF
             LS_CUMSOURCE(Q,1:NS) = ZERO
             DO N = 1, ACLEVEL
               TRANS                = VLIDORT_Out%Main%TS_LOSTRANS(G1,N)
               LS_CUMSOURCE(Q,1:NS) = TRANS * LS_CUMSOURCE(Q,1:NS) + LS_AMSST(Q,N,1:NS)
             ENDDO
             CORRECTED_MS_SURFACEWF(Q,1:NS) = LS_CUMSOURCE(Q,1:NS)
           ENDDO
         ENDIF
      endif

!  Final answers: Add FO Outgoing results to MS fields

      CORRECTED_STOKES(1:NS) = CORRECTED_MS_STOKES(1:NS) + CORRECTED_FO_STOKES(1:NS)
      IF ( DO_PROFILE_LINEARIZATION ) THEN
        DO K = 1, NLAYERS ; IF ( LAYER_VARY_FLAG(K) ) THEN ; DO Q = 1, LAYER_VARY_NUMBER(K)
          CORRECTED_PROFILEWF(Q,K,1:NS) = CORRECTED_MS_PROFILEWF(Q,K,1:NS) + CORRECTED_FO_PROFILEWF(Q,K,1:NS)
        ENDDO ; ENDIF ; ENDDO
      ENDIF
      IF ( DO_SURFACE_LINEARIZATION ) THEN
        DO Q = 1, NSWF
          CORRECTED_SURFACEWF(Q,1:NS) = CORRECTED_MS_SURFACEWF(Q,1:NS) + CORRECTED_FO_SURFACEWF(Q,1:NS)
        ENDDO
      ENDIF

!  3/7/25. Extremely important !!!!!!!!!!!!!!!!!!!!!!!!! 
!    -- Restore the original geometry, otherwise next run will be off!

      VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS = 1
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(:,1:3) = zero
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1:3) = OBSGEOMS_AC(1:3) 

!  Verbose debug

      if ( verbose .and. DO_UPWELLING ) then
         write(242,*) OBSGEOMS_BOA(1:3)
         write(242,*) OBSGEOMS_AC (1:3)
         write(242,*) 'Corrected AC INTENSITY: FO/MS/Total : ',CORRECTED_FO_STOKES(1),CORRECTED_MS_STOKES(1),CORRECTED_STOKES(1)
         write(242,*) 'Standard  AC INTENSITY: FO/MS/Total : ',STANDARD_FO_STOKES(1), STANDARD_MS_STOKES(1), STANDARD_STOKES(1)
      endif

      if ( verbose .and. DO_DNWELLING ) then
         write(242,*) OBSGEOMS_BOA(1:3)
         write(242,*) OBSGEOMS_AC (1:3)
         write(242,*) OBSGEOMS_TOA(1:3)
         write(242,*) 'Corrected AC INTENSITY: FO/MS/Total : ',CORRECTED_FO_STOKES(1),CORRECTED_MS_STOKES(1),CORRECTED_STOKES(1)
         write(242,*) 'Standard  AC INTENSITY: FO/MS/Total : ',STANDARD_FO_STOKES(1), STANDARD_MS_STOKES(1), STANDARD_STOKES(1)
      endif

!  normal Finish

      return
end subroutine MS2pt_AC_LPS_SPHERCORR_V1

!  End module

END MODULE VLIDORT_AC_LPS_SPHERCORR_ROUTINES_m

