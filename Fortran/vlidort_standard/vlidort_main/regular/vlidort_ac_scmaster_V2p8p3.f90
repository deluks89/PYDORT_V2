
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
! #         AC_SPHERICITY_CORRECTION_MASTER (master, public)    #
! #                                                             #
! #         MS2pt_AC_SPHERICITY_CORRECTION_V1                   #
! #                                                             #
! #            AC_to_BOA_geometry_conversion,                   #
! #            AC_to_BOATOA_geometry_conversion                 #
! #                                                             #
! ###############################################################

MODULE VLIDORT_AC_SPHERCORR_ROUTINES_m

!  HERE IS THE HISTORY
!  ===================

!  2-point correction history --
!    V1: 11/20/19. Originally Coded 20-29 November 2019 for VLIDORT Version 2.8.1
!    V2: 12/18/19. Extension to include BOA downwelling situation as well as TOA upwelling.
!    V3: 03/01/20. Renamed, along with addition of 3-point and multi-point Corrections

!  3-point correction history --
!    V1: 03/01/20. New 3pt Correction

!  Multi-point correction history --
!    V1: 03/01/20. New  Correction

!  1/31/21. Developed for Official VLIDORT Version 2.8.3 Package
!     R. Spurr. RT Solutions Inc.

!  11/28/22. Upgrade to Version 2.8.3.
!    -- Gather together sphericity correction code all in one module
!    -- Use single master routine to call required correction (new)

!  2/15/25. Introduce sphericity corrections for airplane viewing
!    -- 2-point, 3-point and Multi-point formalism still in place.
!    -- Aircraft must be at one of the levels in the atmosphere
!    -- Geometry at the aircraft zenith not the same as LOS geometry at BOA or TOA
!         ** necessary to find BOA/TOA geometries.

!  HERE IS THE VLIDORT IMPLEMENTATION
!  ==================================

!  Module files for VLIDORT.

      USE VLIDORT_PARS_m
      USE VLIDORT_IO_DEFS_m
      USE VLIDORT_MASTERS_m
 
!  9/8/22. CODECTOOL Upgrade to Version 2.8.5.
!    -- NEW master routine is public
!    -- Geometry conversion subroutines must be public (needed elsewhere)

!  11/20/23. Version 2.8.3. Sphericity Correction variables are now in type structures
!    -- Remove MSType argument from Master I/O (becomes a local variable)
!    -- 3-point SC now uses weighting factor --> call it Version 2

public  ::  AC_SPHERICITY_CORRECTION_MASTER, &
            AC_to_BOA_geometry_conversion, AC_to_BOATOA_geometry_conversion

!  All three MSType routines are private

private  :: MS2pt_AC_SPHERICITY_CORRECTION_V1

contains

SUBROUTINE AC_SPHERICITY_CORRECTION_MASTER (  &
             do_debug_input, VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,  & ! VLIDORT Inputs, debug flag
             STANDARD_FO_STOKES,  STANDARD_MS_STOKES,  STANDARD_STOKES,  & ! Standard-VLIDORT     outputs
             CORRECTED_FO_STOKES, CORRECTED_MS_STOKES, CORRECTED_STOKES, & ! Sphericity-corrected outputs
             VLIDORT_Out, Twilight_Flag, Fail, Main_Message, Local_Message, Local_Action )  ! VLIDORT TS output, Exceptions

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
      TYPE(VLIDORT_Outputs)        , Intent(inout) :: VLIDORT_Out

!  Main output with explanations

      REAL(fpk), INTENT(OUT) :: STANDARD_FO_STOKES(4)   ! 1. Standard  First-order (SS/DB) STOKES Vector for AC geometry
      REAL(fpk), INTENT(OUT) :: STANDARD_MS_STOKES(4)   ! 2. Standard  Multiple-scattering STOKES Vector for AC geometry
      REAL(fpk), INTENT(OUT) :: STANDARD_STOKES(4)      ! 3. Standard  Total STOKES Vector: 3 = 1 + 2
      REAL(fpk), INTENT(OUT) :: CORRECTED_FO_STOKES(4)  ! 4. Corrected First-order (SS/DB) STOKES Vector for BOA geometry
      REAL(fpk), INTENT(OUT) :: CORRECTED_MS_STOKES(4)  ! 5. Corrected Multiple-scattering STOKES Vector for 2-pt interpolated geometries
      REAL(fpk), INTENT(OUT) :: CORRECTED_STOKES(4)     ! 6. Corrected Total STOKES Vector: 6 = 4 + 5

!  Twilight flag and exception handling

      LOGICAL      , intent(out) :: Twilight_Flag, FAIL
      CHARACTER*(*), intent(out) :: Main_message, Local_Message, Local_Action

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

!  Call to dedicated sphericity correction routine (either 2-point, 3-point or multi-point)
!   -- 11/20/23. Version 2.8.3. 3-point correction has the 3-pt weighting factor argument (Version 2)

      if ( MSType .eq. '2' ) then
        CALL MS2pt_AC_SPHERICITY_CORRECTION_V1 ( &
             do_debug_input, VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,  & ! VLIDORT type structure Inputs
             STANDARD_FO_STOKES,  STANDARD_MS_STOKES,  STANDARD_STOKES,  & ! Standard-VLIDORT      outputs
             CORRECTED_FO_STOKES, CORRECTED_MS_STOKES, CORRECTED_STOKES, & ! Sphericity-corrected outputs
             VLIDORT_Out, Twilight_Flag, Fail, Local_Message, Local_Action ) ! Exceptions and errors
!      else if ( MSType .eq. '3' ) then
!        CALL MS3pt_AC_SPHERICITY_CORRECTION_V1 ( &
!             do_debug_input, VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,  & ! VLIDORT type structure Inputs
!             STANDARD_FO_STOKES,  STANDARD_MS_STOKES,  STANDARD_STOKES,  & ! Standard-VLIDORT      outputs
!             CORRECTED_FO_STOKES, CORRECTED_MS_STOKES, CORRECTED_STOKES, & ! Sphericity-corrected outputs
!             VLIDORT_Out, Twilight_Flag, Fail, Local_Message, Local_Action ) ! Exceptions and errors
!      else if ( MSType .eq. 'M' ) then
!        CALL MSMpt_AC_SPHERICITY_CORRECTION_V1 ( &
!             do_debug_input, VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,  & ! VLIDORT type structure Inputs
!             STANDARD_FO_STOKES,  STANDARD_MS_STOKES,  STANDARD_STOKES,  & ! Standard-VLIDORT      outputs
!             CORRECTED_FO_STOKES, CORRECTED_MS_STOKES, CORRECTED_STOKES, & ! Sphericity-corrected outputs
!             VLIDORT_Out, Twilight_Flag, Fail, Local_Message, Local_Action ) ! Exceptions and errors
      endif

!  Exception handling

      IF ( Fail ) THEN
        Main_message = 'Failure from Call to MS'//MStype//'pt_AC_SPHERICITY_CORRECTION'
      endif

!  return

      return
END SUBROUTINE AC_SPHERICITY_CORRECTION_MASTER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    P R I V A T E    R O U T I N E S
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine MS2pt_AC_SPHERICITY_CORRECTION_V1 ( &
           do_debug_input, VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,    & ! VLIDORT Inputs
           STANDARD_FO_STOKES,  STANDARD_MS_STOKES,  STANDARD_STOKES,    & ! Standard-VLIDORT     outputs
           CORRECTED_FO_STOKES, CORRECTED_MS_STOKES, CORRECTED_STOKES,   & ! Sphericity-corrected outputs
           VLIDORT_OUT, Twilight_Flag, Fail, Local_Message, Local_Action ) ! Exceptions and errors

!  Implicit none

      IMPLICIT NONE

!  Debug input flag

      LOGICAL, intent(in) :: do_debug_input

!  VLIDORT structures (including supplements i/o)

      TYPE(VLIDORT_Fixed_Inputs)   , Intent(in)    :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), Intent(inout) :: VLIDORT_ModIn
      TYPE(VLIDORT_Sup_InOut)      , Intent(inout) :: VLIDORT_Sup
      TYPE(VLIDORT_Outputs)        , Intent(inout) :: VLIDORT_Out

!  Output
!  ------

!  Main output with explanations (pre-initialized)

      REAL(fpk), INTENT(OUT) :: STANDARD_FO_STOKES(4)   ! 1. Standard  First-order (SS/DB) Stokes vector for AC geometry
      REAL(fpk), INTENT(OUT) :: STANDARD_MS_STOKES(4)   ! 2. Standard  Multiple-scattering Stokes vector for AC geometry
      REAL(fpk), INTENT(OUT) :: STANDARD_STOKES(4)      ! 3. Standard  Total Stokes vector: 3 = 1 + 2
      REAL(fpk), INTENT(OUT) :: CORRECTED_FO_STOKES(4)  ! 4. Corrected First-order (SS/DB) Stokes vector for BOA geometry
      REAL(fpk), INTENT(OUT) :: CORRECTED_MS_STOKES(4)  ! 5. Corrected Multiple-scattering Stokes vector for 2-pt interpolated geometries
      REAL(fpk), INTENT(OUT) :: CORRECTED_STOKES(4)     ! 6. Corrected Total Stokes vector: 6 = 4 + 5

!  Twilight flag and exception handling (pre-initialized)

      LOGICAL      , intent(INOUT) :: Twilight_Flag, FAIL
      CHARACTER*(*), intent(INOUT) :: Local_Message, Local_Action

!  Local
!  -----

      LOGICAL   :: DO_DNWELLING, DO_UPWELLING
      INTEGER   :: NLAYERS, ACLEVEL, NGEOMS_LOCAL
      REAL(fpk) :: EARTH_RADIUS, BOA_HEIGHT, AC_HEIGHT, TOA_HEIGHT
      REAL(fpk) :: OBSGEOMS_BOA(3), OBSGEOMS_AC(3), OBSGEOMS_TOA(3)

!  Help Variables for doing the sphericity calculation

      INTEGER   :: N, G, G1, G2, G3, NS, O1
      REAL(fpk) :: MU, MU0, MU1, D01, F0, F1, TRANS, SMSST(4), AMSST(MAXLAYERS,4), COSSCAT

!  Miscellaneous variables

      LOGICAL :: Verbose

!  Start Code
!  ==========

!  Geometry indices
!    G1 = BOAGEOM ; G2 = ACGEOM ; G3 = TOAGEOM ( downwelling only)

      G1 = 1 ; G2 = 2 ; G3 = 3

!  Initialize

      CORRECTED_FO_STOKES = zero
      CORRECTED_MS_STOKES = zero
      CORRECTED_STOKES    = zero
      STANDARD_FO_STOKES  = zero
      STANDARD_MS_STOKES  = zero
      STANDARD_STOKES     = zero

!  Set exceptions

      FAIL          = .false.
      Twilight_Flag = .false.
      Local_Message = ' '
      Local_Action  = ' '

!  Set Proxies (saved values)

      ns      = VLIDORT_FixIn%Cont%TS_nstokes
      nlayers = VLIDORT_FixIn%Cont%TS_nlayers
      DO_DNWELLING = VLIDORT_FixIn%Bool%TS_DO_DNWELLING
      DO_UPWELLING = VLIDORT_FixIn%Bool%TS_DO_UPWELLING
      EARTH_RADIUS = VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS   

!  Set debug output flag

      Verbose   = .false.
!      Verbose   = .true.

!  Set the DO_MSSTS flag for the 2-point sphericity correction
!  There are very strict conditions for this specialist option, as follows :==>
!     1. Either DO_UPWELLING or DO_DNWELLING must be set, Not Both !!!!
!     2. DO_FULLRAD_MODE must be set
!     3. DO_OBSERVATION_GEOMETRY must be set, with N_USER_OBSGEOMS = 1
!     4. DO_FOCORR and DO_FOCORR_OUTGOING must both be set
!     5a. Upwelling   : N_USER_LEVELS = 1, and USER_LEVELS(1) = 0.0            [ TOA upwelling output only ]
!     5b. Downwelling : N_USER_LEVELS = 1, and USER_LEVELS(1) = Real(nlayers)  [ BOA dnwelling output only ]
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
         write(*,*)ACOS(COSSCAT)/DEG_TO_RAD, DO_UPWELLING, DO_DNWELLING
         write(*,*)OBSGEOMS_BOA(1), OBSGEOMS_BOA(2), OBSGEOMS_BOA(3)
         write(*,*)OBSGEOMS_AC(1) , OBSGEOMS_AC(2) , OBSGEOMS_AC(3)
         if ( do_dnwelling ) write(*,*)OBSGEOMS_TOA(1), OBSGEOMS_TOA(2), OBSGEOMS_TOA(3)
      endif

!  Reset the VLIDORT second geometrical input (TOA)

      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(2,1:3) = OBSGEOMS_TOA(1:3)
      VLIDORT_ModIn%MSunrays%TS_SZANGLES(2)                = OBSGEOMS_TOA(1)
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(2)     = OBSGEOMS_TOA(2)
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(2)            = OBSGEOMS_TOA(3)

!  Reset the LIDORT Obervational geometry inputs (BOA and AC). Downwelling, add TOA
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

      CALL VLIDORT_MASTER ( do_debug_input, &
          VLIDORT_FixIn, &
          VLIDORT_ModIn, &
          VLIDORT_Sup,   &
          VLIDORT_Out )

!  Exception handling (simply done here)

      if ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK  .eq. VLIDORT_SERIOUS .or. &
           VLIDORT_Out%Status%TS_STATUS_CALCULATION .eq. VLIDORT_SERIOUS ) then
         Local_message = 'Some errors arising from VLIDORT CALL in MS2pt_AC_SPHERICITY_CORRECTION_V1'
         Local_Action  = 'In Driver, use call to VLIDORT_WRITE_STATUS to examine errors'
         FAIL = .true. ; return
      Endif

!  Store Radiance results for CORRECTED FO LIDORT calculation using BOA Geometry
!    ==> OUTPUT LEVEL IS ALWAYS AT THE AIRCRAFT

      if ( DO_UPWELLING ) then
         CORRECTED_FO_STOKES(1:NS) = VLIDORT_Sup%SS%TS_STOKES_SS(1,G1,1:NS,UPIDX) + VLIDORT_Sup%SS%TS_STOKES_DB(1,G1,1:NS) 
      else if ( DO_DNWELLING ) then
         CORRECTED_FO_STOKES(1:NS)  = VLIDORT_Sup%SS%TS_STOKES_SS(1,G1,1:NS,DNIDX)
      endif

!  Store Radiance results for STANDARD FO+MS LIDORT calculation using AC geometry.
!    ==> OUTPUT LEVEL IS ALWAYS AT THE AIRCRAFT

      if ( DO_UPWELLING ) then
         STANDARD_FO_STOKES(1:NS)  = VLIDORT_Sup%SS%TS_STOKES_SS(1,G2,1:NS,UPIDX) + VLIDORT_Sup%SS%TS_STOKES_DB(1,G2,1:NS) 
         STANDARD_MS_STOKES(1:NS)  = VLIDORT_Out%Main%TS_STOKES(1,G2,1:NS,UPIDX)  - STANDARD_FO_STOKES(1:NS)
         STANDARD_STOKES(1:NS)     = VLIDORT_Out%Main%TS_STOKES(1,G2,1:NS,UPIDX)
      else if ( DO_DNWELLING ) then
         STANDARD_FO_STOKES(1:NS)  = VLIDORT_Sup%SS%TS_STOKES_SS(1,G2,1:NS,DNIDX)
         STANDARD_MS_STOKES(1:NS)  = VLIDORT_Out%Main%TS_STOKES(1,G2,1:NS,DNIDX)  - STANDARD_FO_STOKES(1:NS)
         STANDARD_STOKES(1:NS)     = VLIDORT_Out%Main%TS_STOKES(1,G2,1:NS,DNIDX)
      endif

!  perform 2-POINT SPHERICITY CORRECTION, Radiative Transfer
!  =========================================================

!  initialize

      AMSST = zero ; SMSST = zero

!  set surface MS source term, using BOA geometry, as LOS originates at the surface (UPWELLING ONLY)

      if ( DO_UPWELLING ) then
         SMSST(1:NS)  = VLIDORT_Out%Main%TS_SURF_MSSTS(G1,1:NS)
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
         DO O1 = 1, NS
           AMSST(N,O1) = F0 * VLIDORT_Out%Main%TS_LAYER_MSSTS(G1,O1,N) + F1 * VLIDORT_Out%Main%TS_LAYER_MSSTS(G2,O1,N)
         ENDDO
        ENDDO
      ENDIF

      if ( DO_DNWELLING  ) then
        Mu0 = COS(OBSGEOMS_AC(2) *DEG_TO_RAD)
        Mu1 = COS(OBSGEOMS_TOA(2)*DEG_TO_RAD)
        D01 = ONE / ( Mu1 - Mu0 )
        do N = 1, ACLEVEL
          Mu = 0.5_fpk * ( COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N)) + COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N-1)) )
          F0 = ( Mu1 - Mu ) * D01 ; F1 = ONE - F0
         DO O1 = 1, NS
           AMSST(N,O1) = F0 * VLIDORT_Out%Main%TS_LAYER_MSSTS(G2,O1,N) + F1 * VLIDORT_Out%Main%TS_LAYER_MSSTS(G3,O1,N)
         ENDDO
        ENDDO
      ENDIF

!  Radiative Transfer recursion
!  ============================

!  Perform the radiative transfer recursion for MS_STOKES_SPHERCORR, the spherically-corrected MS field
!   ==> Start with the surface term SMSST, Upwelling scenario

      if ( DO_UPWELLING ) then
         CORRECTED_MS_STOKES(1:NS) = SMSST(1:NS)
         DO N = NLAYERS, ACLEVEL+1, -1
            TRANS = VLIDORT_Out%Main%TS_LOSTRANS(G1,N)
            CORRECTED_MS_STOKES(1:NS) = TRANS * CORRECTED_MS_STOKES(1:NS) + AMSST(N,1:NS)
         ENDDO
      else if ( DO_DNWELLING ) THEN
         CORRECTED_MS_STOKES(1:NS) = ZERO
         DO N = 1, ACLEVEL
            TRANS = VLIDORT_Out%Main%TS_LOSTRANS(G1,N)
            CORRECTED_MS_STOKES(1:NS) = TRANS * CORRECTED_MS_STOKES(1:NS) + AMSST(N,1:NS) 
         ENDDO
      endif

!  Final Radiance answer: Add FO Outgoing result to MS field.

      CORRECTED_STOKES(1:NS) = CORRECTED_MS_STOKES(1:NS) + CORRECTED_FO_STOKES(1:NS)

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
end subroutine MS2pt_AC_SPHERICITY_CORRECTION_V1

!

subroutine AC_to_BOA_geometry_conversion ( &
     Rearth, Hboa, Haircraft, theta_ac, alpha_ac, phi_ac, & ! input  AC/BOA heights, AC Geometry
     cosscat, theta_boa, alpha_boa, phi_boa )               ! output BOA geometry, scattering angle

!  Completely stand-alone 2-point conversion routine.
!   12/18/19. Add direction index for downwelling case

   implicit none

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs/Outputs

   real(fpk), intent(inout)  :: Rearth, Hboa, Haircraft, theta_ac, alpha_ac, phi_ac
   real(fpk), intent(out)    :: theta_boa, alpha_boa, phi_boa, cosscat

!  Local

   real(fpk) :: alpha_ac_R, theta_ac_R, phi_ac_R, dtr, pie, zero, one, small
   real(fpk) :: SolarDirection(3), vsign, vsign_solar, b, px(3)
   real(fpk) :: salpha_ac, calpha_ac, stheta_ac, ctheta_ac, cphi_ac, sphi_ac
   real(fpk) :: salpha_boa, calpha_boa, stheta_boa, ctheta_boa, cphi_boa
   real(fpk) :: term1, term2, rboa, raircraft, cumangle
   real(fpk) :: alpha_boa_R, theta_boa_R, phi_boa_R

!  constants

   zero = 0.0_fpk ; one = 1.0_fpk ; small = 1.0e-12_fpk
   Pie = acos(-1.0_fpk)
   dtr = Pie / 180.0_fpk

   if ( abs(phi_ac)         .lt. small ) phi_boa = zero
   if ( abs(phi_ac-180.0d0) .lt. small ) phi_boa = 180.0d0

!  AC angles in radians (_R) and cos/sin

   alpha_ac_R = alpha_ac * dtr
   salpha_ac  = sin(alpha_ac_R)
   calpha_ac  = cos(alpha_ac_R)

   theta_ac_R    = theta_ac * dtr
   stheta_ac     = sin(theta_ac_R)
   ctheta_ac     = cos(theta_ac_R)
      
   phi_ac_R   = phi_ac * dtr
   cphi_ac    = cos(phi_ac_R)
   sphi_ac    = sin(phi_ac_R)

!  Solar direction (unit vector). UPWELLING, so VSign = -1.0

   vsign       = - one
   vsign_solar = vsign ! vsign_solar = 1.0_fpk
   SolarDirection(1) = - stheta_ac * cphi_ac * vsign_solar
   SolarDirection(2) = - stheta_ac * sphi_ac
   SolarDirection(3) = - ctheta_ac

!  scattering angle

   term1 = salpha_ac * stheta_ac * cphi_ac
   term2 = calpha_ac * ctheta_ac
   cosscat = vsign * term2 + term1 
   if ( cosscat.gt.one)  cosscat = one
   if ( cosscat.lt.-one) cosscat = -one
   if ( abs(cosscat - one) .lt.small ) cosscat = one
   if ( abs(cosscat + one) .lt.small ) cosscat = - one

!  BOA calculation
!  ===============

!  Get Vza at BOA (sine rule)

   rboa = rearth + hboa; raircraft = rearth + haircraft
   salpha_boa = raircraft * salpha_ac / rboa
   calpha_boa = sqrt(one-salpha_boa*salpha_boa)
   alpha_boa_R = asin(salpha_boa)
   cumangle = alpha_boa_R - alpha_ac_R
   alpha_boa = alpha_boa_R / dtr

!  Calculate SZA at BOA, using vector geometry

   px(1) = - rboa * sin(CumAngle)
   px(2) = 0.0_fpk
   px(3) =   rboa * cos(CumAngle)
   b = DOT_PRODUCT(px,SolarDirection)
   ctheta_boa  = - b/rboa
   stheta_boa  = sqrt(one-ctheta_boa*ctheta_boa)
   theta_boa_R = acos(ctheta_boa)
   theta_boa = theta_boa_R / dtr

!  Calculate PHI at BOA using constancy of scattering angle

   if ( phi_ac.eq.zero ) then
     phi_boa = zero
   else if (phi_ac.eq.180.0_fpk) then
     phi_boa = 180.0_fpk
   else
     cphi_boa = (cosscat-vsign*calpha_boa*ctheta_boa)/stheta_boa/salpha_boa
     if ( cphi_boa.gt.one)  cphi_boa = one
     if ( cphi_boa.lt.-one) cphi_boa = -one
     phi_boa_R = acos(cphi_boa)
     if ( phi_ac.gt.180.0_fpk) phi_boa_R = 2.0_fpk * Pie - phi_boa_R
     phi_boa   = phi_boa_R   / dtr
   endif

!write(*,*)acos(cosscat)/dtr, alpha_boa, theta_boa, phi_boa
!write(*,*)acos(cosscat)/dtr, alpha_ac, theta_ac, phi_ac

!  finish

   return
end subroutine AC_to_BOA_geometry_conversion

!

subroutine AC_to_BOATOA_geometry_conversion ( &
     Rearth, Hboa, Haircraft, HToa, theta_ac, alpha_ac, phi_ac,            & ! input  AC/BOA/TOA heights, AC Geometry
     cosscat, theta_boa, alpha_boa, phi_boa, theta_toa, alpha_toa, phi_toa ) ! output BOA/TOA geometry, scattering angle

!  Completely stand-alone geometry conversion routine. Start with AC, derive BOA and TOA
!     -- Usage for the 2-point Downwelling

   implicit none

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs/Outputs

   real(fpk), intent(in)  :: Rearth, Hboa, Haircraft, HToa, theta_ac, alpha_ac, phi_ac
   real(fpk), intent(out) :: theta_boa, alpha_boa, phi_boa, theta_toa, alpha_toa, phi_toa, cosscat

!  Local

   real(fpk) :: alpha_ac_R, theta_ac_R, phi_ac_R, dtr, pie, zero, one, small
   real(fpk) :: SolarDirection(3), vsign, vsign_solar, b, px(3)
   real(fpk) :: salpha_ac, calpha_ac, stheta_ac, ctheta_ac, cphi_ac, sphi_ac
   real(fpk) :: salpha_boa, calpha_boa, stheta_boa, ctheta_boa, cphi_boa
   real(fpk) :: salpha_toa, calpha_toa, stheta_toa, ctheta_toa, cphi_toa
   real(fpk) :: term1, term2, rboa, raircraft, rtoa, cumangle
   real(fpk) :: alpha_boa_R, theta_boa_R, phi_boa_R
   real(fpk) :: alpha_toa_R, theta_toa_R, phi_toa_R

!  constants

   zero = 0.0_fpk ; one = 1.0_fpk ; small = 1.0e-12_fpk
   Pie = acos(-1.0_fpk)
   dtr = Pie / 180.0_fpk

!  AC angles in radians (_R) and cos/sin

   alpha_ac_R = alpha_ac * dtr
   salpha_ac  = sin(alpha_ac_R)
   calpha_ac  = cos(alpha_ac_R)

   theta_ac_R    = theta_ac * dtr
   stheta_ac     = sin(theta_ac_R)
   ctheta_ac     = cos(theta_ac_R)
      
   phi_ac_R   = phi_ac * dtr
   cphi_ac    = cos(phi_ac_R)
   sphi_ac    = sin(phi_ac_R)

!  Solar direction (unit vector). downwelling, so VSign = 1

   vsign       = One
   vsign_solar = - vsign ! vsign_solar = 1.0_fpk
   SolarDirection(1) = - stheta_ac * cphi_ac * vsign_solar
   SolarDirection(2) = - stheta_ac * sphi_ac
   SolarDirection(3) = - ctheta_ac

!  scattering angle

   term1 = salpha_ac * stheta_ac * cphi_ac
   term2 = calpha_ac * ctheta_ac
   cosscat = -vsign * term2 + term1 
   if ( cosscat.gt.one)  cosscat = one
   if ( cosscat.lt.-one) cosscat = -one
   if ( abs(cosscat - one) .lt.small ) cosscat = one
   if ( abs(cosscat + one) .lt.small ) cosscat = - one

!  BOA calculation
!  ===============

!  Get Vza at BOA (sine rule)

   rboa = rearth + hboa; raircraft = rearth + haircraft
   salpha_boa = raircraft * salpha_ac / rboa
   calpha_boa = sqrt(one-salpha_boa*salpha_boa)
   alpha_boa_R = asin(salpha_boa)
   alpha_boa = alpha_boa_R / dtr
   cumangle = alpha_boa_R - alpha_ac_R

!  Calculate SZA at BOA, using vector geometry

   px(1) = - rboa * sin(CumAngle)
   px(2) = 0.0_fpk
   px(3) =  rboa * cos(CumAngle)
   b = DOT_PRODUCT(px,SolarDirection)
   ctheta_boa  = - b/rboa
   stheta_boa  = sqrt(one-ctheta_boa*ctheta_boa)
   theta_boa_R = acos(ctheta_boa)
   theta_boa = theta_boa_R / dtr

!  Calculate PHI at BOA using constancy of scattering angle
!  Calculate PHI at BOA using constancy of scattering angle

   if ( phi_ac.eq.zero ) then
     phi_boa = zero
   else if (phi_ac.eq.180.0_fpk) then
     phi_boa = 180.0_fpk
   else
     cphi_boa = (cosscat+vsign*calpha_boa*ctheta_boa)/stheta_boa/salpha_boa
     if ( cphi_boa.gt.one)  cphi_boa = one
     if ( cphi_boa.lt.-one) cphi_boa = -one
     phi_boa_R = acos(cphi_boa)
     if ( phi_ac.gt.180.0_fpk) phi_boa_R = 2.0_fpk * Pie - phi_boa_R
     phi_boa   = phi_boa_R   / dtr
   endif


!  TOA calculation
!  ===============

!  Get Vza at TOA (sine rule)

   rtoa = rearth + htoa
   salpha_toa = raircraft * salpha_ac / rtoa
   calpha_toa = sqrt(one-salpha_toa*salpha_toa)
   alpha_toa_R = asin(salpha_toa)
   alpha_toa = alpha_toa_R / dtr
   cumangle = alpha_ac_R - alpha_toa_R

!  Calculate SZA at TOA, using vector geometry

   px(1) = rtoa * sin(CumAngle)
   px(2) = 0.0_fpk
   px(3) = rtoa * cos(CumAngle)
   b = DOT_PRODUCT(px,SolarDirection)
   ctheta_toa  = - b/rtoa

   stheta_toa  = sqrt(one-ctheta_toa*ctheta_toa)
   theta_toa_R = acos(ctheta_toa)
   theta_toa = theta_toa_R / dtr

!  Calculate PHI at TOA using constancy of scattering angle

   if ( phi_ac.eq.zero ) then
     phi_toa = zero
   else if (phi_ac.eq.180.0_fpk) then
     phi_toa = 180.0_fpk
   else
     cphi_toa = (cosscat+vsign*calpha_toa*ctheta_toa)/stheta_toa/salpha_toa
     if ( cphi_toa.gt.one)  cphi_toa = one
     if ( cphi_toa.lt.-one) cphi_toa = -one
     phi_toa_R = acos(cphi_toa)
     if ( phi_ac.gt.180.0_fpk) phi_toa_R = 2.0_fpk * Pie - phi_toa_R
     phi_toa   = phi_toa_R   / dtr
   endif

!write(*,*)acos(cosscat)/dtr, alpha_boa, theta_boa, phi_boa
!write(*,*)acos(cosscat)/dtr, alpha_ac, theta_ac, phi_ac
!write(*,*)acos(cosscat)/dtr, alpha_toa, theta_toa, phi_toa

!  finish

   return
end subroutine AC_to_BOATOA_geometry_conversion

!  end module

END MODULE VLIDORT_AC_SPHERCORR_ROUTINES_m
