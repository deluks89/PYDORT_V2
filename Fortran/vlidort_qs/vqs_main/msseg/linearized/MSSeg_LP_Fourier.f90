
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

!  11 April 2023. Revision

      MODULE MSSeg_LP_Fourier_m

      public  MSSEG_LP_FOURIER

      contains

      SUBROUTINE MSSEG_LP_FOURIER ( &
           DO_SURFACE_MSST, DO_UPWELLING, DO_PLANE_PARALLEL, DO_DIRECT_BEAM,        &
           DO_MULTIGEOM, DO_LOCAL_FLUXES, DO_FOURIER0_NSTOKES2, DO_INCLUDE_SURFACE, &
           DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,          &
           DO_PROFILE_LINEARIZATION, LVARY_FLAG, LVARY_NUMBER,                  & ! Linearization control
           FOURIER, NSTOKES, NSTREAMS, NLAYERS, NGEOMS, NCHOICE,                &
           NMOMENTS, NSTREAMS_2, NTOTAL, MUELLER_INDEX, DMAT, DFLUX,            &
           N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2,                  &
           DO_LAYER_SCATTERING, FLUX_FACTOR, FLUXVEC, COS_SZANGLES,             &
           USER_STREAMS, QUAD_STREAMS, QUAD_STRMWTS, QUAD_HALFWTS,              &
           DELTAU_SLANT, DELTAU_VERT, SFI_OPDEPS_SCALED, OMEGA_GREEK,           & ! Input Optical        
           ALBEDO, BRDF_F, BRDF_F_0, USER_BRDF_F, SLTERM_ISOTROPIC, SLTERM_F_0, & ! Input surface 
           BEAM_CUTOFF, TRANS_SOLAR_BEAM, INITIAL_TRANS, AVERAGE_SECANT,        & ! Input Beam param
           T_DELT_MUBAR, LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,  & ! Input Beam param
           L_DELTAU_VERT, L_SFI_OPDEPS_SCALED, L_OMEGA_GREEK,                   & ! Input Lin. Optical
           HGSF, HSOL, HBVP, L_HSOL,                                            & ! Input TS
           MSST_F, SURFACE_MSST_F, LOCAL_FLUXES, PGSF, PSOL, PBVP,              & ! regular outputs and TS
           LP_MSST_F, LP_SURFACE_MSST_F, LOCAL_LP_FLUXES, L_PSOL, L_PBVP,       & ! Linearized outputs and TS
           STATUS, MESSAGE, TRACE_1, TRACE_2 )                                    ! Exception handling

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- ASYMTX_TOLERANCE now introduced and passed into the QHOM call to ASYMTX
!   -- Flux output (LOCAL_FLUXES) controlled by DO_LOCAL_FLUXES input flag
!   -- DO_FOURIER0_NSTOKES2 performance flag introduced.
!   -- Call subroutine specially written to get the discrete-ordinate field for Fluxes

!  Mk2a. 9/21/21. New routine
!  Mk3. Experimental. 8/27/21 and 1/31/22.
!    -- Use internal Type structures. This is the Particular Fourier.
!    -- April 2023. Introduce DO_MULTIGEOM

      USE MSSeg_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES_SQ, MAXSTREAMS, MAXFOURIER, &
                               MAXBEAMS, MAX_SZANGLES, MAX_USER_STREAMS, MAX_GEOMETRIES, MAXEVALUES, &
                               MAX_ATMOSWFS, ZERO, ONE, TWO, PI2, VQS_SUCCESS, VQS_SERIOUS

      use MSSeg_HPTypes_Def_m
      use MSSeg_L_HPTypes_Def_m

      use MSSeg_Partic_Solutions_m
      use MSSeg_Partic_BVProblem_m , Only : MSSEG_BVPSOLUTION_MASTER
      use MSSeg_Partic_PostProcess_m

!  Mk2a. 9/21/21. Use LAPC modules and selected subroutines
!  Mk3.  1/31/22. USe all subroutines as indicated

      use MSSeg_LP_solutions_m
      use MSSeg_LP_BVProblem_m, Only : MSSEG_LP_BVPSOLUTION_MASTER
      use MSSeg_LP_PostProcess_m

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::          DO_SURFACE_MSST
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DIRECT_BEAM
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::          DO_SURFACE_LEAVING
      LOGICAL, INTENT (IN) ::          DO_SL_ISOTROPIC

!  Mk2, 6/21/21. Add the following three inputs
!   -- DO_LOCAL_FLUXES input flag for controlling flux calculations
!   -- DO_FOURIER0_NSTOKES2 performance flag introduced.

      LOGICAL, INTENT (IN) ::          DO_LOCAL_FLUXES 
      LOGICAL, INTENT (IN) ::          DO_FOURIER0_NSTOKES2

!  04/07/23. Mk3. Add the MultiGoem variable
 
      LOGICAL, INTENT (IN) ::          DO_MULTIGEOM ( MAX_GEOMETRIES, 0:MAXFOURIER )

!  Linearization control
!   -- Mk2a. 9/21/21. LC input now

      LOGICAL, INTENT (IN) ::          DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT (IN) ::          LVARY_FLAG   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LVARY_NUMBER ( MAXLAYERS )

!  Integers
!  --------

!  Control integers (basic)

      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NGEOMS
      INTEGER, INTENT (IN) ::          NCHOICE

!  Control integers (derived)

      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2

!  Derived book-keeping quantities
!  -------------------------------

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: DMAT  ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: DFLUX ( MAXSTOKES )

      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING  ( 0:MAXMOMENTS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: FLUX_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: FLUXVEC ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: USER_STREAMS  ( MAX_USER_STREAMS )

!  Quadrature
!  ----------

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_HALFWTS ( MAXSTREAMS )

!  Optical inputs
!  --------------

!  surface reflection
!   -- Mk3. Reduce dimensioning, drop 0:MAXMOMENTS
!   -- 9/21/23. Add BRDF_F argument

      DOUBLE PRECISION, INTENT (IN) :: ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F      ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F_0    ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS   )
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )

!  surface leaving
!   -- Mk3. Reduce dimensioning, drop 0:MAXMOMENTS

      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_ISOTROPIC ( MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_F_0       ( MAXSTOKES, MAXSTREAMS, MAXBEAMS )

!  From the Deltam-scaling setup
!   -- 9/21/23. Added DELTAU_SLANT input for linearized BVP

      DOUBLE PRECISION, INTENT (INOUT) :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: SFI_OPDEPS_SCALED  ( MAX_GEOMETRIES )

!  From the SSALB_INIT setup

      DOUBLE PRECISION, INTENT (INOUT) :: OMEGA_GREEK  ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  From the QSPREP (average secant) setup

      INTEGER, INTENT (INOUT) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: TRANS_SOLAR_BEAM ( MAXBEAMS )

      DOUBLE PRECISION, INTENT (INOUT) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  From the Linearized Deltam-scaling setup

      DOUBLE PRECISION :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION :: L_SFI_OPDEPS_SCALED  ( MAX_GEOMETRIES, MAX_ATMOSWFS )

!  From the Linearized SSALB_INIT setup

      DOUBLE PRECISION :: L_OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  From the Linearized QSPREP (average secant) setup
! -- Mk2a. 9/21/21. LC input now. 9/21/23. Add LP_T_DELT_MUBAR for the Linearized BVP Call

      DOUBLE PRECISION :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Type structures input
!  ---------------------

      type(Homog_GSF)        , intent(in) :: HGSF
      type(Homog_Solutions)  , intent(in) :: HSOL
      type(Homog_BVPSetup)   , intent(in) :: HBVP
      type(L_Homog_Solutions), intent(in) :: L_HSOL

!  Pure output
!  -----------

!  Mk2,  6/21/21. Add the Local Fluxes output (Intent(InOut))
!  Mk2a. 9/21/21. LC output now

      DOUBLE PRECISION, INTENT (OUT)   ::  MSST_F         ( MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT)   ::  SURFACE_MSST_F ( MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) ::  LOCAL_FLUXES   ( MAX_GEOMETRIES, MAXSTOKES, 2, 2 )

      DOUBLE PRECISION, INTENT (OUT)   ::  LP_MSST_F         ( MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT)   ::  LP_SURFACE_MSST_F ( MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) ::  LOCAL_LP_FLUXES   ( MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES, MAXSTOKES, 2, 2 )

!  Type structure output
!  ---------------------

      type(Partic_GSF)      , intent(inout) :: PGSF
      type(Partic_Solutions), intent(inout) :: PSOL
      type(Partic_BVPSolve) , intent(inout) :: PBVP

      type(L_Partic_Solutions), intent(inout) :: L_PSOL
      type(L_Partic_BVPSolve) , intent(inout) :: L_PBVP

!  Bookkeeping and exception handling
!  ----------------------------------

      LOGICAL, INTENT (OUT) ::           DO_INCLUDE_SURFACE
      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (OUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (OUT) :: TRACE_1
      CHARACTER (LEN=*), INTENT (OUT) :: TRACE_2

!  Local variables
!  ---------------

!  downwelling (opposite of upwelling)

      LOGICAL          :: DO_DNWELLING

!  Local transmittances

      DOUBLE PRECISION :: T_UTUP_EIGEN ( MAXEVALUES )
      DOUBLE PRECISION :: T_UTDN_EIGEN ( MAXEVALUES )
      DOUBLE PRECISION :: T_UTDN_MUBAR, TAU_DN, TAU_UP

!  Linearized Local transmittances

      DOUBLE PRECISION :: L_T_UTUP_EIGEN  ( MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UTDN_EIGEN  ( MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_T_UTDN_MUBAR  ( MAX_ATMOSWFS )

!  Local surface variables

      LOGICAL          :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )
      DOUBLE PRECISION :: DIRECT_BEAM  ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION :: ATMOS_ATTN   ( MAXBEAMS )
      DOUBLE PRECISION :: STOKES_DOWNSURF ( MAXSTREAMS, MAXSTOKES )

!  Mk2, 6/21/21. Local variables for flux computations

      DOUBLE PRECISION :: SUM_ActFlux, SUM_RegFlux
      DOUBLE PRECISION :: QSTOKES_F    ( MAXSTREAMS, MAXSTOKES, 2 )
      DOUBLE PRECISION :: QJACOBIANS_F ( MAX_ATMOSWFS, MAXSTREAMS, MAXSTOKES, 2 )

!  Mk2,  6/21/21. Based on VLIDORT Version 2.8.3. Local NSTOKES and related control integers
!  Mk2a. 9/21/21. renamed LC --> LL

      INTEGER :: LL_NSTK, LL_NSTKNSTRM, LL_NSTKNSTRM_2, LL_NTT, LL_NSUBD, LL_NSUPD

!  Flux multiplier and Fourier component numbers

      DOUBLE PRECISION ::  DELTA_FACTOR
      DOUBLE PRECISION ::  SURFACE_FACTOR

!  indices

      INTEGER          :: LAYER, IGEOM, O1, K, Q, KV, QV, I
      LOGICAL          :: DO_RTSOL_VARY
      INTEGER          :: NPARAMS_VARY, LAYER_TO_VARY

!  error tracing

      INTEGER ::           STATUS_SUB
      CHARACTER (LEN=2) :: CF
      CHARACTER (LEN=3) :: CKV

!  ##############
!  initialization
!  ##############

!  module status and message initialization

      STATUS  = VQS_SUCCESS
      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '

!  downwelling

      DO_DNWELLING = .not. DO_UPWELLING

!  Zero output
!  Mk2.  6/21/21. Zero the LOCAL_FLUXES for Fourier = 0
!  Mk2a. 9/21/21. LP output

      MSST_F         = zero ; LP_MSST_F         = zero
      SURFACE_MSST_F = zero ; LP_SURFACE_MSST_F = zero
      if ( FOURIER .EQ. 0 ) THEN
         LOCAL_FLUXES   = zero ; LOCAL_LP_FLUXES = zero
      endif

!  surface reflectance factors

      IF ( FOURIER .EQ. 0 ) THEN
        SURFACE_FACTOR = TWO
        DELTA_FACTOR   = ONE
      ELSE
        SURFACE_FACTOR = ONE
        DELTA_FACTOR   = TWO
      ENDIF

!  Direct beam flag (only if above albedo flag has been set)
!   (DO_REFLECTED_DIRECTBEAM is a stored Bookkeeping variable)

      IF ( DO_DIRECT_BEAM  ) THEN
        IF ( DO_INCLUDE_SURFACE ) THEN
          DO IGEOM = 1, NGEOMS
            DO_REFLECTED_DIRECTBEAM(IGEOM) = .TRUE.
          ENDDO
        ELSE
          DO IGEOM = 1, NGEOMS
            DO_REFLECTED_DIRECTBEAM(IGEOM) = .FALSE.
          ENDDO
        ENDIF
      ENDIF

!  Surface direct beam

      CALL MSSEG_DIRECTBEAM ( &
        DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, &
        FOURIER, NSTOKES, NSTREAMS, NGEOMS, DELTA_FACTOR,           &
        FLUX_FACTOR, FLUXVEC, MUELLER_INDEX, COS_SZANGLES,          &
        ALBEDO, BRDF_F_0, SLTERM_ISOTROPIC, SLTERM_F_0,             &
        TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM, &
        ATMOS_ATTN, DIRECT_BEAM )

!  ###################
!  Spherical functions
!  ###################

!  Get Pi Matrices for this Fourier component
!    Including all beam angles !!

      CALL MSSEG_PARTIC_PIMATRIX_SETUP ( &
        FOURIER, NSTOKES, NLAYERS, NGEOMS, NMOMENTS, &
        COS_SZANGLES, USER_STREAMS, DMAT, HGSF, PGSF )

!  ###################
!  LOCAL NSTOKES USAGE
!  ###################

!  Mk2, 6/21/21. Based on VLIDORT Version 2.8.3. 
!    -- Introduce control for Fast execution of Fourier-zero with LOCAL_NSTOKES = 2

      IF ( DO_Fourier0_NSTOKES2 .and. NSTOKES.eq.3 .and. FOURIER.eq.0 .and. DO_LAMBERTIAN_SURFACE ) THEN
        LL_NSTK = 2 ; LL_NSTKNSTRM = NSTREAMS * LL_NSTK ; LL_NSTKNSTRM_2 = LL_NSTKNSTRM * 2
        LL_NTT  = LL_NSTKNSTRM_2 * NLAYERS
        IF ( NLAYERS .EQ. 1 ) THEN
          LL_NSUBD = 2*LL_NSTKNSTRM - 1 ; LL_NSUPD = LL_NSUBD
        ELSE
          LL_NSUBD = 3*LL_NSTKNSTRM - 1 ; LL_NSUPD = LL_NSUBD
        ENDIF
      ELSE
        LL_NSTK = NSTOKES ; LL_NSTKNSTRM   = NSTKS_NSTRMS ; LL_NSTKNSTRM_2 = NSTKS_NSTRMS_2
        LL_NTT  = NTOTAL  ; LL_NSUBD       = N_SUBDIAG    ; LL_NSUPD       = N_SUPDIAG
      ENDIF

!  ########################################
!  RT differential equation Eigensolutions
!  ########################################

!  Get Post-processing ("user") solutions for Chosen layer
!    MK2, 6/21/21.  Use local variable LL_NSTK

      CALL  MSSEG_UHOM_SOLUTION ( &
        NCHOICE, FOURIER, LL_NSTK, NSTREAMS, NMOMENTS, NGEOMS, &
        DO_UPWELLING, DO_DNWELLING, DO_LAYER_SCATTERING, &
        QUAD_HALFWTS, OMEGA_GREEK, HGSF, PGSF, HSOL, PSOL )

!  Get linearized post-processing solutions for Chosen layer
!  Mk2,  6/21/21.  Use local variable LL_NSTK
!  Mk2a. 9/21/21. Column Linearization

      IF ( DO_PROFILE_LINEARIZATION ) THEN
         DO_RTSOL_VARY = LVARY_FLAG(NCHOICE)
         NPARAMS_VARY  = LVARY_NUMBER(NCHOICE)
         CALL MSSEG_LAP_UHOM_SOLUTION ( &
            NCHOICE, FOURIER, DO_RTSOL_VARY, NPARAMS_VARY,   & ! Input control
            LL_NSTK, NSTREAMS, NMOMENTS, NGEOMS,             & ! Input control
            DO_UPWELLING, DO_DNWELLING, QUAD_HALFWTS,        & ! Input bookkeeping
            DO_LAYER_SCATTERING, OMEGA_GREEK, L_OMEGA_GREEK, & ! Input optical
            HGSF, PGSF, HSOL, PSOL, L_HSOL, L_PSOL )           ! TS Input and Output
      ENDIF

!  Start loop over desired  solutions

      DO IGEOM = 1, NGEOMS
        IF ( DO_MULTIGEOM(IGEOM,FOURIER) ) THEN

!  Solar beam Particular solutions (classical method)
!  --------------------------------------------------

          DO LAYER = 1, NLAYERS

!  Parameter control for the linearization

            IF ( DO_PROFILE_LINEARIZATION ) THEN
              DO_RTSOL_VARY = LVARY_FLAG(LAYER)
              NPARAMS_VARY  = LVARY_NUMBER(LAYER)
            ENDIF

!  get the solution
!    MK2, 6/21/21.  Use local bookkeeping variables LL_NSTK etc.....

            CALL MSSEG_QBEAM_SOLUTION ( &
              LAYER, FOURIER, IGEOM, LL_NSTK, NSTREAMS, NMOMENTS,            &
              NSTREAMS_2, LL_NSTKNSTRM, DO_LAYER_SCATTERING, DFLUX,          &
              FLUX_FACTOR, QUAD_STREAMS, OMEGA_GREEK, BEAM_CUTOFF,           &
              T_DELT_MUBAR, INITIAL_TRANS, AVERAGE_SECANT, PGSF, HGSF, HSOL, &
              PSOL, STATUS_SUB, MESSAGE, TRACE_1 )

!  error flag

            IF ( STATUS_SUB .NE. VQS_SUCCESS ) THEN
              write(CF,'(I2)')FOURIER
              TRACE_2 = 'Error return from MSSEG_QBEAM_SOLUTION, '// &
                        'Called in MSSEG_LP_FOURIER, Fourier # '//CF
              STATUS = VQS_SERIOUS ; RETURN
            ENDIF

!  Linearized solution and error handling
!  Mk2,  6/21/21. Use local bookkeeping variables LL_NSTK etc.....
!  Mk2a. 9/21/21. Column linearization

            IF ( DO_PROFILE_LINEARIZATION ) THEN

              CALL MSSEG_LP_QBEAM_SOLUTION ( &
                LAYER, FOURIER, IGEOM, DO_RTSOL_VARY, NPARAMS_VARY, DO_PLANE_PARALLEL, & ! Input control
                LVARY_FLAG, LVARY_NUMBER, LL_NSTK, NSTREAMS, FLUX_FACTOR,              & ! Input control
                NMOMENTS, NSTREAMS_2, LL_NSTKNSTRM, DFLUX,          & ! Input Bookkeeping
                DO_LAYER_SCATTERING, QUAD_STREAMS, L_OMEGA_GREEK,   & ! Input optical
                BEAM_CUTOFF, AVERAGE_SECANT, LP_AVERAGE_SECANT,     & ! Input average secant
                HGSF, PGSF, HSOL, PSOL, L_HSOL,                     & ! Input TS
                L_PSOL, STATUS_SUB, MESSAGE, TRACE_1 )                ! Output

              IF ( STATUS_SUB .NE. VQS_SUCCESS ) THEN
                write(CF,'(I2)')FOURIER
                TRACE_2 = 'Error from MSSEG_LP_QBEAM_SOLUTION, '// &
                          'Called in MSSEG_LP_FOURIER, Fourier # '//CF
                STATUS = VQS_SERIOUS ; RETURN
              ENDIF

            ENDIF

!  checks out, 09 October 2023
!if ( fourier.eq.0)write(5001,*)nchoice,igeom,layer,L_PSOL%LP_BVEC(1,1,layer,20,1:2),LP_AVERAGE_SECANT(layer,20,igeom,1:2)

!  end layer loop

          END DO

!  checks out, 09 October 2023
!stop '5001,igeom=1,nchoice=1'

!  user solutions for Layer of choice
!    MK2, 6/21/21.  Use local bookkeeping variables LL_NSTK etc.....

!mick fix 7/13/2015 - re-added HELPSTOKES_BEAM to output

          CALL MSSEG_UBEAM_SOLUTION ( &
            NCHOICE, FOURIER, IGEOM,                         &
            DO_UPWELLING, DO_DNWELLING, DO_LAYER_SCATTERING, &
            LL_NSTK, NSTREAMS, NMOMENTS, FLUX_FACTOR,        &
            BEAM_CUTOFF, QUAD_HALFWTS, OMEGA_GREEK,          &
            HGSF, PGSF, PSOL )

!  Linearized User Solutions for Layer of choice
!  Mk2,  6/21/21.  Use local bookkeeping variables LL_NSTK etc.....
!  Mk2a. 9/21/21. Column linearization

          IF ( DO_PROFILE_LINEARIZATION ) THEN
            DO_RTSOL_VARY = LVARY_FLAG(NCHOICE)
            NPARAMS_VARY  = LVARY_NUMBER(NCHOICE)
            CALL MSSEG_LP_UBEAM_SOLUTION ( &
              NCHOICE, FOURIER, IGEOM, DO_RTSOL_VARY, NPARAMS_VARY, LL_NSTK, NSTREAMS, & ! Input control
              DO_PLANE_PARALLEL, LVARY_FLAG, LVARY_NUMBER, DO_UPWELLING, DO_DNWELLING, & ! Input control
              FLUX_FACTOR, DO_LAYER_SCATTERING, NMOMENTS, BEAM_CUTOFF, QUAD_HALFWTS,   & ! Input bookkeeping
              OMEGA_GREEK, L_OMEGA_GREEK, HGSF, PGSF, PSOL, L_PSOL )                     ! Input Optical and TS I/O
          ENDIF

!  #############################################
!         RADIANCE SOLUTIONS
!  #############################################

!  Solve boundary value problem
!  ============================

!  standard case using compressed-band matrices, etc..
!    MK2, 6/21/21.  Use local bookkeeping variables LC_NSTK etc.....
!    Mk2. 8/3/21. Add NCHOICE input (for debug purposes)

          CALL MSSEG_BVPSOLUTION_MASTER ( &
            DO_INCLUDE_SURFACE, DO_REFLECTED_DIRECTBEAM(IGEOM), DO_LAMBERTIAN_SURFACE, &
            IGEOM, FOURIER, NCHOICE, LL_NSTK, NSTREAMS, NLAYERS, SURFACE_FACTOR,       &
            NSTREAMS_2, LL_NTT, LL_NSUBD, LL_NSUPD, LL_NSTKNSTRM, LL_NSTKNSTRM_2,      &
            ALBEDO, QUAD_STRMWTS, MUELLER_INDEX, DIRECT_BEAM,                          &
            HSOL, PSOL, HBVP, PBVP, STATUS_SUB, MESSAGE, TRACE_1 )

!  Exception handling

          IF ( STATUS_SUB .NE. VQS_SUCCESS ) THEN
            write(CF,'(I2)')FOURIER
            TRACE_2 = 'Error return from MSSEG_BVPSOLUTION_MASTER, '// &
                      'Called in MSSEG_LP_FOURIER, Fourier # '//CF
            STATUS = VQS_SERIOUS ; RETURN
          ENDIF

!  Stokes Vector Post Processing
!  =============================

!  Surface term for Regime 1:
!   -- 8/26/24. Now add the possibility of a diffuse flux at the surface

           IF ( Do_Surface_msst ) then

             CALL  MSSEG_POSTPROCESS_SURFACEMSST ( &
               DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,   & ! Input
               IGEOM, FOURIER, LL_NSTK, NSTREAMS, NLAYERS,  & ! Input
               MUELLER_INDEX, SURFACE_FACTOR, QUAD_STRMWTS, & ! Input
               ALBEDO, USER_BRDF_F, HSOL, PSOL, PBVP,       & ! Input 
               STOKES_DOWNSURF, SURFACE_MSST_F )              ! OUTPUT

!  otherwise, layer terms (all regimes) and optional fluxes
!   10/4/23. Ordering of output in MSSEG_POSTPROCESS_LAYERMSST was wrong (TAU_DN before TAU_UP)

           ELSE

             CALL MSSEG_POSTPROCESS_LAYERMSST ( &
               DO_UPWELLING, DO_DNWELLING, IGEOM, FOURIER, NCHOICE, LL_NSTK, & ! Input
               DELTAU_VERT, SFI_OPDEPS_SCALED, BEAM_CUTOFF,                  & ! Input
               AVERAGE_SECANT, INITIAL_TRANS, HSOL, PSOL, PBVP,              & ! Input
               TAU_UP, TAU_DN, T_UTUP_EIGEN, T_UTDN_EIGEN, T_UTDN_MUBAR, MSST_F )  ! OUTPUT

!  End Layer MSST choice, surface or atmosphere

           ENDIF

!  8/26/24. Local Flux calculation is equally valid for the layer or surface MSST option.
!    -- In the surface regime 1 case, this output will give the diffuse fluxes for the surface layer
!    -- In the surface regime 1 case, this works even if the camera is at the ground (that is NO SEGMENTS)

           IF ( DO_LOCAL_FLUXES ) THEN

!  Call new subroutine specially written to get the discrete-ordinate field
!    MK2, 6/21/21.  Use local bookkeeping variable LL_NSTK

             CALL GET_QuadIntens_Levels & 
               ( DO_UPWELLING, NCHOICE, LL_NSTK, NSTREAMS, & ! Input numbers and Flag
                 NLAYERS, HSOL, PSOL, PBVP,                & ! Solutions Input
                 QSTOKES_F )                                 ! OUTPUT

!  Set the diffuse flux outputs.  MK2, 6/21/21.  Use local bookkeeping variable LL_NSTK

             DO O1 = 1, LL_NSTK
               do k = 1, 2
                 SUM_ActFlux = DOT_PRODUCT(QUAD_HALFWTS(1:NSTREAMS),QSTOKES_F(1:NSTREAMS,O1,K))
                 SUM_RegFlux = DOT_PRODUCT(QUAD_STRMWTS(1:NSTREAMS),QSTOKES_F(1:NSTREAMS,O1,K))
                 LOCAL_FLUXES(IGEOM,O1,K,1) = SUM_ActFlux
                 LOCAL_FLUXES(IGEOM,O1,K,2) = SUM_RegFlux * PI2
               ENDDO
             ENDDO

!  End fluxes option

           ENDIF

!  #############################################
!         JACOBIAN SOLUTIONS
!  #############################################

           IF ( DO_PROFILE_LINEARIZATION ) THEN

!  Start over layers that will contain variations
!    - Set flag for variation of layer, and number of variations

             DO LAYER_TO_VARY = 1, NLAYERS
               KV = LAYER_TO_VARY
               IF ( LVARY_FLAG(LAYER_TO_VARY) ) THEN
                 QV = LVARY_NUMBER(LAYER_TO_VARY)

!  Solve the linearized boundary value problem
!  ===========================================

!    MK2, 6/21/21.  Use local bookkeeping variables LC_NSTK etc.....

                 CALL MSSEG_LP_BVPSOLUTION_MASTER ( &
                   DO_INCLUDE_SURFACE, DO_REFLECTED_DIRECTBEAM(IGEOM), DO_LAMBERTIAN_SURFACE,   &
                   DO_PLANE_PARALLEL, FOURIER, IGEOM, LL_NSTK, NSTREAMS, NLAYERS, KV, QV,       &
                   NSTREAMS_2, LL_NTT, LL_NSUBD, LL_NSUPD, LL_NSTKNSTRM, LL_NSTKNSTRM_2,        &
                   SURFACE_FACTOR, ALBEDO, BRDF_F, QUAD_STRMWTS, DIRECT_BEAM,                   &
                   DO_LAYER_SCATTERING, MUELLER_INDEX, DELTAU_SLANT, L_DELTAU_VERT,             &
                   BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, LP_INITIAL_TRANS, LP_T_DELT_MUBAR, &
                   HSOL, HBVP, PSOL, PBVP, L_HSOL, L_PSOL, L_PBVP, STATUS_SUB, MESSAGE, TRACE_1 )

!  Exception handling

                 IF ( STATUS_SUB .NE. VQS_SUCCESS ) THEN
                   write(CF,'(I2)')FOURIER ; WRITE(CKV,'(I3.3)')KV
                   TRACE_2 = 'Error return from MSSEG_LP_BVP_SOLUTION_MASTER, '// &
                             'Called in MSSEG_LP_FOURIER, Fourier # '//CF//', VARYLAYER # '//CKV
                   STATUS = VQS_SERIOUS ; RETURN
                 ENDIF

!  Checks out 10/10/23
!if ( kv.eq.1 .and. fourier.eq.0.and.igeom.eq.2.and.nchoice.eq.1)then
!  do i = 1, nlayers
!    write(5001,*)i,L_PBVP%NCON(3,i,1:2),L_PBVP%PCON(3,i,1:2)
!  enddo
!stop 'after BVP'
!endif

!  Stokes Vector Post Processing
!  =============================

!  Surface term for Regime 1:
!   -- 8/26/24. Now add the possibility of a diffuse flux at the surface

                 IF ( Do_Surface_msst ) then

!  For Jacobians, Start over layers that will contain variations
!    - Set flag for variation of layer, and number of variations
  
                   CALL MSSEG_LP_POSTPROCESS_SURFACEMSST ( &
                     DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, IGEOM, FOURIER, & ! Input control
                     LL_NSTK, NSTREAMS, NLAYERS, KV, QV, MUELLER_INDEX,         & ! Input control
                     SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, USER_BRDF_F,         & ! Input surface stuff
                     HSOL, PBVP, L_HSOL, L_PSOL, L_PBVP, LP_SURFACE_MSST_F )      ! Input TS, Output MSST

!  otherwise, layer terms (all regimes) and optional fluxes
!   10/4/23. Ordering of output in MSSEG_POSTPROCESS_LAYERMSST was wrong (TAU_DN before TAU_UP)
!   10/4/23. Ordering of TAU_UP and TAU_DN changed in MSSEG_LP_POSTPROCESS_LAYERMSST ( for consistency)

                 ELSE

                   CALL MSSEG_LP_POSTPROCESS_LAYERMSST ( &
                     DO_UPWELLING, DO_DNWELLING, IGEOM, FOURIER, NCHOICE, LL_NSTK,     & ! Input Control
                     KV, QV, DELTAU_VERT, L_DELTAU_VERT, L_SFI_OPDEPS_SCALED,          & ! Input control/Optical
                     BEAM_CUTOFF, AVERAGE_SECANT, LP_AVERAGE_SECANT, LP_INITIAL_TRANS, & ! Input Beam stuff
                     HSOL, PSOL, PBVP, L_HSOL, L_PSOL, L_PBVP,                         & ! Input TS
                     TAU_UP, TAU_DN, T_UTUP_EIGEN, T_UTDN_EIGEN, T_UTDN_MUBAR,         & ! Input partial quantities
                     L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, LP_T_UTDN_MUBAR, LP_MSST_F )

!  End Layer MSST choice (surface vs. atmospheric)

                ENDIF

!if ( fourier.eq.0.and.igeom.eq.1.and.nchoice.eq.21)then
!  do i = 1, nlayers
!    write(6001,*)i,LP_MSST_F(i,1:2,igeom,1)
!  enddo
!stop 'after PP'
!endif

!  8/26/24. Local Flux calculation is equally valid for the layer or surface MSST option.
!    -- In the surface regime 1 case, this output will give the diffuse fluxes for the surface layer
!    -- In the surface regime 1 case, this works even if the camera is at the ground (that is NO SEGMENTS)

                IF ( DO_LOCAL_FLUXES ) THEN

                  CALL GET_LP_QuadJacs_Levels & 
                     ( DO_UPWELLING, NCHOICE, LL_NSTK, NSTREAMS, NLAYERS, KV, QV, & ! Input numbers and Flag
                       HSOL, PBVP, L_HSOL, L_PSOL, L_PBVP, QJACOBIANS_F )           ! OUTPUT
 
!  Set the diffuse flux outputs
!    MK2, 6/21/21.  Use local bookkeeping variable LL_NSTK

                   DO Q = 1, QV
                     DO O1 = 1, LL_NSTK
                       do k = 1, 2
                         SUM_ActFlux = DOT_PRODUCT(QUAD_HALFWTS(1:NSTREAMS),QJACOBIANS_F(Q,1:NSTREAMS,O1,K))
                         SUM_RegFlux = DOT_PRODUCT(QUAD_STRMWTS(1:NSTREAMS),QJACOBIANS_F(Q,1:NSTREAMS,O1,K))
                         LOCAL_LP_FLUXES(KV,Q,IGEOM,O1,K,1) = SUM_ActFlux
                         LOCAL_LP_FLUXES(KV,Q,IGEOM,O1,K,2) = SUM_RegFlux
                       ENDDO
                     ENDDO
                   ENDDO

!  End fluxes option

                   ENDIF

!  End profile Jacobian do-loop

              ENDIF
            ENDDO
          ENDIF

!  End loop over Geometrical solutions

        ENDIF
      END DO

!  finish

RETURN
END SUBROUTINE MSSEG_LP_FOURIER

!  End module

END MODULE MSSeg_LP_Fourier_m

