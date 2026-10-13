
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

!  05 April 2023. Revision

      MODULE MSSeg_Partic_Fourier_m

      public ::   MSSEG_PARTIC_FOURIER

      contains

      SUBROUTINE MSSEG_PARTIC_FOURIER ( &
           DO_SURFACE_MSST, DO_UPWELLING, DO_DIRECT_BEAM, DO_MULTIGEOM,    &
           DO_LOCAL_FLUXES, DO_FOURIER0_NSTOKES2, DO_INCLUDE_SURFACE,      &
           DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,     &
           FOURIER, NSTOKES, NSTREAMS, NLAYERS, NGEOMS, NCHOICE,           &
           NMOMENTS, NSTREAMS_2, NTOTAL, MUELLER_INDEX, DMAT, DFLUX,       &
           N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2,             &
           DO_LAYER_SCATTERING, FLUX_FACTOR, FLUXVEC, COS_SZANGLES,        &
           USER_STREAMS, QUAD_STREAMS, QUAD_STRMWTS, QUAD_HALFWTS,         &
           LAMBERTIAN_ALBEDO, DELTAU_VERT, SFI_OPDEPS_SCALED, OMEGA_GREEK, &           
           BRDF_F_0, USER_BRDF_F, SLTERM_ISOTROPIC, SLTERM_F_0,            & 
           BEAM_CUTOFF, TRANS_SOLAR_BEAM, INITIAL_TRANS,                   &
           AVERAGE_SECANT, T_DELT_MUBAR, HGSF, HSOL, HBVP,                 &
           MSST_F, SURFACE_MSST_F, LOCAL_FLUXES, PGSF, PSOL, PBVP,         &
           STATUS, MESSAGE, TRACE_1, TRACE_2 )

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- ASYMTX_TOLERANCE now introduced and passed into the QHOM call to ASYMTX
!   -- Flux output (LOCAL_FLUXES) controlled by DO_LOCAL_FLUXES input flag
!   -- DO_FOURIER0_NSTOKES2 performance flag introduced.
!   -- Call subroutine specially written to get the discrete-ordinate field for Fluxes

!  Mk3. Experimental. 8/27/21.
!    -- Use internal Type structures. This is the Particular Fourier.
!    -- April 2023. Introduce DO_MULTIGEOM

      USE MSSeg_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES_SQ, MAXSTREAMS, MAXFOURIER, &
                               MAXBEAMS, MAX_SZANGLES, MAX_USER_STREAMS, MAX_GEOMETRIES, MAXEVALUES, &
                               ZERO, ONE, TWO, PI2, VQS_SUCCESS, VQS_SERIOUS

      use MSSeg_HPTypes_Def_m

      use MSSeg_Partic_Solutions_m
      use MSSeg_Partic_BVProblem_m
      use MSSeg_Partic_PostProcess_m

      IMPLICIT NONE

!  INPUTS
!  ======

!  Input flags
!  -----------

!  Input flags

      LOGICAL, INTENT (IN) ::          DO_SURFACE_MSST
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DIRECT_BEAM

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::          DO_SURFACE_LEAVING
      LOGICAL, INTENT (IN) ::          DO_SL_ISOTROPIC

!  Mk2, 6/21/21. Add the following inputs
!   -- DO_LOCAL_FLUXES input flag for controlling flux calculations
!   -- DO_FOURIER0_NSTOKES2 performance flag introduced.

      LOGICAL, INTENT (IN) ::          DO_LOCAL_FLUXES 
      LOGICAL, INTENT (IN) ::          DO_FOURIER0_NSTOKES2


!  04/07/23. Mk3. Add the MultiGoem variable
 
      LOGICAL, INTENT (IN) ::          DO_MULTIGEOM ( MAX_GEOMETRIES, 0:MAXFOURIER )

!  Integers
!  --------

!  Control integers (basic)

      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NGEOMS
      INTEGER, INTENT (IN) ::          NCHOICE

!  Control integers (derived). These are local values

      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::          NMOMENTS

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

!  From the QSPREP (average secant) setup

      INTEGER, INTENT (INOUT) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: TRANS_SOLAR_BEAM ( MAXBEAMS )

      DOUBLE PRECISION, INTENT (INOUT) :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  surface reflection
!   -- Mk3. Reduce dimensioning, drop 0:MAXMOMENTS

      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F_0    ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS   )
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )

!  surface leaving
!   -- Mk3. Reduce dimensioning, drop 0:MAXMOMENTS

      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_ISOTROPIC ( MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_F_0 ( MAXSTOKES, MAXSTREAMS, MAXBEAMS )

!  From the Deltam-scaling setup

      DOUBLE PRECISION, INTENT (INOUT) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: SFI_OPDEPS_SCALED  ( MAX_GEOMETRIES )

!  From the SSALB_INIT setup

      DOUBLE PRECISION, INTENT (INOUT) :: OMEGA_GREEK  ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  Type structure input
!  --------------------

      type(Homog_GSF)      , intent(in) :: HGSF
      type(Homog_Solutions), intent(in) :: HSOL
      type(Homog_BVPSetup) , intent(in) :: HBVP

!  OUTPUTS
!  =======

!  Mk2, 6/21/21. Add the Local Fluxes output (Intent(InOut))

      DOUBLE PRECISION, INTENT (OUT)   ::  MSST_F         ( MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT)   ::  SURFACE_MSST_F ( MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) ::  LOCAL_FLUXES   ( MAX_GEOMETRIES, MAXSTOKES, 2, 2 )

!  Type structure output
!  ---------------------

      type(Partic_GSF)      , intent(inout) :: PGSF
      type(Partic_Solutions), intent(inout) :: PSOL
      type(Partic_BVPSolve) , intent(inout) :: PBVP

!  Bookkeeping and exception handling
!  ----------------------------------
!mick fix 6/6/2023 - moved DO_INCLUDE_SURFACE to inputs above and defined as INTENT(IN)

      !LOGICAL, INTENT (OUT) ::           DO_INCLUDE_SURFACE
      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (OUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (OUT) :: TRACE_1
      CHARACTER (LEN=*), INTENT (OUT) :: TRACE_2

!  Local variables
!  ===============

!  downwelling (opposite of upwelling)

      LOGICAL          :: DO_DNWELLING

!  Local transmittances

      DOUBLE PRECISION :: T_UTUP_EIGEN ( MAXEVALUES )
      DOUBLE PRECISION :: T_UTDN_EIGEN ( MAXEVALUES )
      DOUBLE PRECISION :: T_UTDN_MUBAR, TAU_DN, TAU_UP

!  Local surface variables

      LOGICAL          :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )
      DOUBLE PRECISION :: DIRECT_BEAM  ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION :: ATMOS_ATTN   ( MAXBEAMS )
      DOUBLE PRECISION :: STOKES_DOWNSURF ( MAXSTREAMS, MAXSTOKES )

!  Mk2, 6/21/21. Local variables for flux computations

      DOUBLE PRECISION :: SUM_ActFlux, SUM_RegFlux
      DOUBLE PRECISION :: QSTOKES_F ( MAXSTREAMS, MAXSTOKES, 2 )

!  Mk2, 6/21/21. Based on VLIDORT Version 2.8.3. Local NSTOKES and related control integers

      INTEGER :: LC_NSTK, LC_NSTKNSTRM, LC_NSTKNSTRM_2, LC_NTT, LC_NSUBD, LC_NSUPD

!  Flux multiplier and Fourier component numbers

      DOUBLE PRECISION ::  DELTA_FACTOR
      DOUBLE PRECISION ::  SURFACE_FACTOR

!  indices

      INTEGER          :: LAYER, IGEOM, O1, K, N

!  error tracing

      INTEGER ::           STATUS_SUB
      CHARACTER (LEN=2) :: CF

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

      MSST_F         = zero
      SURFACE_MSST_F = zero
      if ( FOURIER .EQ. 0 ) THEN
         LOCAL_FLUXES   = zero
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
        DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,  &
        FOURIER, NSTOKES, NSTREAMS, NGEOMS, DELTA_FACTOR,            &
        FLUX_FACTOR, FLUXVEC, MUELLER_INDEX, COS_SZANGLES,           &
        LAMBERTIAN_ALBEDO, BRDF_F_0, SLTERM_ISOTROPIC, SLTERM_F_0,   &
        TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM, &
        ATMOS_ATTN, DIRECT_BEAM )

!  ###################
!  Spherical functions
!  ###################

!  Get Pi Matrices for this Fourier component

      CALL MSSEG_PARTIC_PIMATRIX_SETUP ( &
        FOURIER, NSTOKES, NLAYERS, NGEOMS, NMOMENTS, &
        COS_SZANGLES, USER_STREAMS, DMAT, HGSF, PGSF )

!  4/7/23. debug for PIMM
!if ( ( Fourier.eq.3.or.Fourier.eq.2 ) .and.nchoice.eq.75 ) then
!      write(*,*)'right after MSSEG_PARTIC_PIMATRIX_SETUP',FOURIER,PGSF%PIMM_11(FOURIER,7)
!endif

!  ###################
!  LOCAL NSTOKES USAGE
!  ###################

!  Mk2, 6/21/21. Based on VLIDORT Version 2.8.3. 
!    -- Introduce control for Fast execution of Fourier-zero with LOCAL_NSTOKES = 2

      IF ( DO_Fourier0_NSTOKES2 .and. NSTOKES.eq.3 .and. FOURIER.eq.0 .and. DO_LAMBERTIAN_SURFACE ) THEN
        LC_NSTK = 2 ; LC_NSTKNSTRM = NSTREAMS * LC_NSTK ; LC_NSTKNSTRM_2 = LC_NSTKNSTRM * 2
        LC_NTT  = LC_NSTKNSTRM_2 * NLAYERS
        IF ( NLAYERS .EQ. 1 ) THEN
          LC_NSUBD = 2*LC_NSTKNSTRM - 1 ; LC_NSUPD = LC_NSUBD
        ELSE
          LC_NSUBD = 3*LC_NSTKNSTRM - 1 ; LC_NSUPD = LC_NSUBD
        ENDIF
      ELSE
        LC_NSTK = NSTOKES ; LC_NSTKNSTRM   = NSTKS_NSTRMS ; LC_NSTKNSTRM_2 = NSTKS_NSTRMS_2
        LC_NTT  = NTOTAL  ; LC_NSUBD       = N_SUBDIAG    ; LC_NSUPD       = N_SUPDIAG
      ENDIF

!  ########################################
!  RT differential equation Eigensolutions
!  ########################################

!   Get user homog solutions for layer NCHOICE
!    MK2, 6/21/21. Use local variables LC_NSTK, LC_NSTKNSTRM

      CALL MSSEG_UHOM_SOLUTION ( &
        NCHOICE, FOURIER, LC_NSTK, NSTREAMS, NMOMENTS, NGEOMS, &
        DO_UPWELLING, DO_DNWELLING, DO_LAYER_SCATTERING, &
        QUAD_HALFWTS, OMEGA_GREEK, HGSF, PGSF, HSOL, PSOL )

!  Start loop over desired solutions

      DO IGEOM = 1, NGEOMS
        IF ( DO_MULTIGEOM(IGEOM,FOURIER) ) THEN

!  Solar beam Particular solutions (classical method)
!  --------------------------------------------------

          DO LAYER = 1, NLAYERS

!  get the solution
!    MK2, 6/21/21.  Use local bookkeeping variables LC_NSTK etc.....

            CALL MSSEG_QBEAM_SOLUTION ( &
              LAYER, FOURIER, IGEOM, LC_NSTK, NSTREAMS, NMOMENTS,            &
              NSTREAMS_2, LC_NSTKNSTRM, DO_LAYER_SCATTERING, DFLUX,          &
              FLUX_FACTOR, QUAD_STREAMS, OMEGA_GREEK, BEAM_CUTOFF,           &
              T_DELT_MUBAR, INITIAL_TRANS, AVERAGE_SECANT, PGSF, HGSF, HSOL, &
              PSOL, STATUS_SUB, MESSAGE, TRACE_1 )

!  debug
!        if ( fourier .eq. 0 ) then
!            if (ngeoms.gt.5) write(*,*)igeom,layer,PSOL%WUPPER(1:4,1,layer)
!            if ( layer.eq.nlayers ) stop 'hhhhhhh'
!        endif

!if ( fourier.gt.0)write(76+fourier,*)NCHOICE,IGEOM,LAYER,PSOL%WUPPER(1:4,1,LAYER)

!  error flag

            IF ( STATUS_SUB .NE. VQS_SUCCESS ) THEN
              write(CF,'(I2)')FOURIER
              TRACE_2 = 'Error return from MSSEG_QBEAM_SOLUTION, '// &
                        'Called in MSSEG_PARTIC_FOURIER, Fourier # '//CF
              STATUS = VQS_SERIOUS ; RETURN
            ENDIF

!  end layer loop

          END DO

!  user solutions for Layer of choice
!    MK2, 6/21/21.  Use local bookkeeping variables LC_NSTK etc.....

          CALL MSSEG_UBEAM_SOLUTION ( &
            NCHOICE, FOURIER, IGEOM,                         &
            DO_UPWELLING, DO_DNWELLING, DO_LAYER_SCATTERING, &
            LC_NSTK, NSTREAMS, NMOMENTS, FLUX_FACTOR,        &
            BEAM_CUTOFF, QUAD_HALFWTS, OMEGA_GREEK,          &
            HGSF, PGSF, PSOL )

!  seems to work
!  if (fourier.eq.1.and.igeom.eq.6)WRITE(901,*)nchoice,PSOL%UPAR_UP_2(1)

!  Solve boundary value problem
!  ----------------------------

!  standard case using compressed-band matrices, etc..
!    MK2, 6/21/21.  Use local bookkeeping variables LC_NSTK etc.....
!    Mk2. 8/3/21. Add NCHOICE input (for debug purposes)

          CALL MSSEG_BVPSOLUTION_MASTER ( &
            DO_INCLUDE_SURFACE, DO_REFLECTED_DIRECTBEAM(IGEOM), DO_LAMBERTIAN_SURFACE, &
            IGEOM, FOURIER, NCHOICE, LC_NSTK, NSTREAMS, NLAYERS, SURFACE_FACTOR,       &
            NSTREAMS_2, LC_NTT, LC_NSUBD, LC_NSUPD, LC_NSTKNSTRM, LC_NSTKNSTRM_2,      &
            LAMBERTIAN_ALBEDO, QUAD_STRMWTS, MUELLER_INDEX, DIRECT_BEAM,               &
            HSOL, PSOL, HBVP, PBVP, STATUS_SUB, MESSAGE, TRACE_1 )

!  Exception handling

          IF ( STATUS_SUB .NE. VQS_SUCCESS ) THEN
            write(CF,'(I2)')FOURIER
            TRACE_2 = 'Error return from MSSEG_BVPSOLUTION_MASTER, '// &
                      'Called in MSSEG_PARTIC_FOURIER, Fourier # '//CF
            STATUS = VQS_SERIOUS ; RETURN
          ENDIF

!if ( fourier.eq.0.and.igeom.eq.1.and.nchoice.eq.nlayers.and..not.Do_Surface_msst ) then
!   do n = 1, LC_NTT
!      write(911,*)n,PBVP%COL2(n,igeom)
!   enddo
!endif

!  Stokes Vector Post Processing
!  =============================
!mick note 6/6/2023 - Although computed, the vars
!                       * MSSEG_POSTPROCESS_SURFACEMSST: STOKES_DOWNSURF
!                       * MSSEG_POSTPROCESS_LAYERMSST  : TAU_DN, TAU_UP, T_UTUP_EIGEN, T_UTDN_EIGEN,
!                                                        T_UTDN_MUBAR
!                     do not appear to be used currently - only the MSSTs from these subs are used.

!  Surface term for Regime 1:
!   -- 8/26/24. Now add the possibility of a diffuse flux at the surface

          IF ( DO_SURFACE_MSST ) then

            CALL MSSEG_POSTPROCESS_SURFACEMSST ( &
              DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,        & ! Input
              IGEOM, FOURIER, LC_NSTK, NSTREAMS, NLAYERS,       & ! Input
              MUELLER_INDEX, SURFACE_FACTOR, QUAD_STRMWTS,      & ! Input
              LAMBERTIAN_ALBEDO, USER_BRDF_F, HSOL, PSOL, PBVP, & ! Input 
              STOKES_DOWNSURF, SURFACE_MSST_F )                   ! OUTPUT


!  Otherwise, layer terms (all regimes) and optional fluxes
!   10/4/23. Ordering of output in MSSEG_POSTPROCESS_LAYERMSST was wrong (TAU_DN before TAU_UP)

          ELSE

            CALL MSSEG_POSTPROCESS_LAYERMSST ( &
              DO_UPWELLING, DO_DNWELLING, IGEOM, FOURIER, NCHOICE, LC_NSTK,   &  ! Input
              DELTAU_VERT, SFI_OPDEPS_SCALED, BEAM_CUTOFF,                    &  ! Input
              AVERAGE_SECANT, INITIAL_TRANS, HSOL, PSOL, PBVP,                &  ! Input
              TAU_UP, TAU_DN, T_UTUP_EIGEN, T_UTDN_EIGEN, T_UTDN_MUBAR, MSST_F ) ! OUTPUT

!  debug
!if ( fourier.eq.0.and.igeom.eq.1 ) then
!   write(951,*)NCHOICE, MSST_F(IGEOM,1)
!endif

!  End Layer MSST choice

          ENDIF

!  8/26/24. Local Flux calculation is equally valid for the layer or surface MSST option.
!    -- In the surface regime 1 case, this output will give the diffuse fluxes for the surface layer
!    -- In the surface regime 1 case, this works even if the camera is at the ground (that is NO SEGMENTS)

          IF ( DO_LOCAL_FLUXES ) THEN

!  Call new subroutine specially written to get the discrete-ordinate field
!    MK2, 6/21/21.  Use local bookkeeping variable LC_NSTK

             CALL GET_QuadIntens_Levels & 
               ( DO_UPWELLING, NCHOICE, LC_NSTK, NSTREAMS, & ! Input numbers and Flag
                 NLAYERS, HSOL, PSOL, PBVP,                & ! Solutions Input
                 QSTOKES_F )                                 ! OUTPUT

!  Set the diffuse flux outputs.  MK2, 6/21/21.  Use local bookkeeping variable LC_NSTK

             DO O1 = 1, LC_NSTK
                do k = 1, 2
                   SUM_ActFlux = DOT_PRODUCT(QUAD_HALFWTS(1:NSTREAMS),QSTOKES_F(1:NSTREAMS,O1,K))
                   SUM_RegFlux = DOT_PRODUCT(QUAD_STRMWTS(1:NSTREAMS),QSTOKES_F(1:NSTREAMS,O1,K))
                   LOCAL_FLUXES(IGEOM,O1,K,1) = SUM_ActFlux
                   LOCAL_FLUXES(IGEOM,O1,K,2) = SUM_RegFlux * PI2
                ENDDO
             ENDDO

!  End fluxes option

          ENDIF

!  End loop over solutions

        ENDIF
      END DO

!  finish

RETURN
END SUBROUTINE MSSEG_PARTIC_FOURIER

!  End module

END MODULE MSSeg_Partic_Fourier_m

