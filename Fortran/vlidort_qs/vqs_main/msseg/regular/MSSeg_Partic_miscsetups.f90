
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

!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension, Only use quantities you need

      MODULE msseg_partic_miscsetups_m

      PRIVATE :: MSSEG_PARTIC_CHECKINPUT, MSSEG_PARTIC_DERIVEINPUT,  &
                 MSSEG_PARTIC_DELTAMSCALE,  MSSEG_CHAPMANFACTORS, MSSEG_QSPREP
      PUBLIC  :: MSSEG_PARTIC_MISCSETUPS
                 
      CONTAINS

      SUBROUTINE MSSEG_PARTIC_MISCSETUPS ( &
        DO_DELTAM_SCALING, DO_PLANE_PARALLEL, NSTOKES, NLAYERS,      & ! Input
        EARTH_RADIUS, N_SZANGLES, N_USER_RELAZMS, N_USER_VZANGLES,   & ! Input
        SZANGLES, USER_RELAZMS, USER_VZANGLES, HEIGHT_GRID,          & ! Input
        DELTAU_VERT_INPUT, DELTAU_VERT, TAUGRID, SFI_OPDEPS, FAC1,   & ! Input
        DO_NO_AZIMUTH, DO_DIRECT_BEAM, COS_SZANGLES, USER_STREAMS,   & ! Output
        NBEAMS, CHAPMAN_FACTORS, DELTAU_SLANT, SFI_OPDEPS_SCALED,    & ! Output
        LAYER_PIS_CUTOFF, TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM, & ! Output
        INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                 & ! Output
        STATUS, NMESSAGES, MESSAGES, ACTIONS )                         ! Output

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need

      USE MSSeg_PARS_m, Only : MAXSTREAMS, MAXLAYERS, MAXBEAMS, MAX_MESSAGES, MAX_GEOMETRIES,  &
                               MAX_SZANGLES, MAX_USER_RELAZMS, MAX_USER_VZANGLES, VQS_SERIOUS


      IMPLICIT NONE

!  Input arguments
!  ===============

!  flags

      LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

!  Numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS

!  Earth radius (km), Coarse grids of heights

      DOUBLE PRECISION, INTENT (IN) :: EARTH_RADIUS
      DOUBLE PRECISION, INTENT (IN) :: HEIGHT_GRID (0:MAXLAYERS)

!  Optical

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SFI_OPDEPS  ( MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (IN) :: FAC1 ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: TAUGRID ( 0:MAXLAYERS )

!  Check variables

      INTEGER, INTENT (INOUT) ::          N_SZANGLES
      INTEGER, INTENT (INOUT) ::          N_USER_RELAZMS
      INTEGER, INTENT (IN) ::             N_USER_VZANGLES
      DOUBLE PRECISION, INTENT (IN) ::    SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::    USER_RELAZMS  ( MAX_USER_RELAZMS )
      DOUBLE PRECISION, INTENT (INOUT) :: USER_VZANGLES ( MAX_USER_VZANGLES )

!  output
!  ======

!  Derive output

      LOGICAL, INTENT (OUT) ::            DO_DIRECT_BEAM
      LOGICAL, INTENT (INOUT) ::          DO_NO_AZIMUTH
      INTEGER, INTENT (OUT) ::            NBEAMS
      DOUBLE PRECISION, INTENT (OUT) ::   COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (OUT) ::   USER_STREAMS ( MAX_USER_VZANGLES )

!  deltam scaling output

      DOUBLE PRECISION, INTENT (OUT) :: SFI_OPDEPS_SCALED  ( MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Chapman factor output

      DOUBLE PRECISION, INTENT (OUT) :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Beam parameterization output

      INTEGER, INTENT (OUT) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: TRANS_SOLAR_BEAM  ( MAXBEAMS )
      LOGICAL, INTENT (OUT) ::          DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Output exception handling

      INTEGER, INTENT (OUT) ::            STATUS
      INTEGER, INTENT (OUT) ::            NMESSAGES
      CHARACTER (LEN=*), INTENT (OUT) ::  MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=*), INTENT (OUT) ::  ACTIONS  ( 0:MAX_MESSAGES )

!  Local

      DOUBLE PRECISION               :: LOCAL_CSZA ( 0:MAXLAYERS, MAXBEAMS )

!  Check

      CALL MSSEG_PARTIC_CHECKINPUT ( &
        N_SZANGLES, N_USER_RELAZMS, N_USER_VZANGLES, & ! Input
        SZANGLES, USER_RELAZMS, USER_VZANGLES,       & ! Input
        DO_DIRECT_BEAM, STATUS, NMESSAGES, MESSAGES, ACTIONS ) ! output

      IF ( STATUS.eq. VQS_SERIOUS ) RETURN

!  Derive

      CALL MSSEG_PARTIC_DERIVEINPUT ( &
        NSTOKES, N_SZANGLES, SZANGLES, N_USER_VZANGLES, USER_VZANGLES,  & ! Input
        NBEAMS, COS_SZANGLES, DO_NO_AZIMUTH, USER_STREAMS )               ! output

!  Chapman factors

      CALL MSSEG_CHAPMANFACTORS ( &
        DO_PLANE_PARALLEL, N_SZANGLES, NLAYERS, &
        SZANGLES, EARTH_RADIUS, HEIGHT_GRID, CHAPMAN_FACTORS )

!  Delta-m scaling of input quantities

      CALL MSSEG_PARTIC_DELTAMSCALE ( &
        DO_DELTAM_SCALING, NLAYERS, NBEAMS, CHAPMAN_FACTORS, & ! Input
        DELTAU_VERT_INPUT, SFI_OPDEPS, FAC1,                 & ! Input
        DELTAU_SLANT, SFI_OPDEPS_SCALED )                      ! output

!  Prepare quasi-spherical attenuation

      CALL MSSEG_QSPREP ( &
        DO_PLANE_PARALLEL, NLAYERS, NBEAMS, COS_SZANGLES, & ! Input
        DELTAU_VERT, TAUGRID, DELTAU_SLANT,               & ! Input
        LAYER_PIS_CUTOFF,  DO_REFLECTED_DIRECTBEAM,       & ! Output
        INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,      & ! Output
        TRANS_SOLAR_BEAM, LOCAL_CSZA )                      ! output

!  Finish

      RETURN
      END SUBROUTINE MSSEG_PARTIC_MISCSETUPS

!

      SUBROUTINE MSSEG_PARTIC_DELTAMSCALE ( &
        DO_DELTAM_SCALING, NLAYERS, NBEAMS, CHAPMAN_FACTORS, & ! Input
        DELTAU_VERT_INPUT, SFI_OPDEPS, FAC1,                 & ! Input
        DELTAU_SLANT, SFI_OPDEPS_SCALED )

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need

      USE MSSeg_PARS_m, Only : MAXLAYERS, MAXBEAMS, MAX_GEOMETRIES, ZERO

      IMPLICIT NONE

!  Input

      LOGICAL, INTENT (IN) ::           DO_DELTAM_SCALING
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (INOUT) ::        NBEAMS

      DOUBLE PRECISION, INTENT (IN)    :: FAC1 ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN)    :: SFI_OPDEPS  ( MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (INOUT) :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN)    :: DELTAU_VERT_INPUT ( MAXLAYERS )

!  Output

      DOUBLE PRECISION, INTENT (OUT)   :: SFI_OPDEPS_SCALED  ( MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT)   :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Local variables

      DOUBLE PRECISION :: DELS
      INTEGER          :: K, N, IB

!  First do the slant input
!  ------------------------

!  Slant optical thickness values

!mick fix 7/23/2014 - initialized for packing
!  Initialize the active workspace only; unused capacity is never consumed.
      DELTAU_SLANT(1:NLAYERS,1:NLAYERS,1:NBEAMS) = ZERO

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          DO K = 1, N
            DELS = CHAPMAN_FACTORS(N,K,IB)
            DELTAU_SLANT(N,K,IB) = DELTAU_VERT_INPUT(K) * DELS
          ENDDO
        ENDDO
      ENDDO

!  DELTAM SCALING
!  ==============

      IF ( DO_DELTAM_SCALING ) THEN

!  Scale partial layer stuff. No !!!!!!!!!!
!    OPDEPS from FO code are Already scaled! Bug 18 November 2015
!    Just need to copy the input to output

        DO IB = 1, NBEAMS
           SFI_OPDEPS_SCALED(IB) = SFI_OPDEPS(IB)
!          SFI_OPDEPS_SCALED(IB) = FAC1(NCHOICE) * SFI_OPDEPS(IB)
        ENDDO

!  Scale layer path thickness values

        DO N = 1, NLAYERS
          DO K = 1, N
            DO IB = 1, NBEAMS
              DELTAU_SLANT(N,K,IB) = DELTAU_SLANT(N,K,IB) * FAC1(K)
            ENDDO
          ENDDO
        ENDDO

!  NO DELTA-M SCALING
!  ==================

!  move input geophysical variables to Workspace quantities

      ELSE

!  Scale partial layer stuff

        DO IB = 1, NBEAMS
           SFI_OPDEPS_SCALED(IB) = SFI_OPDEPS(IB)
        ENDDO

      ENDIF

!  Finish module

      RETURN
      END SUBROUTINE MSSEG_PARTIC_DELTAMSCALE

!

      SUBROUTINE MSSEG_QSPREP ( &
        DO_PLANE_PARALLEL, NLAYERS, NBEAMS, COS_SZANGLES, & ! Input
        DELTAU_VERT, TAUGRID, DELTAU_SLANT,               & ! Input
        LAYER_PIS_CUTOFF,  DO_REFLECTED_DIRECTBEAM,       & ! Output
        INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,      & ! Output
        TRANS_SOLAR_BEAM, LOCAL_CSZA )

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need

      USE MSSeg_PARS_m, Only : MAXLAYERS, MAX_SZANGLES, MAXBEAMS, ZERO, ONE, MAX_TAU_SPATH

      IMPLICIT NONE

      LOGICAL, INTENT (IN)    ::           DO_PLANE_PARALLEL
      INTEGER, INTENT (IN)    ::           NLAYERS
      INTEGER, INTENT (INOUT) ::           NBEAMS
      DOUBLE PRECISION, INTENT (INOUT) ::  COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN)    ::  DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN)    ::  TAUGRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

      INTEGER, INTENT (OUT) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: TRANS_SOLAR_BEAM ( MAXBEAMS )
      LOGICAL, INTENT (OUT) :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: LOCAL_CSZA  ( 0:MAXLAYERS, MAXBEAMS )

!  Local variables
!  ---------------

      INTEGER          :: N, K, IB
      DOUBLE PRECISION :: S_T_0, S_T_1, SEC0, TAU, TAU_SOLAR(MAXBEAMS)
      DOUBLE PRECISION :: TAUSLANT ( 0:MAXLAYERS, MAXBEAMS ), SPHER

!  Initialize

      DO IB = 1, NBEAMS
        LAYER_PIS_CUTOFF(IB) = NLAYERS
      ENDDO

!  TOA  cosines

      DO IB = 1, NBEAMS
        LOCAL_CSZA(0,IB) = COS_SZANGLES(IB)
      ENDDO

!  plane-parallel case
!  -------------------

      IF ( DO_PLANE_PARALLEL ) THEN

       DO IB = 1, NBEAMS
        SEC0 = ONE / COS_SZANGLES(IB)
!        LAYER_PIS_CUTOFF(IB) = NLAYERS
        DO N = 1, NLAYERS
          TAUSLANT(N,IB) = TAUGRID(N) * SEC0
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            IF ( TAUSLANT(N,IB) .GT. MAX_TAU_SPATH ) THEN
              LAYER_PIS_CUTOFF(IB) = N
            ENDIF
            AVERAGE_SECANT(N,IB) = SEC0
            INITIAL_TRANS(N,IB)  = DEXP ( - TAUGRID(N-1) * SEC0 )
            LOCAL_CSZA(N,IB)      = COS_SZANGLES(IB)
          ELSE
            AVERAGE_SECANT(N,IB) = ZERO
            INITIAL_TRANS(N,IB)  = ZERO
            LOCAL_CSZA(N,IB)      = ZERO
          ENDIF
        ENDDO
        TAU_SOLAR(IB) = TAUSLANT(NLAYERS,IB)
       ENDDO

      ELSE

!  pseudo-spherical case
!  ---------------------

       DO IB = 1, NBEAMS

!  Get the total spherical attenuation from layer thickness sums

        TAUSLANT(0,IB) = ZERO
        DO N = 1, NLAYERS
          TAU = ZERO
          DO K = 1, N
            TAU = TAU + DELTAU_SLANT(N,K,IB)
          ENDDO
          TAUSLANT(N,IB) = TAU
        ENDDO
        TAU_SOLAR(IB) = TAUSLANT(NLAYERS,IB)

!  set up the average secant formulation

        S_T_1 = ZERO
        S_T_0 = ONE
        LAYER_PIS_CUTOFF(IB) = NLAYERS
        DO N = 1, NLAYERS
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            IF ( TAUSLANT(N,IB) .GT. MAX_TAU_SPATH ) THEN
              LAYER_PIS_CUTOFF(IB) = N
            ELSE
              S_T_1 = DEXP ( - TAUSLANT(N,IB) )
            ENDIF
            AVERAGE_SECANT(N,IB) = (TAUSLANT(N,IB)-TAUSLANT(N-1,IB)) &
                                     / DELTAU_VERT(N)
            INITIAL_TRANS(N,IB)  = S_T_0
            S_T_0             = S_T_1
          ELSE
            AVERAGE_SECANT(N,IB) = ZERO
            INITIAL_TRANS(N,IB)  = ZERO
          ENDIF
        ENDDO

!  Set the Local solar zenith cosines. Non-refractive-only.

        LOCAL_CSZA(:,IB) = ZERO
        DO N = 1, NLAYERS
            IF(N.LE.LAYER_PIS_CUTOFF(IB)) LOCAL_CSZA(N,IB) = COS_SZANGLES(IB)
        ENDDO

       ENDDO
      ENDIF

!  Set Direct Beam Flag and solar beam total attenuation to surface

      DO IB = 1, NBEAMS
        TRANS_SOLAR_BEAM(IB) = ZERO
        DO_REFLECTED_DIRECTBEAM(IB) = .FALSE.
        IF ( TAU_SOLAR(IB) .LT. MAX_TAU_SPATH ) THEN
           TRANS_SOLAR_BEAM(IB) = DEXP( - TAU_SOLAR(IB) )
           DO_REFLECTED_DIRECTBEAM(IB) = .TRUE.
        ENDIF
      ENDDO

!  layer transmittances

      T_DELT_MUBAR = zero
      DO IB = 1, NBEAMS
       DO N = 1, NLAYERS
        IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
          SPHER = DELTAU_VERT(N) * AVERAGE_SECANT(N,IB)
          IF ( SPHER .LE. MAX_TAU_SPATH ) THEN
            T_DELT_MUBAR(N,IB) = DEXP ( - SPHER )
          ENDIF
        ENDIF
       ENDDO
      ENDDO

! debug
!      do n = 1, nlayers
!        write(*,*)AVERAGE_SECANT(N,1),  INITIAL_TRANS(N,1)
!      enddo
!      pause

!  finish

      RETURN
      END SUBROUTINE MSSEG_QSPREP

!

      SUBROUTINE MSSEG_CHAPMANFACTORS ( &
        DO_PLANE_PARALLEL, N_SZANGLES, NLAYERS, &
        SZANGLES, EARTH_RADIUS, HEIGHT_GRID, CHAPMAN_FACTORS )

!  parameter file

      USE MSSeg_PARS_m, only : MAX_SZANGLES, MAXBEAMS, MAXLAYERS, &
                                 ZERO, ONE, DEG_TO_RAD

      IMPLICIT NONE

!  Input arguments
!  ===============

!  flag for plane parallel case

      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

!  number of coarse layers, beams

      INTEGER, INTENT (IN) ::          NLAYERS, N_SZANGLES

!  True solar zenith angle (degrees)

      DOUBLE PRECISION, INTENT (IN) :: SZANGLES ( MAX_SZANGLES )

!  Earth radius (km)

      DOUBLE PRECISION, INTENT (IN) :: EARTH_RADIUS


!  Coarse grids of heights

      DOUBLE PRECISION, INTENT (IN) :: HEIGHT_GRID (0:MAXLAYERS)

!  Output arguments
!  ================

!  Path segments distance factors

      DOUBLE PRECISION, INTENT (OUT) :: CHAPMAN_FACTORS (MAXLAYERS, MAXLAYERS, MAXBEAMS)

!  Local variables
!  ===============

!  local height arrays

      DOUBLE PRECISION :: H(0:MAXLAYERS)
      DOUBLE PRECISION :: DELZ(MAXLAYERS)

!  help variables

      INTEGER          :: N, K, IB
      DOUBLE PRECISION :: GM_TOA, TH_TOA, MU_TOA

      DOUBLE PRECISION :: STH1, SINTH1, STH2, SINTH2, SECBAR, &
            PHI, SINPHI, RE_LOWER, RE_UPPER, DIST, STH2D

!  Some setup operations
!  =====================

!  initialise output

!  Initialize the active workspace only; unused capacity is never consumed.
      CHAPMAN_FACTORS(1:NLAYERS,1:NLAYERS,1:N_SZANGLES) = ZERO

!  earth radii and heights differences

      DO N = 0, NLAYERS
        H(N) = HEIGHT_GRID(N) + EARTH_RADIUS
      ENDDO
      DO N = 1, NLAYERS
        DELZ(N) = HEIGHT_GRID(N-1)-HEIGHT_GRID(N)
      ENDDO

!  Start Beam loop

      DO IB = 1, N_SZANGLES

!  TOA values

        TH_TOA = SZANGLES(IB) * DEG_TO_RAD
        MU_TOA = COS(TH_TOA)
        GM_TOA = SQRT ( ONE - MU_TOA * MU_TOA )
        STH2D  = ZERO

!  Compute

        IF ( DO_PLANE_PARALLEL ) THEN
          SECBAR = ONE / MU_TOA
          DO N = 1, NLAYERS
            CHAPMAN_FACTORS(N,1:N,IB) = SECBAR
          ENDDO
        ELSE
          IF ( GM_TOA .eq. ZERO ) then
            DO N = 1, NLAYERS
              CHAPMAN_FACTORS(N,1:N,IB) = ZERO
            ENDDO
          ELSE
            DO N = 1, NLAYERS
              SINTH1 = GM_TOA * H(N) / H(0)
              STH1   = ASIN(SINTH1)
              RE_UPPER = H(0)
              DO K = 1, N
                RE_LOWER = RE_UPPER - DELZ(K)
                SINTH2 = RE_UPPER * SINTH1 / RE_LOWER
                STH2   = ASIN(SINTH2)
                PHI    = STH2 - STH1
                SINPHI = SIN(PHI)
                DIST = RE_UPPER * SINPHI / SINTH2
                CHAPMAN_FACTORS(N,K,IB) = DIST / DELZ(K)
                RE_UPPER = RE_LOWER
                SINTH1 = SINTH2
                STH1   = STH2
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!  finish beam loop

      ENDDO

!  end of routine

      RETURN
      END SUBROUTINE MSSEG_CHAPMANFACTORS

!

      SUBROUTINE MSSEG_PARTIC_CHECKINPUT ( &
        N_SZANGLES, N_USER_RELAZMS, N_USER_VZANGLES, &
        SZANGLES, USER_RELAZMS, USER_VZANGLES,       &
        DO_DIRECT_BEAM, STATUS, NMESSAGES, MESSAGES, ACTIONS ) ! output

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need

      USE MSSeg_PARS_m, Only : MAX_SZANGLES, MAX_USER_RELAZMS, MAX_USER_VZANGLES, MAX_MESSAGES, &
                                 ZERO, VQS_SUCCESS, VQS_SERIOUS

      IMPLICIT NONE

!  inputs

      INTEGER, INTENT (INOUT) ::          N_SZANGLES
      INTEGER, INTENT (INOUT) ::          N_USER_RELAZMS
      INTEGER, INTENT (IN) ::             N_USER_VZANGLES

      DOUBLE PRECISION, INTENT (IN) ::    SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::    USER_RELAZMS  ( MAX_USER_RELAZMS )
      DOUBLE PRECISION, INTENT (INOUT) :: USER_VZANGLES ( MAX_USER_VZANGLES )

!  Output

      LOGICAL, INTENT (OUT) ::            DO_DIRECT_BEAM
      INTEGER, INTENT (OUT) ::            STATUS
      INTEGER, INTENT (OUT) ::            NMESSAGES
      CHARACTER (LEN=*), INTENT (OUT) ::  MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=*), INTENT (OUT) ::  ACTIONS  ( 0:MAX_MESSAGES )

!  local variables
!  ---------------

      INTEGER ::           I, NM
      CHARACTER (LEN=2) :: C2
      LOGICAL ::           LOOP

!  Initialize output status

      STATUS = VQS_SUCCESS
      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Check of VQS Basic Input'
      ACTIONS(0)      = 'No Action required for this Task'

      NM     = NMESSAGES
      DO_DIRECT_BEAM = .true.

!  All numbers are now checked against maximum dimensions
!  ------------------------------------------------------

!  check dimensioning for SZA/AZM/VZA/LEVELS

      IF ( N_SZANGLES .GT. MAX_SZANGLES ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
               'Entry under "Number of solar zenith angles" >'// &
               ' allowed Maximum dimension'
        ACTIONS(NM)  = &
            'Re-set input value or increase MAX_SZANGLES dimension '// &
            'in MSSEG_PARS'
        STATUS = VQS_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

      IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
             'Entry under "Number of ..... azimuth angles" >'// &
             ' allowed Maximum dimension'
        ACTIONS(NM)  = &
          'Re-set input value or increase MAX_USER_RELAZMS dimension '// &
          'in MSSEG_PARS'
        STATUS       = VQS_SERIOUS
        NMESSAGES    = NM
        RETURN
      ENDIF

      IF ( N_USER_VZANGLES .GT. MAX_USER_VZANGLES ) THEN
          NM = NM + 1
          MESSAGES(NM) = &
             'Entry under "Number of .....viewing zenith angles" >'// &
             ' allowed Maximum dimension'
          ACTIONS(NM)  = &
          'Re-set input value or increase MAX_USER_VZANGLES dimension '// &
          'in MSSEG_PARS'
          STATUS = VQS_SERIOUS
          NMESSAGES = NM
          RETURN
      ENDIF

!  check viewing geometry input
!  ============================

!  Check solar zenith angle input

      DO I = 1, N_SZANGLES
        IF ( SZANGLES(I) .LT. ZERO .OR. SZANGLES(I).GE.90.0D0 ) THEN
          WRITE(C2,'(I2)')I
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: out-of-range solar angle, no. '//C2
          ACTIONS(NM)  = 'Look at SZANGLES input, should be < 90 & > 0'
          STATUS = VQS_SERIOUS
        ENDIF
      ENDDO

!  Check relative azimuths

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.N_USER_RELAZMS)
        I = I + 1
        IF ( USER_RELAZMS(I) .GT. 360.0D0   .OR. &
             USER_RELAZMS(I) .LT. ZERO ) THEN
          WRITE(C2,'(I2)')I
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: out-of-range azimuth angle, no. '//C2
          ACTIONS(NM)  = 'Look at azimuth angle input, should be in [0,360]'
          LOOP = .FALSE.
          STATUS = VQS_SERIOUS
        ENDIF
      ENDDO

!  check user-defined stream angles (should always be [0,90])

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.N_USER_VZANGLES)
         I = I + 1
         IF ( USER_VZANGLES(I) .GT. 90.0   .OR. &
              USER_VZANGLES(I) .LT. ZERO ) THEN
            WRITE(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: out-of-range user stream, no. '//C2
            ACTIONS(NM)  = 'Look at user viewing zenith angle input'
            LOOP = .FALSE.
            STATUS = VQS_SERIOUS
         ENDIF
      ENDDO

!mick fix

      NMESSAGES = NM

!  Finish

      RETURN
      END SUBROUTINE MSSEG_PARTIC_CHECKINPUT


      SUBROUTINE MSSEG_PARTIC_DERIVEINPUT ( &
        NSTOKES, N_SZANGLES, SZANGLES, N_USER_VZANGLES, USER_VZANGLES,  &
        NBEAMS, COS_SZANGLES, DO_NO_AZIMUTH, USER_STREAMS )

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need

      USE MSSeg_PARS_m, Only : MAX_SZANGLES, MAX_USER_VZANGLES, ZERO, DEG_TO_RAD

      IMPLICIT NONE

!  Inputs

      INTEGER, INTENT (IN) ::             NSTOKES
      INTEGER, INTENT (IN) ::             N_SZANGLES
      DOUBLE PRECISION, INTENT (IN) ::    SZANGLES ( MAX_SZANGLES )

      INTEGER, INTENT (IN) ::             N_USER_VZANGLES
      DOUBLE PRECISION, INTENT (IN) ::    USER_VZANGLES ( MAX_USER_VZANGLES )

!  Modified inputs

      LOGICAL, INTENT (INOUT) ::          DO_NO_AZIMUTH

!  Pure output

      INTEGER, INTENT (OUT) ::            NBEAMS
      DOUBLE PRECISION, INTENT (OUT) ::   COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (OUT) ::   USER_STREAMS ( MAX_USER_VZANGLES )

!  local variables
!  ---------------

      INTEGER ::          I
      LOGICAL ::          LOCAL_NADIR_ONLY

!  New section. Setting the DO_NO_AZIMUTH flag
!    Rt Solutions. 17 January 2006. R. Spurr and V. Natraj.

!  DO_NO_AZIMUTH should be set internally only when:
!     (a) NSTOKES = 1 and nadir view only

      LOCAL_NADIR_ONLY = .FALSE.
      IF ( N_USER_VZANGLES.EQ.1.and.USER_VZANGLES(1).EQ.ZERO ) THEN
        LOCAL_NADIR_ONLY = .TRUE.
      ENDIF
      DO_NO_AZIMUTH = .FALSE.
      IF ( NSTOKES.EQ.1.AND.LOCAL_NADIR_ONLY ) THEN
        DO_NO_AZIMUTH = .TRUE.
      ENDIF

!  Beam cosines

      NBEAMS = N_SZANGLES
      DO I = 1, N_SZANGLES
        COS_SZANGLES(I)  = DCOS ( SZANGLES(I) * DEG_TO_RAD )
      ENDDO

!  User stream cosines and secants

      DO I = 1, N_USER_VZANGLES
         USER_STREAMS(I) = DCOS(DEG_TO_RAD*USER_VZANGLES(I))
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE MSSEG_PARTIC_DERIVEINPUT

      END MODULE msseg_partic_miscsetups_m

