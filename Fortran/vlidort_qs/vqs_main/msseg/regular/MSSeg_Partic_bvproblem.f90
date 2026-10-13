
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

!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- ASYMTX_TOLERANCE now introduced and passed into the QHOM call to ASYMTX
!   -- Flux output (LOCAL_FLUXES) controlled by DO_LOCAL_FLUXES input flag
!   -- DO_FOURIER0_NSTOKES2 performance flag introduced.
!   -- Call subroutine specially written to get the discrete-ordinate field for Fluxes

!  05 April 2023. Revision

      MODULE MSSeg_Partic_BVProblem_m

      PRIVATE :: BVP_BACKSUB, BVP_SURFACE_SETUP_BEAM, BVP_COLUMN_SETUP
      PUBLIC  :: MSSEG_BVPSOLUTION_MASTER

      CONTAINS


      SUBROUTINE MSSEG_BVPSOLUTION_MASTER ( &
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, DO_LAMBERTIAN_SURFACE,       &
        IBEAM, FOURIER, NCHOICE, NSTOKES, NSTREAMS, NLAYERS, SURFACE_FACTOR,    &
        NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        LAMBERTIAN_ALBEDO, QUAD_STRMWTS, MUELLER_INDEX, DIRECT_BEAM,            &
        HSOL, PSOL, HBVP, PBVP, STATUS, MESSAGE, TRACE )

!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- Add NCHOICE (debug only)

      USE MSSeg_PARS_m, only : MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAXLAYERS, MAX_SZANGLES,           &
                                 MAXSTOKES_SQ, MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, MAXSTRMSTKS_2, &
                                 MAXBANDTOTAL, MAXTOTAL, VQS_SUCCESS, VQS_SERIOUS, ZERO, ONE

      use MSSeg_HPTypes_Def_m

      IMPLICIT NONE

!  INPUTS
!  ======

!  Flags

      LOGICAL, INTENT (IN) ::            DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::            DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::            DO_LAMBERTIAN_SURFACE

!  Numbers

      INTEGER, INTENT (IN) ::            FOURIER
      INTEGER, INTENT (IN) ::            IBEAM
      INTEGER, INTENT (IN) ::            NCHOICE
      INTEGER, INTENT (IN) ::            NSTOKES
      INTEGER, INTENT (IN) ::            NSTREAMS
      INTEGER, INTENT (IN) ::            NLAYERS

!  Bookkeeping

      INTEGER, INTENT (IN) ::            NSTREAMS_2
      INTEGER, INTENT (IN) ::            NTOTAL
      INTEGER, INTENT (IN) ::            NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::            NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::            N_SUBDIAG
      INTEGER, INTENT (IN) ::            N_SUPDIAG

!  surface inputs

      DOUBLE PRECISION, INTENT (IN) ::   SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::   LAMBERTIAN_ALBEDO

      DOUBLE PRECISION, INTENT (IN) ::   QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::            MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

!  direct beam

      DOUBLE PRECISION, INTENT (IN) ::   DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )

!  input solution type structures

      type(Homog_Solutions) , intent(in) :: HSOL
      type(Partic_Solutions), intent(in) :: PSOL

!  input BVProblem type structures

      type(Homog_BVPSetup)  , intent(in) :: HBVP

!  output BVProblem solution type structure

      type(Partic_BVPSolve)  , intent(inout) :: PBVP

!  Exception handling

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

      INTEGER ::          STATUS_SUB

!  Intialize Exception handling

      STATUS  = VQS_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  --Additional setups for the albedo layer

      CALL BVP_SURFACE_SETUP_BEAM ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, FOURIER,          &
        SURFACE_FACTOR, NSTOKES, NSTREAMS, NLAYERS, QUAD_STRMWTS,    &
        LAMBERTIAN_ALBEDO, HBVP%AXBID_F, MUELLER_INDEX, PSOL%WLOWER, &
        PBVP%R2_BEAM )

!  --set up Column for solution vector (the "B" as in AX=B)
!   -- 8/3/21. Mk2. Introduce FOURIER and NCHOICE, useful debugs if needed

      CALL BVP_COLUMN_SETUP ( &
        DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFACE,           &
        FOURIER, IBEAM, NCHOICE, NSTOKES, NSTREAMS, NLAYERS, &
        NSTREAMS_2, NTOTAL, NSTKS_NSTRMS, NSTKS_NSTRMS_2,    &
        PBVP%R2_BEAM, DIRECT_BEAM, PSOL%WUPPER, PSOL%WLOWER, &
        PBVP%COL2, PBVP%SCOL2 )

!  --Solve the boundary problem for this Fourier component (back substitution)

      CALL BVP_BACKSUB ( &
        FOURIER, IBEAM, NLAYERS, NTOTAL,                     &
        N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2,  &
        HSOL%K_REAL, HSOL%K_COMPLEX, HBVP%BANDMAT2,          &
        HBVP%IPIVOT, HBVP%SMAT2, HBVP%SIPIVOT,               &
        PBVP%COL2, PBVP%SCOL2, PBVP%LCON, PBVP%MCON,         &
        STATUS_SUB, MESSAGE, TRACE )

!  error tracing

      IF ( STATUS_SUB .NE. VQS_SUCCESS ) THEN
        STATUS = VQS_SERIOUS
        RETURN
      ENDIF

!  finish

      RETURN
      END SUBROUTINE MSSEG_BVPSOLUTION_MASTER

!

      SUBROUTINE BVP_SURFACE_SETUP_BEAM ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, FOURIER_COMPONENT, &
        SURFACE_FACTOR, NSTOKES, NSTREAMS, NLAYERS, QUAD_STRMWTS,     &
        LAMBERTIAN_ALBEDO, AXBID_F, MUELLER_INDEX, WLOWER,            &
        R2_BEAM )

!  This is the Lambertian  surface routine
!     Reflected beam solutions

!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension, Only use quantities you need

      USE MSSeg_PARS_m, only : MAXSTOKES, MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS, MAXSTOKES_SQ, ZERO

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE

      INTEGER, INTENT (IN) ::           FOURIER_COMPONENT
      DOUBLE PRECISION, INTENT (IN) ::  SURFACE_FACTOR

      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS

      DOUBLE PRECISION, INTENT (IN) ::  LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  AXBID_F ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: R2_BEAM ( MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      DOUBLE PRECISION :: REFL_B, REFL_B_S, FACTOR_KERNEL, AXJ
      INTEGER          :: I, J, O1, O2, OM, NL

!  Initialization
!  ==============

!  Zero total reflected contributions

      DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES
          R2_BEAM(I,O1)    = ZERO
        ENDDO
      ENDDO

!  Return with Zeroed values if albedo flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Return if not the Fourier m = 0 component (Lambertian surface)

      IF ( DO_LAMBERTIAN_SURFACE .and. &
             FOURIER_COMPONENT .GT. 0 ) RETURN

!  Skip if non-Lambertian surface

      if ( .not. DO_LAMBERTIAN_SURFACE ) GOTO 878

!  Lambertian surface
!  ==================

!  Integrate Downward streams of particular solutions

      REFL_B = ZERO
      DO J = 1, NSTREAMS
        AXJ = QUAD_STRMWTS(J)
        REFL_B = REFL_B + AXJ*WLOWER(J,1,NLAYERS)
      ENDDO

!  set reflectance

      REFL_B = REFL_B * SURFACE_FACTOR * LAMBERTIAN_ALBEDO
      DO I = 1, NSTREAMS
        R2_BEAM(I,1) = REFL_B
      ENDDO

!  debug

!      DO I = 1, NSTREAMS
!        DO O1 = 1, NSTOKES
!          write(*,'(2i3,1p2e14.6)')I,O1,R2_BEAM(I,O1)
!        ENDDO
!      ENDDO
!      pause'r2beam'

!  Return for Lambertian case

      RETURN

!  BRDF surface
!  ============

!  Continuation point

878   CONTINUE

!  Last layer

      NL = NLAYERS

!  Integrate Downward streams of particular solutions

      FACTOR_KERNEL = SURFACE_FACTOR
      DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES
          REFL_B = ZERO
          DO J = 1, NSTREAMS
            REFL_B_S = ZERO
            DO O2 = 1, NSTOKES
              OM = MUELLER_INDEX(O1,O2)
              REFL_B_S = REFL_B_S + &
                       AXBID_F(I,J,OM) * WLOWER(J,O2,NL)
            ENDDO
            REFL_B = REFL_B + REFL_B_S
          ENDDO
         R2_BEAM(I,O1) = FACTOR_KERNEL * REFL_B
        ENDDO
      ENDDO

!  debug

!      DO I = 1, NSTREAMS
!        DO O1 = 1, NSTOKES
!          write(*,'(2i3,1p2e14.6)')I,O1,R2_BEAM(I,O1)
!        ENDDO
!      ENDDO
!      pause'r2beam'

!  Finish

      RETURN
      END SUBROUTINE BVP_SURFACE_SETUP_BEAM

!

      SUBROUTINE BVP_COLUMN_SETUP ( &
        DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFACE,           &
        FOURIER, IBEAM, NCHOICE, NSTOKES, NSTREAMS, NLAYERS, &
        NSTREAMS_2, NTOTAL, NSTKS_NSTRMS, NSTKS_NSTRMS_2,    &
        R2_BEAM, DIRECT_BEAM, WUPPER, WLOWER,                &
        COL2, SCOL2 )

!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- introduce FOURIER and NCHOICE, useful debugs if needed

      USE MSSeg_PARS_m, only : MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAXBEAMS, &
                                 MAXSTREAMS_2, MAXSTRMSTKS_2, MAXTOTAL, ZERO

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE

      INTEGER, INTENT (IN) ::           IBEAM, FOURIER, NCHOICE
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS

      INTEGER, INTENT (IN) ::           NSTREAMS_2
      INTEGER, INTENT (IN) ::           NTOTAL
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS_2

      DOUBLE PRECISION, INTENT (IN) ::  R2_BEAM ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: COL2 ( MAXTOTAL, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: SCOL2 ( MAXSTRMSTKS_2, MAXBEAMS )

!  Local variables
!  ---------------

      INTEGER          :: I,I1,LAY,LAY1,C0,CM,O1,IR,IROW

!  If Nlayers = 1, go to special case

      IF ( NLAYERS .EQ. 1 ) GO TO 345

!  zero column vector

      DO I = 1, NTOTAL
        COL2(I,IBEAM) = ZERO
      ENDDO

!  Upper boundary for layer 1: no downward diffuse radiation
!  ---------------------------------------------------------

      LAY = 1
      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          COL2(IROW,IBEAM)   = - WUPPER(I,O1,LAY)
        ENDDO
      ENDDO

!  intermediate layer boundaries (will not be done if NLAYERS = 1 )
!  -----------------------------

      DO LAY = 2, NLAYERS

        LAY1 = LAY - 1
        C0 = LAY1*NSTKS_NSTRMS_2 - NSTKS_NSTRMS
        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            COL2(CM,IBEAM) = WUPPER(I,O1,LAY) - WLOWER(I,O1,LAY1)
          ENDDO
        ENDDO

      ENDDO

!  lowest (surface) boundary with albedo (diffuse radiation terms only)
!  -------------------------------------

      LAY = NLAYERS
      C0 = (LAY-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS

!  with non-zero albedo, include integrated downward reflectances

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            COL2(CM,IBEAM) = - WLOWER(I1,O1,LAY) + R2_BEAM(I,O1)
          ENDDO
        ENDDO

!  no albedo, similar code excluding integrated reflectance

      ELSE

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            COL2(CM,IBEAM) = - WLOWER(I1,O1,LAY)
          ENDDO
        ENDDO

      ENDIF

!  Add direct beam solution (only to final level)
!  ----------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              COL2(CM,IBEAM) = COL2(CM,IBEAM) + DIRECT_BEAM(I,IBEAM,O1)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

! debug
!if ( fourier.eq.2 .and.nchoice.eq.75 ) then
!      do i = 1, ntotal
!        write(100+IBEAM,*)i,COL2(i,IBEAM)
!      enddo
!     stop 'hello'
!endif

!  normal completion

      RETURN

!  special case

345   CONTINUE

!  zero column vector

      DO I = 1, NTOTAL
        SCOL2(I,IBEAM) = ZERO
      ENDDO

!  Upper boundary for layer 1: no downward diffuse radiation
!  ---------------------------------------------------------

      LAY = 1
      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          SCOL2(IROW,IBEAM)   = - WUPPER(I,O1,LAY)
        ENDDO
      ENDDO

!  lowest (surface) boundary with albedo (diffuse radiation terms only)
!  -------------------------------------------------------------------

!  with non-zero albedo, include integrated downward reflectances
!  no albedo, similar code excluding integrated reflectance

      C0 = NSTKS_NSTRMS
      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            SCOL2(CM,IBEAM) = - WLOWER(I1,O1,LAY) + R2_BEAM(I,O1)
          ENDDO
        ENDDO
      ELSE
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            SCOL2(CM,IBEAM) = - WLOWER(I1,O1,LAY)
          ENDDO
        ENDDO
      ENDIF

!  Add direct beam solution (only to final level)
!  ----------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              SCOL2(CM,IBEAM) = SCOL2(CM,IBEAM)+DIRECT_BEAM(I,IBEAM,O1)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  finish

      RETURN
      END SUBROUTINE BVP_COLUMN_SETUP

!

      SUBROUTINE BVP_BACKSUB ( &
        FOURIER, IBEAM, NLAYERS, NTOTAL,                     &
        N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2,  &
        K_REAL, K_COMPLEX, BANDMAT2, IPIVOT, SMAT2, SIPIVOT, &
        COL2, SCOL2, LCON, MCON,                             &
        STATUS, MESSAGE, TRACE )

!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension, Only use quantities you need

      USE MSSeg_PARS_m, only : MAXLAYERS, MAXBEAMS, MAXSTRMSTKS, MAXSTRMSTKS_2, &
                                 MAXBANDTOTAL, MAXTOTAL, VQS_SUCCESS, VQS_SERIOUS
      USE LAPACK_TOOLS_m, only : DGBTRS, DGETRS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::             FOURIER
      INTEGER, INTENT (IN) ::             IBEAM
      INTEGER, INTENT (IN) ::             NLAYERS
      INTEGER, INTENT (IN) ::             NTOTAL

      INTEGER, INTENT (IN) ::             N_SUBDIAG
      INTEGER, INTENT (IN) ::             N_SUPDIAG
      INTEGER, INTENT (IN) ::             NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::             NSTKS_NSTRMS_2

      INTEGER, INTENT (IN) ::             K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::             K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::    BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::             IPIVOT ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) ::    SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::             SIPIVOT ( MAXSTRMSTKS_2 )

      DOUBLE PRECISION, INTENT (INOUT) :: COL2 ( MAXTOTAL, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: SCOL2 ( MAXSTRMSTKS_2, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

      INTEGER, INTENT (OUT) ::            STATUS
      CHARACTER (LEN=*), INTENT (INOUT) ::  MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) ::  TRACE

!  local variables
!  ---------------

      INTEGER ::           C0, LAY, K, KO1, K0, K1, K2
      INTEGER ::           IROW, IROW1, IROW_S, IROW1_S
      INTEGER ::           INFO, M
      CHARACTER (LEN=3) :: CI
      CHARACTER (LEN=2) :: CB

!  mick fix 6/23/2021. Local arrays defined
      DOUBLE PRECISION  :: LOC_COL2  ( MAXTOTAL, 1 )
      DOUBLE PRECISION  :: LOC_SCOL2 ( MAXSTRMSTKS_2, 1 )

!  Intialize Exception handling

      STATUS  = VQS_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '
      M = FOURIER

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2
!mick fix 6/23/2021 - pass RHS b vector of Ax = b ("COL2") to a local array with
!                     a 2nd dim of size one to use the LAPACK routine DGBTRS
!                     again as intended

        LOC_COL2(1:NTOTAL,1) = COL2(1:NTOTAL,IBEAM)
        CALL DGBTRS &
           ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, &
              LOC_COL2, MAXTOTAL, INFO )
        COL2(1:NTOTAL,IBEAM) = LOC_COL2(1:NTOTAL,1)

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CB, '(I2)' ) IBEAM
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call in BVP_BACKSUB, Beam # '//CB
          STATUS  = VQS_SERIOUS
          RETURN
        ENDIF

!  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

        DO LAY = 1, NLAYERS
          C0 = (LAY-1)*NSTKS_NSTRMS_2
          KO1 = K_REAL(LAY) + 1
          DO K = 1, K_REAL(LAY)
            IROW = K
            IROW1 = IROW + NSTKS_NSTRMS
            LCON(K,LAY) = COL2(C0+IROW,IBEAM)
            MCON(K,LAY) = COL2(C0+IROW1,IBEAM)
          ENDDO
          DO K = 1, K_COMPLEX(LAY)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(LAY)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(LAY)
            IROW1_S = IROW_S + NSTKS_NSTRMS
            LCON(K1,LAY) = COL2(C0+IROW,   IBEAM)
            LCON(K2,LAY) = COL2(C0+IROW_S, IBEAM)
            MCON(K1,LAY) = COL2(C0+IROW1,  IBEAM)
            MCON(K2,LAY) = COL2(C0+IROW1_S,IBEAM)
          ENDDO
        ENDDO

!  debug-----------------------------------------------

!        if ( do_debug_write) then
!         if (fourier.eq.1 ) then
!         DO LAY = 1, NLAYERS
!          KO1 = K_REAL(LAY) + 1
!          DO K = 1, K_REAL(LAY)
!           write(50,'(4i3,1p2e20.10)')FOURIER,IBEAM,K,LAY,
!     &                LCON(K,LAY), MCON(K,LAY)
!          ENDDO
!          DO K = 1, K_COMPLEX(LAY)
!           K0 = 2*K - 2
!           K1 = KO1 + K0
!           K2 = K1  + 1
!           write(50,'(4i3,1p4e20.10)')FOURIER,IBEAM,K,LAY,
!     &                LCON(K1,LAY), MCON(K1,LAY),
!     &                LCON(K2,LAY), MCON(K2,LAY)
!          ENDDO
!         ENDDO
!        endif
!        ENDIF

!        DO LAY = 20, NLAYERS
!          DO K = 1, K_REAL(LAY)
!             write(97,'(3i5,1p2e15.7)')FOURIER,LAY,K,
!     &                LCON(K,LAY),MCON(K,LAY)
!            ENDDO
!           ENDDO
!          IF ( FOURIER.EQ.0)STOP
!        write(*,*)'reg',fourier,ibeam,lcon(1,7)

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2

!mick fix 6/23/2021 - pass RHS b vector of Ax = b ("SCOL2") to a local array with
!                     a 2nd dim of size one to use the LAPACK routine DGETRS
!                     again as intended
        LOC_SCOL2(1:NSTKS_NSTRMS_2,1) = SCOL2(1:NSTKS_NSTRMS_2,IBEAM)
        CALL DGETRS &
           ( 'N', MAXSTRMSTKS_2, 1, SMAT2, MAXSTRMSTKS_2, SIPIVOT, &
              LOC_SCOL2, MAXSTRMSTKS_2, INFO )
        SCOL2(1:NSTKS_NSTRMS_2,IBEAM) = LOC_SCOL2(1:NSTKS_NSTRMS_2,1)

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CB, '(I2)' ) IBEAM
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRS call (nlayers=1)in BVP_BACKSUB, Beam # '//CB
          STATUS  = VQS_SERIOUS
          RETURN
        ENDIF

!  Set integration constants LCON and MCON for -/+ eigensolutions, all l

        LAY = 1
        KO1 = K_REAL(LAY) + 1
        DO K = 1, K_REAL(LAY)
          IROW = K
          IROW1 = IROW + NSTKS_NSTRMS
          LCON(K,LAY) = SCOL2(IROW,IBEAM)
          MCON(K,LAY) = SCOL2(IROW1,IBEAM)
        ENDDO
        DO K = 1, K_COMPLEX(LAY)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW = K + K_REAL(LAY)
          IROW1 = IROW + NSTKS_NSTRMS
          IROW_S = K + K_REAL(LAY) + K_COMPLEX(LAY)
          IROW1_S = IROW_S + NSTKS_NSTRMS
          LCON(K1,LAY) = SCOL2(IROW,    IBEAM)
          LCON(K2,LAY) = SCOL2(IROW_S,  IBEAM)
          MCON(K1,LAY) = SCOL2(IROW1,   IBEAM)
          MCON(K2,LAY) = SCOL2(IROW1_S, IBEAM)
        ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BVP_BACKSUB


      END MODULE MSSeg_Partic_BVProblem_m

