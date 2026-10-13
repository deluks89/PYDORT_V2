
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

      MODULE Msseg_Homog_Miscsetups_m

!  05 April 2023. Revision

      PUBLIC  :: MSSEG_HOMOG_MISCSETUPS
      PRIVATE :: MSSEG_HOMOG_DERIVEINPUT, MSSEG_HOMOG_CHECKINPUT, MSSEG_HOMOG_DELTAMSCALE

      CONTAINS

      SUBROUTINE MSSEG_HOMOG_MISCSETUPS ( & 
        DO_RAYLEIGH_ONLY, DO_DELTAM_SCALING,                             & ! Input
        NSTOKES, NLAYERS, NSTREAMS, NGREEK_MOMENTS_INPUT, HEIGHT_GRID,   & ! Input
        DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,      & ! Input
        NMOMENTS, NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2,              & ! Output
        NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTOKES_SQ,                        & ! Output
        QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS,          & ! Output
        FLUXVEC, MUELLER_INDEX, DMAT, DFLUX, DO_REAL_EIGENSOLVER,        & ! Output
        DO_LAYER_SCATTERING, DELTAU_VERT, OMEGA_TOTAL, GREEKMAT_TOTAL,   & ! Output
        TAUGRID_INPUT, TAUGRID, TRUNC_FACTOR, FAC1, OMEGA_GREEK,         & ! Output
        STATUS, NMESSAGES, MESSAGES, ACTIONS )                             ! Output

!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- Use public routine GETQUAD2 instead of GAULEG

!  Mk3, 8/29/21. Stripped out beam stuff
!   Greekmats defined with the MAXMOMENTS dimension

      USE MSSEG_PARS_m, Only : MAXSTOKES, MAXSTREAMS, MAXLAYERS, VQS_SERIOUS, &
                               MAXMOMENTS, MAXSTOKES_SQ, MAX_MESSAGES


      IMPLICIT NONE

      LOGICAL, INTENT (INOUT) ::       DO_RAYLEIGH_ONLY
      LOGICAL, INTENT (INOUT) ::       DO_DELTAM_SCALING

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NSTREAMS

!    Mk3, 1/24/22. NGreek_moments_Input does not need to be checked and re-set. BUT IS AN INPUT HERE
!    Mk3, 1/24/22. Dimensioning is always with MAXMOMENTS for Greekmats.

      INTEGER, INTENT (IN) ::          NGREEK_MOMENTS_INPUT
      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

      DOUBLE PRECISION, INTENT (IN) :: HEIGHT_GRID ( 0:MAXLAYERS )

!  derived output (mostly deltam-scaled variables)
!    Mk3, 1/24/22. Dimensioning is always with MAXMOMENTS for Greekmats.

      DOUBLE PRECISION, INTENT (OUT) :: OMEGA_TOTAL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: GREEKMAT_TOTAL ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

      DOUBLE PRECISION, INTENT (OUT) :: TAUGRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: TRUNC_FACTOR ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: FAC1 ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: TAUGRID_INPUT ( 0:MAXLAYERS )

      DOUBLE PRECISION, INTENT (OUT) ::   QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (OUT) ::   QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (OUT) ::   QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (OUT) ::   QUAD_HALFWTS ( MAXSTREAMS )

      INTEGER, INTENT (OUT) ::            NMOMENTS, NSTOKES_SQ
      INTEGER, INTENT (OUT) ::            NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2
      INTEGER, INTENT (OUT) ::            NTOTAL, N_SUBDIAG, N_SUPDIAG

      INTEGER, INTENT (OUT) ::            MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) ::   DMAT ( MAXSTOKES, MAXSTOKES ), DFLUX ( MAXSTOKES)
      DOUBLE PRECISION, INTENT (OUT) ::   FLUXVEC ( MAXSTOKES )

      LOGICAL, INTENT (OUT) :: DO_REAL_EIGENSOLVER ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (OUT) :: DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  Exception handling output

      INTEGER, INTENT (OUT) ::            STATUS
      INTEGER, INTENT (OUT) ::            NMESSAGES
      CHARACTER (LEN=*), INTENT (OUT) ::  MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=*), INTENT (OUT) ::  ACTIONS  ( 0:MAX_MESSAGES )

!  Check input

      CALL MSSEG_HOMOG_CHECKINPUT ( &
        DO_RAYLEIGH_ONLY, DO_DELTAM_SCALING, NSTOKES, NSTREAMS, NLAYERS,  & ! Input
        HEIGHT_GRID, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,             & ! Input
        STATUS, NMESSAGES, MESSAGES, ACTIONS )                  ! output

      IF ( STATUS .eq. VQS_SERIOUS ) RETURN

!  derive input
!    Mk3, 1/24/22. NGreek_moments_Input does not need to be checked and re-set. BUT IS AN INPUT HERE
!    Dimensioning is always with MAXMOMENTS for Greekmats.

      CALL MSSEG_HOMOG_DERIVEINPUT ( &
        DO_RAYLEIGH_ONLY, NSTOKES, NSTREAMS, NLAYERS,          & ! Input
        NGREEK_MOMENTS_INPUT, GREEKMAT_TOTAL_INPUT,            & ! Input
        NMOMENTS, NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,                  & ! Output
        NSTOKES_SQ, NSTKS_NSTRMS, NSTKS_NSTRMS_2, FLUXVEC,                   & ! Output
        QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS,              & ! Output
        MUELLER_INDEX, DMAT, DFLUX, DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING ) ! Output

!  Delta-m scaling of input quantities

      CALL MSSEG_HOMOG_DELTAMSCALE ( &
        DO_DELTAM_SCALING, NSTOKES, NLAYERS, NMOMENTS, MUELLER_INDEX, &
        OMEGA_TOTAL_INPUT, DELTAU_VERT_INPUT, GREEKMAT_TOTAL_INPUT,   &
        TAUGRID_INPUT, OMEGA_TOTAL, DELTAU_VERT, GREEKMAT_TOTAL,      &
        TAUGRID, OMEGA_GREEK, TRUNC_FACTOR, FAC1 )


!  Finish

      RETURN
      END SUBROUTINE MSSEG_HOMOG_MISCSETUPS

!

      SUBROUTINE MSSEG_HOMOG_DELTAMSCALE ( &
        DO_DELTAM_SCALING, NSTOKES, NLAYERS, NMOMENTS, MUELLER_INDEX, &
        OMEGA_TOTAL_INPUT, DELTAU_VERT_INPUT, GREEKMAT_TOTAL_INPUT,   &
        TAUGRID_INPUT, OMEGA_TOTAL, DELTAU_VERT, GREEKMAT_TOTAL,      &
        TAUGRID, OMEGA_GREEK, TRUNC_FACTOR, FAC1 )

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need

      USE MSSeg_PARS_m, Only : MAXSTOKES, MAXLAYERS, MAXMOMENTS, MAXSTOKES_SQ, ZERO, ONE

      IMPLICIT NONE

!  Control inptuts

      LOGICAL, INTENT (IN) ::           DO_DELTAM_SCALING
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

!  optical inputs

      DOUBLE PRECISION, INTENT (IN)  :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN)  :: DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN)  :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )


!  outputs

      DOUBLE PRECISION, INTENT (OUT) :: OMEGA_TOTAL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: GREEKMAT_TOTAL ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

      DOUBLE PRECISION, INTENT (OUT) :: TAUGRID_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: TAUGRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: TRUNC_FACTOR ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: FAC1 ( MAXLAYERS )

!  local variables

      DOUBLE PRECISION :: FDEL, FAC2, DNL1, FDNL1, DNM1, SUMG
      INTEGER          :: K, K1, K2, N, N1, L, NM1, O1, O2, OM

!  Indexing

      INTEGER, DIMENSION(4) :: KTYPE1, KTYPE2

!mick - singularity buster output
!      INTEGER   :: I
      LOGICAL   :: SBUST(6)

!  Setup the two ktype arrays

      KTYPE1 = (/ 1, 6, 11, 16 /)
      KTYPE2 = (/ 2, 5, 12, 15 /)

! taugird input

      TAUGRID_INPUT(0) = ZERO
      DO N = 1, NLAYERS
         TAUGRID_INPUT(N) = TAUGRID_INPUT(N-1) + DELTAU_VERT_INPUT(N)
      ENDDO

!mick fix - initialise

      GREEKMAT_TOTAL = ZERO

!  DELTAM SCALING
!  ==============

      IF ( DO_DELTAM_SCALING ) THEN

!  New section by R. Spurr, RT Solutions Inc.
!   Based in part on the code in VDISORT.

        TAUGRID(0) = ZERO
        NM1  = NMOMENTS+1
        DNM1 = DBLE(2*NM1+1)

!  Scaling for layer input
!  -----------------------

        DO N = 1, NLAYERS

          N1 = N - 1

!  overall truncation factor

          FDEL = GREEKMAT_TOTAL_INPUT(NM1,N,1) / DNM1
          FAC2 = ONE - FDEL
          FAC1(N)         = ONE - FDEL * OMEGA_TOTAL_INPUT(N)
          TRUNC_FACTOR(N) = FDEL

!  Scale Greek Matrix entries

          DO L = 0, NMOMENTS
            DNL1  = DBLE(2*L + 1 )
            FDNL1 = FDEL * DNL1
            IF ( NSTOKES .GT. 1 ) THEN
              DO K = 1, 4
                K1 = KTYPE1(K)
                K2 = KTYPE2(K)
                GREEKMAT_TOTAL(L,N,K1) = ( GREEKMAT_TOTAL_INPUT(L,N,K1) - FDNL1 ) / FAC2
                GREEKMAT_TOTAL(L,N,K2) = GREEKMAT_TOTAL_INPUT(L,N,K2) / FAC2
              ENDDO
            ELSE
              GREEKMAT_TOTAL(L,N,1) = ( GREEKMAT_TOTAL_INPUT(L,N,1) - FDNL1 ) / FAC2
            ENDIF
          ENDDO

!  Maintain phase function normalization

          GREEKMAT_TOTAL(0,N,1) = ONE

!  scale optical depth grid and single scatter albedo

          DELTAU_VERT(N) = DELTAU_VERT_INPUT(N) * FAC1(N)
          OMEGA_TOTAL(N) = OMEGA_TOTAL_INPUT(N) * FAC2 / FAC1(N)
          TAUGRID(N)     = TAUGRID(N1) + DELTAU_VERT(N)

!  end layer loop

        ENDDO

!  NO DELTAM SCALING
!  =================

!  move input geophysical variables to Workspace quantities

      ELSE

        TAUGRID(0) = ZERO
        DO N = 1, NLAYERS
          FAC1(N)        = ONE
          OMEGA_TOTAL(N) = OMEGA_TOTAL_INPUT(N)
          TAUGRID(N)     = TAUGRID_INPUT(N)
          DELTAU_VERT(N) = DELTAU_VERT_INPUT(N)
          DO L = 0, NMOMENTS
            IF ( NSTOKES .GT. 1 ) THEN
              DO K = 1, 4
                K1 = KTYPE1(K)
                K2 = KTYPE2(K)
                GREEKMAT_TOTAL(L,N,K1) = GREEKMAT_TOTAL_INPUT(L,N,K1)
                GREEKMAT_TOTAL(L,N,K2) = GREEKMAT_TOTAL_INPUT(L,N,K2)
              ENDDO
            ELSE
              GREEKMAT_TOTAL(L,N,1) = GREEKMAT_TOTAL_INPUT(L,N,1)
            ENDIF
          ENDDO
        ENDDO

      ENDIF

!mick fix 2/13/2012 - singularity busters added

!  Note: If running a case close to optical property numerical limits,
!        delta-m scaling may modify omega and/or g in such a way as to make
!        them unphysical or introduce instability; therefore, we recheck
!        omega and g AFTER delta-m scaling and slightly adjust them if
!        necessary

      DO N = 1, NLAYERS
        SBUST = .false.

        !Singularity buster for single scatter albedo
        IF (OMEGA_TOTAL(N) > 0.999999999D0) THEN
          OMEGA_TOTAL(N) = 0.999999999D0
          SBUST(1) = .true.
        ELSE IF (OMEGA_TOTAL(N) < 1.0D-9) THEN
          OMEGA_TOTAL(N) = 1.0D-9
          SBUST(2) = .true.
        END IF

        !Singularity buster for asymmetry parameter
        !(1) Divide by 2L+1 where L = 1 to get the asym par
        GREEKMAT_TOTAL(1,N,1) = GREEKMAT_TOTAL(1,N,1)/3.0D0
        !(2) Modify the asym par if necessary
        IF (GREEKMAT_TOTAL(1,N,1) > 0.999999999D0) THEN
          GREEKMAT_TOTAL(1,N,1) = 0.999999999D0
          SBUST(3) = .true.
        ELSE IF (GREEKMAT_TOTAL(1,N,1) < -0.999999999D0) THEN
          GREEKMAT_TOTAL(1,N,1) = -0.999999999D0
          SBUST(4) = .true.
        ELSE IF ((GREEKMAT_TOTAL(1,N,1) >= 0.0D0) .AND. &
                 (GREEKMAT_TOTAL(1,N,1) < 1.0D-9)) THEN
          GREEKMAT_TOTAL(1,N,1) = 1.0D-9
          SBUST(5) = .true.
        ELSE IF ((GREEKMAT_TOTAL(1,N,1) < 0.0D0) .AND. &
                 (GREEKMAT_TOTAL(1,N,1) > -1.0D-9)) THEN
          GREEKMAT_TOTAL(1,N,1) = -1.0D-9
          SBUST(6) = .true.
        END IF
        !(3) Reconstruct the 1st-order phase func moment
        GREEKMAT_TOTAL(1,N,1) = 3.0D0*GREEKMAT_TOTAL(1,N,1)

        !WRITE(*,*)
        !WRITE(*,'(A,I2)') 'FOR LAYER: ',N
        !DO I=1,6
        !  WRITE(*,'(A,I1,A,L1)') '  SBUST(',I,') = ',SBUST(I)
        !ENDDO

      ENDDO
!READ(*,*)

!  phase matrix-weighted OMEGA

      DO N = 1, NLAYERS
        DO L = 0, NMOMENTS
          DO O1 = 1, NSTOKES
            DO O2 = 1, NSTOKES
              OM = MUELLER_INDEX(O1,O2)
              SUMG = OMEGA_TOTAL(N)*GREEKMAT_TOTAL(L,N,OM)
              OMEGA_GREEK(L,N,O1,O2)   = SUMG
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  Finish module

      RETURN
      END SUBROUTINE MSSEG_HOMOG_DELTAMSCALE

!

      SUBROUTINE MSSEG_HOMOG_CHECKINPUT ( &
        DO_RAYLEIGH_ONLY, DO_DELTAM_SCALING, NSTOKES, NSTREAMS, NLAYERS,  & ! Input
        HEIGHT_GRID, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,             & ! Input
        STATUS, NMESSAGES, MESSAGES, ACTIONS )                             ! output

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need

!  Mk3, 1/24/22.
!    NGreek_moments_Input does not need to be checked and re-set. We are MS code only now.
!    Dimensioning is always with MAXMOMENTS for Greekmats.

      USE MSSeg_PARS_m, Only :  MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAXMOMENTS, MAXSTOKES_SQ, MAX_MESSAGES,&
                                ZERO, ONE, VQS_WARNING, VQS_SUCCESS, VQS_SERIOUS, OMEGA_SMALLNUM

      IMPLICIT NONE

!  input

      LOGICAL, INTENT (INOUT) ::          DO_RAYLEIGH_ONLY
      LOGICAL, INTENT (INOUT) ::          DO_DELTAM_SCALING


!     fixed input

      INTEGER, INTENT (IN) ::             NSTOKES
      INTEGER, INTENT (IN) ::             NSTREAMS
      INTEGER, INTENT (IN) ::             NLAYERS

      DOUBLE PRECISION, INTENT (IN) ::    HEIGHT_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::    OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::    GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )


!  Output

      INTEGER, INTENT (OUT) ::            STATUS
      INTEGER, INTENT (OUT) ::            NMESSAGES
      CHARACTER (LEN=*), INTENT (OUT) ::  MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=*), INTENT (OUT) ::  ACTIONS ( 0:MAX_MESSAGES )

!  local variables
!  ---------------

      INTEGER ::           I, L, NM
      CHARACTER (LEN=3) :: C3
      CHARACTER (LEN=2) :: C2
      LOGICAL ::           LOOP

!  Initialize output status

      STATUS = VQS_SUCCESS
      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Check of MSSeg asic HOMOG Input'
      ACTIONS(0)      = 'No Action required for this Task'

      NM     = NMESSAGES

!  All numbers are now checked against maximum dimensions
!  ------------------------------------------------------

!  New Exception handling code, 13 OCTOBER 2010

      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
               'Entry under "Number of half-space streams" >'// ' allowed Maximum dimension'
        ACTIONS(NM)  = &
             'Re-set input value or increase MAXSTREAMS dimension '// 'in MSSEG_PARS'
        STATUS = VQS_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

      IF ( NLAYERS .GT. MAXLAYERS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
               'Entry under "Number of atmospheric layers" >'// ' allowed Maximum dimension'
        ACTIONS(NM)  = &
             'Re-set input value or increase MAXLAYERS dimension '// 'in MSSEG_PARS'
        STATUS = VQS_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

      IF ( NSTOKES .GT. MAXSTOKES ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
         'Entry under'//' "Number of Stokes parameters" >'// ' allowed Maximum dimension'
        ACTIONS(NM)  = 'Re-set input value to 4 or less'
        STATUS = VQS_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

!  check inputs (both file-read and derived)
!  -----------------------------------------

!  no Delta-M scaling with Rayleigh only
!   ---Warning. Turn off delta-M scaling.

      IF ( DO_RAYLEIGH_ONLY ) THEN
        IF ( DO_DELTAM_SCALING ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: No delta-M scaling with Rayleigh-only'
          ACTIONS(NM)  = 'Warning: DO_DELTAM_SCALING turned off internally'
          STATUS = VQS_WARNING
          DO_DELTAM_SCALING = .FALSE.
        ENDIF
      ENDIF

!  check number of input expansion coefficient moments (non-Rayleigh)

!  Mk3, 1/24/22.  NGreek_moments_Input does not need to be checked and re-set. 
!      IF ( .NOT.DO_RAYLEIGH_ONLY  ) THEN
!        IF ( DO_DELTAM_SCALING ) THEN
!          IF ( NGREEK_MOMENTS_INPUT.LT.2*NSTREAMS ) THEN
!            NM = NM + 1
!            MESSAGES(NM) = 'Bad input: Fewer than 2N expansion moments with delta-M'
!            ACTIONS(NM)  = 'Warning: Re-set NGREEK_MOMENTS_INPUT to 2N internally'
!            STATUS = VQS_WARNING
!            NGREEK_MOMENTS_INPUT = 2*NSTREAMS
!          ENDIF
!        ENDIF
!      ELSE IF ( DO_RAYLEIGH_ONLY ) THEN
!        IF ( NGREEK_MOMENTS_INPUT.NE.2 ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = 'Bad input: Rayleigh-only, expansion momemts NOT = 2'
!          ACTIONS(NM)  = 'Warning: Set NGREEK_MOMENTS_INPUT = 2 internally'
!          STATUS = VQS_WARNING
!          NGREEK_MOMENTS_INPUT = 2
!        ENDIF
!      ENDIF

!  Check height grid input (Chapman function only)

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.NLAYERS)
         I = I + 1
         IF ( HEIGHT_GRID(I-1).LE.HEIGHT_GRID(I) ) THEN
            WRITE(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM)   = 'Bad input: Height-grid not monotonically decreasing; Layer '//C2
            ACTIONS(NM) = 'Look at Height-grid input'
            LOOP = .FALSE.
            STATUS = VQS_SERIOUS
         ENDIF
      ENDDO

!  check single scatter albedos (too large)

     DO L = 1, NLAYERS
        IF ( OMEGA_TOTAL_INPUT(L).GT.ONE-OMEGA_SMALLNUM ) THEN
          NM = NM + 1
          WRITE(C3,'(I3)')L
          MESSAGES(NM) = 'Bad input: SS-albedo too close to 1, layer '//C3
          ACTIONS(NM)  = 'Check SS-albedo input'
          STATUS = VQS_SERIOUS
        ELSE IF ( OMEGA_TOTAL_INPUT(L).LT.OMEGA_SMALLNUM ) THEN
          NM = NM + 1
          WRITE(C3,'(I3)')L
          MESSAGES(NM) = 'Bad input: SS-albedo too close to 0, layer '//C3
          ACTIONS(NM)  = 'Check SS-albedo input'
          STATUS = VQS_SERIOUS
        ENDIF
     ENDDO

!  check single scatter albedos (too small)

     DO L = 1, NLAYERS
         IF ( OMEGA_TOTAL_INPUT(L).LT.OMEGA_SMALLNUM ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: SS-albedo too close to 0, layer '//C3
          ACTIONS(NM)  = 'Check SS-albedo input'
          STATUS = VQS_SERIOUS
         ENDIF
     ENDDO

!  Check first phase function moments (normalization)

     DO L = 1, NLAYERS
        IF ( GREEKMAT_TOTAL_INPUT(0,L,1).NE.ONE ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGES(NM) = 'First phase moment (GREEK_11) not 1 for layer '//C3
          ACTIONS(NM)  = 'Check First phase function moment'
          STATUS = VQS_SERIOUS
        ENDIF
     ENDDO

!mick fix

      NMESSAGES = NM

!  Finish

      RETURN
      END SUBROUTINE MSSEG_HOMOG_CHECKINPUT


      SUBROUTINE MSSEG_HOMOG_DERIVEINPUT ( &
        DO_RAYLEIGH_ONLY, NSTOKES, NSTREAMS, NLAYERS, & ! Input
        NGREEK_MOMENTS_INPUT, GREEKMAT_TOTAL_INPUT,   & ! Input
        NMOMENTS, NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,     & ! Output
        NSTOKES_SQ, NSTKS_NSTRMS, NSTKS_NSTRMS_2, FLUXVEC,      & ! Output
        QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS, & ! Output
        MUELLER_INDEX, DMAT, DFLUX, DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING ) ! Output

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- Update use of auxiliary routine GETQUAD2

!  Mk3, 1/24/22.
!    NGreek_moments_Input does not need to be checked and re-set. BUT IS AN INPUT HERE
!    Dimensioning is always with MAXMOMENTS for Greekmats.

      USE MSSeg_PARS_m, Only : MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXLAYERS,  &
                                 MAX_PSOLS, MAXSTOKES_SQ, ZERO, ONE, HALF

      USE MSSeg_AUX_m, Only : GETQUAD2

      IMPLICIT NONE

!  Inputs

      LOGICAL, INTENT (IN) ::             DO_RAYLEIGH_ONLY
      INTEGER, INTENT (IN) ::             NSTOKES
      INTEGER, INTENT (IN) ::             NSTREAMS
      INTEGER, INTENT (IN) ::             NLAYERS
      INTEGER, INTENT (IN) ::             NGREEK_MOMENTS_INPUT

      DOUBLE PRECISION, INTENT (IN) ::    GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

!  Pure output

      DOUBLE PRECISION, INTENT (OUT) ::   QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (OUT) ::   QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (OUT) ::   QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (OUT) ::   QUAD_HALFWTS ( MAXSTREAMS )

      INTEGER, INTENT (OUT) ::            NMOMENTS
      INTEGER, INTENT (OUT) ::            NSTREAMS_2
      INTEGER, INTENT (OUT) ::            NTOTAL
      INTEGER, INTENT (OUT) ::            N_SUBDIAG
      INTEGER, INTENT (OUT) ::            N_SUPDIAG
      INTEGER, INTENT (OUT) ::            NSTOKES_SQ
      INTEGER, INTENT (OUT) ::            NSTKS_NSTRMS
      INTEGER, INTENT (OUT) ::            NSTKS_NSTRMS_2

      INTEGER, INTENT (OUT) ::            MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) ::   DMAT ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) ::   FLUXVEC ( MAXSTOKES ), DFLUX(MAXSTOKES)

      LOGICAL, INTENT (OUT) :: DO_REAL_EIGENSOLVER ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (OUT) :: DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  local variables
!  ---------------

      INTEGER ::          I, N, M
      INTEGER ::          O1, O2, O1_S, MS2, P, Q, QC, L, EC

      DOUBLE PRECISION :: DMAT_PSOLS ( MAXSTOKES, MAXSTOKES, MAX_PSOLS )
      DOUBLE PRECISION :: DMAT_PSOLS_FLUXVEC ( MAXSTOKES,MAX_PSOLS )

!  Flux vector set unity input.

      FLUXVEC(1) = ONE
      FLUXVEC(2) = ZERO
      FLUXVEC(3) = ZERO
      FLUXVEC(4) = ZERO

!  Number of moments. Isotropic option removed 17 January 2006.

      IF ( DO_RAYLEIGH_ONLY ) THEN
        NMOMENTS = 2
      ENDIF
      IF ( .NOT.DO_RAYLEIGH_ONLY ) THEN
        NMOMENTS = MIN ( 2 * NSTREAMS - 1, NGREEK_MOMENTS_INPUT )
      ENDIF

!  total quadratures (up and down)

      NSTREAMS_2 = 2*NSTREAMS

!  Additional quantities (Stokes)

      NSTOKES_SQ     = NSTOKES * NSTOKES
      NSTKS_NSTRMS   = NSTOKES * NSTREAMS
      NSTKS_NSTRMS_2 = NSTOKES * NSTREAMS_2

!  Mueller index

      DO O1 = 1, MAXSTOKES
        O1_S = MAXSTOKES*(O1 - 1)
        DO O2 = 1, MAXSTOKES
          MUELLER_INDEX(O1,O2) = O1_S + O2
        ENDDO
      ENDDO

!  D matrices, and multiply by Fluxvector

      DO O1 = 1, MAXSTOKES
        DFLUX(O1) = ZERO
        DO O2 = 1, MAXSTOKES
          DMAT(O1,O2) = ZERO
          DO P = 1, 1
            DMAT_PSOLS(O1,O2,P) = ZERO
          ENDDO
        ENDDO
      ENDDO
      MS2 = MAXSTOKES/2
      DO O1 = 1, MS2
        O2 = O1 + MS2
        DMAT(O1,O1) = ONE
        DMAT(O2,O2) = -ONE
        DMAT_PSOLS(O1,O1,1) = ONE
        DMAT_PSOLS_FLUXVEC(O1,1) = FLUXVEC(O1)
        DFLUX(O1) = FLUXVEC(O1)
      ENDDO

!  Set Quadrature abscissae and weights
!  Mk2, 6/21/21. Update use of auxiliary routine GETQUAD2

      CALL GETQUAD2 ( ZERO, ONE, NSTREAMS, QUAD_STREAMS, QUAD_WEIGHTS )

!  set auxiliary quantities

      DO I = 1, NSTREAMS
        QUAD_STRMWTS(I) = QUAD_STREAMS(I)*QUAD_WEIGHTS(I)
        QUAD_HALFWTS(I) = HALF * QUAD_WEIGHTS(I)
      ENDDO

!  size of boundary value problem matrices and vectors

      NTOTAL = NLAYERS*NSTKS_NSTRMS_2

!  number of sub and super diagonals in band matrix (boundary value problem)

      IF ( NLAYERS .EQ. 1 ) THEN
        N_SUBDIAG = 2*NSTKS_NSTRMS - 1
        N_SUPDIAG = 2*NSTKS_NSTRMS - 1
      ELSE
        N_SUBDIAG = 3*NSTKS_NSTRMS - 1
        N_SUPDIAG = 3*NSTKS_NSTRMS - 1
      ENDIF

!  initialise (M = Fourier index) to normal mode.

      DO M = 0, NMOMENTS
        DO N = 1, NLAYERS
          DO_LAYER_SCATTERING(M,N) = .TRUE.
        ENDDO
      ENDDO

      DO M = 3, NMOMENTS
        QC = NMOMENTS - M + 1
        DO N = 1, NLAYERS
          Q = 0
          DO L = M, NMOMENTS
            IF(GREEKMAT_TOTAL_INPUT(L,N,1).EQ.ZERO)Q=Q+1
          ENDDO
          DO_LAYER_SCATTERING(M,N) = (Q.LT.QC)
        ENDDO
      ENDDO

!  Set the Eigensolver Mask. New section 21 December 2005.
!   1. Always set the Real Eigensolver for Rayleigh-only or Scalar-only
!   2. Otherwise, examine the occurrence of the "epsilon" Greek constant

      IF ( DO_RAYLEIGH_ONLY .OR. NSTOKES.EQ.1 ) THEN
       DO M = 0, NMOMENTS
        DO N = 1, NLAYERS
         DO_REAL_EIGENSOLVER(M,N) = .TRUE.
        ENDDO
       ENDDO
      ELSE
       DO M = 0, NMOMENTS
        DO N = 1, NLAYERS
         IF ( DO_LAYER_SCATTERING(M,N) ) THEN
           EC = 0
           DO L = M, NMOMENTS
            IF (GREEKMAT_TOTAL_INPUT(L,N,12).NE.ZERO) EC = EC + 1
           ENDDO
           DO_REAL_EIGENSOLVER(M,N) = (EC.EQ.ZERO)
!mick fix 2/17/11 - added ELSE
         ELSE
           DO_REAL_EIGENSOLVER(M,N) = .TRUE.
         ENDIF
        ENDDO
       ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE MSSEG_HOMOG_DERIVEINPUT

      END MODULE Msseg_Homog_Miscsetups_m

