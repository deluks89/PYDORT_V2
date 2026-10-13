
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

      MODULE MSSeg_Partic_Solutions_m

      PRIVATE
      PUBLIC :: MSSEG_UHOM_SOLUTION,  &
                MSSEG_QBEAM_SOLUTION, &
                MSSEG_UBEAM_SOLUTION, &
                MSSEG_DIRECTBEAM,     &
                MSSEG_PARTIC_PIMATRIX_SETUP

      CONTAINS

!

      SUBROUTINE MSSEG_UHOM_SOLUTION ( &
        GIVEN_LAYER, FOURIER, NSTOKES, NSTREAMS, NMOMENTS, NGEOMS, &
        DO_UPWELLING, DO_DNWELLING, DO_LAYER_SCATTERING, &
        QUAD_HALFWTS, OMEGA_GREEK, HGSF, PGSF, HSOL, PSOL )

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension
!   -- Use module statements, only using relevant quantities

      USE MSSeg_PARS_m, only : MAXSTOKES, MAXLAYERS, MAXSTREAMS, MAX_USER_STREAMS, MAXMOMENTS, &
                               MAXSTREAMS_2, MAXEVALUES, ZERO, ONE, TWO                   

      use MSSeg_HPTypes_Def_m

      IMPLICIT NONE

!  INPUTS
!  ======

!  control and bookkeeping

      INTEGER, INTENT (IN) :: GIVEN_LAYER
      INTEGER, INTENT (IN) :: FOURIER
      INTEGER, INTENT (IN) :: NSTOKES
      INTEGER, INTENT (IN) :: NSTREAMS
      INTEGER, INTENT (IN) :: NMOMENTS
      INTEGER, INTENT (IN) :: NGEOMS

      LOGICAL, INTENT (IN) :: DO_UPWELLING
      LOGICAL, INTENT (IN) :: DO_DNWELLING
      LOGICAL, INTENT (IN) :: DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  Quadrature and optical

      DOUBLE PRECISION, INTENT (IN) ::  QUAD_HALFWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  input type structures (Generalized spherical functions)

      type(Homog_GSF) , intent(in) :: HGSF
      type(Partic_GSF), intent(in) :: PGSF

!  input type structure (Homogeneous solutions)

      type(Homog_solutions), intent(in) :: HSOL

!  output typs structure (particular solutions)

      type(Partic_solutions), intent(inout) :: PSOL

!  Local variables
!  ---------------

      INTEGER ::          UM, J, L, N, M, K, O1, O2, KO1, K0, K1, K2
      DOUBLE PRECISION :: GAUX_R  ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: GAUX_CR ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: GAUX_CI ( 0:MAXMOMENTS, MAXSTOKES )

!  for real solutions

      DOUBLE PRECISION :: SM, SN, SP, SNS, SPS, A5, SUMG

!  for complex solutions

!      DOUBLE PRECISION :: A_USER, A_HELP, B_USER, B_HELP, HPC, HMC
      DOUBLE PRECISION :: SMCR, SNCR, SPCR, SNSCR, SPSCR, SGCR, HPC, HMC
      DOUBLE PRECISION :: SMCI, SNCI, SPCI, SNSCI, SPSCI, SGCI

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

!  If there is no scattering avoid these solutions
!  Warning - Solutions are zeroed in the scalar code before returning
!  @@@ Rob Fix. Zero them here too.

!mick fix 2/1/2011 - modified IF to initialize solutions
!      IF ( .NOT. DO_LAYER_SCATTERING(M,N) ) RETURN

      IF ( .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
       IF ( DO_DNWELLING ) THEN
          PSOL%UHOM_DNDN(:,:,:) = ZERO
          PSOL%UHOM_DNUP(:,:,:) = ZERO
       ENDIF
       IF ( DO_UPWELLING ) THEN
          PSOL%UHOM_UPDN(:,:,:) = ZERO
          PSOL%UHOM_UPUP(:,:,:) = ZERO
       ENDIF
       PSOL%HELPSTOKES(M,:,:)  = ZERO
       RETURN
      ENDIF

!  User defined solutions (Real)
!  -----------------------------

      DO K = 1, HSOL%K_REAL(N)

!  For each moment, do inner sum over computational angles
!  for the positive and negative eigenvectors.
!    Gvec = STOKES. Eq 34. GAUx = Greekmat.Gvec

       DO L = M, NMOMENTS
        DO O1 = 1, NSTOKES
         SPS = ZERO
         DO J = 1, NSTREAMS
          SP = ZERO
          SM = ZERO
          A5 = QUAD_HALFWTS(J)
          DO O2 = 1, NSTOKES
           SP = SP + HGSF%PI_XQP    (L,J,O1,O2)*A5*HSOL%SOLA_XPOS(J,O2,K,N)
           SM = SM + HGSF%PI_XQM_PRE(L,J,O1,O2)*A5*HSOL%SOLB_XNEG(J,O2,K,N)
          ENDDO
          SPS = SPS + SP + SM
         ENDDO
         PSOL%HELPSTOKES(L,K,O1) = SPS
        ENDDO
        DO O1 = 1, NSTOKES
         SUMG = ZERO
         DO O2 = 1, NSTOKES
          SUMG = SUMG + OMEGA_GREEK(L,N,O1,O2) * PSOL%HELPSTOKES(L,K,O2)
         ENDDO
         GAUX_R(L,O1) = SUMG
        ENDDO
       ENDDO

!  Now sum over all harmonic contributions (downwelling)

       IF ( DO_DNWELLING ) THEN
        DO UM = 1, NGEOMS
         DO O1 = 1, NSTOKES
          SP = ZERO
          SN = ZERO
          DO L = M, NMOMENTS
           SPS = ZERO
           SNS = ZERO
           DO O2 = 1, NSTOKES
            SPS = SPS + PGSF%PI_XUP(L,UM,O1,O2)      * GAUX_R(L,O2)
            SNS = SNS + PGSF%PI_XUM_POST(L,UM,O1,O2) * GAUX_R(L,O2)
           ENDDO
           SP = SP + SPS
           SN = SN + SNS
          ENDDO
          PSOL%UHOM_DNDN(UM,O1,K) = SP
          PSOL%UHOM_DNUP(UM,O1,K) = SN
         ENDDO
        ENDDO
       ENDIF

!  Now sum over all harmonic contributions (Upwelling)

       IF ( DO_UPWELLING ) THEN
        DO UM = 1, NGEOMS
         DO O1 = 1, NSTOKES
          SP = ZERO
          SN = ZERO
          DO L = M, NMOMENTS
           SPS = ZERO
           SNS = ZERO
           DO O2 = 1, NSTOKES
            SPS = SPS + PGSF%PI_XUP_PRE(L,UM,O1,O2)  * GAUX_R(L,O2)
            SNS = SNS + PGSF%PI_XUM(L,UM,O1,O2)      * GAUX_R(L,O2)
           ENDDO
           SP = SP + SPS
           SN = SN + SNS
          ENDDO
          PSOL%UHOM_UPDN(UM,O1,K) = SN
          PSOL%UHOM_UPUP(UM,O1,K) = SP
         ENDDO
        ENDDO
       ENDIF

!  Debug (important)

!        IF ( DO_DEBUG_WRITE .AND. DO_FDTEST) THEN
!         DO UM = LOCAL_UM_START, N_USER_STREAMS
!          DO O1 = 1, NSTOKES
!            WRITE(94,'(4I3,1P4E20.10)')M,K,UM,O1, &
!               UHOM_UPUP(UM,O1,K), &
!               UHOM_UPDN(UM,O1,K), &
!               UHOM_DNUP(UM,O1,K), &
!               UHOM_DNDN(UM,O1,K)
!           ENDDO
!         ENDDO
!        ENDIF

!  end loop over real eigensolutions

      ENDDO

!  User defined solutions (Complex)
!  -------------------------------

!  offset

      KO1 = HSOL%K_REAL(N) + 1

      DO K = 1, HSOL%K_COMPLEX(N)

       K0 = 2 * K - 2
       K1 = KO1 + K0
       K2 = K1  + 1

!  For each moment, do inner sum over computational angles
!  for the positive and negative eigenvectors.
!    Gvec = STOKES. Eq 34. GAUx = Greekmat.Gvec

       DO L = M, NMOMENTS
        DO O1 = 1, NSTOKES
         SPSCR = ZERO
         SPSCI = ZERO
         DO J = 1, NSTREAMS
          SPCR = ZERO
          SPCI = ZERO
          SMCR = ZERO
          SMCI = ZERO
          DO O2 = 1, NSTOKES
           HPC  =  HGSF%PI_XQP(L,J,O1,O2)*QUAD_HALFWTS(J)
           SPCR = SPCR + HPC * HSOL%SOLA_XPOS(J,O2,K1,N)
           SPCI = SPCI + HPC * HSOL%SOLA_XPOS(J,O2,K2,N)
           HMC  =  HGSF%PI_XQM_PRE(L,J,O1,O2)*QUAD_HALFWTS(J)
           SMCR = SMCR + HMC * HSOL%SOLB_XNEG(J,O2,K1,N)
           SMCI = SMCI + HMC * HSOL%SOLB_XNEG(J,O2,K2,N)
          ENDDO
          SPSCR = SPSCR + SPCR + SMCR
          SPSCI = SPSCI + SPCI + SMCI
         ENDDO
         PSOL%HELPSTOKES(L,K1,O1) = SPSCR
         PSOL%HELPSTOKES(L,K2,O1) = SPSCI
        ENDDO
        DO O1 = 1, NSTOKES
         SGCR = ZERO
         SGCI = ZERO
         DO O2 = 1, NSTOKES
          SGCR = SGCR + OMEGA_GREEK(L,N,O1,O2) * PSOL%HELPSTOKES(L,K1,O2)
          SGCI = SGCI + OMEGA_GREEK(L,N,O1,O2) * PSOL%HELPSTOKES(L,K2,O2)
         ENDDO
         GAUX_CR(L,O1) = SGCR
         GAUX_CI(L,O1) = SGCI
        ENDDO
       ENDDO

!  Now sum over all harmonic contributions (downwelling)

       IF ( DO_DNWELLING ) THEN
        DO UM = 1, NGEOMS
         DO O1 = 1, NSTOKES
          SPCR = ZERO
          SNCR = ZERO
          SPCI = ZERO
          SNCI = ZERO
          DO L = M, NMOMENTS
           SPSCR = ZERO
           SNSCR = ZERO
           SPSCI = ZERO
           SNSCI = ZERO
           DO O2 = 1, NSTOKES
            HPC = PGSF%PI_XUP(L,UM,O1,O2)
            SPSCR = SPSCR + HPC * GAUX_CR(L,O2)
            SPSCI = SPSCI + HPC * GAUX_CI(L,O2)
            HMC = PGSF%PI_XUM_POST(L,UM,O1,O2)
            SNSCR = SNSCR + HMC * GAUX_CR(L,O2)
            SNSCI = SNSCI + HMC * GAUX_CI(L,O2)
           ENDDO
           SPCR = SPCR + SPSCR
           SPCI = SPCI + SPSCI
           SNCR = SNCR + SNSCR
           SNCI = SNCI + SNSCI
          ENDDO
          PSOL%UHOM_DNDN(UM,O1,K1) = SPCR
          PSOL%UHOM_DNUP(UM,O1,K1) = SNCR
          PSOL%UHOM_DNDN(UM,O1,K2) = SPCI
          PSOL%UHOM_DNUP(UM,O1,K2) = SNCI
         ENDDO
        ENDDO
       ENDIF

!  Now sum over all harmonic contributions (upwelling)

       IF ( DO_UPWELLING ) THEN
        DO UM = 1, NGEOMS
         DO O1 = 1, NSTOKES
          SPCR = ZERO
          SNCR = ZERO
          SPCI = ZERO
          SNCI = ZERO
          DO L = M, NMOMENTS
           SPSCR = ZERO
           SNSCR = ZERO
           SPSCI = ZERO
           SNSCI = ZERO
           DO O2 = 1, NSTOKES
            HPC = PGSF%PI_XUP_PRE(L,UM,O1,O2)
            SPSCR = SPSCR + HPC * GAUX_CR(L,O2)
            SPSCI = SPSCI + HPC * GAUX_CI(L,O2)
            HMC = PGSF%PI_XUM(L,UM,O1,O2)
            SNSCR = SNSCR + HMC * GAUX_CR(L,O2)
            SNSCI = SNSCI + HMC * GAUX_CI(L,O2)
           ENDDO
           SPCR = SPCR + SPSCR
           SPCI = SPCI + SPSCI
           SNCR = SNCR + SNSCR
           SNCI = SNCI + SNSCI
          ENDDO
          PSOL%UHOM_UPDN(UM,O1,K1) = SNCR
          PSOL%UHOM_UPUP(UM,O1,K1) = SPCR
          PSOL%UHOM_UPDN(UM,O1,K2) = SNCI
          PSOL%UHOM_UPUP(UM,O1,K2) = SPCI
         ENDDO
        ENDDO
       ENDIF

!  end loop over complex eigensolutions

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE MSSEG_UHOM_SOLUTION

!

      SUBROUTINE MSSEG_QBEAM_SOLUTION ( &
        GIVEN_LAYER, FOURIER, IBEAM, NSTOKES, NSTREAMS, NMOMENTS,  &
        NSTREAMS_2, NSTKS_NSTRMS, DO_LAYER_SCATTERING, DFLUX,      &
        FLUX_FACTOR, QUAD_STREAMS, OMEGA_GREEK, LAYER_PIS_CUTOFF,  &
        T_DELT_MUBAR, INITIAL_TRANS, AVERAGE_SECANT, PGSF, HGSF, HSOL,   &
        PSOL, STATUS, MESSAGE, TRACE )

!  This is the classical Chandrasekhar beam solution.
!  ( plane parallel or average secant only)
!  Linear Matrix algebra.

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension
!   -- Use module statements, only using relevant quantities

      USE MSSeg_PARS_m, only : MAXSTOKES, MAXLAYERS, MAXBEAMS, MAXSTREAMS, MAXMOMENTS, MAXSTREAMS_2,  &
                               MAXEVALUES, MAXSTRMSTKS, VQS_SUCCESS, VQS_SERIOUS, ZERO, HALF, PI4

      use MSSeg_HPTypes_Def_m

      USE LAPACK_TOOLS_m, only : DGETRF, DGETRS

      IMPLICIT NONE

!  INPUTS
!  ======

!  Basic control

      INTEGER, INTENT (IN) ::           GIVEN_LAYER
      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NMOMENTS

!  bookkeeping

      INTEGER, INTENT (IN) ::           NSTREAMS_2
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  Quadrature

      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STREAMS ( MAXSTREAMS )

!  optical and flux

      DOUBLE PRECISION, INTENT (IN) ::  DFLUX ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  FLUX_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  Beam parameterization

      INTEGER, INTENT (IN) ::           LAYER_PIS_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  input type structures (Generalized spherical functions)

      type(Partic_GSF), intent(in) :: PGSF
      type(Homog_GSF) , intent(in) :: HGSF

!  input type structures (Homogeneous solutions)

      type(Homog_Solutions) , intent(in) :: HSOL

!  OUTPUTS
!  =======

!  output solution type structure

      type(Partic_Solutions), intent(inout) :: PSOL

!  Exception handling

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  help variables

      INTEGER ::          I, J, I1, L, N, M
      INTEGER ::          O1, O2, O3, IR, JR, IROW, JROW
      DOUBLE PRECISION :: SUM, GSUM, INV_X0SQ, SECBAR, F1
      DOUBLE PRECISION :: WSUM, WDIF, TRANS1, TRANS2, WVEC_D, WVEC_U
      DOUBLE PRECISION :: TPA, TMA, S_TPA, S_TMA, S_HELP, HELP

      DOUBLE PRECISION :: HELP_QFUNC ( MAXLAYERS, 0:MAXMOMENTS, MAXSTOKES )

!  Local error handling

      INTEGER ::           INFO
      CHARACTER (LEN=3) :: CI, C3
      CHARACTER (LEN=2) :: CB

!  Initial section
!  ---------------

!  Initialize status

      STATUS  = VQS_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  initialise indices

      N = GIVEN_LAYER
      M = FOURIER
      F1 = FLUX_FACTOR / PI4

!  No particular solution beyond the cutoff layer
!  OR No scattering in this layer then:
!    ---> Zero the boundary layer values and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) .OR. .NOT. DO_LAYER_SCATTERING(M,N) ) THEN

        PSOL%WUPPER(1:NSTREAMS_2,1:NSTOKES,N) = ZERO
        PSOL%WLOWER(1:NSTREAMS_2,1:NSTOKES,N) = ZERO
        PSOL%BVEC  (1:NSTREAMS_2,1:NSTOKES,N) = ZERO

        PSOL%QSUMVEC_SAVE = ZERO
        PSOL%QDIFVEC_SAVE = ZERO
        PSOL%QVEC_SAVE = ZERO
        PSOL%QDIF_SAVE = ZERO
        PSOL%QMAT_SAVE = ZERO
        PSOL%QPIVOT = ZERO

        RETURN
      ENDIF

!  set local values

      SECBAR   = AVERAGE_SECANT(N,IBEAM)
      INV_X0SQ = SECBAR * SECBAR

!  Initialise matrix and column vector for solution
!  Changed at FMI 28 september 2004

      DO J = 1, NSTKS_NSTRMS
        DO I = 1, NSTKS_NSTRMS
          PSOL%QMAT_SAVE(I,J) = ZERO
        ENDDO
        PSOL%QVEC_SAVE(J) = ZERO
      ENDDO

!  Set up sum and difference vectors for Beam source terms
!  ( sum vector may be required again in linearization )
!  Auxiliary matrix for Q functions

      DO L = M, NMOMENTS
        DO O1 = 1, NSTOKES
          SUM = ZERO
          DO O2 = 1, NSTOKES
            GSUM = ZERO
            DO O3 = 1, NSTOKES
              GSUM = GSUM + PGSF%PI_X0P(L,IBEAM,N,O2,O3) * DFLUX(O3)
            ENDDO
            SUM = SUM + OMEGA_GREEK(L,N,O1,O2) * GSUM
          ENDDO
          HELP_QFUNC(N,L,O1) = SUM
        ENDDO
      ENDDO

      DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES
          S_TPA = ZERO
          S_TMA = ZERO
          DO L = M, NMOMENTS
            TPA = ZERO
            TMA = ZERO
            DO O2 = 1, NSTOKES
              TPA = TPA + HGSF%PI_XQP(L,I,O1,O2)      * HELP_QFUNC(N,L,O2)
              TMA = TMA + HGSF%PI_XQM_POST(L,I,O1,O2) * HELP_QFUNC(N,L,O2)
            ENDDO
            S_TPA = S_TPA + TPA
            S_TMA = S_TMA + TMA
          ENDDO
          PSOL%QSUMVEC_SAVE(I,O1) = F1 * ( S_TPA + S_TMA ) / QUAD_STREAMS(I)
          PSOL%QDIFVEC_SAVE(I,O1) = F1 * ( S_TPA - S_TMA ) / QUAD_STREAMS(I)
        ENDDO
      ENDDO

!  solution matrix for the reduced problem
!  ( matrix should be saved in the LU decomposition form)
!  Changed at FMI 28 september

      DO J = 1, NSTKS_NSTRMS
        DO I = 1, NSTKS_NSTRMS
          PSOL%QMAT_SAVE(I,J) = HSOL%EIGENMAT_SAVE(I,J,N)
        ENDDO
      ENDDO

      DO I = 1, NSTKS_NSTRMS
        PSOL%QMAT_SAVE(I,I) = PSOL%QMAT_SAVE(I,I) - INV_X0SQ
      ENDDO

!  RHS vector for the reduced problem
!  ( this vector will be the answer after the linear algebra solution,
!    and may be needed again if there is linearization )

      DO I = 1, NSTREAMS
        IR = ( I - 1 ) * NSTOKES
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          HELP = ZERO
          DO J = 1, NSTREAMS
            S_HELP = ZERO
            DO O2 = 1, NSTOKES
              S_HELP = S_HELP + HSOL%DAB(I,J,O1,O2,N)*PSOL%QSUMVEC_SAVE(J,O2)
            ENDDO
            HELP = HELP + S_HELP
          ENDDO
          PSOL%QVEC_SAVE(IROW) = HELP + PSOL%QDIFVEC_SAVE(I,O1) * SECBAR
        ENDDO
      ENDDO

!  L-U decomposition of the solution matrix

      CALL DGETRF ( NSTKS_NSTRMS, NSTKS_NSTRMS, &
              PSOL%QMAT_SAVE, MAXSTRMSTKS, PSOL%QPIVOT, INFO )

      IF ( INFO .GT. 0 ) THEN
        WRITE(CI, '(I3)' ) INFO
        WRITE(CB, '(I2)' ) IBEAM
        WRITE(C3, '(I3)' ) N
        MESSAGE = 'argument i illegal value, for i = '//CI
        TRACE   = 'DGETRF call in MSSEG_QBEAM_SOLUTION, Beam/layer'// &
            ' numbers = '//CB//', '//C3
        STATUS  = VQS_SERIOUS
        RETURN
      ENDIF

!  Solution of reduced problem by back-substitution

      CALL DGETRS &
          ( 'N', NSTKS_NSTRMS, 1, &
             PSOL%QMAT_SAVE, MAXSTRMSTKS, PSOL%QPIVOT, &
             PSOL%QVEC_SAVE, MAXSTRMSTKS, INFO )

      IF ( INFO .LT. 0 ) THEN
        WRITE(CI, '(I3)' ) INFO
        WRITE(CB, '(I2)' ) IBEAM
        WRITE(C3, '(I3)' ) N
        MESSAGE = 'argument i illegal value, for i = '//CI
        TRACE   = 'DGETRS call in MSSEG_QBEAM_SOLUTION, Beam/layer'// &
            ' numbers = '//CB//', '//C3
        STATUS  = VQS_SERIOUS
        RETURN
      ENDIF

!  Assigning beam particular integral solution vector

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        IR = ( I - 1 ) * NSTOKES
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          HELP = ZERO
          DO J = 1, NSTREAMS
            S_HELP = ZERO
            JR = ( J - 1 ) * NSTOKES
            DO O2 = 1, NSTOKES
              JROW = JR + O2
              S_HELP = S_HELP + HSOL%SAB(I,J,O1,O2,N)*PSOL%QVEC_SAVE(JROW)
            ENDDO
            HELP = HELP + S_HELP
          ENDDO
          PSOL%QDIF_SAVE(IROW) = ( HELP - PSOL%QSUMVEC_SAVE(I,O1) ) / SECBAR
          WSUM = PSOL%QVEC_SAVE(IROW)
          WDIF = PSOL%QDIF_SAVE(IROW)
          WVEC_D = HALF * ( WSUM + WDIF )
          WVEC_U = HALF * ( WSUM - WDIF )
          PSOL%BVEC(I,O1,N)  = WVEC_D
          PSOL%BVEC(I1,O1,N) = WVEC_U
          IF ( O1 .GT. 2 ) PSOL%BVEC(I1,O1,N) = - WVEC_U
        ENDDO
      ENDDO

!  Values at the layer boundaries
!      (transmittance factors have been determined in SETUPS module)

      TRANS1 = INITIAL_TRANS(N,IBEAM)
      TRANS2 = T_DELT_MUBAR(N,IBEAM) * TRANS1
      DO I = 1, NSTREAMS_2
        DO O1 = 1, NSTOKES
          PSOL%WUPPER(I,O1,N) = PSOL%BVEC(I,O1,N)*TRANS1
          PSOL%WLOWER(I,O1,N) = PSOL%BVEC(I,O1,N)*TRANS2
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE MSSEG_QBEAM_SOLUTION

!

      SUBROUTINE MSSEG_UBEAM_SOLUTION ( &
        GIVEN_LAYER, FOURIER, IBEAM,                     &
        DO_UPWELLING, DO_DNWELLING, DO_LAYER_SCATTERING, &
        NSTOKES, NSTREAMS, NMOMENTS, FLUX_FACTOR,        &
        LAYER_PIS_CUTOFF, QUAD_HALFWTS, OMEGA_GREEK,     &
        HGSF, PGSF, PSOL )

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension
!   -- Use module statements, only using relevant quantities

      USE MSSeg_PARS_m, only : MAXSTOKES, MAXLAYERS, MAXSTREAMS, MAX_USER_STREAMS, &
                                 MAXBEAMS, MAXMOMENTS, MAXSTREAMS_2, ZERO, PI4

      use MSSeg_HPTypes_Def_m

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::           GIVEN_LAYER
      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM

      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_DNWELLING
      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NMOMENTS
      DOUBLE PRECISION, INTENT (IN) ::  FLUX_FACTOR

      INTEGER, INTENT (IN) ::           LAYER_PIS_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_HALFWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  input type structures (Generalized spherical functions)

      type(Homog_GSF) , intent(in) :: HGSF
      type(Partic_GSF), intent(in) :: PGSF

!  input/out type structure (Partic solutions)

      type(Partic_Solutions), intent(inout) :: PSOL

!  Local variables
!  ---------------

      INTEGER ::          L, N, M, O1, O2, J, J1
      DOUBLE PRECISION :: SUM2, T2, A5, F1
      DOUBLE PRECISION :: SP, SM, SPS, SPT
      DOUBLE PRECISION :: HELP_Q2 ( 0:MAXMOMENTS, MAXSTOKES )

      DOUBLE PRECISION, DIMENSION(MAXSTOKES)  :: &
        SGW =  (/ 1.0d0, 1.0d0, -1.0d0, -1.0d0 /)

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER
      F1 = FLUX_FACTOR / PI4

      PSOL%UPAR_UP_2 = ZERO  ; PSOL%UPAR_DN_2 = ZERO ; PSOL%HELPSTOKES_BEAM = zero

!  No particular solution beyond the cutoff layer.
!  OR No scattering in this layer

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) .OR. &
             .NOT. DO_LAYER_SCATTERING(M,N) ) RETURN

!  For each moment, do inner sum over computational angles
!    Equivalent to Gvec = GVECTOR. Eq 34. GAUx = Greekmat.Gvec

      DO L = M, NMOMENTS
        DO O1 = 1, NSTOKES
          SPS = ZERO
          DO J = 1, NSTREAMS
           J1 = J + NSTREAMS
           SP = ZERO
           SM = ZERO
           A5 = QUAD_HALFWTS(J)
           DO O2 = 1, NSTOKES
            SP = SP +         HGSF%PI_XQP    (L,J,O1,O2)*A5*PSOL%BVEC(J,O2,N)
            SM = SM + SGW(O2)*HGSF%PI_XQM_PRE(L,J,O1,O2)*A5*PSOL%BVEC(J1,O2,N)
           ENDDO
           SPS = SPS + SP + SM
          ENDDO
          PSOL%HELPSTOKES_BEAM(L,O1) = SPS
         ENDDO
         DO O1 = 1, NSTOKES
          SPT = ZERO
          DO O2 = 1, NSTOKES
           SPT = SPT + OMEGA_GREEK(L,N,O1,O2) * PSOL%HELPSTOKES_BEAM(L,O2)
          ENDDO
          HELP_Q2(L,O1) = SPT
         ENDDO
       ENDDO

!  Now sum over all harmonic contributions (Upwelling)
!    Direct and integarated contributions

      IF ( DO_UPWELLING ) THEN
        DO O1 = 1, NSTOKES
          T2 = ZERO
          DO L = M, NMOMENTS
            SUM2 = ZERO
            DO O2 = 1, NSTOKES
              SUM2 = SUM2 + HELP_Q2(L,O2)*PGSF%PI_XUM(L,IBEAM,O1,O2)
            ENDDO
            T2 = T2 + SUM2
          ENDDO
          PSOL%UPAR_UP_2(O1) = T2
        ENDDO
      ENDIF

!  Now sum over all harmonic contributions (Downwelling)
!    Direct and integarated contributions

      IF ( DO_DNWELLING ) THEN
        DO O1 = 1, NSTOKES
          T2 = ZERO
          DO L = M, NMOMENTS
            SUM2 = ZERO
            DO O2 = 1, NSTOKES
              SUM2 = SUM2 + HELP_Q2(L,O2)*PGSF%PI_XUP(L,IBEAM,O1,O2)
            ENDDO
            T2 = T2 + SUM2
          ENDDO
          PSOL%UPAR_DN_2(O1) = T2
         ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE MSSEG_UBEAM_SOLUTION

!

      SUBROUTINE MSSEG_PARTIC_PIMATRIX_SETUP ( &
        FOURIER, NSTOKES, NLAYERS, NGEOMS, NMOMENTS, &
        COS_SZANGLES, USER_STREAMS, DMAT, HGSF, PGSF )

!  Notes
!  -----

!  This is equivalent to the Legendre setup modules in LIDORT

!---------old comments ---------------------------------
!  This needs modification for the case of Refractive atmosphere,
!  because the solar zenith angle is not constant.
!-------------------------------------------------------

!  PI-matrix setup, following recipe of Siewert (1982).
!  Single Fourier component only.

!  Tested against benchmark results in Vestrucci & Siewert (1984), Problem

!  original coding, September 2002, R. Spurr SAO
!  Multibeam SZA coding, July 2004, R. Spurr SAO
!  Coding for refractive geometry case, R. Spurr, RT Solutions, May 2005

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need

!  Mk3. 04/07/23. Add Fourier dimension to PIMM arrays

      USE MSSeg_PARS_m, only : MAXSTOKES, MAXLAYERS, MAXSTREAMS, MAX_USER_STREAMS, MAXBEAMS, &
                               MAXMOMENTS, MAX_SZANGLES, MAX_ALLSTRMS_P1, ZERO, HALF, ONE, TWO, THREE

      use MSSeg_HPTypes_Def_m

      IMPLICIT NONE

!  Input control

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           NGEOMS

!  Angles

      DOUBLE PRECISION, INTENT (IN) ::  COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  USER_STREAMS  ( MAX_USER_STREAMS )

!  Indices and Dmatrix from DERIVE_INPUTS (mick check)

      DOUBLE PRECISION, INTENT (IN)  :: DMAT ( MAXSTOKES, MAXSTOKES )

!  input type structure (Generalized spherical functions)

      type(Homog_GSF) , intent(in) :: HGSF

!  Output type structure (Generalized spherical functions)

      type(Partic_GSF), intent(inout) :: PGSF

!  Local Output (generalized spherical functions)

      DOUBLE PRECISION :: PI ( 0:MAXMOMENTS, MAX_ALLSTRMS_P1, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: UMXSQ         ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION :: RT_UMXSQ      ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION :: UPXSQ_D_UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION :: PLEG20        ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION :: M2X_D_UMXSQ   ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION :: XUMXSQ        ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION :: X_RT_UMXSQ    ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION :: X2LP1         ( 0:MAXMOMENTS, MAX_ALLSTRMS_P1 )

!  local help variables

      INTEGER          :: NBEAMS, N_USER_STREAMS
      DOUBLE PRECISION :: XSQ, HX, HY, HZ, H1, H2, RCONST_20, RCONST_21
      DOUBLE PRECISION :: PL20, RL20, PL21, RL21, TL21, XM(MAX_ALLSTRMS_P1)
      INTEGER          :: M, I, L, SK1, SK2, SK3, SK4, N
      INTEGER          :: NS, NSA, IB, IP

!  Set integer M = Fourier number

      M = FOURIER
      NSA = 0

!  04/07/23. Initialize them for the first call

      if ( M.eq.0 ) then
        PGSF%PIMM_11 = 0.0d0
        PGSF%PIMM_KM = 0.0d0
      endif

!  Local NSTOKES

!      LOCAL_NSTOKES = NSTOKES
!      IF ( NSTOKES .GT. 1 ) LOCAL_NSTOKES = 4

!  Local NSTOKES

!      LOCAL_NSTOKES = NSTOKES
!      IF ( NSTOKES .GT. 1 ) LOCAL_NSTOKES = 4

      N_USER_STREAMS = NGEOMS
      NBEAMS         = NGEOMS

!  total number of angles

      NS  = N_USER_STREAMS ! User streams
      NSA = NS + NBEAMS    ! Add beams

!  constants

      RCONST_20 = 0.25D0 * DSQRT(6.0D0)
      RCONST_21 = 2.0D0 * RCONST_20

!  local array of all streams (every Fourier component)

      DO I = 1, N_USER_STREAMS
        XM(I) = USER_STREAMS(I)
      ENDDO
      DO IB = 1, NBEAMS
        IP = IB + N_USER_STREAMS
        XM(IP) = COS_SZANGLES(IB)
      ENDDO

!  factors associated with stream values
!   Special case when XM(I) = ONE, factors are zeroed, and
!    those that give singularities are avoided later (see below)
!   R. Spurr and V. Natraj, 16 january 2006

!mick fix 7/28/2014 - do for each fourier right now to make VLIDORT threadsafe
!  1/27/22. No Saving for Every Fourier.

      DO I = 1, NSA
         XSQ = XM(I) * XM(I)
         PLEG20(I)        = HALF*(THREE*XSQ-ONE)
         UMXSQ(I)         = ONE - XSQ
         XUMXSQ(I)        = UMXSQ(I) * XM(I)
         IF ( XM(I) .EQ. ONE ) THEN
            RT_UMXSQ(I)      = ZERO
            X_RT_UMXSQ(I)    = ZERO
            UPXSQ_D_UMXSQ(I) = ZERO
            M2X_D_UMXSQ(I)   = ZERO
         ELSE
            RT_UMXSQ(I)      = DSQRT(UMXSQ(I))
            X_RT_UMXSQ(I)    = RT_UMXSQ(I) * XM(I)
            UPXSQ_D_UMXSQ(I) = (ONE + XSQ) / UMXSQ(I)
            M2X_D_UMXSQ(I)   = - TWO * XM(I) / UMXSQ(I)
         ENDIF
         DO L = 0, NMOMENTS
            X2LP1(L,I) = XM(I) * HGSF%DF_2LP1(L)
         ENDDO
      ENDDO

!  PI calculation
!  --------------

      PI = ZERO

!  M = 0 component. [Siewert (1982), Eqs (28) to (31)]

      IF ( M .EQ. 0 ) THEN

!  .. L = 0
        DO I = 1, NSA
          PI(0,I,1,1) = PI(0,I,1,1) + ONE
          PI(0,I,4,4) = PI(0,I,4,4) + ONE
        ENDDO
!  .. L = 1
        DO I = 1, NSA
          PI(1,I,1,1) = PI(1,I,1,1) + XM(I)
          PI(1,I,4,4) = PI(1,I,4,4) + XM(I)
        ENDDO
!  .. L = 2
        DO I = 1, NSA
          PL20 = PLEG20(I)
          RL20 = RCONST_20*UMXSQ(I)
          PI(2,I,1,1) = PI(2,I,1,1) + PL20
          PI(2,I,2,2) = PI(2,I,2,2) + RL20
          PI(2,I,3,3) = PI(2,I,2,2)
          PI(2,I,4,4) = PI(2,I,1,1)
        ENDDO
!  .. L > 2
        DO L = 2, NMOMENTS - 1
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(L,I)*PI(L,I,SK1,SK2)
                HY = HGSF%YLM_DIAG(L,SK1)*PI(L-1,I,SK1,SK2)
                PI(L+1,I,SK1,SK2) = ( HX - HY ) * HGSF%XLM_DIAG(L,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  M = 1 component. [Siewert (1982), Eqs (32) to (34), plus (35)]

      ELSE IF ( M .EQ. 1 ) THEN

!  .. L = 1
        DO I = 1, NSA
          PI(1,I,1,1) = PI(1,I,1,1) + RT_UMXSQ(I)
          PI(1,I,4,4) = PI(1,I,4,4) + RT_UMXSQ(I)
        ENDDO
!  .. L = 2
        DO I = 1, NSA
          PL21 =   THREE     * X_RT_UMXSQ(I)
          RL21 = - RCONST_21 * X_RT_UMXSQ(I)
          TL21 = + RCONST_21 * RT_UMXSQ(I)
          PI(2,I,1,1) = PI(2,I,1,1) + PL21
          PI(2,I,2,2) = PI(2,I,2,2) + RL21
          PI(2,I,2,3) = PI(2,I,2,3) + TL21
          PI(2,I,3,3) = PI(2,I,2,2)
          PI(2,I,4,4) = PI(2,I,1,1)
          PI(2,I,3,2) = PI(2,I,2,3)
        ENDDO
!  .. L > 2
        DO L = 2, NMOMENTS - 1
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(L,I)*PI(L,I,SK1,SK2)
                HY = HGSF%YLM_DIAG(L,SK1)*PI(L-1,I,SK1,SK2)
                HZ = ZERO
                DO SK3 = 1, NSTOKES
                  HZ = HZ + HGSF%ZLM(L,SK1,SK3)*PI(L,I,SK3,SK2)
                ENDDO
                PI(L+1,I,SK1,SK2) = ( HX - HY + HZ ) * HGSF%XLM_DIAG(L,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  M > 1 components. [Siewert (1982), Eqs (36) to (38), plus (35)]
!  Limiting case of XM(I) = ONE requires special treatment to avoid NaN
!   R. Spurr and V. Natraj, 16 january 2006

!  04/07/23. Re-set PIMM arrays with Fourier index.
!     -- for M = 2, Set the PIMM arrays with Fourier index 2
!     -- for M > 2, set PIMM arrays with index M interms of PIMM arrays with index M-1

      ELSE

!  .. L = M
!  --------

!  04/07/23. Re-set PIMM arrays with Fourier indices.
!     -- for M = 2, Set the PIMM arrays with Fourier index 2
!     -- for M > 2, set PIMM arrays with index M interms of PIMM arrays with index M-1

        IF ( M .EQ. 2 ) THEN
          DO I = 1, NSA
            PGSF%PIMM_11(M,I) =   THREE * UMXSQ(I)
            IF ( XM(I).EQ.ONE) THEN
              PGSF%PIMM_KM(M,I) =   RCONST_21
            ELSE
              PGSF%PIMM_KM(M,I) =   RCONST_21 * UMXSQ(I)
            ENDIF
          ENDDO
        ELSE
          H1 = HGSF%DF_2LP1(M-1)
          H2 = HGSF%DF_LP1(M-1)*H1/HGSF%DF_RT_LP3XLM1(M-1)
          DO I = 1, NSA
!  Old       PGSF%PIMM_11(I)   = H1 * RT_UMXSQ(I) * PGSF%PIMM_11(I)
            PGSF%PIMM_11(M,I) = H1 * RT_UMXSQ(I) * PGSF%PIMM_11(M-1,I)
            IF ( XM(I).EQ.ONE) THEN
              PGSF%PIMM_KM(M,I) = ZERO
            ELSE
!  Old         PGSF%PIMM_KM(I)   = H2 * RT_UMXSQ(I) * PGSF%PIMM_KM(I)
              PGSF%PIMM_KM(M,I) = H2 * RT_UMXSQ(I) * PGSF%PIMM_KM(M-1,I)
            ENDIF
          ENDDO
        ENDIF

!  04/07/23. Use PIMM arrays with current Fourier index.

        DO I = 1, NSA
          PI(M,I,1,1) = PI(M,I,1,1) + PGSF%PIMM_11(M,I)
          IF ( XM(I).EQ.ONE ) THEN
            PI(M,I,2,2) = PI(M,I,2,2) + PGSF%PIMM_KM(M,I) * TWO
            PI(M,I,2,3) = PI(M,I,2,3) - PGSF%PIMM_KM(M,I) * TWO
          ELSE
            PI(M,I,2,2) = PI(M,I,2,2) + PGSF%PIMM_KM(M,I) * UPXSQ_D_UMXSQ(I)
            PI(M,I,2,3) = PI(M,I,2,3) + PGSF%PIMM_KM(M,I) * M2X_D_UMXSQ(I)
          ENDIF
          PI(M,I,3,3) = PI(M,I,2,2)
          PI(M,I,4,4) = PI(M,I,1,1)
          PI(M,I,3,2) = PI(M,I,2,3)
        ENDDO

!  .. L = M + 1
!  ------------

        IF ( M .LT. NMOMENTS ) THEN
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(M,I)*PI(M,I,SK1,SK2)
                HZ = ZERO
                DO SK3 = 1, NSTOKES
                  HZ = HZ + HGSF%ZLM(M,SK1,SK3)*PI(M,I,SK3,SK2)
                ENDDO
                PI(M+1,I,SK1,SK2) = ( HX + HZ ) * HGSF%XLM_DIAG(M,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!  .. L > M + 1
        DO L = M + 1, NMOMENTS - 1
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(L,I)*PI(L,I,SK1,SK2)
                HY = HGSF%YLM_DIAG(L,SK1)*PI(L-1,I,SK1,SK2)
                HZ = ZERO
                DO SK3 = 1, NSTOKES
                  HZ = HZ + HGSF%ZLM(L,SK1,SK3)*PI(L,I,SK3,SK2)
                ENDDO
                PI(L+1,I,SK1,SK2) = ( HX - HY + HZ ) * HGSF%XLM_DIAG(L,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDIF

!  Normalized output.
!  ------------------

      DO L = M, NMOMENTS

!  .. at positive user_defined angles

        DO I = 1, N_USER_STREAMS
          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              PGSF%PI_XUP(L,I,SK1,SK2) = PI(L,I,SK1,SK2) * HGSF%PINORM(L,M)
            ENDDO
          ENDDO
          DO SK1 = 1, NSTOKES
             DO SK2 = 1, NSTOKES
                H2 = ZERO
                H1 = ZERO
                DO SK3 = 1, NSTOKES
                  H1 = H1 + DMAT(SK1,SK3) * PGSF%PI_XUP(L,I,SK3,SK2)
                  H2 = H2 + DMAT(SK3,SK2) * PGSF%PI_XUP(L,I,SK1,SK3)
                ENDDO
                PGSF%PI_XUP_PRE(L,I,SK1,SK2)  = H1
                PGSF%PI_XUM_POST(L,I,SK1,SK2) = H2 * HGSF%PISIGN(L,M)
                H1 = ZERO
                DO SK3 = 1, NSTOKES
                  H2 = ZERO
                  DO SK4 = 1, NSTOKES
                    H2 = H2 + PGSF%PI_XUP(L,I,SK3,SK4)*DMAT(SK4,SK2)
                  ENDDO
                  H1 = H1 + DMAT(SK1,SK3)*H2
                ENDDO
                PGSF%PI_XUM(L,I,SK1,SK2) = H1 * HGSF%PISIGN(L,M)
             ENDDO
          ENDDO
        ENDDO

!  .. at solar zenith angles

        DO IB = 1, NBEAMS
           IP = IB + N_USER_STREAMS
           DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                PGSF%PI_X0P(L,IB,1,SK1,SK2) = PI(L,IP,SK1,SK2) * HGSF%PINORM(L,M)
                DO N = 2, NLAYERS
                  PGSF%PI_X0P(L,IB,N,SK1,SK2)= PGSF%PI_X0P(L,IB,1,SK1,SK2)
                ENDDO
              ENDDO
           ENDDO
        ENDDO

!  end loop moments

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE MSSEG_PARTIC_PIMATRIX_SETUP

      SUBROUTINE MSSEG_DIRECTBEAM ( &
        DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,  &
        FOURIER, NSTOKES, NSTREAMS, NBEAMS, DELTA_FACTOR,            &
        FLUX_FACTOR, FLUXVEC, MUELLER_INDEX, COS_SZANGLES,           &
        LAMBERTIAN_ALBEDO, BRDF_F_0, SLTERM_ISOTROPIC, SLTERM_F_0,   &
        TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM, &
        ATMOS_ATTN, DIRECT_BEAM )

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need

      USE MSSeg_PARS_m, Only : MAXSTREAMS, MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAXMOMENTS, MAX_SZANGLES, &
                                 MAXSTOKES, MAXSTOKES_SQ, ZERO, ONE, FOUR, PI4, MAX_TAU_QPATH, MAX_TAU_SPATH, MAX_TAU_UPATH


      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::           DO_SURFACE_LEAVING
      LOGICAL, INTENT (IN) ::           DO_SL_ISOTROPIC

!  input control integers

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NBEAMS

!  input control and geometry

      DOUBLE PRECISION, INTENT (IN) ::  DELTA_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  FLUX_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  FLUXVEC ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  COS_SZANGLES ( MAX_SZANGLES )
      INTEGER         , INTENT (IN) ::  MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

!  Input  reflectance

      DOUBLE PRECISION, INTENT (IN) ::  LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  BRDF_F_0 ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS   )

!  Input surface leaving

      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_ISOTROPIC ( MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_F_0 ( MAXSTOKES, MAXSTREAMS, MAXBEAMS )

!  Input attenuation

      DOUBLE PRECISION, INTENT (IN) ::  TRANS_SOLAR_BEAM ( MAXBEAMS )
      LOGICAL, INTENT (IN) ::           DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  output

      DOUBLE PRECISION, INTENT (OUT) :: ATMOS_ATTN ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )

!  Local variables
!  ---------------

      DOUBLE PRECISION :: X0_FLUX, X0_BOA, ATTN, REFLEC, SL, HELP, SUM
      INTEGER          :: I, O1, O2, IB, M, OM

!  Initialize
!  ----------

!   Safety first!  Return if there is no reflection.

      DIRECT_BEAM      = ZERO

!  Fourier component

      M = FOURIER

!  Attenuation of solar beam
!  -------------------------

!  New code to deal with refractive geometry case
!   R. Spurr, 7 May 2005. RT Solutions Inc.

      DO IB = 1, NBEAMS
       IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN

!  There should be no flux factor here.
!    Bug fixed 18 November 2005. Earlier Italian Job!!
!       X0_FLUX        = FOUR * X0_BOA * FLUX_FACTOR / DELTA_FACTOR

        X0_BOA         = COS_SZANGLES(IB)
        X0_FLUX        = FOUR * X0_BOA / DELTA_FACTOR
        X0_FLUX        = FLUX_FACTOR * X0_FLUX / PI4             ! New
        ATTN           = X0_FLUX * TRANS_SOLAR_BEAM(IB)
        ATMOS_ATTN(IB) = ATTN

!  Lambertian case
!  ---------------

!  Set value for Fourier = 0
!  Natural light only, so only the first Stokes component is nonzero

        IF ( DO_LAMBERTIAN_SURFACE ) THEN
          REFLEC = ATTN * LAMBERTIAN_ALBEDO
          IF ( FOURIER .EQ. 0 ) THEN
            DIRECT_BEAM(1:NSTREAMS,IB,1) = REFLEC
          ENDIF
        ENDIF

!  Non-Lambertian case
!  -------------------

        IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN

!  Solar beam reflected into quad directions

          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              SUM = ZERO
              DO O2 = 1, NSTOKES
                OM = MUELLER_INDEX(O1,O2)
                SUM = SUM + FLUXVEC(O2) * BRDF_F_0(OM,I,IB)
              ENDDO
              DIRECT_BEAM(I,IB,O1) = ATTN * SUM
            ENDDO
          ENDDO


        ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ START
!  New Surface-Leaving stuff 17 May 2012
!  Corrected implementation, 30 July 2012
!    Normalized to Flux-factor / DELTA_Factor
!    Delta_Factor = 1.0 for the Isotropic or non-iso Fourier = 0 cases

        IF ( DO_SURFACE_LEAVING ) THEN
          HELP = FLUX_FACTOR / DELTA_FACTOR
          IF ( DO_SL_ISOTROPIC .and. M.EQ.0 ) THEN
            DO O1 = 1, NSTOKES
              SL = SLTERM_ISOTROPIC(O1,IB) * HELP
              DO I = 1, NSTREAMS
                DIRECT_BEAM(I,IB,O1) = DIRECT_BEAM(I,IB,O1) + SL
              ENDDO
            ENDDO
          ELSE
            DO O1 = 1, NSTOKES
              DO I = 1, NSTREAMS
                SL = SLTERM_F_0(O1,I,IB) * HELP
                DIRECT_BEAM(I,IB,O1) = DIRECT_BEAM(I,IB,O1) + SL
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END

!  end direct beam calculation

       ENDIF
      ENDDO

!  finish

      RETURN
      END SUBROUTINE MSSEG_DIRECTBEAM

      END MODULE MSSeg_Partic_Solutions_m

