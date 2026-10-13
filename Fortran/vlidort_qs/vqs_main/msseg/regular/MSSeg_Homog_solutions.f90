
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

      MODULE MSSeg_Homog_Solutions_m

!  05 April 2023. Revision

      PUBLIC   :: MSSEG_QHOM_SOLUTION, MSSEG_HOMOG_PIMATRIX_SETUP

      CONTAINS

      SUBROUTINE MSSEG_QHOM_SOLUTION ( &
        GIVEN_LAYER, FOURIER, NSTOKES, NSTREAMS, NMOMENTS, NSTKS_NSTRMS,      &
        QUAD_STREAMS, QUAD_HALFWTS, DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING, &
        ASYMTX_TOLERANCE, DELTAU_VERT, OMEGA_GREEK, HGSF, &
        HSOL, REAL_KSQ, IMAG_KSQ, LEFT_EVEC, RITE_EVEC,   &
        EIGENDEGEN, EIGENMASK_R, EIGENMASK_C, KEIGEN_CSQ, &
        STATUS, MESSAGE, TRACE )

!  Numerical solution of Eigenproblem.

!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension
!   -- Use module statements, only using relevant quantities
!   -- Add ASYMTX_TOLERANCE to argument list

      USE MSSEG_PARS_m, Only : MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAX_USER_LEVELS, &
                               MAXSTREAMS_2, MAXSTRMSTKS, MAXEVALUES, MAXMOMENTS, &
                               VQS_SERIOUS, VQS_SUCCESS, &
                               ZERO, ONE, HALF, MINUS_ONE, MAX_TAU_QPATH

      use MSSEG_HPTypes_Def_m

      USE MSSEG_AUX_m   , only : ASYMTX 
      USE LAPACK_TOOLS_m, only : DGEEV

      IMPLICIT NONE

!  Inputs
!  ------

!  control numbers

      INTEGER, INTENT (IN) ::           GIVEN_LAYER
      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS

!  Quadrature

      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_HALFWTS ( MAXSTREAMS )

!  Bookkeeping

      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      LOGICAL, INTENT (IN) ::           DO_REAL_EIGENSOLVER ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  ASYMTX_TOLERANCE

!  optical

      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  input type structure (Generalized spherical functions)

      type(Homog_GSF), intent(inout) :: HGSF

!  Output
!  ------

!  output solution type structure

      type(Homog_Solutions), intent(inout) :: HSOL

!  Auxiliary output

      DOUBLE PRECISION, INTENT (INOUT) ::  REAL_KSQ ( MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) ::  IMAG_KSQ ( MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) ::  LEFT_EVEC ( MAXSTRMSTKS, MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) ::  RITE_EVEC ( MAXSTRMSTKS, MAXSTRMSTKS )
      LOGICAL, INTENT (INOUT) ::           EIGENDEGEN ( MAXSTRMSTKS, MAXLAYERS )
      INTEGER, INTENT (INOUT) ::           EIGENMASK_R ( MAXEVALUES )
      INTEGER, INTENT (INOUT) ::           EIGENMASK_C ( MAXEVALUES )
      DOUBLE PRECISION, INTENT (INOUT) ::  KEIGEN_CSQ ( MAXEVALUES )

!  Exception handling

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  local matrix for eigenvalue computation

      DOUBLE PRECISION :: EIGENMAT ( MAXSTRMSTKS, MAXSTRMSTKS )

!  output from Eigenpackage module ASYMTX
!    Following output now stored in commons:
!      DOUBLE PRECISION REAL_KSQ(MAXSTRMSTKS)
!      DOUBLE PRECISION RITE_EVEC(MAXSTRMSTKS,MAXSTRMSTKS)

!  Additional diagnostic output from ASYMTX, and Work space.
!    Dimension of scratch array WK increased to 4 times the basic dimens
!    R. Spurr, RTS, December 2005: More work space needed to avoid seg-f
!    Internal dimensioning in ASYMTX is the same.

!     DOUBLE PRECISION ::  WK(16*MAXSTRMSTKS)
      DOUBLE PRECISION ::  WK(4*MAXSTRMSTKS)
      INTEGER ::           IER
      LOGICAL ::           ASYMTX_FAILURE
      CHARACTER (LEN=3) :: CN, CI

!  output for Eigenpackage DGEEV ( LAPACK)
!    Following now stored in commons...............
!      DOUBLE PRECISION REAL_KSQ(MAXSTRMSTKS),IMAG_KSQ(MAXSTRMSTKS)
!      DOUBLE PRECISION LEFT_EVEC(MAXSTRMSTKS,MAXSTRMSTKS)
!      DOUBLE PRECISION RITE_EVEC(MAXSTRMSTKS,MAXSTRMSTKS)

!  Additional diagnostic output and Work space from DGEEV.
!    Dimension LWORK increased from 4 to 64 times the basic dimension
!    J. Kujanpaa, FMI, August 2005: More work space for DGEEV needed.

      INTEGER, PARAMETER :: LWORK = 64*MAXSTRMSTKS
      DOUBLE PRECISION ::   WORK ( LWORK )
      INTEGER ::            DGEEV_INFO

!  optimization of Sum and Difference matrices SAB and DAB
!  ...requires intermediate storage of these two arrays.
!   Timing tests by J. Kujanpaa, FMI, August 2005.

      DOUBLE PRECISION :: &
          H2PARR ( MAXSTOKES, 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS ), &
          H2MARR ( MAXSTOKES, 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS )

!  Miscellaneous local variables

      INTEGER ::          I, J, I1, L, N,  NMINST
      INTEGER ::          M, AA, AA1, K, K_R, K_C, K_AA, KO1, K0, K1, K2
      INTEGER ::          O1, O2, O3, O4, IR, IROW, JC, JCOL
      DOUBLE PRECISION :: DP, DM, FAC, SC_R, XINV, HELP
      DOUBLE PRECISION :: H1P, H1M, NORM_R
      DOUBLE PRECISION :: FWD_H1, FWD_H2
      DOUBLE PRECISION :: XSQ_R, YSQ_R, RKSQ, IKSQ, KMUT, KVAL
      DOUBLE PRECISION :: CARG, SARG, ARG, MOD_KVAL, H1, H2
      DOUBLE PRECISION :: FWD_R, FWD_I, H2P, H2M, FPD, FND

!  original complex variables were:
!      COMPLEX*16      FWD_H1_C, FWD_H2_C, ADJ_H1_C, ADJ_H2_C
!      COMPLEX*16      HELP_C, SC_C, H1_C, H2_C

!  replacement for complex variables

      DOUBLE PRECISION :: HELP_CR, HELP_CI
      DOUBLE PRECISION :: FWD_H1_CR, FWD_H1_CI
      DOUBLE PRECISION :: FWD_H2_CR, FWD_H2_CI
      DOUBLE PRECISION :: SC_CR, SC_CI, MODKRT
      DOUBLE PRECISION :: FPD1, FPD2, FND1, FND2

!  Code start
!  ----------

!  initialise exception handling

      STATUS  = VQS_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

!  Zeroing some of the type  structure output

      HSOL%FWD_SUMVEC = zero
      HSOL%FWD_DIFVEC = zero
      HSOL%SAB(:,:,:,:,N)       = zero
      HSOL%DAB(:,:,:,:,N)       = zero
      HSOL%SOLA_XPOS(:,:,:,N)   = zero
      HSOL%SOLB_XNEG(:,:,:,N)   = zero

!  Initialization

      KO1 = 0

!  initialise status

      STATUS = VQS_SUCCESS

!  Flipper number

      NMINST = MIN(2,NSTOKES)

!  If there is no scattering in this layer, we can set
!  the solution vectors directly and move on.
!    Just a transmitttance case.
!mick fix 2/12/2019 - initialized  EIGENMASK_R, RITE_EVEC, FWD_SUMVEC, FWD_DIFVEC, SAB, DAB

      IF ( .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
        HSOL%K_COMPLEX(N) = 0
        HSOL%K_REAL(N)    = NSTKS_NSTRMS
        RITE_EVEC = ZERO
        DO AA = 1, NSTKS_NSTRMS
          EIGENMASK_R(AA) = AA
          RITE_EVEC(AA,AA) = ONE
        ENDDO
        DO K = 1, HSOL%K_REAL(N)
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              HSOL%FWD_SUMVEC(I,O1,K) = ONE
            ENDDO
          ENDDO
        ENDDO
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          XINV = ONE / QUAD_STREAMS(I)
          DO O1 = 1, NMINST
            K = NSTOKES*(I-1) + O1
            HSOL%KEIGEN(K,N) = XINV
            HSOL%SOLA_XPOS(I,O1,K,N)  = ONE
            HSOL%SOLB_XNEG(I1,O1,K,N) = ONE
          ENDDO
          DO O1 = 3, NSTOKES
            K = NSTOKES*(I-1) + O1
            HSOL%KEIGEN(K,N) = XINV
            HSOL%SOLA_XPOS(I,O1,K,N)  = MINUS_ONE
            HSOL%SOLB_XNEG(I1,O1,K,N) = MINUS_ONE
          ENDDO
        ENDDO
        GO TO 3456
      ENDIF

!  Scattering solutions
!  ====================

!  Construct Eigenmatrix
!  ---------------------

!  zero the eigensolution matrix
!   Not needed. Redundancy test by J. Kujanpaa, FMI, August 2005.
!      DO I = 1, NSTKS_NSTRMS
!        DO J = 1, NSTKS_NSTRMS
!          EIGENMAT(I,J) = ZERO
!        ENDDO
!      ENDDO

!  Develop Sum and Difference matrices, part I
!   Introduction of holding arrays by J. Kujanpaa, FMI, August 2005

      DO J = 1, NSTREAMS
       DO O2 = 1, NSTOKES
         DP = ZERO
         DM = ZERO
         DO L = M, NMOMENTS
          DO O3 = 1, NSTOKES
           H2P = ZERO
           H2M = ZERO
           DO O4 = 1, NSTOKES
            H2P = H2P + OMEGA_GREEK(L,N,O3,O4)*HGSF%PI_XQP(L,J,O4,O2)
            H2M = H2M + OMEGA_GREEK(L,N,O3,O4)*HGSF%PI_XQM_PRE(L,J,O4,O2)
           ENDDO
           H2PARR(O3,L,O2,J) = H2P
           H2MARR(O3,L,O2,J) = H2M
          ENDDO
         ENDDO
        ENDDO
      ENDDO

!  Develop Sum and Difference matrices, part II of calculation
!   Use of holding arrays by J. Kujanpaa, FMI, August 2005

      DO I = 1, NSTREAMS
       XINV = - ONE / QUAD_STREAMS(I)
       DO J = 1, NSTREAMS
        FAC = XINV * QUAD_HALFWTS(J)
        DO O1 = 1, NSTOKES
         DO O2 = 1, NSTOKES
          DP = ZERO
          DM = ZERO
          DO L = M, NMOMENTS
           H1P = ZERO
           H1M = ZERO
           DO O3 = 1, NSTOKES
            H1P = H1P + HGSF%PI_XQP(L,I,O1,O3)*H2PARR(O3,L,O2,J)
            H1M = H1M + HGSF%PI_XQP(L,I,O1,O3)*H2MARR(O3,L,O2,J)
           ENDDO
           DP = DP + H1P
           DM = DM + H1M
          ENDDO
          HSOL%SAB(I,J,O1,O2,N) = FAC * ( DP + DM )
          HSOL%DAB(I,J,O1,O2,N) = FAC * ( DP - DM )
         ENDDO
        ENDDO
       ENDDO
       DO O1 = 1, NSTOKES
        HSOL%SAB(I,I,O1,O1,N) = HSOL%SAB(I,I,O1,O1,N) - XINV
        HSOL%DAB(I,I,O1,O1,N) = HSOL%DAB(I,I,O1,O1,N) - XINV
       ENDDO
      ENDDO

!  Note. In some circumstances, H2PARR = H2MARR
!        Then DP = DM in the above code, and DAB is Diagonal.
!        Then there are pair-degenerate eigenvalues, equal to the
!        discrete ordinate directions.

!  Compute Eigenmatrix

      DO I = 1, NSTREAMS
       IR = NSTOKES*(I-1)
       DO J = 1, NSTREAMS
        JC = NSTOKES*(J-1)
        DO O1 = 1, NSTOKES
         IROW = IR + O1
         DO O2 = 1, NSTOKES
          JCOL = JC + O2
          H1 = ZERO
          DO K = 1, NSTREAMS
           H2 = ZERO
           DO O3 = 1, NSTOKES
            H2 = H2 + HSOL%DAB(I,K,O1,O3,N) * HSOL%SAB(K,J,O3,O2,N)
           ENDDO
           H1 = H1 + H2
          ENDDO
          EIGENMAT(IROW,JCOL) = H1
         ENDDO
        ENDDO
       ENDDO
      ENDDO

!  save Eigenmatrix (original is destroyed, want to use later)
!    Looping order changed by J. Kujanpaa, FMI, August 2005

      DO J = 1, NSTKS_NSTRMS
        DO I = 1, NSTKS_NSTRMS
         HSOL%EIGENMAT_SAVE(I,J,N) = EIGENMAT(I,J)
        ENDDO
      ENDDO

!      if (n.eq.24) then
!      do j = 1, 9
!      if(nstokes.eq.3)write(33,'(i4,1p12e11.3)')J,(EIGENMAT(i,J),i=1,9)
!      if(nstokes.eq.4)write(34,'(i4,1p12e11.3)')J,(EIGENMAT(i,J),i=1,12
!      enddo
!      pause
!      endif

!  Debug for testing new eigensolver
!  7 August 2006. Result --> NEW LAPACK ROUTINES NOT FASTER

!      IF ( N.EQ.1) then
!        OPEN(88,file='esolver.inp',status='unknown')
!      ENDIF
!      write(88,'(2i5,l2)')n,NSTKS_NSTRMS,DO_REAL_EIGENSOLVER(M,N)
!      DO I = 1, NSTKS_NSTRMS
!        write(88,'(1p32e20.10)')(eigenmat(i,j),j=1,NSTKS_NSTRMS)
!      enddo
!      if ( n.eq.nlayers)then
!        close(88)
!        pause'end debug'
!      ENDIF

!  Eigensolver package
!  -------------------

      IF ( DO_REAL_EIGENSOLVER(M,N) ) THEN

!  Here is the DISORT ASYMTX package, as used in LIDORT scalar codes
!    This is sufficient for Real Symmetric Eigenmatrices.
!  Mk2, 6/21/21.. Add ASYMTX_TOLERANCE to the call

       CALL  ASYMTX ( ASYMTX_TOLERANCE, &
              EIGENMAT, NSTKS_NSTRMS, MAXSTRMSTKS, MAXSTRMSTKS, &
              RITE_EVEC, REAL_KSQ, IER, WK, &
              MESSAGE, ASYMTX_FAILURE )

!  Exception handling 1

       IF ( ASYMTX_FAILURE  ) THEN
         WRITE(CN,'(I3)')N
         TRACE   ='ASYMTX error in MSSEG_QHOM_SOLUTION, Layer='//CN
         STATUS  = VQS_SERIOUS
         RETURN
       ENDIF

!  Exception handling 2

       IF ( IER.GT.0 ) THEN
         WRITE(CI,'(I3)')IER
         WRITE(CN,'(I3)')N
         MESSAGE = 'eigenvalue '//CI//' has not converged'
         TRACE   = 'ASYMTX error in MSSEG_QHOM_SOLUTION, Layer ='//CN
         STATUS  = VQS_SERIOUS
         RETURN
       ENDIF

!  renormalize - this is vital for ASMTX output

       DO AA = 1, NSTKS_NSTRMS
        NORM_R = ZERO
        DO I = 1, NSTKS_NSTRMS
          NORM_R = NORM_R + RITE_EVEC(I,AA)*RITE_EVEC(I,AA)
        ENDDO
        NORM_R = DSQRT(NORM_R)
        DO I = 1, NSTKS_NSTRMS
          RITE_EVEC(I,AA) = RITE_EVEC(I,AA) / NORM_R
        ENDDO
       ENDDO

      ELSE

!  Complex roots : Must use DGEEV (LAPACK module)
!  DGEEV is 1.7-2.0 times slower than ASYMTX because look for complex ro
!  Also first call to DGEEV is 20 times slower.

       CALL DGEEV &
            ( 'V', 'V', NSTKS_NSTRMS, EIGENMAT, &
              MAXSTRMSTKS, REAL_KSQ, IMAG_KSQ, &
              LEFT_EVEC, MAXSTRMSTKS, RITE_EVEC, MAXSTRMSTKS, &
              WORK, LWORK, DGEEV_INFO )

!  Exception handling

       IF ( DGEEV_INFO .NE. 0 ) THEN
        WRITE(CN,'(I3)')N
        TRACE = 'DGEEV call, layer '//CN//' in MSSEG_QHOM_SOLUTION'
        IF ( DGEEV_INFO .LT. 0 ) THEN
         WRITE(CI,'(I3)')IABS(DGEEV_INFO)
         MESSAGE='Argument # '//CI//' had an illegal value'
        ELSE
         MESSAGE='QR algorithm in DGEEV failed to compute all e-values'
        ENDIF
        STATUS = VQS_SERIOUS
        RETURN
       ENDIF

      ENDIF

!  Sort real and complex eigenvalues; get masks
!  --------------------------------------------

!  Initialise counts

       K_R  = 0
       K_C  = 0
       K_AA = 0

!  all eigenvalues are real if ASYMTX is used

      IF ( DO_REAL_EIGENSOLVER(M,N) ) THEN
       DO AA = 1, NSTKS_NSTRMS
        K_R = K_R + 1
        EIGENMASK_R(K_R) = AA
       ENDDO
      ENDIF

!  Degeneracy mask

      IF ( DO_REAL_EIGENSOLVER(M,N) ) THEN
       DO AA = 1, NSTKS_NSTRMS
         KVAL = DSQRT(REAL_KSQ(AA))
         EIGENDEGEN(AA,N) = .FALSE.
         DO K = 1, NSTKS_NSTRMS
           IF (K.NE.AA.AND..NOT.EIGENDEGEN(AA,N)) THEN
             KMUT= DSQRT(REAL_KSQ(K))
             EIGENDEGEN(AA,N) = (DABS(KMUT-KVAL).LT.1.0D-04)
           ENDIF
         ENDDO
       ENDDO
      ENDIF

!  For the DGEEV eigensolver--
!    real eigenvalues if IMAG_KSQ = 0, otherwise complex in pairs

      IF ( .NOT.DO_REAL_EIGENSOLVER(M,N) ) THEN
       DO AA = 1, NSTKS_NSTRMS
        IF ( IMAG_KSQ(AA) .EQ. ZERO ) THEN
          K_R = K_R + 1
          EIGENMASK_R(K_R) = AA
        ELSE
          K_AA = K_AA + 1
          IF ( K_AA .EQ. 1 ) THEN
            K_C = K_C + 1
            EIGENMASK_C(K_C) = AA
          ELSE IF ( K_AA .EQ. 2 ) THEN
            K_AA = 0
          ENDIF
        ENDIF
       ENDDO
      ENDIF

!  save number of eigenvalues

      HSOL%K_REAL(N)    = K_R
      HSOL%K_COMPLEX(N) = K_C

!  Real Solutions
!  --------------

!  start loop over real eigenvalues

      DO K = 1, HSOL%K_REAL(N)

!  get eigenvector entry

        AA = EIGENMASK_R(K)

!  eigenvalue and separation constants (save the eigenvalues)

        HSOL%KEIGEN(K,N) = DSQRT(REAL_KSQ(AA))
        SC_R    = ONE / HSOL%KEIGEN(K,N)

!  set Right eigenvector norms
!    This piece of code should not be necessary as
!    (a) DGEEV norms are all = 1 upon output
!    (b) ASMTYX output has been renormalized to unity

        NORM_R = ZERO
        DO I = 1, NSTKS_NSTRMS
          NORM_R = NORM_R + RITE_EVEC(I,AA)*RITE_EVEC(I,AA)
        ENDDO
        NORM_R = DSQRT(NORM_R)

!  Forward Sum-vector = normalized Right eigenvector

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            HSOL%FWD_SUMVEC(I,O1,K) = RITE_EVEC(IROW,AA)/NORM_R
          ENDDO
        ENDDO

!  Find Forward difference vectors (Siewert's notation)

        DO I = 1, NSTREAMS
         DO O1 = 1, NSTOKES
          FWD_H1 = ZERO
          DO J = 1, NSTREAMS
           FWD_H2 = ZERO
           DO O2 = 1, NSTOKES
            FWD_H2 = FWD_H2 + HSOL%SAB (I,J,O1,O2,N)*HSOL%FWD_SUMVEC(J,O2,K)
           ENDDO
           FWD_H1 = FWD_H1 + FWD_H2
          ENDDO
          HSOL%FWD_DIFVEC(I,O1,K) = FWD_H1 * SC_R
         ENDDO
        ENDDO

!  assign Forward solutions PHI (Siewert's notation)
!    --->  first N are "DOWN", last N are "UP" (streams)
!    --->  Use symmetry properties to set -ve eigensolutions

        DO I = 1, NSTREAMS
          XINV = HALF
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
            FPD = XINV * ( HSOL%FWD_SUMVEC(I,O1,K) + HSOL%FWD_DIFVEC(I,O1,K) )
            FND = XINV * ( HSOL%FWD_SUMVEC(I,O1,K) - HSOL%FWD_DIFVEC(I,O1,K) )
            HSOL%SOLA_XPOS(I,O1,K,N)  = FPD
            HSOL%SOLB_XNEG(I,O1,K,N)  = FND
            HSOL%SOLA_XPOS(I1,O1,K,N) = FND
            HSOL%SOLB_XNEG(I1,O1,K,N) = FPD
            IF ( O1 .GT. 2 ) THEN
              HSOL%SOLA_XPOS(I1,O1,K,N)  = - FND
              HSOL%SOLB_XNEG(I1,O1,K,N)  = - FPD
            ENDIF
          ENDDO
        ENDDO

!  end loop over real eigenvalues

      ENDDO

!  Complex solutions
!  -----------------

!  Offsets

      KO1 = HSOL%K_REAL(N) + 1

!  start loop over complex eigenvalues

      DO K = 1, HSOL%K_COMPLEX(N)

!  Bookkeeping indices

        K0 = 2 * ( K - 1 )
        K1 = KO1 + K0
        K2 = K1  + 1

!  get eigenvector entry

        AA  = EIGENMASK_C(K)
        AA1 = AA + 1

!  set eigenvalue and separation constants for one of pair of Conjugates
!   ---> (the complex conjugate follows from this information)

        RKSQ = REAL_KSQ(AA) * REAL_KSQ(AA)
        IKSQ = IMAG_KSQ(AA) * IMAG_KSQ(AA)
        MOD_KVAL = DSQRT(RKSQ + IKSQ)
        ARG = HALF * DATAN ( IMAG_KSQ(AA) / REAL_KSQ(AA) )
        CARG = DCOS(ARG)
        SARG = DSIN(ARG)

!  original complex --------------------------------------------------
!        KEIGEN_C(K,N) = DSQRT(MOD_KVAL) * CMPLX ( CARG, SARG )
!        SC_C  = ONE_C / KEIGEN_C(K,N)
!  -------------------------------------------------------------------

!  replacement complex

        KEIGEN_CSQ(K) = MOD_KVAL
        MODKRT        = DSQRT(MOD_KVAL)
        HSOL%KEIGEN(K1,N)  = MODKRT * CARG
        HSOL%KEIGEN(K2,N)  = MODKRT * SARG
        SC_CR =   HSOL%KEIGEN(K1,N) / MOD_KVAL
        SC_CI = - HSOL%KEIGEN(K2,N) / MOD_KVAL

!  set Right eigenvector norms
!   ---> Follows convention of way that eigenvectors are fomatted
!        upon output from the LAPACK module DGEEV.

        NORM_R = ZERO
        DO I = 1, NSTKS_NSTRMS
          XSQ_R = RITE_EVEC(I,AA)  * RITE_EVEC(I,AA)
          YSQ_R = RITE_EVEC(I,AA1) * RITE_EVEC(I,AA1)
          NORM_R = NORM_R + XSQ_R + YSQ_R
        ENDDO
        NORM_R = DSQRT(NORM_R)

!  Forward Sum-vector = normalized Right eigenvector

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            FWD_R = RITE_EVEC(IROW,AA) /NORM_R
            FWD_I = RITE_EVEC(IROW,AA1)/NORM_R
            HSOL%FWD_SUMVEC(I,O1,K1) = FWD_R
            HSOL%FWD_SUMVEC(I,O1,K2) = FWD_I
          ENDDO
        ENDDO

!  Find Forward difference vectors (Siewert's notation)

        DO I = 1, NSTREAMS
         DO O1 = 1, NSTOKES
          FWD_H1_CR = ZERO
          FWD_H1_CI = ZERO
          DO J = 1, NSTREAMS
           FWD_H2_CR = ZERO
           FWD_H2_CI = ZERO
           DO O2 = 1, NSTOKES
            FWD_H2_CR = FWD_H2_CR + HSOL%SAB(I,J,O1,O2,N) * HSOL%FWD_SUMVEC(J,O2,K1)
            FWD_H2_CI = FWD_H2_CI + HSOL%SAB(I,J,O1,O2,N) * HSOL%FWD_SUMVEC(J,O2,K2)
           ENDDO
           FWD_H1_CR = FWD_H1_CR + FWD_H2_CR
           FWD_H1_CI = FWD_H1_CI + FWD_H2_CI
          ENDDO
          HSOL%FWD_DIFVEC(I,O1,K1) = FWD_H1_CR*SC_CR - FWD_H1_CI*SC_CI
          HSOL%FWD_DIFVEC(I,O1,K2) = FWD_H1_CR*SC_CI + FWD_H1_CI*SC_CR
         ENDDO
        ENDDO

!  assign Forward solutions PHI (Siewert's notation)
!    --->  first N are "DOWN", last N are "UP" (streams)
!    --->  Use symmetry properties to set -ve eigensolutions

        DO I = 1, NSTREAMS
          XINV = HALF
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
            FPD1 = XINV * ( HSOL%FWD_SUMVEC(I,O1,K1) + HSOL%FWD_DIFVEC(I,O1,K1) )
            FPD2 = XINV * ( HSOL%FWD_SUMVEC(I,O1,K2) + HSOL%FWD_DIFVEC(I,O1,K2) )
            FND1 = XINV * ( HSOL%FWD_SUMVEC(I,O1,K1) - HSOL%FWD_DIFVEC(I,O1,K1) )
            FND2 = XINV * ( HSOL%FWD_SUMVEC(I,O1,K2) - HSOL%FWD_DIFVEC(I,O1,K2) )
            HSOL%SOLA_XPOS(I,O1,K1,N)   = FPD1
            HSOL%SOLB_XNEG(I,O1,K1,N)   = FND1
            HSOL%SOLA_XPOS(I,O1,K2,N)   = FPD2
            HSOL%SOLB_XNEG(I,O1,K2,N)   = FND2
            HSOL%SOLA_XPOS(I1,O1,K1,N)  = FND1
            HSOL%SOLB_XNEG(I1,O1,K1,N)  = FPD1
            HSOL%SOLA_XPOS(I1,O1,K2,N)  = FND2
            HSOL%SOLB_XNEG(I1,O1,K2,N)  = FPD2
            IF ( O1 .GT. 2 ) THEN
              HSOL%SOLA_XPOS(I1,O1,K1,N)  = - FND1
              HSOL%SOLB_XNEG(I1,O1,K1,N)  = - FPD1
              HSOL%SOLA_XPOS(I1,O1,K2,N)  = - FND2
              HSOL%SOLB_XNEG(I1,O1,K2,N)  = - FPD2
            ENDIF
          ENDDO
        ENDDO

!  end eigenvalue loop

      ENDDO

!  control point for avoiding eigensolver
!  ======================================

 3456 CONTINUE

      HSOL%T_DELT_EIGEN(:,N) = ZERO

!  Eigenstream transmittance factors for whole layer
!  -------------------------------------------------

      DO K = 1, HSOL%K_REAL(N)
         HELP = HSOL%KEIGEN(K,N) * DELTAU_VERT(N)
         IF ( HELP .LE. MAX_TAU_QPATH ) HSOL%T_DELT_EIGEN(K,N) = EXP ( - HELP)
      ENDDO
      DO K = 1, HSOL%K_COMPLEX(N)
         K0 = 2 * K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
         HELP = DELTAU_VERT(N) * HSOL%KEIGEN(K1,N)
         IF ( HELP .LE. MAX_TAU_QPATH ) THEN
           HELP_CR = EXP ( - HELP)
           HELP_CI = DELTAU_VERT(N) * HSOL%KEIGEN(K2,N)
           HSOL%T_DELT_EIGEN(K1,N) =   HELP_CR * COS ( HELP_CI )
           HSOL%T_DELT_EIGEN(K2,N) = - HELP_CR * SIN ( HELP_CI )
         ENDIF
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE MSSEG_QHOM_SOLUTION

!

      SUBROUTINE MSSEG_HOMOG_PIMATRIX_SETUP ( &
        FOURIER, NSTOKES, NSTREAMS, NMOMENTS, QUAD_STREAMS, DMAT, HGSF )


!  Notes
!  -----

!  This is equivalent to the Legendre setup modules in LIDORT

!---------old comments ---------------------------------
!  This needs modification for the case of Refractive atmosphere,
!  because the solar zenith angle is not constant.
!-------------------------------------------------------

!  PI-matrix setup, following recipe of Siewert (1982).
!  Single Fourier component only.

!  Tested against benchmark results in Vestrucci & Siewert (1984), Probl

!  original coding, September 2002, R. Spurr SAO
!  Multibeam SZA coding, July 2004, R. Spurr SAO
!  Coding for refractive geometry case, R. Spurr, RT Solutions, May 2005

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need

      USE MSSeg_PARS_m, only : MAXSTOKES, MAXSTREAMS, MAXMOMENTS, ZERO, HALF, ONE, TWO, THREE

!  Mk3. Introduce the type structure

      use MSSeg_HPTypes_Def_m

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NMOMENTS

      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STREAMS ( MAXSTREAMS )

!  Indices and Dmatrix from DERIVE_INPUTS (mick check)

      DOUBLE PRECISION, INTENT (IN)  :: DMAT ( MAXSTOKES, MAXSTOKES )

!  Output GSFs for homogeneous solutions
!  -------------------------------------

!  Must be intent(inout) as some quantities are saved only from the Fourier 0 calculations

      type(Homog_GSF), intent(inout) :: HGSF

!  Local matrices
!  --------------

      DOUBLE PRECISION :: PI ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

!  local variables

!     INTEGER          :: LOCAL_NSTOKES
      DOUBLE PRECISION :: FAC, XSQ, DF_LPM, DF_LP1MM, ZLM_VALUE
      DOUBLE PRECISION :: HX, HY, HZ, H1, H2, RCONST_20, RCONST_21
      DOUBLE PRECISION :: PL20, RL20, PL21, RL21, TL21, XM(MAXSTREAMS)
      INTEGER          :: M, I, L, SK1, SK2, SK3, SK4, NS, M1

!  Set integer M = Fourier number

      M = FOURIER

!  total number of angles

      NS  = NSTREAMS

!  constants

      RCONST_20 = 0.25D0 * DSQRT(6.0D0)
      RCONST_21 = 2.0D0 * RCONST_20

!  Coefficient matrix PINORM (for normalization) = [ (L-m)!/(L+m)! ] ^1/
!  ---------------------------------------------------------------------

!  .. first entry = 1 / (2m)!    ---- be careful with overflow

      HGSF%PHI(M,M) = ONE
      FAC = ONE
      DO L = 1, 2*M
        FAC = FAC * DBLE(L)
      ENDDO
      HGSF%PHI(M,M) = HGSF%PHI(M,M) / FAC

!  .. Other entries by recurrence

      DO L = M + 1, NMOMENTS
        FAC = DBLE(L-M)/DBLE(L+M)
        HGSF%PHI(L,M) = HGSF%PHI(L-1,M)*FAC
      ENDDO

!  .. Square root

      DO L = M , NMOMENTS
        HGSF%PINORM(L,M) = SQRT(HGSF%PHI(L,M))
      ENDDO

!  Additional saved quantities
!  ---------------------------

!  Only for the first fundamental harmonic

      IF ( M .EQ. 0 ) THEN

!  Sign matrix

        DO M1 = 0, NMOMENTS
          DO L = M1, NMOMENTS
            IF (MOD((L-M1),2).EQ.0) THEN
              HGSF%PISIGN(L,M1) = ONE
            ELSE
              HGSF%PISIGN(L,M1) = -ONE
            ENDIF
          ENDDO
        ENDDO

!  floating point integer values

        DO L = 0, NMOMENTS
          HGSF%DF_L(L)    = DBLE(L)
          HGSF%DF_LP1(L)  = DBLE(L+1)
          HGSF%DF_2LP1(L) = DBLE(2*L+1)
        ENDDO

        HGSF%DF_LSQM4(2) = ZERO
        DO L = 3, NMOMENTS
          HGSF%DF_LSQM4(L) = DSQRT(DBLE(L*L-4))
        ENDDO

        DO L = 2, NMOMENTS
          HGSF%DF_RT_LP3XLM1(L) = DSQRT(DBLE((L+3)*(L-1)))
          HGSF%ZHELP(L) = TWO*HGSF%DF_2LP1(L)/HGSF%DF_LP1(L)/HGSF%DF_L(L)
        ENDDO

      ENDIF

!  local array of all streams (every Fourier component)

     DO I = 1, NSTREAMS
       XM(I) = QUAD_STREAMS(I)
     ENDDO

!  factors associated with stream values
!   Special case when XM(I) = ONE, factors are zeroed, and
!    those that give singularities are avoided later (see below)
!   R. Spurr and V. Natraj, 16 january 2006

      IF ( M .EQ. 0 )  THEN
        DO I = 1, NS
          XSQ = XM(I) * XM(I)
          HGSF%PLEG20(I)        = HALF*(THREE*XSQ-ONE)
          HGSF%UMXSQ(I)         = ONE - XSQ
          HGSF%XUMXSQ(I)        = HGSF%UMXSQ(I) * XM(I)
          IF ( XM(I) .EQ. ONE ) THEN
           HGSF%RT_UMXSQ(I)      = ZERO
           HGSF%X_RT_UMXSQ(I)    = ZERO
           HGSF%UPXSQ_D_UMXSQ(I) = ZERO
           HGSF%M2X_D_UMXSQ(I)   = ZERO
          ELSE
           HGSF%RT_UMXSQ(I)      = DSQRT(HGSF%UMXSQ(I))
           HGSF%X_RT_UMXSQ(I)    = HGSF%RT_UMXSQ(I) * XM(I)
           HGSF%UPXSQ_D_UMXSQ(I) = (ONE + XSQ) / HGSF%UMXSQ(I)
           HGSF%M2X_D_UMXSQ(I)   = - TWO * XM(I) / HGSF%UMXSQ(I)
          ENDIF
          DO L = 0, NMOMENTS
            HGSF%X2LP1(L,I) = XM(I) * HGSF%DF_2LP1(L)
          ENDDO
        ENDDO
      ENDIF

!  XYZ matrices
!  ------------

!  Inverse XLM diagonal matrices (Siewert (1982), Eq. 35a)
!  YLM diagonal matrices (Siewert (1982), Eq. 35b)
!  ZLM matrices. (Siewert (1982), Eq. 35c)

!  .. for the azimuth-independent harmonic

      IF ( M .EQ. 0 ) THEN

        DO L = 2, NMOMENTS
          HGSF%XLM_DIAG(L,1) = HGSF%DF_LP1(L)
          HGSF%XLM_DIAG(L,2) = HGSF%DF_RT_LP3XLM1(L)
          HGSF%XLM_DIAG(L,3) = HGSF%XLM_DIAG(L,2)
          HGSF%XLM_DIAG(L,4) = HGSF%XLM_DIAG(L,1)
          DO SK1 = 1, NSTOKES
            HGSF%XLM_DIAG(L,SK1) = ONE/HGSF%XLM_DIAG(L,SK1)
          ENDDO
        ENDDO

        DO L = 2, NMOMENTS
          HGSF%YLM_DIAG(L,1) = HGSF%DF_L(L)
          HGSF%YLM_DIAG(L,2) = HGSF%DF_LSQM4(L)
          HGSF%YLM_DIAG(L,3) = HGSF%YLM_DIAG(L,2)
          HGSF%YLM_DIAG(L,4) = HGSF%YLM_DIAG(L,1)
        ENDDO

!  .. for the other harmonics

      ELSE

        DO L = M, NMOMENTS
          DF_LP1MM = DBLE ( L + 1 - M )
          HGSF%XLM_DIAG(L,1) = DF_LP1MM
          HGSF%XLM_DIAG(L,1) = ONE/HGSF%XLM_DIAG(L,1)
          IF ( L .EQ. 1 ) THEN
            HGSF%XLM_DIAG(L,2) = ZERO
          ELSE
            HGSF%XLM_DIAG(L,2) = HGSF%DF_RT_LP3XLM1(L) * DF_LP1MM / HGSF%DF_LP1(L)
            HGSF%XLM_DIAG(L,2) = ONE / HGSF%XLM_DIAG(L,2)
          ENDIF
          HGSF%XLM_DIAG(L,3) = HGSF%XLM_DIAG(L,2)
          HGSF%XLM_DIAG(L,4) = HGSF%XLM_DIAG(L,1)
        ENDDO

        DO SK1 = 1, NSTOKES
          HGSF%YLM_DIAG(M,SK1) = ZERO
        ENDDO
        DO L = M + 1, NMOMENTS
          DF_LPM = DBLE ( L + M )
          HGSF%YLM_DIAG(L,1) = DF_LPM
          HGSF%YLM_DIAG(L,2) = DF_LPM * HGSF%DF_LSQM4(L) / HGSF%DF_L(L)
          HGSF%YLM_DIAG(L,3) = HGSF%YLM_DIAG(L,2)
          HGSF%YLM_DIAG(L,4) = HGSF%YLM_DIAG(L,1)
        ENDDO

        DO L = 2, NMOMENTS
          ZLM_VALUE = DBLE(M) * HGSF%ZHELP(L)
          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              HGSF%ZLM(L,SK1,SK2) = ZERO
            ENDDO
          ENDDO
          HGSF%ZLM(L,2,3) = ZLM_VALUE
          HGSF%ZLM(L,3,2) = ZLM_VALUE
        ENDDO

      ENDIF

!  PI calculation
!  --------------
!  Initialise all the PI matrices for given harmonic
!mick fix 6/6/2023 - applied alternate PI initialization since Siewert eqs
!                    below automatically define some PI elements > NSTOKES

      !DO I = 1, NS
      !  DO L = M, NMOMENTS
      !    DO SK1 = 1, NSTOKES
      !      DO SK2 = 1, NSTOKES
      !        PI(L,I,SK1,SK2) = ZERO
      !      ENDDO
      !    ENDDO
      !  ENDDO
      !ENDDO
      PI(M:NMOMENTS,1:NS,1:MAXSTOKES,1:MAXSTOKES) = ZERO

!  M = 0 component. [Siewert (1982), Eqs (28) to (31)]

      IF ( M .EQ. 0 ) THEN

!  .. L = 0
        DO I = 1, NS
          PI(0,I,1,1) = PI(0,I,1,1) + ONE
          PI(0,I,4,4) = PI(0,I,4,4) + ONE
        ENDDO
!  .. L = 1
        DO I = 1, NS
          PI(1,I,1,1) = PI(1,I,1,1) + XM(I)
          PI(1,I,4,4) = PI(1,I,4,4) + XM(I)
        ENDDO
!  .. L = 2
        DO I = 1, NS
          PL20 = HGSF%PLEG20(I)
          RL20 = RCONST_20 * HGSF%UMXSQ(I)
          PI(2,I,1,1) = PI(2,I,1,1) + PL20
          PI(2,I,2,2) = PI(2,I,2,2) + RL20
          PI(2,I,3,3) = PI(2,I,2,2)
          PI(2,I,4,4) = PI(2,I,1,1)
        ENDDO
!  .. L > 2
        DO L = 2, NMOMENTS - 1
          DO I = 1, NS
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = HGSF%X2LP1(L,I)*PI(L,I,SK1,SK2)
                HY = HGSF%YLM_DIAG(L,SK1)*PI(L-1,I,SK1,SK2)
                PI(L+1,I,SK1,SK2) = ( HX - HY ) * HGSF%XLM_DIAG(L,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  M = 1 component. [Siewert (1982), Eqs (32) to (34), plus (35)]

      ELSE IF ( M .EQ. 1 ) THEN

!  .. L = 1
        DO I = 1, NS
          PI(1,I,1,1) = PI(1,I,1,1) + HGSF%RT_UMXSQ(I)
          PI(1,I,4,4) = PI(1,I,4,4) + HGSF%RT_UMXSQ(I)
        ENDDO
!  .. L = 2
        DO I = 1, NS
          PL21 =   THREE     * HGSF%X_RT_UMXSQ(I)
          RL21 = - RCONST_21 * HGSF%X_RT_UMXSQ(I)
          TL21 = + RCONST_21 * HGSF%RT_UMXSQ(I)
          PI(2,I,1,1) = PI(2,I,1,1) + PL21
          PI(2,I,2,2) = PI(2,I,2,2) + RL21
          PI(2,I,2,3) = PI(2,I,2,3) + TL21
          PI(2,I,3,3) = PI(2,I,2,2)
          PI(2,I,4,4) = PI(2,I,1,1)
          PI(2,I,3,2) = PI(2,I,2,3)
        ENDDO
!  .. L > 2
        DO L = 2, NMOMENTS - 1
          DO I = 1, NS
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = HGSF%X2LP1(L,I)*PI(L,I,SK1,SK2)
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

      ELSE

!  .. L = M
        IF ( M .EQ. 2 ) THEN
          DO I = 1, NS
            HGSF%PIMM_11(I) =   THREE     * HGSF%UMXSQ(I)
            IF ( XM(I).EQ.ONE) THEN
              HGSF%PIMM_KM(I) =   RCONST_21
            ELSE
              HGSF%PIMM_KM(I) =   RCONST_21 * HGSF%UMXSQ(I)
            ENDIF
          ENDDO
        ELSE
          H1 = HGSF%DF_2LP1(M-1)
          H2 = HGSF%DF_LP1(M-1)*H1/HGSF%DF_RT_LP3XLM1(M-1)
          DO I = 1, NS
            HGSF%PIMM_11(I) = H1 * HGSF%RT_UMXSQ(I) * HGSF%PIMM_11(I)
            IF ( XM(I).EQ.ONE) THEN
              HGSF%PIMM_KM(I) = ZERO
            ELSE
              HGSF%PIMM_KM(I) = H2 * HGSF%RT_UMXSQ(I) * HGSF%PIMM_KM(I)
            ENDIF
          ENDDO
        ENDIF

        DO I = 1, NS
          PI(M,I,1,1) = PI(M,I,1,1) + HGSF%PIMM_11(I)
          IF ( XM(I).EQ.ONE ) THEN
            PI(M,I,2,2) = PI(M,I,2,2) + HGSF%PIMM_KM(I) * TWO
            PI(M,I,2,3) = PI(M,I,2,3) - HGSF%PIMM_KM(I) * TWO
          ELSE
            PI(M,I,2,2) = PI(M,I,2,2) + HGSF%PIMM_KM(I) * HGSF%UPXSQ_D_UMXSQ(I)
            PI(M,I,2,3) = PI(M,I,2,3) + HGSF%PIMM_KM(I) * HGSF%M2X_D_UMXSQ(I)
          ENDIF
          PI(M,I,3,3) = PI(M,I,2,2)
          PI(M,I,4,4) = PI(M,I,1,1)
          PI(M,I,3,2) = PI(M,I,2,3)
        ENDDO
!  .. L = M + 1
        IF ( M .LT. NMOMENTS ) THEN
          DO I = 1, NS
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = HGSF%X2LP1(M,I)*PI(M,I,SK1,SK2)
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
          DO I = 1, NS
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = HGSF%X2LP1(L,I)*PI(L,I,SK1,SK2)
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

!  Normalized output
!  -----------------

      DO L = M, NMOMENTS

!  .. at quadrature streams

        DO I = 1, NSTREAMS

          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              HGSF%PI_XQP(L,I,SK1,SK2) = PI(L,I,SK1,SK2) * HGSF%PINORM(L,M)
            ENDDO
          ENDDO

          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES

              H1 = ZERO
              H2 = ZERO
              DO SK3 = 1, NSTOKES
                H1 = H1 + DMAT(SK1,SK3) * HGSF%PI_XQP(L,I,SK3,SK2)
                H2 = H2 + HGSF%PI_XQP(L,I,SK1,SK3) * DMAT(SK3,SK2)
              ENDDO

              HGSF%PI_XQP_PRE(L,I,SK1,SK2)  = H1
              HGSF%PI_XQM_PRE(L,I,SK1,SK2)  = H1 * HGSF%PISIGN(L,M)
              HGSF%PI_XQM_POST(L,I,SK1,SK2) = H2 * HGSF%PISIGN(L,M)

              H1 = ZERO
              DO SK3 = 1, NSTOKES
                H2 = ZERO
                DO SK4 = 1, NSTOKES
                  H2 = H2 + HGSF%PI_XQP(L,I,SK3,SK4)*DMAT(SK4,SK2)
                ENDDO
                H1 = H1 + DMAT(SK1,SK3)*H2
              ENDDO
              HGSF%PI_XQM(L,I,SK1,SK2)  = H1 * HGSF%PISIGN(L,M)

            ENDDO
          ENDDO

        ENDDO

!  end loop moments

      ENDDO

!  debug

!      write(77,*)'Fourier',M
!      DO L = M, NMOMENTS
!       write(77,*)'Moment ',L
!       DO SK1 = 1, 4
!        WRITE(77,'(I3,1p4e15.6)')SK1,(PI_XUP(L,4,SK1,SK2),SK2=1,4)
!       ENDDO
!      ENDDO

!  Finish

      RETURN
      END SUBROUTINE MSSEG_HOMOG_PIMATRIX_SETUP

!

      END MODULE MSSeg_Homog_Solutions_m

