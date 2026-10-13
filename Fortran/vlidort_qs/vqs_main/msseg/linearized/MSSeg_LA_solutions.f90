
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

!  Mk2a. 9/21/21. Created LC routines

!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- ASYMTX_TOLERANCE now introduced and passed into the QHOM call to ASYMTX
!   -- Flux output (LOCAL_FLUXES) controlled by DO_LOCAL_FLUXES input flag
!   -- DO_FOURIER0_NSTOKES2 performance flag introduced.
!   -- Call subroutine specially written to get the discrete-ordinate field for Fluxes
!   -- Use public routine GETQUAD2 instead of GAULEG

!  Mk3. 1/31/22. Only need LA_SOLUTIONS for this module

      MODULE MSSeg_LA_Solutions_m

      PRIVATE
      PUBLIC :: MSSEG_LA_QHOM_SOLUTION

      CONTAINS

      SUBROUTINE MSSEG_LA_QHOM_SOLUTION ( &
        GIVEN_LAYER, FOURIER, DO_VARY, N_PARAMETERS, NSTOKES,   & ! Input control
        NSTREAMS, NMOMENTS, NSTKS_NSTRMS, NSTKS_NSTRMS_2,       & ! Input Control/bookkeeping
        QUAD_STREAMS, QUAD_HALFWTS, DO_REAL_EIGENSOLVER,        & ! Input bookkeeping
        DELTAU_VERT, OMEGA_GREEK, L_DELTAU_VERT, L_OMEGA_GREEK, & ! Input Optical
        HGSF, HSOL, REAL_KSQ, IMAG_KSQ, LEFT_EVEC, RITE_EVEC,   & ! Input Homog Solutions
        EIGENDEGEN, EIGENMASK_R, EIGENMASK_C, KEIGEN_CSQ,       & ! Input Homog Solutions 
        L_HSOL, STATUS, MESSAGE, TRACE )                          ! Output TS L_Homog, Exceptions 

!  Linearization of the homogeneous solutions.

!  Mk2, 6/21/21. Use module names now have the "_m" extension, Only use quantities you need
!  Mk3. 1/31/22. Use type structure inputs, reaarranged argument list
!  Mk3. 4/11/23. Resurrected and modified for practical use

      USE MSSEG_PARS_m, Only : MAXSTREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS, MAXMOMENTS, &
                               MAXSTREAMS_2, MAXSTRMSTKS, MAXSTRMSTKS_2, MAXSTRMSTKS_21, &
                               MAXSTRMSTKS_P1, MAXEVALUES, &
                               VQS_SUCCESS, VQS_SERIOUS, ZERO, HALF, ONE, TWO

      USE MSSEG_HPTypes_Def_m
      USE MSSEG_L_HPTypes_Def_m

      USE LAPACK_TOOLS_m, Only : DGETRS, DGETRF

      IMPLICIT NONE

!  inputs
!  ------

!  Control

      INTEGER, INTENT (IN) ::           GIVEN_LAYER
      INTEGER, INTENT (IN) ::           FOURIER
      LOGICAL, INTENT (IN) ::           DO_VARY
      INTEGER, INTENT (IN) ::           N_PARAMETERS
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NMOMENTS

!  bookkeeping

      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS_2
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_HALFWTS ( MAXSTREAMS )
      LOGICAL, INTENT (IN) ::           DO_REAL_EIGENSOLVER ( 0:MAXMOMENTS, MAXLAYERS )

!  optical

      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_VERT   ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_GREEK   ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  L_OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  Input type structure (Generalized spherical functions)

      type(Homog_GSF), intent(inout) :: HGSF

!  Input solution type structure

      type(Homog_Solutions), intent(inout) :: HSOL

!  other homogeneous solution inputs

      DOUBLE PRECISION, INTENT (INOUT) :: REAL_KSQ ( MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) :: IMAG_KSQ ( MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) :: LEFT_EVEC ( MAXSTRMSTKS, MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) :: RITE_EVEC ( MAXSTRMSTKS, MAXSTRMSTKS )
      LOGICAL, INTENT (INOUT) ::          EIGENDEGEN ( MAXSTRMSTKS, MAXLAYERS )
      INTEGER, INTENT (INOUT) ::          EIGENMASK_R ( MAXEVALUES )
      INTEGER, INTENT (INOUT) ::          EIGENMASK_C ( MAXEVALUES )
      DOUBLE PRECISION, INTENT (INOUT) :: KEIGEN_CSQ ( MAXEVALUES )

!  Output solution type structure

      type(L_Homog_Solutions), intent(inout) :: L_HSOL

!  Exception handling

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  local matrices for eigenvalue computation
!  -----------------------------------------

!  For the real calculation based on DGEEV

      DOUBLE PRECISION :: HMAT1(MAXSTRMSTKS,MAXSTRMSTKS), HVEC1(MAXSTRMSTKS,MAX_ATMOSWFS)
      DOUBLE PRECISION :: HVEC_AUX(MAXSTRMSTKS,MAX_ATMOSWFS)

!  For the complex calculations based on DGEEV

      DOUBLE PRECISION :: HMAT2   (MAXSTRMSTKS_2,MAXSTRMSTKS_2)  , HVEC2   (MAXSTRMSTKS_2,MAX_ATMOSWFS)
      DOUBLE PRECISION :: HMATRED (MAXSTRMSTKS_21,MAXSTRMSTKS_21), HVECRED (MAXSTRMSTKS_21,MAX_ATMOSWFS)
      DOUBLE PRECISION :: HVEC_CR (MAXSTRMSTKS,MAX_ATMOSWFS)     , HVEC_CI (MAXSTRMSTKS,MAX_ATMOSWFS)
      DOUBLE PRECISION :: ADJN_CR(MAXSTRMSTKS), ADJN_CI(MAXSTRMSTKS)

!  For the real calculation based on ASYMTX

      DOUBLE PRECISION :: HMAT1A(MAXSTRMSTKS_P1,MAXSTRMSTKS_P1), HVEC1A(MAXSTRMSTKS_P1,MAX_ATMOSWFS)

!  Local linearized vectors

      DOUBLE PRECISION :: L_FWD_DIFVEC_R (MAXSTREAMS,MAXSTOKES), L_FWD_SUMVEC_R( MAXSTREAMS,MAXSTOKES,MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_FWD_DIFVEC_CR(MAXSTREAMS,MAXSTOKES), L_FWD_SUMVEC_CR(MAXSTREAMS,MAXSTOKES)
      DOUBLE PRECISION :: L_FWD_DIFVEC_CI(MAXSTREAMS,MAXSTOKES), L_FWD_SUMVEC_CI(MAXSTREAMS,MAXSTOKES)

!  Pivoting arrays for LAPACK linear-algebra solvers

      INTEGER ::          LAPACK_IPIV1(MAXSTRMSTKS)
      INTEGER ::          LAPACK_IPIV21(MAXSTRMSTKS_21)
      INTEGER ::          LAPACK_INFO, NSTKS_NSTRMS_21, NSTMSTKS_P1
      INTEGER ::          LAPACK_IPIV1A(MAXSTRMSTKS_P1)

!  Miscellaneous local variables

      INTEGER ::          I, J, I1, J1, JC, IR, JCOL, IROW, IROW1, ILAST
      INTEGER ::          L, N, M, Q, O1, O2, O3, O4
      INTEGER ::          K, KO1, K0, K1, K2
      DOUBLE PRECISION :: DP, DM, KSQD, KVAL, KVAL2, XINV, FAC
      DOUBLE PRECISION :: H1P, H1M, H2P, H2M, H1, H2, L_OG
      DOUBLE PRECISION :: TBAS, TBAS_CR, TBAS_CI, LTQ, LCR, LCI

!  additional variables

      INTEGER           :: AA, AA1, TC, TR, I_NONTRIVIAL
      DOUBLE PRECISION  :: SR, SCR, SCI, HR, HCR, HCI, MODULUS
      DOUBLE PRECISION  :: KVAL_CR, KVAL_CI, KVAL2_CR, KVAL2_CI
      DOUBLE PRECISION  :: MOD_KVAL, SEPCON_CR, SEPCON_CI
      DOUBLE PRECISION  :: DENOM_CR, DENOM_CI, LFPD, LFND
      DOUBLE PRECISION  :: LFPD1, LFPD2, LFND1, LFND2
      CHARACTER (LEN=3) :: CN, CI

!  Start of code
!  -------------

!  initialise exception handling

      STATUS  = VQS_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Layer and Fourier number

      N = GIVEN_LAYER
      M = FOURIER

!  nothing to do if not varying

      IF ( .NOT. DO_VARY ) RETURN

!  Linearize the Eigenmatrix
!  =========================

!  set up linearizations of SAB and DAB. All real variables.
!   Thoroughly debugged, November 2005

      DO I = 1, NSTREAMS
       XINV = - ONE / QUAD_STREAMS(I)
       DO J = 1, NSTREAMS
        FAC =  QUAD_HALFWTS(J) * XINV
        DO Q = 1, N_PARAMETERS
         DO O1 = 1, NSTOKES
          DO O2 = 1, NSTOKES
           DP = ZERO
           DM = ZERO
           DO L = M, NMOMENTS
            H1P = ZERO
            H1M = ZERO
            DO O3 = 1, NSTOKES
             H2P = ZERO
             H2M = ZERO
             DO O4 = 1, NSTOKES
              L_OG = L_OMEGA_GREEK(L,N,O3,O4,Q)
              H2P = H2P + L_OG * HGSF%PI_XQP(L,J,O4,O2)
              H2M = H2M + L_OG * HGSF%PI_XQM_PRE(L,J,O4,O2)
             ENDDO
             H1P = H1P + HGSF%PI_XQP(L,I,O1,O3)*H2P
             H1M = H1M + HGSF%PI_XQP(L,I,O1,O3)*H2M
            ENDDO
            DP = DP + H1P
            DM = DM + H1M
           ENDDO
           L_HSOL%L_SAB(I,J,O1,O2,N,Q) = FAC * ( DP + DM )
           L_HSOL%L_DAB(I,J,O1,O2,N,Q) = FAC * ( DP - DM )
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

!  set up linearized eigenmatrices, saved for all layers. Real variables
!   Thoroughly debugged, November 2005

      DO Q = 1, N_PARAMETERS
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
             H2 = H2 + L_HSOL%L_DAB(I,K,O1,O3,N,Q) *   HSOL%SAB  (K,J,O3,O2,N) &
                     +   HSOL%DAB  (I,K,O1,O3,N)   * L_HSOL%L_SAB(K,J,O3,O2,N,Q)
            ENDDO
            H1 = H1 + H2
           ENDDO
           L_HSOL%L_EIGENMAT(IROW,JCOL,N,Q) = H1
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Linearization of Real Eigensolutions
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  start eigenvalue loop

      DO K = 1, HSOL%K_REAL(N)

!  some basic quantities

        AA    = EIGENMASK_R(K)
        KVAL  = HSOL%KEIGEN(K,N)
        KVAL2 = 2.0D0 * KVAL
        KSQD  = KVAL * KVAL

!   REAL SOLUTIONS ONLY (ASYMTX linearization)
!   ==========================================

       IF ( DO_REAL_EIGENSOLVER(M,N) ) THEN

!  Count the number of entries in the eigenvector ( = TC)
!    If the eigenvector is trivial (one entry, TC = 1).....
!    ......Zero the results and skip linearization for this eigenvalue

        TC = 0
        DO JCOL = 1,  NSTKS_NSTRMS
          if ( RITE_EVEC(JCOL,AA).NE.ZERO ) TC = TC + 1
        ENDDO
!Rob fix 4/12/12
        IF ( TC.EQ.1 .and. NSTKS_NSTRMS.gt.1 ) THEN
          DO Q = 1, N_PARAMETERS
            HVEC1A(1,Q) = ZERO
            DO I = 1, NSTKS_NSTRMS
              HVEC1A(I+1,Q) = ZERO
            ENDDO
          ENDDO
          GO TO 4555
        ENDIF

!  Nothing to do also if there is degeneracy
!    NEW CODE by R. Spurr 22 August 2008.

        IF ( EIGENDEGEN(K,N) ) THEN
          DO Q = 1, N_PARAMETERS
            HVEC1A(1,Q) = ZERO
            DO I = 1, NSTKS_NSTRMS
              HVEC1A(I+1,Q) = ZERO
            ENDDO
          ENDDO
          GO TO 4555
        ENDIF

!  initialise solution matrix HMAT1A (important to do this!)
!    HMAT1A is destroyed upon passage through the LAPACK linear modules

        DO I = 1, MAXSTRMSTKS_P1
!          LAPACK_IPIV1A(I) = 0
          DO J = 1, MAXSTRMSTKS_P1
            HMAT1A(I,J) = ZERO
          ENDDO
        ENDDO

!  Determine solution matrix HMAT (A in the matrix equation A.x = B)

        NSTMSTKS_P1 = NSTKS_NSTRMS + 1
        DO I = 1, NSTKS_NSTRMS
          DO J = 1, NSTKS_NSTRMS
            HMAT1A(I,J+1) = - HSOL%EIGENMAT_SAVE(I,J,N)
          ENDDO
          HMAT1A(I,I+1) = HMAT1A(I,I+1) + KSQD
          HMAT1A(I,1)   = TWO * KVAL * RITE_EVEC(I,AA)
          HMAT1A(NSTMSTKS_P1,I+1) = RITE_EVEC(I,AA)
        ENDDO
        HMAT1A(NSTMSTKS_P1,1) = ZERO

!        HMAT1A(NSTMSTKS_P1,1) = 1.879d-15   ! shoudl technically be zer

!  solution column vectors (B in the matrix equation A.x = B).
!    (One for each parameter to be varied)

        DO Q = 1, N_PARAMETERS
          DO I = 1, NSTKS_NSTRMS
            H1 = ZERO
            DO J = 1, NSTKS_NSTRMS
              H1 = H1 + L_HSOL%L_EIGENMAT(I,J,N,Q) * RITE_EVEC(J,AA)
            ENDDO
            HVEC1A(I,Q) = H1
          ENDDO
          HVEC1A(NSTMSTKS_P1,Q) = ZERO
        ENDDO

!  Find any degenerate columns J, and replace the (J,J) entry with
!   a small number. This is a quick fix that ensures that the LAPACK
!   routine does not fall over because of degeneracy.
!    --This problem is only present in the vector case.

!  Fiddle code

!        DO J = 1, NSTKS_NSTRMS
!          NC = 0
!          DO I = 1,  NSTKS_NSTRMS
!           IF ( HMAT1A(I,J+1).EQ.ZERO) NC = NC + 1
!          ENDDO
!          IF ( NC .EQ. NSTKS_NSTRMS ) THEN
!            do I = 1, NSTKS_NSTRMS
!               HMAT1A(J,I+1) = 1.237d-15*DBLE(i)
!            enddo
!            if ( HMAT1A(J+1,J+1).EQ.ZERO)
!     &    HMAT1A(J+1,J+1) = 1.056d-15
!          ENDIF
!        ENDDO

!  ultimate fiddle

!        DO J = 1, NSTMSTKS_P1
!           Do I = 1, NSTMSTKS_P1
!             IF (HMAT1A(I,J).EQ.ZERO )HMAT1A(I,J) = 1.0d-15*(I+J)
!           ENDDO
!        ENDDO

!  solve System for the linearized sum-vector + eigenvalue, together
!     Solve matrix system  using LAPACK modules DGETRF and DEGTRS

!        IF (m.eq.2.and.n.eq.18) THEN
!           WRITE(45, *) m, k, N, NSTMSTKS_P1, LAPACK_INFO
!           WRITE(45, '(17D24.12)')
!     &          ((HMAT1A(j, i), i = 1, NSTMSTKS_P1), j = 1, NSTMSTKS_P1
!c           if(lapack_info.ne.0)STOP'2,18,9'
!        ENDIF

!#######################################################################

!  Resume normal code

        CALL DGETRF  ( NSTMSTKS_P1, NSTMSTKS_P1, HMAT1A, &
                       MAXSTRMSTKS_P1,  LAPACK_IPIV1A, LAPACK_INFO )

!  Exception handling 1

        IF ( LAPACK_INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) LAPACK_INFO
          WRITE(CN, '(I3)' ) N
          MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE   = 'DGETRF call # 1 in MSSEG_LA_QHOM_SOLUTION,'// &
                      ' layer # '//CN
          STATUS  = VQS_SERIOUS
          RETURN
        ELSE IF ( LAPACK_INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) LAPACK_INFO
          WRITE(CN, '(I3)' ) N
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRF call # 1 in MSSEG_LA_QHOM_SOLUTION,'// &
                      ' layer # '//CN
          STATUS  = VQS_SERIOUS
          RETURN
        ENDIF

        CALL DGETRS ('N', NSTMSTKS_P1, N_PARAMETERS, HMAT1A, &
                     MAXSTRMSTKS_P1, LAPACK_IPIV1A,  HVEC1A, &
                     MAXSTRMSTKS_P1, LAPACK_INFO )

!  Exception handling 1

        IF ( LAPACK_INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) LAPACK_INFO
          WRITE(CN, '(I3)' ) N
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRS call # 1 in MSSEG_LA_QHOM_SOLUTION,'// &
                      ' layer # '//CN
          STATUS  = VQS_SERIOUS
          RETURN
        ENDIF

!  Continuation point if linearized skipped

 4555   CONTINUE

!  Start loop over varying parameters, to assign answers

        DO Q = 1, N_PARAMETERS

          L_HSOL%L_KEIGEN(K,N,Q) = HVEC1A(1,Q)
          DO I = 1, NSTREAMS
           IR = NSTOKES*(I-1)
           DO O1 = 1, NSTOKES
            IROW = IR + O1
            L_FWD_SUMVEC_R(I,O1,Q) = HVEC1A(IROW+1,Q)
           ENDDO
          ENDDO

!  Debug (important)
!   Tested 28 December 2005.

!         IF ( DO_DEBUG_WRITE ) THEN ; IF (Q.EQ.1.and.m.eq.2.and.n.eq.1) THEN
!           write(37,'(2i3,L3,1p2e20.12)') M,K,EIGENDEGEN(K,N),KEIGEN(K,N),L_KEIGEN(K,N,Q)
!           DO I = 1, NSTREAMS ; DO O1 = 1, NSTOKES
!              write(37,'(4i3,1p2e20.12)') m,k,i,o1,RITE_EVEC(I,k ),L_FWD_SUMVEC_R(I,O1,Q)
!           ENDDO ; ENDDO
!         ENDIF ; ENDIF

        ENDDO

!   REAL SOLUTIONS (DGEEV linearization)
!   ====================================

       ELSE

!  initialise

       I_NONTRIVIAL = 0

!  Count the number of entries in the eigenvector ( = TC)
!    If the eigenvector is trivial (one entry, TC = 1).....
!    ......Zero the results and skip linearization for this eigenvalue

        TC = 0
        DO JCOL = 1,  NSTKS_NSTRMS
          if ( RITE_EVEC(JCOL,AA).NE.ZERO ) TC = TC + 1
        ENDDO
        IF ( TC .EQ. 1 ) THEN
          DO Q = 1, N_PARAMETERS
            L_HSOL%L_KEIGEN(K,N,Q)  = ZERO
            DO I = 1, NSTKS_NSTRMS
              HVEC1(I,Q) = ZERO
            ENDDO
          ENDDO
          GO TO 4556
        ENDIF

!  Establish the first non-trivial row
!      This is necessary to avoid underflow problems

        TR = 0
        DO JCOL = 1,  NSTKS_NSTRMS
          if ( DABS(RITE_EVEC(JCOL,AA)).GT.1.0D-10 ) THEN
            TR = JCOL
            GO TO 3442
          ENDIF
        ENDDO
 3442   CONTINUE

!  Use the adjoint method to get L(k) = <X^,L(G)X> / 2k <X^,X>

        SR  = ZERO
        DO JCOL = 1,  NSTKS_NSTRMS
          SR = SR + LEFT_EVEC(JCOL,AA) * RITE_EVEC(JCOL,AA)
        ENDDO
        MODULUS = ONE / KVAL2 / SR
        DO Q = 1, N_PARAMETERS
          DO IROW = 1, NSTKS_NSTRMS
            HR = ZERO
            DO JCOL = 1, NSTKS_NSTRMS
              HR = HR + L_HSOL%L_EIGENMAT(IROW,JCOL,N,Q) * RITE_EVEC(JCOL,AA)
            ENDDO
            HVEC_AUX(IROW,Q) = HR
          ENDDO
          HR = ZERO
          DO JCOL = 1,  NSTKS_NSTRMS
            HR = HR + LEFT_EVEC(JCOL,AA) * HVEC_AUX(JCOL,Q)
          ENDDO
          L_HSOL%L_KEIGEN(K,N,Q) = HR * MODULUS
        ENDDO

!  Now solve for the linearized eigenvector L(X) from A.L(X) = B
!     Determine matrix HMAT (=A) and vectors HVEC (B)

!  Initialise first, as HMAT is destroyed after LAPACK routines

        DO IROW = 1, NSTKS_NSTRMS
          DO JCOL = 1, NSTKS_NSTRMS
            HMAT1(IROW,JCOL) = ZERO
          ENDDO
        ENDDO

!  set up the matrix HMAT
!    Use the linearization of the Eigen-equation for all but one
!    of the rows. For the FIRST row, use the linearization of the
!    normalization condition on the eigenvectors.

        DO IROW = 1, NSTKS_NSTRMS
          IF ( IROW.NE.TR ) THEN
           DO JCOL = 1, NSTKS_NSTRMS
            HMAT1(IROW,JCOL) = - HSOL%EIGENMAT_SAVE(IROW,JCOL,N)
           ENDDO
           HMAT1(IROW,IROW) = HMAT1(IROW,IROW) + KSQD
          ELSE
           DO JCOL = 1, NSTKS_NSTRMS
            HMAT1(IROW,JCOL) = RITE_EVEC(JCOL,AA)
           ENDDO
          ENDIF
        ENDDO

!  set up the HVEC:
!    Use the linearization of the Eigen-equation for all but one
!    of the rows. For the FIRST row, use the linearization of the
!    normalization condition on the eigenvectors.

        DO Q = 1, N_PARAMETERS
          I_NONTRIVIAL = 0
          DO IROW = 1, NSTKS_NSTRMS
           IF ( IROW.NE.TR ) THEN
            HVEC1(IROW,Q) = HVEC_AUX(IROW,Q) - &
               KVAL2 * L_HSOL%L_KEIGEN(K,N,Q) * RITE_EVEC(IROW,AA)
            IF ( DABS(HVEC1(IROW,Q)).LT.1.0D-12)THEN
              HVEC1(IROW,Q) = ZERO
            ELSE
              I_NONTRIVIAL = I_NONTRIVIAL + 1
            ENDIF
           ELSE
            HVEC1(IROW,Q) = ZERO
           ENDIF
          ENDDO
        ENDDO

!  Solve matrix system  using LAPACK modules DGETRF and DEGTRS
!  Only do if vector has non-trivial entries

        IF ( I_NONTRIVIAL.GT.0) THEN

          CALL DGETRF  ( NSTKS_NSTRMS, NSTKS_NSTRMS, HMAT1, &
                       MAXSTRMSTKS,  LAPACK_IPIV1, LAPACK_INFO )

!  Exception handling 2

          IF ( LAPACK_INFO .GT. 0 ) THEN
            WRITE(CI, '(I3)' ) LAPACK_INFO
            WRITE(CN, '(I3)' ) N
            MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
            TRACE   = 'DGETRF call # 2 in MSSEG_LA_QHOM_SOLUTION,'// &
                      ' layer # '//CN
            STATUS  = VQS_SERIOUS
            RETURN
          ELSE IF ( LAPACK_INFO .LT. 0 ) THEN
            WRITE(CI, '(I3)' ) LAPACK_INFO
            WRITE(CN, '(I3)' ) N
            MESSAGE = 'argument i illegal value, for i = '//CI
            TRACE   = 'DGETRF call # 2 in MSSEG_LA_QHOM_SOLUTION,'// &
                      ' layer # '//CN
            STATUS  = VQS_SERIOUS
            RETURN
          ENDIF

          CALL DGETRS ('N', NSTKS_NSTRMS, N_PARAMETERS, HMAT1, &
                     MAXSTRMSTKS, LAPACK_IPIV1,       HVEC1, &
                     MAXSTRMSTKS, LAPACK_INFO )

!  Exception handling 2

          IF ( LAPACK_INFO .LT. 0 ) THEN
            WRITE(CI, '(I3)' ) LAPACK_INFO
            WRITE(CN, '(I3)' ) N
            MESSAGE = 'argument i illegal value, for i = '//CI
            TRACE   = 'DGETRS call # 2 in MSSEG_LA_QHOM_SOLUTION,'// &
                      ' layer # '//CN
            STATUS  = VQS_SERIOUS
            RETURN
          ENDIF

        ENDIF

!  Continuation point if linearization is skipped

 4556   CONTINUE

!  Start loop over varying parameters, and assign
!    linearized sum vector = linearized eigensolution

        DO Q = 1, N_PARAMETERS

         DO I = 1, NSTREAMS
           IR = NSTOKES*(I-1)
           DO O1 = 1, NSTOKES
            IROW = IR + O1
            L_FWD_SUMVEC_R(I,O1,Q) = HVEC1(IROW,Q)
           ENDDO
         ENDDO

!  Debug (important)
!   Tested 28 December 2005.

!         IF ( DO_DEBUG_WRITE ) THEN
!          IF (Q.EQ.1.and.m.eq.0.and.n.eq.1) THEN
!           write(37,'(2i3,1p2e20.12)') M,K,KEIGEN(K,N),L_KEIGEN(K,N,Q)
!           DO I = 1, NSTREAMS ; DO O1 = 1, NSTOKES
!                 write(37,'(4i3,1p2e20.12)')
!     &        m,k,i,o1,FWD_SUMVEC(I,O1,K),L_FWD_SUMVEC_R(I,O1,Q)
!            ENDDO
!           ENDDO
!          ENDIF
!         ENDIF

        ENDDO

!  End Clause ASYMTX vs. DGEEV

       ENDIF

!  Resume General Code
!  ===================

!  Start loop over varying parameters

       DO Q = 1, N_PARAMETERS

!  linearized difference vector

         DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           SR = ZERO
           DO J = 1, NSTREAMS
            JC = NSTOKES*(J-1)
            DO O2 = 1, NSTOKES
             JCOL = JC + O2
             SR = SR - L_HSOL%L_SAB(I,J,O1,O2,N,Q) *   HSOL%FWD_SUMVEC    (J,O2,K) &
                     -   HSOL%SAB  (I,J,O1,O2,N)   * L_FWD_SUMVEC_R(J,O2,Q)
            ENDDO
           ENDDO
           HR = ( - SR - L_HSOL%L_KEIGEN(K,N,Q) * HSOL%FWD_DIFVEC(I,O1,K) )
           L_FWD_DIFVEC_R(I,O1) = HR / KVAL
          ENDDO
         ENDDO

!  assign Forward solutions PHI (Siewert's notation)
!    --->  first N are "DOWN", last N are "UP" (streams)
!    --->  Use symmetry properties to set -ve eigensolutions

         DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
           LFPD = HALF * ( L_FWD_SUMVEC_R(I,O1,Q)+L_FWD_DIFVEC_R(I,O1) )
           LFND = HALF * ( L_FWD_SUMVEC_R(I,O1,Q)-L_FWD_DIFVEC_R(I,O1) )
           L_HSOL%L_SOLA_XPOS(I,O1,K,N,Q)  = LFPD
           L_HSOL%L_SOLB_XNEG(I,O1,K,N,Q)  = LFND
           L_HSOL%L_SOLA_XPOS(I1,O1,K,N,Q)  = LFND
           L_HSOL%L_SOLB_XNEG(I1,O1,K,N,Q)  = LFPD
           IF ( O1 .GT. 2 ) THEN
            L_HSOL%L_SOLA_XPOS(I1,O1,K,N,Q)  = - LFND
            L_HSOL%L_SOLB_XNEG(I1,O1,K,N,Q)  = - LFPD
           ENDIF
          ENDDO
         ENDDO

!  End parameter loop

       ENDDO

!  End real eigenvalue loop

      ENDDO

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Linearization of COMPLEX Eigensolutions
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Only if flagged

      IF ( .NOT. DO_REAL_EIGENSOLVER(M,N) ) THEN

!  Offset

       KO1 = HSOL%K_REAL(N) + 1

!  start loop over eigensolutions

       DO K = 1, HSOL%K_COMPLEX(N)

!  Bookkeeping indices

        K0 = 2 * ( K - 1 )
        K1 = KO1 + K0
        K2 = K1  + 1

!  Basic quantities

        AA  = EIGENMASK_C(K)
        AA1 = AA + 1
        KVAL_CR  = HSOL%KEIGEN(K1,N)
        KVAL_CI  = HSOL%KEIGEN(K2,N)
        KVAL2_CR = TWO * KVAL_CR
        KVAL2_CI = TWO * KVAL_CI
        MOD_KVAL = KEIGEN_CSQ(K)
        SEPCON_CR =   HSOL%KEIGEN(K1,N) / MOD_KVAL
        SEPCON_CI = - HSOL%KEIGEN(K2,N) / MOD_KVAL

!  Establish the zero imaginary part of one eigenvector component

        TR = 0
        DO JCOL = 1,  NSTKS_NSTRMS
          if ( RITE_EVEC(JCOL,AA1).EQ.ZERO ) THEN
            TR = JCOL + NSTKS_NSTRMS
            GO TO 3443
          ENDIF
        ENDDO
 3443   CONTINUE

!  Use the adjoint theory to get L(k) from 2kL(k) = <X^,L(G)X> / <X^,X>
!      MUST USE THE COMPLEX CONJUGATE OF THE LEFT EIGENVECTOR.

!  --------------- Adjoint  norm

        SCR = ZERO
        SCI = ZERO
        DO JCOL = 1,  NSTKS_NSTRMS
          HCR = + LEFT_EVEC(JCOL,AA)  * RITE_EVEC(JCOL,AA) &
                + LEFT_EVEC(JCOL,AA1) * RITE_EVEC(JCOL,AA1)
          SCR = SCR + HCR
          HCI = - LEFT_EVEC(JCOL,AA1) * RITE_EVEC(JCOL,AA) &
                + LEFT_EVEC(JCOL,AA)  * RITE_EVEC(JCOL,AA1)
          SCI = SCI + HCI
          ADJN_CR(JCOL) = + KVAL2_CR * RITE_EVEC(JCOL,AA) &
                          - KVAL2_CI * RITE_EVEC(JCOL,AA1)
          ADJN_CI(JCOL) = + KVAL2_CI * RITE_EVEC(JCOL,AA) &
                          + KVAL2_CR * RITE_EVEC(JCOL,AA1)
        ENDDO
        DENOM_CR = KVAL2_CR * SCR - KVAL2_CI * SCI
        DENOM_CI = KVAL2_CR * SCI + KVAL2_CI * SCR
        MODULUS  = ONE / ( DENOM_CR*DENOM_CR + DENOM_CI * DENOM_CI )

!  ---------------- Adjoint numerator, and final computation

        DO Q = 1, N_PARAMETERS
         DO IROW = 1, NSTKS_NSTRMS
          HCR = ZERO
          HCI = ZERO
          DO JCOL = 1, NSTKS_NSTRMS
           HCR = HCR + L_HSOL%L_EIGENMAT(IROW,JCOL,N,Q) * RITE_EVEC(JCOL,AA)
           HCI = HCI + L_HSOL%L_EIGENMAT(IROW,JCOL,N,Q) * RITE_EVEC(JCOL,AA1)
          ENDDO
          HVEC_CR(IROW,Q) = HCR
          HVEC_CI(IROW,Q) = HCI
         ENDDO
         HCR = ZERO
         HCI = ZERO
         DO JCOL = 1,  NSTKS_NSTRMS
          HCR = HCR + LEFT_EVEC(JCOL,AA)  * HVEC_CR(JCOL,Q) &
                    + LEFT_EVEC(JCOL,AA1) * HVEC_CI(JCOL,Q)
          HCI = HCI - LEFT_EVEC(JCOL,AA1) * HVEC_CR(JCOL,Q) &
                    + LEFT_EVEC(JCOL,AA)  * HVEC_CI(JCOL,Q)
         ENDDO
         L_HSOL%L_KEIGEN(K1,N,Q) = ( HCR*DENOM_CR + HCI*DENOM_CI )*MODULUS
         L_HSOL%L_KEIGEN(K2,N,Q) = ( HCI*DENOM_CR - HCR*DENOM_CI )*MODULUS
        ENDDO

!  Now solve for the linearized eigenvector L(X) from A.L(X) = B
!     Determine matrix HMAT (=A) and vectors HVEC (B)

!  ---------- Initialise first, as HMAT is destroyed after LAPACK routin

        DO IROW = 1, NSTKS_NSTRMS_2
          DO JCOL = 1, NSTKS_NSTRMS_2
            HMAT2(IROW,JCOL) = ZERO
          ENDDO
        ENDDO

!  ---------- set up matrix for Full solution

        DO I = 1, NSTKS_NSTRMS
          I1 = I + NSTKS_NSTRMS
          DO J = 1, NSTKS_NSTRMS
            J1 = J + NSTKS_NSTRMS
            HMAT2(I,J)   = - HSOL%EIGENMAT_SAVE(I,J,N)
            HMAT2(I1,J1) = - HSOL%EIGENMAT_SAVE(I,J,N)
          ENDDO
          HMAT2(I,I)   = HMAT2(I,I)   + REAL_KSQ(AA)
          HMAT2(I1,I1) = HMAT2(I1,I1) + REAL_KSQ(AA)
          HMAT2(I,I1)  = - IMAG_KSQ(AA)
          HMAT2(I1,I)  = + IMAG_KSQ(AA)
        ENDDO

!  ---------- set up the vectors for full solution

        DO Q = 1, N_PARAMETERS
          DO I = 1, NSTKS_NSTRMS
           I1 = I + NSTKS_NSTRMS
           HCR = + L_HSOL%L_KEIGEN(K1,N,Q) * ADJN_CR(I) &
                 - L_HSOL%L_KEIGEN(K2,N,Q) * ADJN_CI(I)
           HCI = + L_HSOL%L_KEIGEN(K1,N,Q) * ADJN_CI(I) &
                 + L_HSOL%L_KEIGEN(K2,N,Q) * ADJN_CR(I)
           HVEC2(I,Q)  = HVEC_CR(I,Q) - HCR
           HVEC2(I1,Q) = HVEC_CI(I,Q) - HCI
          ENDDO
        ENDDO

!  ---------- Strike out TR row/column to get reduced matrix and vectors

!mick mod
        !NSTKS_NSTRMS_2  = 2 * NSTKS_NSTRMS
        NSTKS_NSTRMS_21 = 2 * NSTKS_NSTRMS - 1
        DO I = 1, TR - 1
          DO J = 1, TR - 1
            HMATRED(I,J) = HMAT2(I,J)
          ENDDO
          DO J = TR + 1, NSTKS_NSTRMS_2
            HMATRED(I,J-1) = HMAT2(I,J)
          ENDDO
          DO Q = 1, N_PARAMETERS
            HVECRED(I,Q) = HVEC2(I,Q)
          ENDDO
        ENDDO
        DO I = TR + 1, NSTKS_NSTRMS_2
          DO J = 1, TR - 1
            HMATRED(I-1,J) = HMAT2(I,J)
          ENDDO
          DO J = TR + 1, NSTKS_NSTRMS_2
            HMATRED(I-1,J-1) = HMAT2(I,J)
          ENDDO
          DO Q = 1, N_PARAMETERS
            HVECRED(I-1,Q) = HVEC2(I,Q)
          ENDDO
        ENDDO

!  ----------  Replace the last row by linearized normalization constrai
!                ( IF (TR.LT.NSTKS_NSTRMS_2) THEN  possible exception...
!                Should not be necessary to flag this case )
!
        ILAST =  NSTKS_NSTRMS_21
        DO J = 1, NSTKS_NSTRMS
          HMATRED(ILAST,J) = RITE_EVEC(J,AA)
        ENDDO
        DO J = NSTKS_NSTRMS+1, TR - 1
          HMATRED(ILAST,J) = RITE_EVEC(J-NSTKS_NSTRMS,AA1)
        ENDDO
        DO J = TR + 1, NSTKS_NSTRMS_2
          HMATRED(ILAST,J-1) = RITE_EVEC(J-NSTKS_NSTRMS,AA1)
        ENDDO
        DO Q = 1, N_PARAMETERS
          HVECRED(ILAST,Q) = ZERO
        ENDDO

!  ---------- Solve matrix system using LAPACK modules DGETRF and DEGTRS

        CALL DGETRF  ( NSTKS_NSTRMS_21, NSTKS_NSTRMS_21, HMATRED, &
                       MAXSTRMSTKS_21,  LAPACK_IPIV21,   LAPACK_INFO )

!  Exception handling 3

        IF ( LAPACK_INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) LAPACK_INFO
          WRITE(CN, '(I3)' ) N
          MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE   = 'DGETRF call # 3 in MSSEG_LA_QHOM_SOLUTION,'// &
                      ' layer # '//CN
          STATUS  = VQS_SERIOUS
          RETURN
        ELSE IF ( LAPACK_INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) LAPACK_INFO
          WRITE(CN, '(I3)' ) N
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRF call # 3 in MSSEG_LA_QHOM_SOLUTION,'// &
                      ' layer # '//CN
          STATUS  = VQS_SERIOUS
          RETURN
        ENDIF

        CALL DGETRS ('N', NSTKS_NSTRMS_21, N_PARAMETERS, HMATRED, &
                     MAXSTRMSTKS_21, LAPACK_IPIV21,      HVECRED, &
                     MAXSTRMSTKS_21, LAPACK_INFO )

!  Exception handling 3

        IF ( LAPACK_INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) LAPACK_INFO
          WRITE(CN, '(I3)' ) N
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRS call # 3 in MSSEG_LA_QHOM_SOLUTION,'// &
                      ' layer # '//CN
          STATUS  = VQS_SERIOUS
          RETURN
        ENDIF

!  Start loop over varying parameters

        DO Q = 1, N_PARAMETERS

!  linearized sum vector = linearized eigensolution

         DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
           IROW  = IR + O1
           IROW1 = IROW + NSTKS_NSTRMS
           L_FWD_SUMVEC_CR(I,O1) = HVECRED(IROW,Q)
           IF ( IROW1.LT.TR)THEN
             L_FWD_SUMVEC_CI(I,O1) = HVECRED(IROW1,Q)
           ELSE IF (IROW1.EQ.TR) THEN
             L_FWD_SUMVEC_CI(I,O1) = ZERO
           ELSE
             L_FWD_SUMVEC_CI(I,O1) = HVECRED(IROW1-1,Q)
           ENDIF
          ENDDO
         ENDDO

!  Debug code (important)

!     IF (DO_DEBUG_WRITE) THEN ; IF ( Q.EQ.1.AND.M.EQ.0.AND.N.EQ.1) THEN
!       WRITE(91,'(2I3,1P4E20.10)') M,K,KEIGEN(K1,N),KEIGEN(K2,N),L_KEIGEN(K1,N,Q),L_KEIGEN(K2,N,Q)
!       DO I = 1, NSTREAMS ; DO O1 = 1, NSTOKES
!         WRITE(91,'(4I3,1P4E20.10)')M,K,I,O1,FWD_SUMVEC(I,O1,K1),FWD_SUMVEC(I,O1,K2),L_FWD_SUMVEC_CR(I,O1),L_FWD_SUMVEC_CI(I,O1)
!       ENDDO ; ENDDO
!     ENDIF ; ENDIF

!  Determine the linearized difference vector

         DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           SCR = ZERO
           SCI = ZERO
           DO J = 1, NSTREAMS
            JC = NSTOKES*(J-1)
            DO O2 = 1, NSTOKES
             JCOL = JC + O2
             SCR = SCR - L_HSOL%L_SAB(I,J,O1,O2,N,Q) *   HSOL%FWD_SUMVEC     (J,O2,K1) &
                       -   HSOL%SAB  (I,J,O1,O2,N)   * L_FWD_SUMVEC_CR(J,O2)
             SCI = SCI - L_HSOL%L_SAB(I,J,O1,O2,N,Q) *   HSOL%FWD_SUMVEC     (J,O2,K2) &
                       -   HSOL%SAB  (I,J,O1,O2,N)   * L_FWD_SUMVEC_CI(J,O2)
            ENDDO
           ENDDO
           HCR = - SCR - L_HSOL%L_KEIGEN(K1,N,Q) * HSOL%FWD_DIFVEC(I,O1,K1) &
                       + L_HSOL%L_KEIGEN(K2,N,Q) * HSOL%FWD_DIFVEC(I,O1,K2)
           HCI = - SCI - L_HSOL%L_KEIGEN(K2,N,Q) * HSOL%FWD_DIFVEC(I,O1,K1) &
                       - L_HSOL%L_KEIGEN(K1,N,Q) * HSOL%FWD_DIFVEC(I,O1,K2)
           L_FWD_DIFVEC_CR(I,O1) = HCR*SEPCON_CR - HCI*SEPCON_CI
           L_FWD_DIFVEC_CI(I,O1) = HCI*SEPCON_CR + HCR*SEPCON_CI
          ENDDO
         ENDDO

!  assign Forward solutions PHI (Siewert's notation)
!    --->  first N are "DOWN", last N are "UP" (streams)
!    --->  Use symmetry properties to set -ve eigensolutions

         DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
           LFND1 = HALF * (L_FWD_SUMVEC_CR(I,O1)-L_FWD_DIFVEC_CR(I,O1))
           LFPD1 = HALF * (L_FWD_SUMVEC_CR(I,O1)+L_FWD_DIFVEC_CR(I,O1))
           LFND2 = HALF * (L_FWD_SUMVEC_CI(I,O1)-L_FWD_DIFVEC_CI(I,O1))
           LFPD2 = HALF * (L_FWD_SUMVEC_CI(I,O1)+L_FWD_DIFVEC_CI(I,O1))
           L_HSOL%L_SOLA_XPOS(I,O1,K1,N,Q)  = LFPD1
           L_HSOL%L_SOLB_XNEG(I,O1,K1,N,Q)  = LFND1
           L_HSOL%L_SOLA_XPOS(I,O1,K2,N,Q)  = LFPD2
           L_HSOL%L_SOLB_XNEG(I,O1,K2,N,Q)  = LFND2
           L_HSOL%L_SOLA_XPOS(I1,O1,K1,N,Q)  = LFND1
           L_HSOL%L_SOLB_XNEG(I1,O1,K1,N,Q)  = LFPD1
           L_HSOL%L_SOLA_XPOS(I1,O1,K2,N,Q)  = LFND2
           L_HSOL%L_SOLB_XNEG(I1,O1,K2,N,Q)  = LFPD2
           IF ( O1 .GT. 2 ) THEN
            L_HSOL%L_SOLA_XPOS(I1,O1,K1,N,Q)  = - LFND1
            L_HSOL%L_SOLB_XNEG(I1,O1,K1,N,Q)  = - LFPD1
            L_HSOL%L_SOLA_XPOS(I1,O1,K2,N,Q)  = - LFND2
            L_HSOL%L_SOLB_XNEG(I1,O1,K2,N,Q)  = - LFPD2
           ENDIF
          ENDDO
         ENDDO

!  End parameter loop

        ENDDO

!  End complex eigenvalue loop

       ENDDO

!  End clause for not using the real eigensolver

      ENDIF

!  Linearized Eigenstream transmittance factors for whole layer
!  ------------------------------------------------------------

!  Real

      DO K = 1, HSOL%K_REAL(N)
         TBAS = - HSOL%T_DELT_EIGEN(K,N) * DELTAU_VERT(N)
         DO Q = 1, N_PARAMETERS
            LTQ = L_HSOL%L_KEIGEN(K,N,Q) + HSOL%KEIGEN(K,N)*L_DELTAU_VERT(Q,N)
            L_HSOL%L_T_DELT_EIGEN(K,N,Q) = TBAS * LTQ
         ENDDO
      ENDDO

!  Complex

      DO K = 1, HSOL%K_COMPLEX(N)
         K0 = 2 * K - 2
         K1 = KO1 + K0
         K2 = K1  + 1
         TBAS_CR = - HSOL%T_DELT_EIGEN(K1,N) *  DELTAU_VERT(N)
         TBAS_CI = - HSOL%T_DELT_EIGEN(K2,N) *  DELTAU_VERT(N)
         DO Q = 1, N_PARAMETERS
            LCR = L_HSOL%L_KEIGEN(K1,N,Q) + HSOL%KEIGEN(K1,N) * L_DELTAU_VERT(Q,N)
            LCI = L_HSOL%L_KEIGEN(K2,N,Q) + HSOL%KEIGEN(K2,N) * L_DELTAU_VERT(Q,N)
            L_HSOL%L_T_DELT_EIGEN(K1,N,Q) = TBAS_CR * LCR - TBAS_CI * LCI
            L_HSOL%L_T_DELT_EIGEN(K2,N,Q) = TBAS_CI * LCR + TBAS_CR * LCI
         ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE MSSEG_LA_QHOM_SOLUTION

!  Finish module

      END MODULE MSSeg_LA_Solutions_m
