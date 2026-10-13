
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
!   -- Use public routine GETQUAD2 instead of GAULEG

!  Mk2a. 9/21/21. Created LC routines
!  Mk3.  1/31/22. Rearranged arguments lists and used type structures, LC routines
!  Mk3.  4/11/23. Readied for use. Some renaming

      MODULE MSSeg_LC_Solutions_m

      PRIVATE
      PUBLIC :: MSSEG_LAC_UHOM_SOLUTION, &
                MSSEG_LC_QBEAM_SOLUTION, MSSEG_LC_UBEAM_SOLUTION

      CONTAINS

!

      SUBROUTINE MSSEG_LAC_UHOM_SOLUTION ( &
        GIVEN_LAYER, FOURIER, DO_VARY, N_PARAMETERS,    & ! Input control
        NSTOKES, NSTREAMS, NMOMENTS, NGEOMS,            & ! Input control
        DO_UPWELLING, DO_DNWELLING, QUAD_HALFWTS,       & ! Input bookkeeping
        DO_LAYER_SCATTERING,OMEGA_GREEK, L_OMEGA_GREEK, & ! Input optical
        HGSF, PGSF, HSOL, PSOL, L_HSOL, L_PSOL )          ! TS Input and Output

!  Mk2, 6/21/21.  Use module names now have the "_m" extension, Only use quantities you need
!  Mk3. 1/31/22. Rearranged arguments lists and used type structures

      USE MSSeg_PARS_m, Only : MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS, MAX_USER_STREAMS, &
                               MAXMOMENTS, MAXEVALUES, MAXSTREAMS_2, MAXSTRMSTKS, ZERO, ONE, TWO

      use MSSEG_HPTypes_Def_m
      use MSSEG_L_HPTypes_Def_m

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: GIVEN_LAYER
      INTEGER, INTENT (IN) :: FOURIER
      LOGICAL, INTENT (IN) :: DO_VARY
      INTEGER, INTENT (IN) :: N_PARAMETERS

      INTEGER, INTENT (IN) :: NSTOKES
      INTEGER, INTENT (IN) :: NSTREAMS
      INTEGER, INTENT (IN) :: NMOMENTS
      INTEGER, INTENT (IN) :: NGEOMS

!  flags

      LOGICAL, INTENT (IN) :: DO_UPWELLING
      LOGICAL, INTENT (IN) :: DO_DNWELLING

!  bookkeeping

      DOUBLE PRECISION, INTENT (IN) ::  QUAD_HALFWTS ( MAXSTREAMS )
      LOGICAL, INTENT (IN) :: DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  optical

      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_GREEK  ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  input type structures (Generalized spherical functions)

      type(Homog_GSF) , intent(in) :: HGSF
      type(Partic_GSF), intent(in) :: PGSF

!  input type structure (Homogeneous and particular solutions)

      type(Homog_solutions) , intent(in) :: HSOL
      type(Partic_solutions), intent(in) :: PSOL

!  input type structure (Linearized homogeneous solutions)
 
      type(L_Homog_solutions) , intent(in) :: L_HSOL

!  Output type structure (Linearized particular solutions)
 
      type(L_Partic_Solutions) , intent(inout) :: L_PSOL

!  Local variables
!  ---------------

      DOUBLE PRECISION :: L_HELPSTOKES_R(MAXSTOKES)
      DOUBLE PRECISION :: L_HELPSTOKES_CR(MAXSTOKES)
      DOUBLE PRECISION :: L_HELPSTOKES_CI(MAXSTOKES)
      DOUBLE PRECISION :: L_GAUX_R(0:MAXMOMENTS,MAXSTOKES)
      DOUBLE PRECISION :: L_GAUX_CR(0:MAXMOMENTS,MAXSTOKES)
      DOUBLE PRECISION :: L_GAUX_CI(0:MAXMOMENTS,MAXSTOKES)

      INTEGER ::          UM, J, L, N, M, Q, O1, O2
      INTEGER ::          K, KO1, K0, K1, K2
      DOUBLE PRECISION :: SPS, H, SP, SM, SG, SN, SNS, HPC, HMC
      DOUBLE PRECISION :: SPSCR, SPSCI, SPCR, SPCI
      DOUBLE PRECISION :: SNSCR, SNSCI, SNCR, SNCI
      DOUBLE PRECISION :: SR, SI, SMCR, SMCI

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

!  Initialize output

      L_PSOL%L_UHOM_UPUP = zero
      L_PSOL%L_UHOM_UPDN = zero
      L_PSOL%L_UHOM_DNUP = zero
      L_PSOL%L_UHOM_DNDN = zero

!  Eigenvector interpolation to user-defined angles
!  ================================================

!  Only do this if variation flag,
!    Or if there is a solution for scattering in this layer

      IF ( .NOT.DO_VARY .OR. .NOT. DO_LAYER_SCATTERING(M,N) ) RETURN

!  For each parameter

      DO Q = 1, N_PARAMETERS

!  User defined solutions (Real)
!  -----------------------------

       DO K = 1, HSOL%K_REAL(N)

!  For each moment, do inner sum over computational angles
!  for the positive and negative linearized eigenvectors.
!  Linearized the product with OMEGA_GREEK ---> L_GAUX

        DO L = M, NMOMENTS
         DO O1 = 1, NSTOKES
          SPS = ZERO
          DO J = 1, NSTREAMS
           H = QUAD_HALFWTS(J)
           SP = ZERO
           SM = ZERO
           DO O2 = 1, NSTOKES
            SP = SP + H*HGSF%PI_XQP    (L,J,O1,O2)*L_HSOL%L_SOLA_XPOS(J,O2,K,N,Q)
            SM = SM + H*HGSF%PI_XQM_PRE(L,J,O1,O2)*L_HSOL%L_SOLB_XNEG(J,O2,K,N,Q)
           ENDDO
           SPS = SPS + SP + SM
          ENDDO
          L_HELPSTOKES_R(O1) = SPS
         ENDDO
         DO O1 = 1, NSTOKES
          SG = ZERO
          DO O2 = 1, NSTOKES
           SG = SG + OMEGA_GREEK(L,N,O1,O2)     * L_HELPSTOKES_R(O2) &
                   + L_OMEGA_GREEK(L,N,O1,O2,Q) *   PSOL%HELPSTOKES(L,K,O2)
          ENDDO
          L_GAUX_R(L,O1) = SG
         ENDDO
        ENDDO

!  Now sum over all harmonic contributions (downwelling)

        IF ( DO_DNWELLING ) THEN
         DO UM = 1, NGEOMS
          DO O1 = 1, NSTOKES
           SP = ZERO ; SN = ZERO
           DO L = M, NMOMENTS
            SPS = ZERO ; SNS = ZERO
            DO O2 = 1, NSTOKES
             SPS = SPS + PGSF%PI_XUP(L,UM,O1,O2)      * L_GAUX_R(L,O2)
             SNS = SNS + PGSF%PI_XUM_POST(L,UM,O1,O2) * L_GAUX_R(L,O2)
            ENDDO
            SP = SP + SPS
            SN = SN + SNS
           ENDDO
           L_PSOL%L_UHOM_DNDN(UM,O1,K,Q) = SP
           L_PSOL%L_UHOM_DNUP(UM,O1,K,Q) = SN
          ENDDO
         ENDDO
        ENDIF

!  Now sum over all harmonic contributions (Upwelling)

        IF ( DO_UPWELLING ) THEN
         DO UM = 1, NGEOMS
          DO O1 = 1, NSTOKES
           SP = ZERO ; SN = ZERO
           DO L = M, NMOMENTS
            SPS = ZERO ; SNS = ZERO
            DO O2 = 1, NSTOKES
             SPS = SPS + PGSF%PI_XUP_PRE(L,UM,O1,O2)  * L_GAUX_R(L,O2)
             SNS = SNS + PGSF%PI_XUM(L,UM,O1,O2)      * L_GAUX_R(L,O2)
            ENDDO
            SP = SP + SPS
            SN = SN + SNS
           ENDDO
           L_PSOL%L_UHOM_UPDN(UM,O1,K,Q) = SN
           L_PSOL%L_UHOM_UPUP(UM,O1,K,Q) = SP
          ENDDO
         ENDDO
        ENDIF

!  Debug (important)

!        IF ( DO_DEBUG_WRITE.AND.Q.EQ.2) THEN
!         DO UM = 1, NGEOMS
!          DO O1 = 1, NSTOKES
!            WRITE(95,'(4I3,1P8E20.10)')M,K,UM,O1, &
!               L_UHOM_UPUP(UM,O1,K,Q), &
!               L_UHOM_UPDN(UM,O1,K,Q), &
!               L_UHOM_DNUP(UM,O1,K,Q), &
!               L_UHOM_DNDN(UM,O1,K,Q), &
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
            SPCR = SPCR + HPC * L_HSOL%L_SOLA_XPOS(J,O2,K1,N,Q)
            SPCI = SPCI + HPC * L_HSOL%L_SOLA_XPOS(J,O2,K2,N,Q)
            HMC  =  HGSF%PI_XQM_PRE(L,J,O1,O2)*QUAD_HALFWTS(J)
            SMCR = SMCR + HMC * L_HSOL%L_SOLB_XNEG(J,O2,K1,N,Q)
            SMCI = SMCI + HMC * L_HSOL%L_SOLB_XNEG(J,O2,K2,N,Q)
           ENDDO
           SPSCR = SPSCR + SPCR + SMCR
           SPSCI = SPSCI + SPCI + SMCI
          ENDDO
          L_HELPSTOKES_CR(O1) = SPSCR
          L_HELPSTOKES_CI(O1) = SPSCI
         ENDDO
         DO O1 = 1, NSTOKES
          SR = ZERO
          SI = ZERO
          DO O2 = 1, NSTOKES
           SR = SR + OMEGA_GREEK(L,N,O1,O2)     *L_HELPSTOKES_CR(O2) &
                   + L_OMEGA_GREEK(L,N,O1,O2,Q) *  PSOL%HELPSTOKES(L,K1,O2)
           SI = SI + OMEGA_GREEK(L,N,O1,O2)     *L_HELPSTOKES_CI(O2) &
                   + L_OMEGA_GREEK(L,N,O1,O2,Q) *  PSOL%HELPSTOKES(L,K2,O2)
          ENDDO
          L_GAUX_CR(L,O1) = SR
          L_GAUX_CI(L,O1) = SI
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
             SPSCR = SPSCR + HPC * L_GAUX_CR(L,O2)
             SPSCI = SPSCI + HPC * L_GAUX_CI(L,O2)
             HMC = PGSF%PI_XUM_POST(L,UM,O1,O2)
             SNSCR = SNSCR + HMC * L_GAUX_CR(L,O2)
             SNSCI = SNSCI + HMC * L_GAUX_CI(L,O2)
            ENDDO
            SPCR = SPCR + SPSCR
            SPCI = SPCI + SPSCI
            SNCR = SNCR + SNSCR
            SNCI = SNCI + SNSCI
           ENDDO
           L_PSOL%L_UHOM_DNDN(UM,O1,K1,Q) = SPCR
           L_PSOL%L_UHOM_DNUP(UM,O1,K1,Q) = SNCR
           L_PSOL%L_UHOM_DNDN(UM,O1,K2,Q) = SPCI
           L_PSOL%L_UHOM_DNUP(UM,O1,K2,Q) = SNCI
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
             SPSCR = SPSCR + HPC * L_GAUX_CR(L,O2)
             SPSCI = SPSCI + HPC * L_GAUX_CI(L,O2)
             HMC = PGSF%PI_XUM(L,UM,O1,O2)
             SNSCR = SNSCR + HMC * L_GAUX_CR(L,O2)
             SNSCI = SNSCI + HMC * L_GAUX_CI(L,O2)
            ENDDO
            SPCR = SPCR + SPSCR
            SPCI = SPCI + SPSCI
            SNCR = SNCR + SNSCR
            SNCI = SNCI + SNSCI
           ENDDO
           L_PSOL%L_UHOM_UPDN(UM,O1,K1,Q) = SNCR
           L_PSOL%L_UHOM_UPUP(UM,O1,K1,Q) = SPCR
           L_PSOL%L_UHOM_UPDN(UM,O1,K2,Q) = SNCI
           L_PSOL%L_UHOM_UPUP(UM,O1,K2,Q) = SPCI
          ENDDO
         ENDDO
        ENDIF

!  Debug (important)

!        IF ( DO_DEBUG_WRITE.AND.Q.EQ.2) THEN
!         DO UM = 1, NGEOMS
!          DO O1 = 1, NSTOKES
!            WRITE(97,'(4I3,1P8E20.10)')M,K,UM,O1, &
!               L_UHOM_UPUP(UM,O1,K2,Q), &
!               L_UHOM_UPDN(UM,O1,K2,Q), &
!               L_UHOM_DNUP(UM,O1,K2,Q), &
!               L_UHOM_DNDN(UM,O1,K2,Q), &
!               UHOM_UPUP(UM,O1,K2), &
!               UHOM_UPDN(UM,O1,K2), &
!               UHOM_DNUP(UM,O1,K2), &
!               UHOM_DNDN(UM,O1,K2)
!           ENDDO
!         ENDDO
!        ENDIF

!  end loop over complex eigensolutions

       ENDDO

!  End loop over parameters

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE MSSEG_LAC_UHOM_SOLUTION

!

      SUBROUTINE MSSEG_LC_QBEAM_SOLUTION ( &
        GIVEN_LAYER, FOURIER, IBEAM, DOVARY, N_PARAMETERS, & ! Input control
        DO_PLANE_PARALLEL, NSTOKES, NSTREAMS, FLUX_FACTOR, & ! Input control
        NMOMENTS, NSTREAMS_2, NSTKS_NSTRMS, DFLUX,         & ! Input Bookkeeping
        DO_LAYER_SCATTERING, QUAD_STREAMS, L_OMEGA_GREEK,  & ! Input optical
        BEAM_CUTOFF, AVERAGE_SECANT, LC_AVERAGE_SECANT,    & ! Input average secant
        HGSF, PGSF, HSOL, PSOL, L_HSOL,                    & ! Input TS
        L_PSOL, STATUS, MESSAGE, TRACE )                     ! Output

!  This is the classical Chandrasekhar beam solution.
!  ( plane parallel or average secant only)
!  Linear Matrix algebra.

!  Mk 2a. Linearization LC subroutine added

!  Mk2, 6/21/21. Use module names now have the "_m" extension, Only use quantities you need
!  Mk3. 1/31/21. Rearranged arguments lists and used type structures

      USE MSSeg_PARS_m, Only : MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS, MAXMOMENTS, MAXSTREAMS_2, &
                               MAXSTRMSTKS, MAXEVALUES, VQS_SUCCESS, VQS_SERIOUS, ZERO, HALF, TWO, PI4

      use MSSEG_HPTypes_Def_m
      use MSSEG_L_HPTypes_Def_m

      USE LAPACK_TOOLS_m, Only : DGETRS

      IMPLICIT NONE

!  Subroutine Arguments
!  --------------------

!  Basic control

      INTEGER, INTENT (IN) ::           GIVEN_LAYER
      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      LOGICAL, INTENT (IN) ::           DOVARY
      INTEGER, INTENT (IN) ::           N_PARAMETERS

      LOGICAL, INTENT (IN) ::           DO_PLANE_PARALLEL
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      DOUBLE PRECISION, INTENT (IN) ::  FLUX_FACTOR

!  bookkeeping

      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           NSTREAMS_2
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      DOUBLE PRECISION, INTENT (IN) ::  DFLUX ( MAXSTOKES )

      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STREAMS ( MAXSTREAMS )

!  Linearized optical

      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  average secant inputs

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT    ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  input type structures (Generalized spherical functions)

      type(Homog_GSF) , intent(in) :: HGSF
      type(Partic_GSF), intent(in) :: PGSF

!  input type structure (Homogeneous and particular solutions)

      type(Homog_solutions) , intent(in) :: HSOL
      type(Partic_solutions), intent(in) :: PSOL

!  input type structure (Linearized homogeneous solutions)
 
      type(L_Homog_solutions) , intent(in) :: L_HSOL

!  Output type structure (Linearized particular solutions)
 
      type(L_Partic_solutions) , intent(inout) :: L_PSOL

!  Exception handling

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  Saved vector

      DOUBLE PRECISION ::  QSAVE(MAXSTRMSTKS)

!  linearization arrays

      DOUBLE PRECISION ::  L_QVEC(MAXSTRMSTKS,MAX_ATMOSWFS)

      DOUBLE PRECISION ::  L_QSUMVEC(MAXSTREAMS,MAXSTOKES,MAX_ATMOSWFS)
      DOUBLE PRECISION ::  L_QDIFVEC(MAXSTREAMS,MAXSTOKES,MAX_ATMOSWFS)
      DOUBLE PRECISION ::  HELP_TEMPO(0:MAXMOMENTS,MAXSTOKES)

!  help variables

      DOUBLE PRECISION ::  HELP, WSUM, WDIF, SECBAR
      DOUBLE PRECISION ::  TM1, TM2, TM3, L_WVEC_D, L_WVEC_U
      DOUBLE PRECISION ::  SUM, GSUM, S_TPA, S_TMA, TPA, TMA, F1

      INTEGER          ::  I, J, I1, L, N, M, INFO, IB, LVARY, Q
      INTEGER          ::  O1, O2, O3, IR, JC, IROW, JCOL
      LOGICAL          ::  DO_FIRST
      CHARACTER (LEN=3) :: CN, CI

!  Start of code
!  -------------

!  initialise Exception handling

      STATUS  = VQS_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Set layer, fourier, beam index

      N  = GIVEN_LAYER
      M  = FOURIER
      IB = IBEAM
      F1 = FLUX_FACTOR / PI4

!  safety

      SECBAR = ZERO

!  set up driving vector (using saved results)
!  This must be done regardless of whether layer N is varying or not.

      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          QSAVE(IROW) = PSOL%QDIFVEC_SAVE(I,O1) + TWO * AVERAGE_SECANT(N,IB) * PSOL%QVEC_SAVE(IROW)
        ENDDO
      ENDDO

!  Linearization for layer N is in two parts:
!    1A. Linearization due to variations in Layer N itself (Profiles)
!    1B. Linearization due to columns
!    2. Linearization due to variations in layers K < N (Profiles)

!  Part 1.
!  =======

!  Set layer to vary

      LVARY = N

!  Only a solution if the layer is active and not below Cutoff.

      DO_FIRST = ( N .LE. BEAM_CUTOFF(IB) ) .AND. DO_LAYER_SCATTERING(M,N)

!  solar zenith cosine for this layer

      SECBAR = AVERAGE_SECANT(N,IB)

!  No particular solution beyond the cutoff layer.
!  Or no scattering in this layer
!    [ Zero the boundary layer values and start Part 2 )

      IF (.NOT. DOVARY .OR. .NOT. DO_FIRST  ) THEN
        DO I = 1, NSTREAMS_2
         DO O1 = 1, NSTOKES
          DO Q = 1, N_PARAMETERS
            L_PSOL%LC_BVEC(I,O1,N,Q) = ZERO
          ENDDO
         ENDDO
        ENDDO
        RETURN
      ENDIF

!  For each varying parameter

      DO Q = 1, N_PARAMETERS

!  Auxiliary matrix for Q functions

        DO L = M, NMOMENTS
          DO O1 = 1, NSTOKES
            SUM = ZERO
            DO O2 = 1, NSTOKES
              GSUM = ZERO
              DO O3 = 1, NSTOKES
                GSUM = GSUM + PGSF%PI_X0P(L,IBEAM,N,O2,O3) * DFLUX(O3)
              ENDDO
              SUM = SUM + L_OMEGA_GREEK(L,N,O1,O2,Q) * GSUM
            ENDDO
            HELP_TEMPO(L,O1) = SUM
          ENDDO
        ENDDO

!  Set up linearized sum and difference vectors for Beam source terms

        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            S_TPA = ZERO
            S_TMA = ZERO
            DO L = M, NMOMENTS
              TPA = ZERO
              TMA = ZERO
              DO O2 = 1, NSTOKES
                TPA = TPA + HGSF%PI_XQP(L,I,O1,O2)      * HELP_TEMPO(L,O2)
                TMA = TMA + HGSF%PI_XQM_POST(L,I,O1,O2) * HELP_TEMPO(L,O2)
              ENDDO
              S_TPA = S_TPA + TPA
              S_TMA = S_TMA + TMA
            ENDDO
            L_QSUMVEC(I,O1,Q) = F1 * ( S_TPA + S_TMA ) / QUAD_STREAMS(I)
            L_QDIFVEC(I,O1,Q) = F1 * ( S_TPA - S_TMA ) / QUAD_STREAMS(I)
          ENDDO
        ENDDO

!   setup linearized RHS vector
!  ( use results from the original solution )
!    WARNING. BE CAREFUL OF SIGNS on TM1, TM2, TM3
!     slgihtly different from the scalar model case.

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          HELP = ZERO
          DO J = 1, NSTREAMS
           JC = NSTOKES*(J-1)
           DO O2 = 1, NSTOKES
            JCOL = JC + O2
            TM1 = L_HSOL%L_EIGENMAT(IROW,JCOL,N,Q) * PSOL%QVEC_SAVE(JCOL)
            TM2 =   HSOL%DAB  (I,J,O1,O2,N)   *    L_QSUMVEC     (J,O2,Q)
            TM3 = L_HSOL%L_DAB(I,J,O1,O2,N,Q) * PSOL%QSUMVEC_SAVE(J,O2)
            HELP = HELP - TM1 + TM2 + TM3
!            HELP = HELP - TM1 - TM2 - TM3
           ENDDO
          ENDDO
          L_QVEC(IROW,Q)  = HELP + L_QDIFVEC(I,O1,Q) * SECBAR
         ENDDO
        ENDDO

!  end parameter loop

      ENDDO

!  additional terms for the quasi-spherical case
!  ( layers greater than one )

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
        IF ( N.GT.1 ) THEN
          DO Q = 1, N_PARAMETERS
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                HELP = QSAVE(IROW) * LC_AVERAGE_SECANT(N,IB,Q)
                L_QVEC(IROW,Q) = L_QVEC(IROW,Q) + HELP
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Solve problem by back substitution for all of the RHS vectors
!  ( uses L_U decomposition of QMAT_SAVE from original solution )

      CALL DGETRS &
              ('N',NSTKS_NSTRMS,N_PARAMETERS,PSOL%QMAT_SAVE, &
                MAXSTRMSTKS,PSOL%QPIVOT,L_QVEC,MAXSTRMSTKS,INFO)

!  Exception handling 1

      IF ( INFO .LT. 0 ) THEN
        WRITE(CI, '(I3)' ) INFO ; WRITE(CN, '(I3)' ) N
        MESSAGE = 'argument i illegal value, for i = '//CI
        TRACE   = 'DGETRS call # 1 in MSSEG_LC_QBEAM_SOLUTION,'// &
                  ' Beam P.I. linearization (N,N), layer # '//CN
        STATUS  = VQS_SERIOUS ; RETURN
      ENDIF

!  assign solutions for the quasi-spherical case, N > 1
!    Note change of sign with HELP, this is because SAB as
!    defined in VLIDORT is the negative of SAB in the scalar model.

      IF ( .NOT. DO_PLANE_PARALLEL .AND. N.GT.1 ) THEN

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         I1 = I + NSTREAMS
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          TM3 = - PSOL%QDIF_SAVE(IROW) / SECBAR
          DO Q = 1, N_PARAMETERS
           HELP = ZERO
           DO J = 1, NSTREAMS
            JC = NSTOKES*(J-1)
            DO O2 = 1, NSTOKES
             JCOL = JC + O2
             TM1 =   HSOL%SAB  (I,J,O1,O2,N)   *    L_QVEC     (JCOL,Q)
             TM2 = L_HSOL%L_SAB(I,J,O1,O2,N,Q) * PSOL%QVEC_SAVE(JCOL)
             HELP = HELP + TM1 + TM2
            ENDDO
           ENDDO
           WSUM = L_QVEC(IROW,Q)
           TM2 = ( HELP - L_QSUMVEC(I,O1,Q) ) / SECBAR
           WDIF = LC_AVERAGE_SECANT(N,IB,Q) * TM3 + TM2
           L_WVEC_D = HALF * ( WSUM + WDIF )
           L_WVEC_U = HALF * ( WSUM - WDIF )
           L_PSOL%LC_BVEC(I,O1,N,Q)  = L_WVEC_D
           L_PSOL%LC_BVEC(I1,O1,N,Q) = L_WVEC_U
           IF ( O1 .GT. 2 ) L_PSOL%LC_BVEC(I1,O1,N,Q) = - L_WVEC_U

          ENDDO
         ENDDO
        ENDDO

!  assign solutions for plane/parallel & quasi-spherical case N = 1

      ELSE

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         I1 = I + NSTREAMS
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          DO Q = 1, N_PARAMETERS
           HELP = ZERO
           DO J = 1, NSTREAMS
            JC = NSTOKES*(J-1)
            DO O2 = 1, NSTOKES
             JCOL = JC + O2
             TM1 =   HSOL%SAB  (I,J,O1,O2,N)   *    L_QVEC     (JCOL,Q)
             TM2 = L_HSOL%L_SAB(I,J,O1,O2,N,Q) * PSOL%QVEC_SAVE(JCOL)
             HELP = HELP - TM1 - TM2
            ENDDO
           ENDDO
           WSUM = L_QVEC(IROW,Q)
           WDIF = ( - HELP - L_QSUMVEC(I,O1,Q) ) / SECBAR
           L_WVEC_D = HALF * ( WSUM + WDIF )
           L_WVEC_U = HALF * ( WSUM - WDIF )
           L_PSOL%LC_BVEC(I,O1,N,Q)  = L_WVEC_D
           L_PSOL%LC_BVEC(I1,O1,N,Q) = L_WVEC_U
           IF ( O1 .GT. 2 ) L_PSOL%LC_BVEC(I1,O1,N,Q) = - L_WVEC_U
          ENDDO
         ENDDO
        ENDDO

      ENDIF

!  Part 2.
!  =======

!  No part 2 (Profile Jacobians only). 
!  Version 2.8, goto statement removed.

!  Continuation point
! 2222 CONTINUE


!  Finish

      RETURN
      END SUBROUTINE MSSEG_LC_QBEAM_SOLUTION

!

      SUBROUTINE MSSEG_LC_UBEAM_SOLUTION ( &
        GIVEN_LAYER, FOURIER, IBEAM, DO_VARY, N_PARAMETERS,         & ! Input control
        DO_UPWELLING, DO_DNWELLING, NSTOKES, NSTREAMS, FLUX_FACTOR, & ! Input control
        DO_LAYER_SCATTERING, NMOMENTS, BEAM_CUTOFF, QUAD_HALFWTS,   & ! Input bookkeeping
        OMEGA_GREEK, L_OMEGA_GREEK, HGSF, PGSF, PSOL, L_PSOL )        ! Input Optical and TS I/O

!  Mk2.  6/21/21. Use module names now have the "_m" extension, Only use quantities you need
!  Mk2a. 9/21/21. New routine
!  Mk3.  1/31/21. Rearranged arguments lists and used type structures

      USE MSSEG_PARS_m, Only : MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS, &
                               MAX_USER_STREAMS, MAXMOMENTS, MAXSTREAMS_2, MAXSTRMSTKS, PI4, ZERO

      use MSSEG_HPTypes_Def_m
      use MSSEG_L_HPTypes_Def_m

      IMPLICIT NONE

!  Arguments
!  ---------

!  Control

      INTEGER, INTENT (IN) ::           GIVEN_LAYER
      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      LOGICAL, INTENT (IN) ::           DO_VARY
      INTEGER, INTENT (IN) ::           N_PARAMETERS

      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_DNWELLING
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS

      DOUBLE PRECISION, INTENT (IN) ::  FLUX_FACTOR

!  bookkeeping

      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )
      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_HALFWTS ( MAXSTREAMS )

!  input optical

      DOUBLE PRECISION, INTENT (IN) :: OMEGA_GREEK  ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  input type structures (Generalized spherical functions)

      type(Homog_GSF) , intent(in) :: HGSF
      type(Partic_GSF), intent(in) :: PGSF

!  input type structure (particular solutions)

      type(Partic_solutions), intent(in) :: PSOL

!  Output type structure (Linearized particular solutions)
 
      type(L_Partic_solutions) , intent(inout) :: L_PSOL

!  Local variables
!  ---------------

      LOGICAL ::          DO_LAYER, DO_FIRST
      INTEGER ::          IB, L, N, M, O1, O2, J, J1, Q

      DOUBLE PRECISION :: SUM2, T2, H, HH, F1
      DOUBLE PRECISION :: SP, SM, SPS, SPT, SGW(MAXSTOKES)

      DOUBLE PRECISION :: L_HELP_Q2(0:MAXMOMENTS,MAXSTOKES)
      DOUBLE PRECISION :: L_HELPSTOKES_BEAM(MAXSTOKES)

      SGW = (/ 1.0d0, 1.0d0, -1.0d0, -1.0d0 /)

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER
      F1 = FLUX_FACTOR / PI4
      IB = IBEAM
      DO_LAYER = .FALSE.

!  Initialize output

      L_PSOL%LC_UPAR_UP_2 = ZERO  ; L_PSOL%LC_UPAR_DN_2 = ZERO

!  Part 1. Linearizations for quantities varying in Layer N
!  ========================================================

!  Check existence
!  ---------------

!  Only a solution if the layer is active and not below Cutoff.

      DO_FIRST = ( N .LE. BEAM_CUTOFF(IB) ) .AND. DO_LAYER_SCATTERING(M,N)

!  If no solution or no variation, remove GOTO statement and return instead (NO PART 2)

      IF ( .NOT. DO_VARY .OR. .NOT. DO_FIRST ) RETURN

!  start parameter loop

      DO Q = 1, N_PARAMETERS

!  For each moment, do inner sum over computational angles
!   (Linearization of Diffuse scattering of beam)

        DO L = M, NMOMENTS
         DO O1 = 1, NSTOKES
          SPS = ZERO
          DO J = 1, NSTREAMS
           J1 = J + NSTREAMS
           SP = ZERO
           SM = ZERO
           H = QUAD_HALFWTS(J)
           DO O2 = 1, NSTOKES
            HH = H * SGW(O2)
            SP = SP + H  * HGSF%PI_XQP    (L,J,O1,O2)*L_PSOL%LC_BVEC(J,O2,N,Q)
            SM = SM + HH * HGSF%PI_XQM_PRE(L,J,O1,O2)*L_PSOL%LC_BVEC(J1,O2,N,Q)
           ENDDO
           SPS = SPS + SP + SM
          ENDDO
          L_HELPSTOKES_BEAM(O1) = SPS
         ENDDO
         DO O1 = 1, NSTOKES
          SPT = ZERO
          DO O2 = 1, NSTOKES
           SPT = SPT + OMEGA_GREEK(L,N,O1,O2)   * L_HELPSTOKES_BEAM(O2) &
                   + L_OMEGA_GREEK(L,N,O1,O2,Q) * PSOL%HELPSTOKES_BEAM(L,O2)
          ENDDO
          L_HELP_Q2(L,O1) = SPT
         ENDDO
        ENDDO

!  Now sum over all harmonic contributions (Upwelling)

       IF ( DO_UPWELLING  ) THEN
         DO O1 = 1, NSTOKES
          T2 = ZERO
          DO L = M, NMOMENTS
           SUM2 = ZERO
           DO O2 = 1, NSTOKES
            SUM2 = SUM2 + L_HELP_Q2(L,O2)*PGSF%PI_XUM(L,IBEAM,O1,O2)
           ENDDO
           T2 = T2 + SUM2
          ENDDO
          L_PSOL%LC_UPAR_UP_2(O1,Q) = T2
         ENDDO
       ENDIF

!  Now sum over all harmonic contributions (Downwelling)

       IF ( DO_DNWELLING  ) THEN
         DO O1 = 1, NSTOKES
          T2 = ZERO
          DO L = M, NMOMENTS
           SUM2 = ZERO
           DO O2 = 1, NSTOKES
            SUM2 = SUM2 + L_HELP_Q2(L,O2)*PGSF%PI_XUP(L,IBEAM,O1,O2)
           ENDDO
           T2 = T2 + SUM2
          ENDDO
          L_PSOL%LC_UPAR_DN_2(O1,Q) = T2
         ENDDO
       ENDIF

!  end parameter loop

      ENDDO

!  Part 2. Version 2.8, removed.
!  Continuation point
! 2222 CONTINUE

!  Finish

      RETURN
      END SUBROUTINE MSSEG_LC_UBEAM_SOLUTION

!  Finish module

      END MODULE MSSeg_LC_Solutions_m
