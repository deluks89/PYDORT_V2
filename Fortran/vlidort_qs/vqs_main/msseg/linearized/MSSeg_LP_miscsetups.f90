
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

! ###############################################################
! #                                                             #
! # MSSEG_Adapted Subroutines in this Module                    #
! #                                                             #
! #            MSSEG_LP_MISCSETUPS (calling:--)                 #
! #                                                             #
! #              MSSEG_LP_DELTAMSCALE     (private)             #
! #              MSSEG_LP_QSPREP          (private)             #
! #                                                             #
! ###############################################################

!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- ASYMTX_TOLERANCE now introduced and passed into the QHOM call to ASYMTX
!   -- Flux output (LOCAL_FLUXES) controlled by DO_LOCAL_FLUXES input flag
!   -- DO_FOURIER0_NSTOKES2 performance flag introduced.
!   -- Call subroutine specially written to get the discrete-ordinate field for Fluxes
!   -- Use public routine GETQUAD2 instead of GAULEG

!  Mk2a. LP Subroutines created

!  Mk3. 4/11/23. (1/31/22 created). Modified for MK3 use

      MODULE MSSeg_LP_Miscsetups_m

      PRIVATE :: MSSEG_LP_DELTAMSCALE, MSSEG_LP_QSPREP
      PUBLIC  :: MSSEG_LP_MISCSETUPS

      CONTAINS


      SUBROUTINE MSSEG_LP_MISCSETUPS( &
        DO_DELTAM_SCALING, DO_PLANE_PARALLEL, NLAYERS, NBEAMS, & ! Input
        NCHOICE, LVARY_FLAG, LVARY_NUMBER, OMEGA_TOTAL_INPUT,  & ! Input
        DELTAU_VERT, DELTAU_SLANT, TRUNC_FACTOR, FAC1,         & ! Input
        BEAM_CUTOFF, AVERAGE_SECANT, T_DELT_MUBAR,             & ! Input
        L_SFI_OPDEPS, L_TRUNC_FACTOR, L_DELTAU_VERT,           & ! Input
        L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT,              & ! Input
        L_DELTAU_SLANT, L_SFI_OPDEPS_SCALED,                   & ! Output
        LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR )   ! Output

!  Mk2a, 9/21/21. LAP top-level subroutine created by R. Spurr
!   -- Modeled after LAP subroutine of the same name

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- Remove SFI_OPDEPS array from argument list.

!  Mk3. 1/31/22. Re-engineered software.

      USE MSSeg_PARS_m, Only : MAXSTOKES, MAXLAYERS, MAXSTREAMS, MAXBEAMS, MAX_ATMOSWFS, MAX_GEOMETRIES, &
                               MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, ZERO

      IMPLICIT NONE

!  Inputs
!  ------

      LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          NCHOICE

      INTEGER, INTENT (IN) ::          LVARY_NUMBER ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          LVARY_FLAG   ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_TOTAL_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT_INPUT ( MAX_ATMOSWFS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT    ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_SFI_OPDEPS  ( MAX_GEOMETRIES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_TRUNC_FACTOR ( MAX_ATMOSWFS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: TRUNC_FACTOR ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: FAC1 ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  output

      DOUBLE PRECISION, INTENT (OUT) :: L_DELTAU_SLANT   ( MAX_ATMOSWFS, MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: L_SFI_OPDEPS_SCALED  ( MAX_GEOMETRIES, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  miscellaneous setup operations, Master routine

!  Delta-m scaling of input quantities
!   Mk2, 6/21/21. remove SFI_OPDEPS array from argument list
!   Mk3, 1/31/22. Only for the beam stuff

      CALL MSSEG_LP_DELTAMSCALE ( &
        DO_DELTAM_SCALING, NLAYERS, NBEAMS, NCHOICE, LVARY_FLAG, LVARY_NUMBER,  &
        OMEGA_TOTAL_INPUT, DELTAU_SLANT, TRUNC_FACTOR, FAC1,                    &
        L_OMEGA_TOTAL_INPUT, L_DELTAU_VERT_INPUT, L_TRUNC_FACTOR, L_SFI_OPDEPS, &
        L_DELTAU_SLANT, L_SFI_OPDEPS_SCALED )

!  Prepare quasi-spherical attenuation
!   Mk2a. 9/21/21. Only this routine is new
!   Mk3.  1/31/22. Only for the beam stuff

      CALL MSSEG_LP_QSPREP ( &
        DO_PLANE_PARALLEL, NLAYERS, NBEAMS, LVARY_NUMBER, &
        DELTAU_VERT, DELTAU_SLANT, L_DELTAU_VERT,         &
        BEAM_CUTOFF,  T_DELT_MUBAR, AVERAGE_SECANT,       &
        LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR )

!  Finish

      RETURN
      END SUBROUTINE MSSEG_LP_MISCSETUPS

!

      SUBROUTINE MSSEG_LP_DELTAMSCALE ( &
        DO_DELTAM_SCALING, NLAYERS, NBEAMS, NCHOICE, LVARY_FLAG, LVARY_NUMBER,  &
        OMEGA_TOTAL_INPUT, DELTAU_SLANT, TRUNC_FACTOR, FAC1,                    &
        L_OMEGA_TOTAL_INPUT, L_DELTAU_VERT_INPUT, L_TRUNC_FACTOR, L_SFI_OPDEPS, &
        L_DELTAU_SLANT, L_SFI_OPDEPS_SCALED )

!  module, dimensions and numbers

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- remove SFI_OPDEPS array from argument list

      USE MSSeg_PARS_m, Only : MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES, MAXBEAMS, zero

      IMPLICIT NONE

!  inputs

      LOGICAL, INTENT (IN) ::           DO_DELTAM_SCALING
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           NCHOICE

      LOGICAL, INTENT (IN) ::           LVARY_FLAG   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           LVARY_NUMBER ( MAXLAYERS )

!  optical

      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: TRUNC_FACTOR ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: FAC1 ( MAXLAYERS )

!  Linearized optical

      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_TOTAL_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_TRUNC_FACTOR ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_SFI_OPDEPS  ( MAX_GEOMETRIES, MAX_ATMOSWFS )
!  output

      DOUBLE PRECISION, INTENT (OUT) :: L_DELTAU_SLANT ( MAX_ATMOSWFS, MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: L_SFI_OPDEPS_SCALED  ( MAX_GEOMETRIES, MAX_ATMOSWFS )

!  local variables

      DOUBLE PRECISION :: L_FAC1, DELS
      INTEGER          :: N, Q, K, IB

!  slant optical thickness values
!    (Revert to the input values of deltau_slant)

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          DO K = 1, N
            DELS = DELTAU_SLANT(N,K,IB)/FAC1(K)
            IF ( LVARY_FLAG(K) ) THEN
              DO Q = 1, LVARY_NUMBER(K)
                L_DELTAU_SLANT(Q,N,K,IB) = L_DELTAU_VERT_INPUT(Q,K) * DELS
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!  DELTAM SCALING
!  ==============

!  Rob fix, 10/17/21. Must initialize

      L_SFI_OPDEPS_SCALED = ZERO

!  New section added 21 December 2005 by R. Spurr
!   Linearization of the scaled delta-M inputs

      IF ( DO_DELTAM_SCALING ) THEN

!  Scale linearized slant path optical thickness values

         DO N = 1, NLAYERS
            DO K = 1, N
              IF ( LVARY_FLAG(K) ) THEN
                DO Q = 1, LVARY_NUMBER(K)
                  L_FAC1 = L_TRUNC_FACTOR(Q,K) &
                           + TRUNC_FACTOR(K)   * L_OMEGA_TOTAL_INPUT(Q,K)
                  L_FAC1 = L_FAC1 *  OMEGA_TOTAL_INPUT(K)
                  DO IB = 1, NBEAMS
                    DELS = DELTAU_SLANT(N,K,IB)/FAC1(K)
                    L_DELTAU_SLANT(Q,N,K,IB) = - L_FAC1 * DELS + &
                                 FAC1(K) * L_DELTAU_SLANT(Q,N,K,IB)
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
         ENDDO

!  Scale partial layer stuff. No !!!!!!!!!!
!    OPDEPS from FO code are Already scaled! Bug 18 November 2015
!    Just need to copy the input to output

         DO IB = 1, NBEAMS
           K = NCHOICE
           IF ( LVARY_FLAG(K) ) THEN
             DO Q = 1, LVARY_NUMBER(K)
!               L_SFI_OPDEPS_SCALED(IB,Q) = L_SFI_OPDEPS(IB,Q) * FAC1(K) &
!                                        +   SFI_OPDEPS(IB)   * L_TRUNC_FACTOR(Q,K)
               L_SFI_OPDEPS_SCALED(IB,Q) = L_SFI_OPDEPS(IB,Q)
             ENDDO
           ENDIF
         ENDDO

!  NO DELTAM SCALING
!  =================

!  move input geophysical variables to Workspace quantities

      ELSE

!  Copy partial layer stuff

         DO IB = 1, NBEAMS
           K = NCHOICE
           IF ( LVARY_FLAG(K) ) THEN
             DO Q = 1, LVARY_NUMBER(K)
               L_SFI_OPDEPS_SCALED(IB,Q) = L_SFI_OPDEPS(IB,Q) 
             ENDDO
           ENDIF
         ENDDO

      ENDIF

!  Finish module

      RETURN
      END SUBROUTINE MSSEG_LP_DELTAMSCALE

!

      SUBROUTINE MSSEG_LP_QSPREP ( &
        DO_PLANE_PARALLEL, NLAYERS, NBEAMS, LVARY_NUMBER, &
        DELTAU_VERT, DELTAU_SLANT, L_DELTAU_VERT,         &
        BEAM_CUTOFF, T_DELT_MUBAR, AVERAGE_SECANT,        &
        LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR )

!  Mk2a. 9/21/21. LC subroutine created by R. Spurr
!     -- Modeled after the LP subroutine of the same name

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need

      USE MSSeg_PARS_m, Only : MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, ZERO

      IMPLICIT NONE

!  Input control

      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          LVARY_NUMBER ( MAXLAYERS )

!  input otpical

      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  input beam parameterization

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  output linearized parameterization

      DOUBLE PRECISION, INTENT (OUT) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          N, Q, K, IB
      DOUBLE PRECISION :: WDEL, VAR, RHO, FAC, DELT, LAMDA

!  Initialize outputs

      LP_T_DELT_MUBAR   = ZERO
      LP_INITIAL_TRANS  = ZERO
      LP_AVERAGE_SECANT = ZERO

!  linearization of Initial transmittances (Profile linearization)
!  =======================================

!   Bug fixed, 12 August 2005 for linearization of INITIAL_TRANS
!         Use Logarithmic derivative !!!!
!         Reason: avoids exceptions if INITIAL_TRANS underflows

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          IF ( N .GT. 1 ) THEN
            IF  (N.LE.BEAM_CUTOFF(IB) ) THEN
              DO K = 1, N-1
                DO Q = 1, LVARY_NUMBER(N)
                  LP_INITIAL_TRANS(N,K,IB,Q) = - L_DELTAU_VERT(Q,K) * DELTAU_SLANT(N-1,K,IB)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO

!  linearization of average secants for pseudo-spherical case
!  ==========================================================

!   (average secant = 1/mu-0 = constant for plane parallel)

       IF ( .NOT. DO_PLANE_PARALLEL ) THEN
        DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            IF ( N .GT. 1 ) THEN
              IF  (N.LE.BEAM_CUTOFF(IB) ) THEN
                DELT  = DELTAU_VERT(N)
                LAMDA = AVERAGE_SECANT(N,IB)
                FAC   = ( DELTAU_SLANT(N,N,IB) / DELT ) - LAMDA
                DO Q = 1, LVARY_NUMBER(N)
                  LP_AVERAGE_SECANT(N,N,IB,Q) = L_DELTAU_VERT(Q,N) * FAC
                ENDDO
                DO K = 1, N-1
                  FAC = ( DELTAU_SLANT(N,K,IB) - DELTAU_SLANT(N-1,K,IB) ) / DELT
                  DO Q = 1, LVARY_NUMBER(K)
                    LP_AVERAGE_SECANT(N,K,IB,Q) = L_DELTAU_VERT(Q,K) * FAC
                  ENDDO
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!  Profile Linearization of Whole layer Transmittance factors
!  ==========================================================

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS

         WDEL  = T_DELT_MUBAR(N,IB)
         VAR   = - DELTAU_VERT(N) * WDEL
         LAMDA = AVERAGE_SECANT(N,IB)
         FAC   = VAR * AVERAGE_SECANT(N,IB)

!  Pseudo-spherical

          IF ( .NOT. DO_PLANE_PARALLEL ) THEN

            IF ( N .EQ. 1 ) THEN
              DO Q = 1, LVARY_NUMBER(N)
                LP_T_DELT_MUBAR(N,N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
              ENDDO
            ELSE
              IF  ( N.LE.BEAM_CUTOFF(IB) ) THEN
                DO Q = 1, LVARY_NUMBER(N)
                  RHO = LP_AVERAGE_SECANT(N,N,IB,Q)
                  LP_T_DELT_MUBAR(N,N,IB,Q) = L_DELTAU_VERT(Q,N) * FAC + VAR * RHO
                ENDDO
                DO K = 1, N-1
                  DO Q = 1, LVARY_NUMBER(K)
                    RHO = LP_AVERAGE_SECANT(N,K,IB,Q)
                    LP_T_DELT_MUBAR(N,K,IB,Q) = VAR * RHO
                  ENDDO
                ENDDO
              ENDIF
            ENDIF

!  Plane-parallel

          ELSE IF ( DO_PLANE_PARALLEL ) THEN

            IF  (N.LE.BEAM_CUTOFF(IB) ) THEN
              DO Q = 1, LVARY_NUMBER(N)
                LP_T_DELT_MUBAR(N,N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
              ENDDO
            ENDIF

          ENDIF

!  end layer and beam loops

        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE MSSEG_LP_QSPREP

!  Finish Module

END MODULE MSSeg_LP_Miscsetups_m

