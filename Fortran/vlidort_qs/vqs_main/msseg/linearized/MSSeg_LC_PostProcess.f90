
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

!  Mk3. 1/31/22. Completely new module, with 2 subroutines taking over usual VLIDORT Jacobian routines
!  Mk3. 4/11/23. Modified for MK3 use

      MODULE MSSeg_LC_PostProcess_m

      PRIVATE
      PUBLIC :: MSSEG_LC_POSTPROCESS_LAYERMSST, MSSEG_LC_POSTPROCESS_SURFACEMSST, GET_LC_QuadJacs_Levels

      CONTAINS


      SUBROUTINE MSSEG_LC_POSTPROCESS_LAYERMSST ( &
        DO_UPWELLING, DO_DNWELLING, IGEOM, FOURIER, NCHOICE, NSTOKES,     & ! Input Control
        N_COLUMN_WFS, DELTAU_VERT, L_DELTAU_VERT, L_SFI_OPDEPS_SCALED,    & ! Input control/Optical
        BEAM_CUTOFF, AVERAGE_SECANT, LC_AVERAGE_SECANT, LC_INITIAL_TRANS, & ! Input Beam stuff
        HSOL, PSOL, PBVP, L_HSOL, L_PSOL, L_PBVP,                         & ! Input TS
        TAU_UP, TAU_DN, T_UTUP_EIGEN, T_UTDN_EIGEN, T_UTDN_MUBAR,         & ! Input partial quantities
        L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, LC_T_UTDN_MUBAR, LC_MSST_F )

!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- Add NCHOICE (debug only)

!  Mk3, 1/31/22. New subroutine made from existing material.

!  10/4/23. Ordering of TAU_UP and TAU_DN changed for consistency

      USE MSSeg_PARS_m, only : MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAXLAYERS, MAX_GEOMETRIES, &
                               MAXEVALUES, MAX_ATMOSWFS, MAX_TAU_SPATH, MAX_TAU_QPATH, ZERO

      use MSSEG_HPTypes_Def_m
      use MSSEG_L_HPTypes_Def_m

      IMPLICIT NONE

!  INPUTS
!  ======

!  Flags

      LOGICAL, INTENT (IN) ::            DO_UPWELLING
      LOGICAL, INTENT (IN) ::            DO_DNWELLING

!  Numbers

      INTEGER, INTENT (IN) ::            FOURIER
      INTEGER, INTENT (IN) ::            IGEOM
      INTEGER, INTENT (IN) ::            NCHOICE
      INTEGER, INTENT (IN) ::            NSTOKES
      INTEGER, INTENT (IN) ::            N_COLUMN_WFS

!  Optical inputs

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT   ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_SFI_OPDEPS_SCALED  ( MAX_GEOMETRIES, MAX_ATMOSWFS )

!  Beam inputs

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  input solution type structures

      type(Homog_Solutions) , intent(in) :: HSOL
      type(Partic_Solutions), intent(in) :: PSOL

!  Input BVProblem solution type structure

      type(Partic_BVPSolve)  , intent(in) :: PBVP

!  input type structure (Linearized homogeneous solutions)
 
      type(L_Homog_Solutions) , intent(in) :: L_HSOL

!  Input type structures (Linearized particular solutions and BVP)
 
      type(L_Partic_Solutions) , intent(in) :: L_PSOL
      type(L_Partic_BVPSolve)  , intent(in) :: L_PBVP

!  input from regular postprocess subroutine

      DOUBLE PRECISION, INTENT (IN)  :: TAU_DN, TAU_UP
      DOUBLE PRECISION, INTENT (IN)  :: T_UTUP_EIGEN ( MAXEVALUES )
      DOUBLE PRECISION, INTENT (IN)  :: T_UTDN_EIGEN ( MAXEVALUES )
      DOUBLE PRECISION, INTENT (IN)  :: T_UTDN_MUBAR

!  OUTPUTS
!  =======

      DOUBLE PRECISION, INTENT (OUT)  :: LC_MSST_F      ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT)  :: L_T_UTUP_EIGEN ( MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT)  :: L_T_UTDN_EIGEN ( MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT)  :: LC_T_UTDN_MUBAR ( MAX_ATMOSWFS )

!  Local variables
!  ---------------

      DOUBLE PRECISION :: L_TRANSH(MAX_ATMOSWFS), L_TAU_DN(MAX_ATMOSWFS), L_TAU_UP(MAX_ATMOSWFS)
      DOUBLE PRECISION :: LUXR,  MUXR,  LLUXR,  MLUXR,  NUXR,  PUXR, H1, H2, SHOM_R
      DOUBLE PRECISION :: LUXR1, MUXR1, LLUXR1, MLUXR1, NUXR1, PUXR1, SHOM_CR
      DOUBLE PRECISION :: LUXR2, MUXR2, LLUXR2, MLUXR2, NUXR2, PUXR2, SPHER, L_SPHER, SPAR
      DOUBLE PRECISION :: HELP, HELP_CR, HELP_CI, L_HELP, L_HELP_CR, L_HELP_CI, COS_CI, SIN_CI

      INTEGER          :: O1, K, K0, K1, K2, KO1, N, Q

!  Optical depths. Must used scaled quantities.
!    Rob Fix, 5/30/18. L_DELTAU_VERT is double normalized

      DO Q = 1, N_COLUMN_WFS
         L_TAU_DN(Q) = L_SFI_OPDEPS_SCALED(IGEOM,Q)
         L_TAU_UP(Q) = L_DELTAU_VERT(Q,NCHOICE) * DELTAU_VERT(NCHOICE) - L_SFI_OPDEPS_SCALED(IGEOM,Q)
      ENDDO

!  This checks out, 10/4/23
!if ( Fourier.eq.1)write(1002,'(2i4,1p4e18.8)')NCHOICE,IGEOM,L_DELTAU_VERT(1:2,NCHOICE),L_SFI_OPDEPS_SCALED(IGEOM,1:2)

!  Linearized Transmittances
!  -------------------------

      N = NCHOICE
      KO1 = HSOL%K_REAL(N) + 1

!  Initialize them

      L_T_UTUP_EIGEN = zero ; L_T_UTDN_EIGEN = zero ; LC_T_UTDN_MUBAR = zero
      DO Q = 1, N_COLUMN_WFS
         L_TRANSH(Q) = LC_INITIAL_TRANS(NCHOICE,IGEOM,Q)
      ENDDO

!  checks out 10/4/23
!if ( fourier.eq.0)write(1000,'(2i4,1p2e18.8)')NCHOICE,IGEOM,L_TRANSH(1:2)

!  Real eigenvalues

      DO K = 1, HSOL%K_REAL(N)
         DO Q = 1, N_COLUMN_WFS
            L_HELP = L_HSOL%L_KEIGEN(K,N,Q) * TAU_DN + HSOL%KEIGEN(K,N) * L_TAU_DN(Q)
            L_T_UTDN_EIGEN(K,Q) = - L_HELP * T_UTDN_EIGEN(K)
            L_HELP = L_HSOL%L_KEIGEN(K,N,Q) * TAU_UP + HSOL%KEIGEN(K,N) * L_TAU_UP(Q)
            L_T_UTUP_EIGEN(K,Q) = - L_HELP * T_UTUP_EIGEN(K)
         ENDDO
!if ( fourier.eq.0.and.igeom.eq.1)write(1001,'(3i4,1p2e18.8)')NCHOICE,IGEOM,K,L_T_UTUP_EIGEN(K,1:2)
!if ( fourier.eq.0.and.igeom.eq.1)write(1001,'(3i4,1p4e18.8)')NCHOICE,IGEOM,K,HSOL%KEIGEN(K,N), L_TAU_UP(1:2),T_UTUP_EIGEN(K)
!if ( fourier.eq.0.and.igeom.eq.1)write(1001,'(3i4,1p4e18.8)')NCHOICE,IGEOM,K,L_HSOL%L_KEIGEN(K,N,1:2),L_TAU_UP(1:2)
      ENDDO

!  Complex eigenvalues

      DO K = 1, HSOL%K_COMPLEX(N)
         K0 = 2*K-2 ; K1 = KO1 + K0 ;  K2 = K1  + 1
         HELP = TAU_DN * HSOL%KEIGEN(K1,N)
         IF ( HELP .LE. MAX_TAU_QPATH ) THEN
            HELP_CR = DEXP ( - HELP )
            HELP_CI = TAU_DN * HSOL%KEIGEN(K2,N)
            COS_CI = DCOS ( HELP_CI )
            SIN_CI = DSIN ( HELP_CI )
            DO Q = 1, N_COLUMN_WFS
               L_HELP    = L_HSOL%L_KEIGEN(K1,N,Q) * TAU_DN + HSOL%KEIGEN(K1,N) * L_TAU_DN(Q)
               L_HELP_CR = L_HELP * HELP_CR
               L_HELP_CI = L_HSOL%L_KEIGEN(K2,N,Q) * TAU_DN + HSOL%KEIGEN(K2,N) * L_TAU_DN(Q)
               L_T_UTDN_EIGEN(K1,Q) =  - HELP_CR * SIN_CI * L_HELP_CI + L_HELP_CR * COS_CI
               L_T_UTDN_EIGEN(K2,Q) =  - HELP_CR * COS_CI * L_HELP_CI - L_HELP_CR * SIN_CI
            ENDDO
         ENDIF
         HELP = TAU_UP * HSOL%KEIGEN(K1,N)
         IF ( HELP .LE. MAX_TAU_QPATH ) THEN
            HELP_CR = DEXP ( - HELP )
            HELP_CI = TAU_UP * HSOL%KEIGEN(K2,N)
            COS_CI = DCOS ( HELP_CI )
            SIN_CI = DSIN ( HELP_CI )
            DO Q = 1, N_COLUMN_WFS 
               L_HELP    = L_HSOL%L_KEIGEN(K1,N,Q) * TAU_UP + HSOL%KEIGEN(K1,N) * L_TAU_UP(Q)
               L_HELP_CR = L_HELP * HELP_CR
               L_HELP_CI = L_HSOL%L_KEIGEN(K2,N,Q) * TAU_UP + HSOL%KEIGEN(K2,N) * L_TAU_UP(Q)
               L_T_UTUP_EIGEN(K1,Q) =  - HELP_CR * SIN_CI * L_HELP_CI + L_HELP_CR * COS_CI
               L_T_UTUP_EIGEN(K2,Q) =  - HELP_CR * COS_CI * L_HELP_CI - L_HELP_CR * SIN_CI
            ENDDO
         ENDIF
      ENDDO

!  Particular integral term

      IF ( N .LE. BEAM_CUTOFF(IGEOM) ) THEN
         SPHER = TAU_DN * AVERAGE_SECANT(N,IGEOM)
         DO Q = 1, N_COLUMN_WFS
            L_SPHER = L_TAU_DN(Q) *    AVERAGE_SECANT(N,IGEOM) + &
                        TAU_DN    * LC_AVERAGE_SECANT(N,IGEOM,Q)
            IF ( SPHER .LE. MAX_TAU_SPATH ) THEN        
               LC_T_UTDN_MUBAR(Q) = ( L_TRANSH(Q) - L_SPHER ) * T_UTDN_MUBAR
            ENDIF
         ENDDO

!  Checks out, 10/4/23
!if(fourier.eq.0)write(1003,'(2i4,1p7e17.7)')nchoice,igeom,T_UTDN_MUBAR,L_TRANSH(1),TAU_DN,L_TAU_DN(1),&
!                                AVERAGE_SECANT(N,IGEOM) ,LC_AVERAGE_SECANT(N,IGEOM,1),LC_T_UTDN_MUBAR(1)

      ENDIF

!  Upwelling Jacobian Contributions
!  --------------------------------

      IF ( DO_UPWELLING ) THEN

!  Homogeneous solutions
!  =====================

!    MK2, 6/21/21.  Use local bookkeeping variable LL_NSTK

!  Loops

         DO O1 = 1, NSTOKES
            DO Q = 1, N_COLUMN_WFS

!  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, HSOL%K_REAL(N)
               LUXR  = PBVP%LCON(K,N)   *   PSOL%UHOM_UPDN(IGEOM,O1,K)
               MUXR  = PBVP%MCON(K,N)   *   PSOL%UHOM_UPUP(IGEOM,O1,K)
               LLUXR = PBVP%LCON(K,N)   * L_PSOL%L_UHOM_UPDN(IGEOM,O1,K,Q)
               MLUXR = PBVP%MCON(K,N)   * L_PSOL%L_UHOM_UPUP(IGEOM,O1,K,Q)
                  NUXR  = L_PBVP%NCON(K,N,Q) *   PSOL%UHOM_UPDN(IGEOM,O1,K)
                  PUXR  = L_PBVP%PCON(K,N,Q) *   PSOL%UHOM_UPUP(IGEOM,O1,K)
                  H1 = ( NUXR + LLUXR ) * T_UTDN_EIGEN(K) + LUXR * L_T_UTDN_EIGEN(K,Q)
                  H2 = ( PUXR + MLUXR ) * T_UTUP_EIGEN(K) + MUXR * L_T_UTUP_EIGEN(K,Q)
                  SHOM_R = SHOM_R + H1 + H2
               ENDDO

!  Complex homogeneous solutions

               SHOM_CR = ZERO
               DO K = 1, HSOL%K_COMPLEX(N)
                 K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
 
                 LUXR1 =  PBVP%LCON(K1,N)   *   PSOL%UHOM_UPDN(IGEOM,O1,K1) - &
                          PBVP%LCON(K2,N)   *   PSOL%UHOM_UPDN(IGEOM,O1,K2)
                 LUXR2 =  PBVP%LCON(K1,N)   *   PSOL%UHOM_UPDN(IGEOM,O1,K2) + &
                          PBVP%LCON(K2,N)   *   PSOL%UHOM_UPDN(IGEOM,O1,K1)
                 MUXR1 =  PBVP%MCON(K1,N)   *   PSOL%UHOM_UPUP(IGEOM,O1,K1) - &
                          PBVP%MCON(K2,N)   *   PSOL%UHOM_UPUP(IGEOM,O1,K2)
                 MUXR2 =  PBVP%MCON(K1,N)   *   PSOL%UHOM_UPUP(IGEOM,O1,K2) + &
                          PBVP%MCON(K2,N)   *   PSOL%UHOM_UPUP(IGEOM,O1,K1)
 
                 LLUXR1 = PBVP%LCON(K1,N)   * L_PSOL%L_UHOM_UPDN(IGEOM,O1,K1,Q) - &
                          PBVP%LCON(K2,N)   * L_PSOL%L_UHOM_UPDN(IGEOM,O1,K2,Q)
                 LLUXR2 = PBVP%LCON(K1,N)   * L_PSOL%L_UHOM_UPDN(IGEOM,O1,K2,Q) + &
                          PBVP%LCON(K2,N)   * L_PSOL%L_UHOM_UPDN(IGEOM,O1,K1,Q)
                 MLUXR1 = PBVP%MCON(K1,N)   * L_PSOL%L_UHOM_UPUP(IGEOM,O1,K1,Q) - &
                          PBVP%MCON(K2,N)   * L_PSOL%L_UHOM_UPUP(IGEOM,O1,K2,Q)
                 MLUXR2 = PBVP%MCON(K1,N)   * L_PSOL%L_UHOM_UPUP(IGEOM,O1,K2,Q) + &
                          PBVP%MCON(K2,N)   * L_PSOL%L_UHOM_UPUP(IGEOM,O1,K1,Q)
 
                 NUXR1 =  L_PBVP%NCON(K1,N,Q) *   PSOL%UHOM_UPDN(IGEOM,O1,K1) - &
                          L_PBVP%NCON(K2,N,Q) *   PSOL%UHOM_UPDN(IGEOM,O1,K2)
                 NUXR2 =  L_PBVP%NCON(K1,N,Q) *   PSOL%UHOM_UPDN(IGEOM,O1,K2) + &
                          L_PBVP%NCON(K2,N,Q) *   PSOL%UHOM_UPDN(IGEOM,O1,K1)
                 PUXR1 =  L_PBVP%PCON(K1,N,Q) *   PSOL%UHOM_UPUP(IGEOM,O1,K1) - &
                          L_PBVP%PCON(K2,N,Q) *   PSOL%UHOM_UPUP(IGEOM,O1,K2)
                 PUXR2 =  L_PBVP%PCON(K1,N,Q) *   PSOL%UHOM_UPUP(IGEOM,O1,K2) + &
                          L_PBVP%PCON(K2,N,Q) *   PSOL%UHOM_UPUP(IGEOM,O1,K1)

                 H1 =    ( NUXR1 + LLUXR1 ) * T_UTDN_EIGEN(K1) &
                       - ( NUXR2 + LLUXR2 ) * T_UTDN_EIGEN(K2) &
                       +  LUXR1       * L_T_UTDN_EIGEN(K1,Q) &
                       -  LUXR2       * L_T_UTDN_EIGEN(K2,Q)
                 H2 =    ( PUXR1 + MLUXR1 ) * T_UTUP_EIGEN(K1) &
                       - ( PUXR2 + MLUXR2 ) * T_UTUP_EIGEN(K2) &
                       +  MUXR1       * L_T_UTUP_EIGEN(K1,Q) &
                       -  MUXR2       * L_T_UTUP_EIGEN(K2,Q)
 
                 SHOM_CR = SHOM_CR + H1 + H2

               ENDDO

!  Final contribution and end loops

               LC_MSST_F(Q,IGEOM,O1) = SHOM_R + SHOM_CR
            ENDDO
         ENDDO

!  Add particular integral contribution
!  ====================================

!    MK2, 6/21/21.  Use local bookkeeping variable LL_NSTK

         DO O1 = 1, NSTOKES
            DO Q = 1, N_COLUMN_WFS
               SPAR =  L_PSOL%LC_UPAR_UP_2(O1,Q) *    T_UTDN_MUBAR + &
                         PSOL%UPAR_UP_2(O1)      * LC_T_UTDN_MUBAR(Q)
               LC_MSST_F(Q,IGEOM,O1) = LC_MSST_F(Q,IGEOM,O1) + SPAR
            ENDDO

         ENDDO

!  Checks out, 10/4/23
!if ( fourier.eq.0)write(1003,'(2i4,1p4e18.8)')nchoice,igeom,L_PSOL%LC_UPAR_UP_2(1,1:2),LC_T_UTDN_MUBAR(1:2)

!  End upwelling

      ENDIF

!  Downwelling Jacobian Contributions
!  ----------------------------------

      IF ( DO_DNWELLING ) THEN

!  Homogeneous solutions
!  =====================

!  Loops
!    MK2, 6/21/21.  Use local bookkeeping variable LL_NSTK

         DO O1 = 1, NSTOKES
            DO Q = 1, N_COLUMN_WFS

!  Real homogeneous solutions
!   10/14/21. Oh Boy bug, with H1 and H2. Wrong use of UTUP/UTDN transmittances

               SHOM_R = ZERO
               DO K = 1, HSOL%K_REAL(N)
                  LUXR  = PBVP%LCON(K,N)   *   PSOL%UHOM_DNDN(IGEOM,O1,K)
                  MUXR  = PBVP%MCON(K,N)   *   PSOL%UHOM_DNUP(IGEOM,O1,K)
                  LLUXR = PBVP%LCON(K,N)   * L_PSOL%L_UHOM_DNDN(IGEOM,O1,K,Q)
                  MLUXR = PBVP%MCON(K,N)   * L_PSOL%L_UHOM_DNUP(IGEOM,O1,K,Q)
                  NUXR  = L_PBVP%NCON(K,N,Q) *   PSOL%UHOM_DNDN(IGEOM,O1,K)
                  PUXR  = L_PBVP%PCON(K,N,Q) *   PSOL%UHOM_DNUP(IGEOM,O1,K)
!                  H1 = ( NUXR + LLUXR ) *   T_UTUP_EIGEN(K) + LUXR    * L_T_UTUP_EIGEN(K,Q)
!                  H2 = ( PUXR + MLUXR ) *   T_UTDN_EIGEN(K) + MUXR    * L_T_UTDN_EIGEN(K,Q)
                  H1 = ( NUXR + LLUXR ) * T_UTDN_EIGEN(K) + LUXR * L_T_UTDN_EIGEN(K,Q)
                  H2 = ( PUXR + MLUXR ) * T_UTUP_EIGEN(K) + MUXR * L_T_UTUP_EIGEN(K,Q)
                  SHOM_R = SHOM_R + H1 + H2
               ENDDO

!  Complex homogeneous solutions
!   10/14/21. Oh Boy bug, with H1 and H2. Wrong use of UTUP/UTDN transmittances

               SHOM_CR = ZERO
               DO K = 1, HSOL%K_COMPLEX(N)
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
 
                  LUXR1 =  PBVP%LCON(K1,N)   *   PSOL%UHOM_DNDN(IGEOM,O1,K1) - &
                           PBVP%LCON(K2,N)   *   PSOL%UHOM_DNDN(IGEOM,O1,K2)
                  LUXR2 =  PBVP%LCON(K1,N)   *   PSOL%UHOM_DNDN(IGEOM,O1,K2) + &
                           PBVP%LCON(K2,N)   *   PSOL%UHOM_DNDN(IGEOM,O1,K1)
                  MUXR1 =  PBVP%MCON(K1,N)   *   PSOL%UHOM_DNUP(IGEOM,O1,K1) - &
                           PBVP%MCON(K2,N)   *   PSOL%UHOM_DNUP(IGEOM,O1,K2)
                  MUXR2 =  PBVP%MCON(K1,N)   *   PSOL%UHOM_DNUP(IGEOM,O1,K2) + &
                           PBVP%MCON(K2,N)   *   PSOL%UHOM_DNUP(IGEOM,O1,K1)

                  LLUXR1 = PBVP%LCON(K1,N)   * L_PSOL%L_UHOM_DNDN(IGEOM,O1,K1,Q) - &
                           PBVP%LCON(K2,N)   * L_PSOL%L_UHOM_DNDN(IGEOM,O1,K2,Q)
                  LLUXR2 = PBVP%LCON(K1,N)   * L_PSOL%L_UHOM_DNDN(IGEOM,O1,K2,Q) + &
                           PBVP%LCON(K2,N)   * L_PSOL%L_UHOM_DNDN(IGEOM,O1,K1,Q)
                  MLUXR1 = PBVP%MCON(K1,N)   * L_PSOL%L_UHOM_DNUP(IGEOM,O1,K1,Q) - &
                           PBVP%MCON(K2,N)   * L_PSOL%L_UHOM_DNUP(IGEOM,O1,K2,Q)
                  MLUXR2 = PBVP%MCON(K1,N)   * L_PSOL%L_UHOM_DNUP(IGEOM,O1,K2,Q) + &
                           PBVP%MCON(K2,N)   * L_PSOL%L_UHOM_DNUP(IGEOM,O1,K1,Q)
 
                  NUXR1 =  L_PBVP%NCON(K1,N,Q) *   PSOL%UHOM_DNDN(IGEOM,O1,K1) - &
                           L_PBVP%NCON(K2,N,Q) *   PSOL%UHOM_DNDN(IGEOM,O1,K2)
                  NUXR2 =  L_PBVP%NCON(K1,N,Q) *   PSOL%UHOM_DNDN(IGEOM,O1,K2) + &
                           L_PBVP%NCON(K2,N,Q) *   PSOL%UHOM_DNDN(IGEOM,O1,K1)
                  PUXR1 =  L_PBVP%PCON(K1,N,Q) *   PSOL%UHOM_DNUP(IGEOM,O1,K1) - &
                           L_PBVP%PCON(K2,N,Q) *   PSOL%UHOM_DNUP(IGEOM,O1,K2)
                  PUXR2 =  L_PBVP%PCON(K1,N,Q) *   PSOL%UHOM_DNUP(IGEOM,O1,K2) + &
                           L_PBVP%PCON(K2,N,Q) *   PSOL%UHOM_DNUP(IGEOM,O1,K1)

!                  H1 =    ( NUXR1 + LLUXR1 ) * T_UTUP_EIGEN(K1)     - ( NUXR2 + LLUXR2 ) * T_UTUP_EIGEN(K2) &
!                                  +  LUXR1   * L_T_UTUP_EIGEN(K1,Q) -           LUXR2    * L_T_UTUP_EIGEN(K2,Q)
!                  H2 =    ( PUXR1 + MLUXR1 ) * T_UTDN_EIGEN(K1)     - ( PUXR2 + MLUXR2 ) * T_UTDN_EIGEN(K2) &
!                                  +  MUXR1   * L_T_UTDN_EIGEN(K1,Q) -           MUXR2    * L_T_UTDN_EIGEN(K2,Q)
 
                  H1 =    ( NUXR1 + LLUXR1 ) * T_UTDN_EIGEN(K1)     - ( NUXR2 + LLUXR2 ) * T_UTDN_EIGEN(K2) &
                                  +  LUXR1   * L_T_UTDN_EIGEN(K1,Q) -           LUXR2    * L_T_UTDN_EIGEN(K2,Q)
                  H2 =    ( PUXR1 + MLUXR1 ) * T_UTUP_EIGEN(K1)     - ( PUXR2 + MLUXR2 ) * T_UTUP_EIGEN(K2) &
                                  +  MUXR1   * L_T_UTUP_EIGEN(K1,Q) -           MUXR2    * L_T_UTUP_EIGEN(K2,Q)
 
                  SHOM_CR = SHOM_CR + H1 + H2

               ENDDO

!  Final contribution and end loops

               LC_MSST_F(Q,IGEOM,O1) = SHOM_R + SHOM_CR

            ENDDO 
         ENDDO

!  Add particular integral contribution
!  ====================================

!  Rob debug, 5/31/18
!    MK2, 6/21/21.  Use local bookkeeping variable LL_NSTK

         DO O1 = 1, NSTOKES
            DO Q = 1, N_COLUMN_WFS
               SPAR =  L_PSOL%LC_UPAR_DN_2(O1,Q) *    T_UTDN_MUBAR + &
                         PSOL%UPAR_DN_2(O1)      * LC_T_UTDN_MUBAR(Q)
               LC_MSST_F(Q,IGEOM,O1) = LC_MSST_F(Q,IGEOM,O1) + SPAR
            ENDDO
         ENDDO

!  End downwelling

      ENDIF

!  finish

      RETURN
      END SUBROUTINE MSSEG_LC_POSTPROCESS_LAYERMSST

!

       SUBROUTINE MSSEG_LC_POSTPROCESS_SURFACEMSST ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, IGEOM, FOURIER,    & ! Input control
        NSTOKES, NSTREAMS, NLAYERS, N_COLUMN_WFS, MUELLER_INDEX,      & ! Input control
        SURFACE_FACTOR, QUAD_STRMWTS, LAMBERTIAN_ALBEDO, USER_BRDF_F, & ! Input surface stuff
        HSOL, PBVP, L_HSOL, L_PSOL, L_PBVP, LC_SURFACE_MSST_F )         ! Input TS, Output MSST


!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- Add NCHOICE (debug only)

      USE MSSeg_PARS_m, only : MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAXLAYERS, MAX_GEOMETRIES, &
                               MAXEVALUES, MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_ATMOSWFS, ZERO

      use MSSEG_HPTypes_Def_m
      use MSSEG_L_HPTypes_Def_m

      IMPLICIT NONE

!  INPUTS
!  ======

!  Flags

      LOGICAL, INTENT (IN) ::            DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::            DO_LAMBERTIAN_SURFACE

!  Numbers

      INTEGER, INTENT (IN) ::            FOURIER
      INTEGER, INTENT (IN) ::            IGEOM
      INTEGER, INTENT (IN) ::            NSTOKES
      INTEGER, INTENT (IN) ::            NLAYERS
      INTEGER, INTENT (IN) ::            NSTREAMS
      INTEGER, INTENT (IN) ::            N_COLUMN_WFS

!  bookkeeping

      INTEGER         , INTENT (IN) :: MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  Surface inputs

      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )

!  input solution type structures

      type(Homog_Solutions) , intent(in) :: HSOL

!  Input BVProblem solution type structure

      type(Partic_BVPSolve)  , intent(in) :: PBVP

!  input type structure (Linearized homogeneous solutions)
 
      type(L_Homog_solutions) , intent(in) :: L_HSOL

!  Input type structures (Linearized particular solutions and BVP)
 
      type(L_Partic_solutions) , intent(in) :: L_PSOL
      type(L_Partic_BVPSolve)  , intent(in) :: L_PBVP

!  OUTPUTS
!  =======

      DOUBLE PRECISION, INTENT (INOUT)  :: LC_SURFACE_MSST_F ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXSTOKES )

!  Local variables
!  ---------------

      DOUBLE PRECISION :: LXR,  LLXR,  MLXR,  NXR,  PXR,  SHOM_R, HOM1, HOM2, SPAR
      DOUBLE PRECISION :: LXR1, LLXR1, MLXR1, NXR1, PXR1, SHOM
      DOUBLE PRECISION :: LXR2, LLXR2, NXR2, L_DOWN ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: SHOM_CR, HOM1CR, HOM2CR, HOM3CR, REFLEC, S_REFLEC, KMULT

      INTEGER          :: I, O1, O2, OM, K, K0, K1, K2, KO1, N, J, Q

!  Surface term
!  ============

!   Lowest layer

      N = NLAYERS
      KO1 = HSOL%K_REAL(N) + 1

!  Linearized Downward Stokes vector at surface at computational angles
!    reflectance integrand  a(j).x(j).I(-j). Two cases:

!    MK2, 6/21/21.  Use local bookkeeping variable LL_NSTK

      IF ( DO_INCLUDE_SURFACE ) THEN
         DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
               DO Q = 1, N_COLUMN_WFS

!  Real homogeneous solutions

                  SHOM_R = ZERO
                  DO K = 1, HSOL%K_REAL(N)
                     LXR  = PBVP%LCON(K,N)   *   HSOL%SOLA_XPOS(I,O1,K,N)
                     LLXR = PBVP%LCON(K,N)   * L_HSOL%L_SOLA_XPOS(I,O1,K,N,Q)
                     MLXR = PBVP%MCON(K,N)   * L_HSOL%L_SOLB_XNEG(I,O1,K,N,Q)
                     NXR  = L_PBVP%NCON(K,N,Q) *   HSOL%SOLA_XPOS(I,O1,K,N)
                     PXR  = L_PBVP%PCON(K,N,Q) *   HSOL%SOLB_XNEG(I,O1,K,N)
                     HOM1 = ( NXR + LLXR ) *   HSOL%T_DELT_EIGEN(K,N) &
                                  +  LXR   * L_HSOL%L_T_DELT_EIGEN(K,N,Q)
                     HOM2 = PXR + MLXR
                     SHOM_R = SHOM_R + HOM1 + HOM2
                  ENDDO

!  Complex homogeneous solutions

                  SHOM_CR = ZERO
                  DO K = 1, HSOL%K_COMPLEX(N)
                     K0 = 2 * K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1

                     NXR1  =   L_PBVP%NCON(K1,N,Q) *   HSOL%SOLA_XPOS(I,O1,K1,N) &
                             - L_PBVP%NCON(K2,N,Q) *   HSOL%SOLA_XPOS(I,O1,K2,N)
                     NXR2  =   L_PBVP%NCON(K1,N,Q) *   HSOL%SOLA_XPOS(I,O1,K2,N) &
                             + L_PBVP%NCON(K2,N,Q) *   HSOL%SOLA_XPOS(I,O1,K1,N)
                     PXR1  =   L_PBVP%PCON(K1,N,Q) *   HSOL%SOLB_XNEG(I,O1,K1,N) &
                             - L_PBVP%PCON(K2,N,Q) *   HSOL%SOLB_XNEG(I,O1,K2,N)
  
                     LXR1  =   PBVP%LCON(K1,N) *   HSOL%SOLA_XPOS(I,O1,K1,N) &
                             - PBVP%LCON(K2,N) *   HSOL%SOLA_XPOS(I,O1,K2,N)
                     LXR2  =   PBVP%LCON(K1,N) *   HSOL%SOLA_XPOS(I,O1,K2,N) &
                             + PBVP%LCON(K2,N) *   HSOL%SOLA_XPOS(I,O1,K1,N)
  
                     LLXR1  =   PBVP%LCON(K1,N) * L_HSOL%L_SOLA_XPOS(I,O1,K1,N,Q) &
                              - PBVP%LCON(K2,N) * L_HSOL%L_SOLA_XPOS(I,O1,K2,N,Q)
                     LLXR2  =   PBVP%LCON(K1,N) * L_HSOL%L_SOLA_XPOS(I,O1,K2,N,Q) &
                              + PBVP%LCON(K2,N) * L_HSOL%L_SOLA_XPOS(I,O1,K1,N,Q)
                     MLXR1  =   PBVP%MCON(K1,N) * L_HSOL%L_SOLB_XNEG(I,O1,K1,N,Q) &
                              - PBVP%MCON(K2,N) * L_HSOL%L_SOLB_XNEG(I,O1,K2,N,Q)
  
                     HOM1CR =   ( NXR1 + LLXR1 ) *   HSOL%T_DELT_EIGEN(K1,N) &
                              - ( NXR2 + LLXR2 ) *   HSOL%T_DELT_EIGEN(K2,N)
                     HOM2CR =             LXR1   * L_HSOL%L_T_DELT_EIGEN(K1,N,Q) &
                                     - LXR2   * L_HSOL%L_T_DELT_EIGEN(K2,N,Q)
                     HOM3CR = PXR1 + MLXR1
                     SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR
 
                  ENDDO

!  real part and add particular solution

                  SHOM = SHOM_R + SHOM_CR

!  Add Particular integral linearization

                  SPAR = L_PSOL%L_WLOWER(I,O1,N,Q)

!  Final result
!mick fix 2/12/2019 - applied "QUAD_STRMWTS" to SPAR + SHOM

                  L_DOWN(I,O1,Q) = QUAD_STRMWTS(I) * ( SPAR + SHOM )

!  End loops over Q, O1 and I

               ENDDO
            ENDDO
         ENDDO

!  End clause

      ENDIF

!  set the output
!   Validated, 02 May 2012 --> Same as BOA_DIFFUSE_SOURCE(1,1) Lambertian only.
!mick fix 7/13/2015 - added Q loops in IF block below
!mick fix 2/12/2019 - added or adjusted initializing of REFLEC in IF block below
!    MK2, 6/21/21.  Use local bookkeeping variable LL_NSTK (BRDF cases)

      IF ( DO_INCLUDE_SURFACE ) THEN
         IF ( DO_LAMBERTIAN_SURFACE .and. FOURIER .EQ. 0 ) THEN
            KMULT = SURFACE_FACTOR * LAMBERTIAN_ALBEDO
            O1 = 1
            DO Q = 1, N_COLUMN_WFS
               REFLEC = ZERO
               DO J = 1, NSTREAMS
                   REFLEC = REFLEC + L_DOWN(J,O1,Q)
               ENDDO
               LC_SURFACE_MSST_F(Q,IGEOM,O1) = KMULT * REFLEC
            ENDDO
         ELSE IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN
            KMULT = SURFACE_FACTOR
            DO O1 = 1, NSTOKES
               DO Q = 1, N_COLUMN_WFS
                  REFLEC = ZERO
                  DO J = 1, NSTREAMS
                     S_REFLEC = ZERO
                     DO O2 = 1, NSTOKES
                        OM = MUELLER_INDEX(O1,O2)
                        S_REFLEC = S_REFLEC + L_DOWN(J,O2,Q) * USER_BRDF_F(OM,IGEOM,J)
                     ENDDO
                     REFLEC = REFLEC + S_REFLEC
                  ENDDO
                  LC_SURFACE_MSST_F(Q,IGEOM,O1) = KMULT * REFLEC
               ENDDO
            ENDDO
         ENDIF
      ENDIF

!  finish

RETURN
END SUBROUTINE MSSEG_LC_POSTPROCESS_SURFACEMSST

!

SUBROUTINE GET_LC_QuadJacs_Levels & 
           ( DO_UPWELLING, NCHOICE, NSTOKES, NSTREAMS, NLAYERS, N_COLUMN_WFS, & ! Input numbers and Flag
             HSOL, PBVP, L_HSOL, L_PSOL, L_PBVP, QATMOSWF_F )                   ! OUTPUT

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- new subroutine specially written to get the discrete-ordinate Jacobian field

!  Mk2a. 9/21/21. New routine
!  Mk3.  1/31/22. Use type structures and rearrange argument lists

      USE MSSeg_PARS_m, Only : MAXSTREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS, &
                               MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO

      use MSSEG_HPTypes_Def_m
      use MSSEG_L_HPTypes_Def_m

      IMPLICIT NONE

!  Flag and numbers

      LOGICAL, INTENT (IN) ::          DO_UPWELLING      
      INTEGER, INTENT (IN) ::          NCHOICE
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_COLUMN_WFS

!  input solution type structures

      type(Homog_Solutions) , intent(in) :: HSOL

!  Input BVProblem solution type structure

      type(Partic_BVPSolve)  , intent(in) :: PBVP

!  input type structure (Linearized homogeneous solutions)
 
      type(L_Homog_solutions) , intent(in) :: L_HSOL

!  Input type structures (Linearized particular solutions and BVP)
 
      type(L_Partic_solutions) , intent(in) :: L_PSOL
      type(L_Partic_BVPSolve)  , intent(in) :: L_PBVP

!  output
!  ======

      DOUBLE PRECISION, INTENT (INOUT) :: QATMOSWF_F ( MAX_ATMOSWFS, MAXSTREAMS, MAXSTOKES, 2 )

!  local variables
!  ---------------

      INTEGER ::          L, L1, N, NL, I, I1, O1, K, KO1, K0, K1, K2, Q
      DOUBLE PRECISION :: SPAR, SHOM_R, SHOM_CR, SHOM
      DOUBLE PRECISION :: HOM1, HOM2, HOM1CR, HOM2CR, HOM3CR
      DOUBLE PRECISION :: NXR, PXR, NXR1, NXR2, PXR1, PXR2
      DOUBLE PRECISION :: LXR, MXR, LXR1, MXR1, LXR2, MXR2
      DOUBLE PRECISION :: LLXR, MLXR, LLXR1, MLXR1, LLXR2, MLXR2

!  Upwelling field
!  ---------------

      IF ( DO_UPWELLING ) THEN
        DO L = -1, 0
          NL = NCHOICE + L ; L1 = L + 2 ; N = NL + 1

!  Lowest level

          IF ( NL .EQ. NLAYERS ) THEN

!  Stokes and streams Parameter loops

            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              DO O1 = 1, NSTOKES
                DO Q = 1, N_COLUMN_WFS

!  Real homogeneous solutions

                   SHOM_R = ZERO
                   DO K = 1, HSOL%K_REAL(NL)
                      LXR  = PBVP%LCON(K,NL)   *   HSOL%SOLA_XPOS(I1,O1,K,NL)
                      LLXR = PBVP%LCON(K,NL)   * L_HSOL%L_SOLA_XPOS(I1,O1,K,NL,Q)
                      MLXR = PBVP%MCON(K,NL)   * L_HSOL%L_SOLB_XNEG(I1,O1,K,NL,Q)
                      NXR  = L_PBVP%NCON(K,NL,Q) *   HSOL%SOLA_XPOS(I1,O1,K,NL)
                      PXR  = L_PBVP%PCON(K,NL,Q) *   HSOL%SOLB_XNEG(I1,O1,K,NL)
                      HOM1 = ( LLXR + NXR ) * HSOL%T_DELT_EIGEN(K,NL) + LXR  * L_HSOL%L_T_DELT_EIGEN(K,NL,Q)
                      HOM2 = PXR + MLXR
                      SHOM_R = SHOM_R + HOM1 + HOM2
                   ENDDO

!  Complex homogeneous solutions

                   SHOM_CR = ZERO
                   KO1 = HSOL%K_REAL(NL) + 1
                   DO K = 1, HSOL%K_COMPLEX(NL)
                      K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                      NXR1  =   L_PBVP%NCON(K1,NL,Q) *   HSOL%SOLA_XPOS(I1,O1,K1,NL) &
                              - L_PBVP%NCON(K2,NL,Q) *   HSOL%SOLA_XPOS(I1,O1,K2,NL)
                      NXR2  =   L_PBVP%NCON(K1,NL,Q) *   HSOL%SOLA_XPOS(I1,O1,K2,NL) &
                              + L_PBVP%NCON(K2,NL,Q) *   HSOL%SOLA_XPOS(I1,O1,K1,NL)
                      PXR1  =   L_PBVP%PCON(K1,NL,Q) *   HSOL%SOLB_XNEG(I1,O1,K1,NL) &
                              - L_PBVP%PCON(K2,NL,Q) *   HSOL%SOLB_XNEG(I1,O1,K2,NL)
                      LXR1  =   PBVP%LCON(K1,NL) *   HSOL%SOLA_XPOS(I1,O1,K1,NL) &
                              - PBVP%LCON(K2,NL) *   HSOL%SOLA_XPOS(I1,O1,K2,NL)
                      LXR2  =   PBVP%LCON(K1,NL) *   HSOL%SOLA_XPOS(I1,O1,K2,NL) &
                              + PBVP%LCON(K2,NL) *   HSOL%SOLA_XPOS(I1,O1,K1,NL)
                      LLXR1  =   PBVP%LCON(K1,NL) * L_HSOL%L_SOLA_XPOS(I1,O1,K1,NL,Q) &
                               - PBVP%LCON(K2,NL) * L_HSOL%L_SOLA_XPOS(I1,O1,K2,NL,Q)
                      LLXR2  =   PBVP%LCON(K1,NL) * L_HSOL%L_SOLA_XPOS(I1,O1,K2,NL,Q) &
                               + PBVP%LCON(K2,NL) * L_HSOL%L_SOLA_XPOS(I1,O1,K1,NL,Q)
                      MLXR1  =   PBVP%MCON(K1,NL) * L_HSOL%L_SOLB_XNEG(I1,O1,K1,NL,Q) &
                               - PBVP%MCON(K2,NL) * L_HSOL%L_SOLB_XNEG(I1,O1,K2,NL,Q)
                      HOM1CR =   ( NXR1 + LLXR1 ) *   HSOL%T_DELT_EIGEN(K1,NL) &
                               - ( NXR2 + LLXR2 ) *   HSOL%T_DELT_EIGEN(K2,NL)
                      HOM2CR =             LXR1   * L_HSOL%L_T_DELT_EIGEN(K1,NL,Q) &
                                         - LXR2   * L_HSOL%L_T_DELT_EIGEN(K2,NL,Q)
                      HOM3CR = PXR1 + MLXR1
                      SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR
                   ENDDO

!  real part, add particular solution, complete result

                   SHOM = SHOM_R + SHOM_CR
                   SPAR = L_PSOL%L_WLOWER(I1,O1,NL,Q)
                   QATMOSWF_F(Q,I,O1,L1) = SPAR + SHOM

!  Finish Q, I and O1 loops

                ENDDO
              ENDDO
            ENDDO

!  End lowest level clause

          ENDIF

!  For other levels in the atmosphere
!  ==================================

          IF ( NL .NE. NLAYERS ) THEN

!  Stokes and streams and parameter loops

            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              DO O1 = 1, NSTOKES
                DO Q = 1, N_COLUMN_WFS

!  Real homogeneous solutions

                   SHOM_R = ZERO
                   DO K = 1, HSOL%K_REAL(N)
                      LXR  = PBVP%LCON(K,N)   *   HSOL%SOLA_XPOS(I1,O1,K,N)
                      MXR  = PBVP%MCON(K,N)   *   HSOL%SOLB_XNEG(I1,O1,K,N)
                      LLXR = PBVP%LCON(K,N)   * L_HSOL%L_SOLA_XPOS(I1,O1,K,N,Q)
                      MLXR = PBVP%MCON(K,N)   * L_HSOL%L_SOLB_XNEG(I1,O1,K,N,Q)
                      NXR  = L_PBVP%NCON(K,N,Q) *   HSOL%SOLA_XPOS(I1,O1,K,N)
                      PXR  = L_PBVP%PCON(K,N,Q) *   HSOL%SOLB_XNEG(I1,O1,K,N)
                      HOM1 =   NXR + LLXR
                      HOM2 = ( PXR + MLXR ) *   HSOL%T_DELT_EIGEN(K,N) &
                                   + MXR    * L_HSOL%L_T_DELT_EIGEN(K,N,Q)
                      SHOM_R = SHOM_R + HOM1 + HOM2
                   ENDDO

!  Complex homogeneous solutions

                   SHOM_CR = ZERO
                   KO1 = HSOL%K_REAL(N) + 1
                   DO K = 1, HSOL%K_COMPLEX(N)
                      K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                      NXR1  =   L_PBVP%NCON(K1,N,Q) *   HSOL%SOLA_XPOS(I1,O1,K1,N) &
                              - L_PBVP%NCON(K2,N,Q) *   HSOL%SOLA_XPOS(I1,O1,K2,N)
                      PXR1  =   L_PBVP%PCON(K1,N,Q) *   HSOL%SOLB_XNEG(I1,O1,K1,N) &
                              - L_PBVP%PCON(K2,N,Q) *   HSOL%SOLB_XNEG(I1,O1,K2,N)
                      PXR2  =   L_PBVP%PCON(K1,N,Q) *   HSOL%SOLB_XNEG(I1,O1,K2,N) &
                              + L_PBVP%PCON(K2,N,Q) *   HSOL%SOLB_XNEG(I1,O1,K1,N)
                      MXR1  =   PBVP%MCON(K1,N) *   HSOL%SOLB_XNEG(I1,O1,K1,N) &
                              - PBVP%MCON(K2,N) *   HSOL%SOLB_XNEG(I1,O1,K2,N)
                      MXR2  =   PBVP%MCON(K1,N) *   HSOL%SOLB_XNEG(I1,O1,K2,N) &
                              + PBVP%MCON(K2,N) *   HSOL%SOLB_XNEG(I1,O1,K1,N)
                      LLXR1  =   PBVP%LCON(K1,N) * L_HSOL%L_SOLA_XPOS(I1,O1,K1,N,Q) &
                               - PBVP%LCON(K2,N) * L_HSOL%L_SOLA_XPOS(I1,O1,K2,N,Q)
                      MLXR1  =   PBVP%MCON(K1,N) * L_HSOL%L_SOLB_XNEG(I1,O1,K1,N,Q) &
                               - PBVP%MCON(K2,N) * L_HSOL%L_SOLB_XNEG(I1,O1,K2,N,Q)
                      MLXR2  =   PBVP%MCON(K1,N) * L_HSOL%L_SOLB_XNEG(I1,O1,K2,N,Q) &
                               + PBVP%MCON(K2,N) * L_HSOL%L_SOLB_XNEG(I1,O1,K1,N,Q)
                      HOM1CR =   ( PXR1 + MLXR1 ) *   HSOL%T_DELT_EIGEN(K1,N) &
                               - ( PXR2 + MLXR2 ) *   HSOL%T_DELT_EIGEN(K2,N)
                      HOM2CR =             MXR1   * L_HSOL%L_T_DELT_EIGEN(K1,N,Q) &
                                         - MXR2   * L_HSOL%L_T_DELT_EIGEN(K2,N,Q)
                      HOM3CR = NXR1 + LLXR1
                      SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR
                   ENDDO

!  real part, add particular solution, complete result

                   SHOM = SHOM_R + SHOM_CR
                   SPAR = L_PSOL%L_WUPPER(I1,O1,N,Q)
                   QATMOSWF_F(Q,I,O1,L1) = SPAR + SHOM

!  Finish Q, I and O1 loops

                ENDDO
              ENDDO
            ENDDO

          ENDIF

!  end upwelling calculation

        ENDDO
      ENDIF

!  Downwelling field
!  -----------------

      IF ( .not.DO_UPWELLING ) THEN
        DO L = -1, 0
          NL = NCHOICE + L ; L1 = L + 2 ; N = NL

!  Highest level

          IF ( NL .EQ. 0 ) THEN
            DO I = 1, NSTREAMS
              DO Q = 1, N_COLUMN_WFS
                DO O1 = 1, NSTOKES
                  QATMOSWF_F(Q,I,O1,L1) = ZERO
                ENDDO
              ENDDO
            ENDDO

!  other levels

          ELSE IF ( NL.NE.0 ) THEN

!  Stokes and streams and parameter loops

            DO I = 1, NSTREAMS
              DO O1 = 1, NSTOKES
                DO Q = 1, N_COLUMN_WFS

!  Real homogeneous solutions

                   SHOM_R = ZERO
                   DO K = 1, HSOL%K_REAL(N)
                      LXR  = PBVP%LCON(K,N)   *   HSOL%SOLA_XPOS(I,O1,K,N)
                      LLXR = PBVP%LCON(K,N)   * L_HSOL%L_SOLA_XPOS(I,O1,K,N,Q)
                      MLXR = PBVP%MCON(K,N)   * L_HSOL%L_SOLB_XNEG(I,O1,K,N,Q)
                      NXR  = L_PBVP%NCON(K,N,Q) *   HSOL%SOLA_XPOS(I,O1,K,N)
                      PXR  = L_PBVP%PCON(K,N,Q) *   HSOL%SOLB_XNEG(I,O1,K,N)
                      HOM1 = ( NXR + LLXR ) *   HSOL%T_DELT_EIGEN(K,N) &
                                   +  LXR   * L_HSOL%L_T_DELT_EIGEN(K,N,Q)
                      HOM2 = PXR + MLXR
                      SHOM_R = SHOM_R + HOM1 + HOM2
                   ENDDO

!  Complex homogeneous solutions

                   SHOM_CR = ZERO
                   KO1 = HSOL%K_REAL(N) + 1
                   DO K = 1, HSOL%K_COMPLEX(N)
                      K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                      NXR1  =   L_PBVP%NCON(K1,N,Q) *   HSOL%SOLA_XPOS(I,O1,K1,N) &
                              - L_PBVP%NCON(K2,N,Q) *   HSOL%SOLA_XPOS(I,O1,K2,N)
                      NXR2  =   L_PBVP%NCON(K1,N,Q) *   HSOL%SOLA_XPOS(I,O1,K2,N) &
                              + L_PBVP%NCON(K2,N,Q) *   HSOL%SOLA_XPOS(I,O1,K1,N)
                      PXR1  =   L_PBVP%PCON(K1,N,Q) *   HSOL%SOLB_XNEG(I,O1,K1,N) &
                              - L_PBVP%PCON(K2,N,Q) *   HSOL%SOLB_XNEG(I,O1,K2,N)
                      LXR1  =   PBVP%LCON(K1,N) *   HSOL%SOLA_XPOS(I,O1,K1,N) &
                              - PBVP%LCON(K2,N) *   HSOL%SOLA_XPOS(I,O1,K2,N)
                      LXR2  =   PBVP%LCON(K1,N) *   HSOL%SOLA_XPOS(I,O1,K2,N) &
                              + PBVP%LCON(K2,N) *   HSOL%SOLA_XPOS(I,O1,K1,N)
                      LLXR1  =   PBVP%LCON(K1,N) * L_HSOL%L_SOLA_XPOS(I,O1,K1,N,Q) &
                               - PBVP%LCON(K2,N) * L_HSOL%L_SOLA_XPOS(I,O1,K2,N,Q)
                      LLXR2  =   PBVP%LCON(K1,N) * L_HSOL%L_SOLA_XPOS(I,O1,K2,N,Q) &
                               + PBVP%LCON(K2,N) * L_HSOL%L_SOLA_XPOS(I,O1,K1,N,Q)
                      MLXR1  =   PBVP%MCON(K1,N) * L_HSOL%L_SOLB_XNEG(I,O1,K1,N,Q) &
                               - PBVP%MCON(K2,N) * L_HSOL%L_SOLB_XNEG(I,O1,K2,N,Q)
                      HOM1CR =   ( NXR1 + LLXR1 ) *   HSOL%T_DELT_EIGEN(K1,N) &
                               - ( NXR2 + LLXR2 ) *   HSOL%T_DELT_EIGEN(K2,N)
                      HOM2CR =             LXR1   * L_HSOL%L_T_DELT_EIGEN(K1,N,Q) &
                                         - LXR2   * L_HSOL%L_T_DELT_EIGEN(K2,N,Q)
                      HOM3CR = PXR1 + MLXR1
                      SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR
                   ENDDO

!  real part, add particular solution, complete result

                   SHOM = SHOM_R + SHOM_CR
                   SPAR = L_PSOL%L_WLOWER(I,O1,N,Q)
                   QATMOSWF_F(Q,I,O1,L1) = SPAR + SHOM

!  Finish Q, I and O1 loops

                ENDDO
              ENDDO
            ENDDO

          ENDIF

!  End Downwelling calculation

        ENDDO
      ENDIF

!  done

RETURN
END SUBROUTINE GET_LC_QuadJacs_Levels

!  End module

END MODULE MSSeg_LC_PostProcess_m

