
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

!  Mk3. 1/31/21. Completely new module, with 2 subroutines taking over the function
!  05 April 2023. Revision

      MODULE MSSeg_Partic_PostProcess_m

      PRIVATE
      PUBLIC :: MSSEG_POSTPROCESS_LAYERMSST, MSSEG_POSTPROCESS_SURFACEMSST, GET_QuadIntens_Levels

      CONTAINS


      SUBROUTINE MSSEG_POSTPROCESS_LAYERMSST ( &
        DO_UPWELLING, DO_DNWELLING, IGEOM, FOURIER, NCHOICE, NSTOKES, DELTAU_VERT,       &
        SFI_OPDEPS_SCALED, BEAM_CUTOFF, AVERAGE_SECANT, INITIAL_TRANS, HSOL, PSOL, PBVP, & 
        TAU_UP, TAU_DN, T_UTUP_EIGEN, T_UTDN_EIGEN, T_UTDN_MUBAR, MSST_F )

!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- Add NCHOICE (debug only)

!  Mk3, 1/31/22. New subroutine made from existing material.

      USE MSSeg_PARS_m, only : MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAXLAYERS, MAX_GEOMETRIES, &
                               MAXEVALUES, MAX_TAU_SPATH, MAX_TAU_QPATH, ZERO

      use MSSeg_HPTypes_Def_m

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

!  Optical inputs

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SFI_OPDEPS_SCALED  ( MAX_GEOMETRIES )

!  Beam inputs

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF    ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )

!  input solution type structures

      type(Homog_Solutions) , intent(in) :: HSOL
      type(Partic_Solutions), intent(in) :: PSOL

!  Input BVProblem solution type structure

      type(Partic_BVPSolve)  , intent(in) :: PBVP

!  OUTPUTS
!  =======
!mick fix 6/6/2023 - changed intent of MSST_F from OUT to INOUT

      DOUBLE PRECISION, INTENT (OUT)  :: TAU_DN, TAU_UP
      DOUBLE PRECISION, INTENT (OUT)  :: T_UTUP_EIGEN ( MAXEVALUES )
      DOUBLE PRECISION, INTENT (OUT)  :: T_UTDN_EIGEN ( MAXEVALUES )
      DOUBLE PRECISION, INTENT (OUT)  :: T_UTDN_MUBAR

      DOUBLE PRECISION, INTENT (INOUT) :: MSST_F ( MAX_GEOMETRIES, MAXSTOKES )

!  Local variables
!  ---------------

      DOUBLE PRECISION :: TRANSH, HELP, HELP_CR, HELP_CI, SPHER
      DOUBLE PRECISION :: LUXR, MUXR, H_R, LUX_CR, LUX_CI, MUX_CR, MUX_CI, H_CR, SFOR2

      INTEGER          :: O1, K, K0, K1, K2, KO1, N

!  Optical depths. Must used scaled quantities.
!    ALREADY SCALED

      TAU_DN = SFI_OPDEPS_SCALED(IGEOM)
      TAU_UP = DELTAU_VERT(NCHOICE) - SFI_OPDEPS_SCALED(IGEOM)

!  Transmittances
!  --------------

      N = NCHOICE
      KO1 = HSOL%K_REAL(N) + 1

!  Initialize them

      T_UTUP_EIGEN = zero ; T_UTDN_EIGEN = zero ; T_UTDN_MUBAR = zero

!  Real eigenvalues

      DO K = 1, HSOL%K_REAL(N)
         HELP = HSOL%KEIGEN(K,N) * TAU_DN
         IF ( HELP .LE. MAX_TAU_QPATH ) T_UTDN_EIGEN(K) = DEXP(-HELP)
         HELP = HSOL%KEIGEN(K,N) * TAU_UP
         IF ( HELP .LE. MAX_TAU_QPATH ) T_UTUP_EIGEN(K) = DEXP(-HELP)
      ENDDO

!  Complex eigenvalues

      DO K = 1, HSOL%K_COMPLEX(N)
         K0 = 2*K-2 ; K1 = KO1 + K0 ;  K2 = K1  + 1
         HELP = TAU_DN * HSOL%KEIGEN(K1,N)
         IF ( HELP .LE. MAX_TAU_QPATH ) THEN
            HELP_CR = DEXP ( - HELP )
            HELP_CI = TAU_DN * HSOL%KEIGEN(K2,N)
            T_UTDN_EIGEN(K1) =   HELP_CR * DCOS ( HELP_CI )
            T_UTDN_EIGEN(K2) = - HELP_CR * DSIN ( HELP_CI )
         ENDIF
         HELP = TAU_UP * HSOL%KEIGEN(K1,N)
         IF ( HELP .LE. MAX_TAU_QPATH ) THEN
            HELP_CR = DEXP ( - HELP )
            HELP_CI = TAU_UP * HSOL%KEIGEN(K2,N)
            T_UTUP_EIGEN(K1) =   HELP_CR * DCOS ( HELP_CI )
            T_UTUP_EIGEN(K2) = - HELP_CR * DSIN ( HELP_CI )
         ENDIF
      ENDDO

!  Particular integral term

      IF ( N .LE. BEAM_CUTOFF(IGEOM) ) THEN
         TRANSH = INITIAL_TRANS(NCHOICE, IGEOM)
         SPHER = TAU_DN * AVERAGE_SECANT(NCHOICE,IGEOM)
         IF ( SPHER .LE. MAX_TAU_SPATH ) T_UTDN_MUBAR = TRANSH * DEXP ( - SPHER )
      ENDIF

!  Upwelling Radiance Contribution
!  -------------------------------

      IF ( DO_UPWELLING ) THEN

!  Homogeneous term
!    MK2, 6/21/21.  Use local bookkeeping variable LC_NSTK

         DO O1 = 1, NSTOKES
           H_R = ZERO
           DO K = 1, HSOL%K_REAL(N)
             LUXR = PBVP%LCON(K,N) * PSOL%UHOM_UPDN(IGEOM,O1,K)
             MUXR = PBVP%MCON(K,N) * PSOL%UHOM_UPUP(IGEOM,O1,K)
             H_R = H_R + LUXR * T_UTDN_EIGEN(K) +  MUXR * T_UTUP_EIGEN(K)
           ENDDO
           H_CR = ZERO
           DO K = 1, HSOL%K_COMPLEX(N)
             K0 = 2*K-2 ; K1 = KO1 + K0 ;  K2 = K1  + 1
             LUX_CR = PBVP%LCON(K1,N) * PSOL%UHOM_UPDN(IGEOM,O1,K1) - PBVP%LCON(K2,N) * PSOL%UHOM_UPDN(IGEOM,O1,K2)
             LUX_CI = PBVP%LCON(K1,N) * PSOL%UHOM_UPDN(IGEOM,O1,K2) + PBVP%LCON(K2,N) * PSOL%UHOM_UPDN(IGEOM,O1,K1)
             MUX_CR = PBVP%MCON(K1,N) * PSOL%UHOM_UPUP(IGEOM,O1,K1) - PBVP%MCON(K2,N) * PSOL%UHOM_UPUP(IGEOM,O1,K2)
             MUX_CI = PBVP%MCON(K1,N) * PSOL%UHOM_UPUP(IGEOM,O1,K2) + PBVP%MCON(K2,N) * PSOL%UHOM_UPUP(IGEOM,O1,K1)
             H_CR = H_CR + LUX_CR * T_UTDN_EIGEN(K1) - LUX_CI * T_UTDN_EIGEN(K2) &
                         + MUX_CR * T_UTUP_EIGEN(K1) - MUX_CI * T_UTUP_EIGEN(K2)
           ENDDO
           MSST_F(IGEOM,O1) = H_R + H_CR
         ENDDO

!  Add particular integral contribution
!    MK2, 6/21/21.  Use local bookkeeping variable LC_NSTK

         DO O1 = 1, NSTOKES
            SFOR2 = T_UTDN_MUBAR * PSOL%UPAR_UP_2(O1)
            MSST_F(IGEOM,O1) = MSST_F(IGEOM,O1) + SFOR2
         ENDDO

      ENDIF

!  Downwelling Radiance Contribution
!  ---------------------------------

      IF ( DO_DNWELLING ) THEN

!  homogeneous Term

         DO O1 = 1, NSTOKES
           H_R = ZERO
           DO K = 1, HSOL%K_REAL(N)
             LUXR = PBVP%LCON(K,N) * PSOL%UHOM_DNDN(IGEOM,O1,K)
             MUXR = PBVP%MCON(K,N) * PSOL%UHOM_DNUP(IGEOM,O1,K)
             H_R = H_R + LUXR * T_UTDN_EIGEN(K) + MUXR * T_UTUP_EIGEN(K)
           ENDDO
           H_CR = ZERO
           KO1 = HSOL%K_REAL(N) + 1
           DO K = 1, HSOL%K_COMPLEX(N)
             K0 = 2*K-2 ; K1 = KO1 + K0 ;  K2 = K1  + 1
             LUX_CR = PBVP%LCON(K1,N) * PSOL%UHOM_DNDN(IGEOM,O1,K1) - PBVP%LCON(K2,N) * PSOL%UHOM_DNDN(IGEOM,O1,K2)
             LUX_CI = PBVP%LCON(K1,N) * PSOL%UHOM_DNDN(IGEOM,O1,K2) + PBVP%LCON(K2,N) * PSOL%UHOM_DNDN(IGEOM,O1,K1)
             MUX_CR = PBVP%MCON(K1,N) * PSOL%UHOM_DNUP(IGEOM,O1,K1) - PBVP%MCON(K2,N) * PSOL%UHOM_DNUP(IGEOM,O1,K2)
             MUX_CI = PBVP%MCON(K1,N) * PSOL%UHOM_DNUP(IGEOM,O1,K2) + PBVP%MCON(K2,N) * PSOL%UHOM_DNUP(IGEOM,O1,K1)
             H_CR = H_CR + LUX_CR * T_UTDN_EIGEN(K1) - LUX_CI * T_UTDN_EIGEN(K2)  &
                         + MUX_CR * T_UTUP_EIGEN(K1) - MUX_CI * T_UTUP_EIGEN(K2)
           ENDDO
           MSST_F(IGEOM,O1) = H_R + H_CR
         ENDDO

!  Add particular integral contribution
!    MK2, 6/21/21.  Use local bookkeeping variable LC_NSTK

         DO O1 = 1, NSTOKES
            SFOR2 = T_UTDN_MUBAR * PSOL%UPAR_DN_2(O1)
            MSST_F(IGEOM,O1) = MSST_F(IGEOM,O1) + SFOR2
         ENDDO

!  Debug check, 19 November 2015
!    if ( ngeoms.gt.5 .and.fourier.eq.0) write(*,*)'sfor2',IGEOM,MSST_F(IGEOM,1),SFOR2

      ENDIF

!  finish

      RETURN
      END SUBROUTINE MSSEG_POSTPROCESS_LAYERMSST

!

      SUBROUTINE MSSEG_POSTPROCESS_SURFACEMSST ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, IGEOM, FOURIER, NSTOKES, NSTREAMS, NLAYERS, &
        MUELLER_INDEX, SURFACE_FACTOR, QUAD_STRMWTS, LAMBERTIAN_ALBEDO, USER_BRDF_F, HSOL, PSOL, PBVP, & 
        STOKES_DOWNSURF, SURFACE_MSST_F )

!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- Add NCHOICE (debug only)

      USE MSSeg_PARS_m, only : MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAXLAYERS, MAX_GEOMETRIES, &
                               MAXEVALUES, MAXSTOKES_SQ, MAX_USER_STREAMS, ZERO

      use MSSeg_HPTypes_Def_m

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

!  bookkeeping

      INTEGER         , INTENT (IN) :: MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  Surface inputs

      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )

!  input solution type structures

      type(Homog_Solutions) , intent(in) :: HSOL
      type(Partic_Solutions), intent(in) :: PSOL

!  Input BVProblem solution type structure

      type(Partic_BVPSolve)  , intent(in) :: PBVP

!  OUTPUTS
!  =======

      DOUBLE PRECISION, INTENT (INOUT)  :: SURFACE_MSST_F ( MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT)  :: STOKES_DOWNSURF ( MAXSTREAMS, MAXSTOKES )

!  Local variables
!  ---------------

      DOUBLE PRECISION :: LXR, MXR, LXR_CR, LXR_CI, MXR_CR, STOTAL, SPAR
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, HOM1, HOM1CR, HOM2CR, REFLEC, S_REFLEC, KMULT

      INTEGER          :: I, O1, O2, OM, K, K0, K1, K2, KO1, N, J

!  Surface term
!  ============

!   Lowest layer

      N = NLAYERS
      KO1 = HSOL%K_REAL(N) + 1

!  Downward Stokes vector at surface at computational angles
!    reflectance integrand  a(j).x(j).I(-j)
!    MK2, 6/21/21.  Use local bookkeeping variable LC_NSTK

      IF ( DO_INCLUDE_SURFACE ) THEN
         DO I = 1, NSTREAMS
           DO O1 = 1, NSTOKES
             SPAR = PSOL%WLOWER(I,O1,N)
             SHOM_R = ZERO
             DO K = 1, HSOL%K_REAL(N)
               LXR  = PBVP%LCON(K,N)*HSOL%SOLA_XPOS(I,O1,K,N)
               MXR  = PBVP%MCON(K,N)*HSOL%SOLB_XNEG(I,O1,K,N)
               HOM1 = LXR * HSOL%T_DELT_EIGEN(K,N) + MXR
               SHOM_R = SHOM_R + HOM1
             ENDDO
             SHOM_CR = ZERO
             DO K = 1, HSOL%K_COMPLEX(N)
               K0 = 2*K-2 ; K1 = KO1 + K0 ;  K2 = K1  + 1
               LXR_CR  = PBVP%LCON(K1,N) * HSOL%SOLA_XPOS(I,O1,K1,N) - PBVP%LCON(K2,N) * HSOL%SOLA_XPOS(I,O1,K2,N)
               LXR_CI  = PBVP%LCON(K1,N) * HSOL%SOLA_XPOS(I,O1,K2,N) + PBVP%LCON(K2,N) * HSOL%SOLA_XPOS(I,O1,K1,N)
               MXR_CR  = PBVP%MCON(K1,N) * HSOL%SOLB_XNEG(I,O1,K1,N) - PBVP%MCON(K2,N) * HSOL%SOLB_XNEG(I,O1,K2,N)
               HOM1CR  = LXR_CR*HSOL%T_DELT_EIGEN(K1,N) - LXR_CI*HSOL%T_DELT_EIGEN(K2,N)
               HOM2CR  = MXR_CR
               SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
             ENDDO
             STOTAL = SPAR + SHOM_R + SHOM_CR
             STOKES_DOWNSURF(I,O1) =  QUAD_STRMWTS(I) * STOTAL
           ENDDO
         ENDDO
      ENDIF

!  set the output
!   Validated, 02 May 2012 --> Same as BOA_DIFFUSE_SOURCE(1,1). Lambertian only.
!    MK2, 6/21/21.  Use local bookkeeping variable LC_NSTK

      IF ( DO_INCLUDE_SURFACE ) THEN
        IF ( DO_LAMBERTIAN_SURFACE .and. FOURIER .EQ. 0 ) THEN
          KMULT = SURFACE_FACTOR * LAMBERTIAN_ALBEDO
          O1 = 1 ; REFLEC = ZERO
          DO J = 1, NSTREAMS
            REFLEC = REFLEC + STOKES_DOWNSURF(J,O1)
          ENDDO
          SURFACE_MSST_F(IGEOM,O1) = KMULT * REFLEC
        ELSE IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN
           KMULT = SURFACE_FACTOR
           DO O1 = 1, NSTOKES
              DO J = 1, NSTREAMS
                S_REFLEC = ZERO
                DO O2 = 1, NSTOKES
                  OM = MUELLER_INDEX(O1,O2)
                  S_REFLEC = S_REFLEC + STOKES_DOWNSURF(J,O2) * USER_BRDF_F(OM,IGEOM,J)
                ENDDO
                REFLEC = REFLEC + S_REFLEC
              ENDDO
              SURFACE_MSST_F(IGEOM,O1) = KMULT * REFLEC
           ENDDO
        ENDIF
      ENDIF

!  finish

      RETURN
      END SUBROUTINE MSSEG_POSTPROCESS_SURFACEMSST

!

      SUBROUTINE GET_QuadIntens_Levels & 
           ( DO_UPWELLING, NCHOICE, NSTOKES, NSTREAMS, & ! Flag, numbers
             NLAYERS, HSOL, PSOL, PBVP,                & ! Solutions
             QSTOKES_F )                                 ! OUTPUT

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- new subroutine specially written to get the discrete-ordinate field

      USE MSSeg_PARS_m, Only : MAXSTREAMS, MAXSTOKES, MAXLAYERS,  &
                               MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO

      use MSSeg_HPTypes_Def_m

      IMPLICIT NONE

!  Flag and numbers

      LOGICAL, INTENT (IN) ::          DO_UPWELLING      
      INTEGER, INTENT (IN) ::          NCHOICE
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS

!  input solution type structures

      type(Homog_Solutions) , intent(in) :: HSOL
      type(Partic_Solutions), intent(in) :: PSOL

!  Input BVProblem solution type structure

      type(Partic_BVPSolve) , intent(in) :: PBVP

!  output

      DOUBLE PRECISION, INTENT (INOUT) :: QSTOKES_F ( MAXSTREAMS, MAXSTOKES, 2 )

!  local variables
!  ---------------

      INTEGER ::          L, L1, N, NL, I, I1, O1, K, KO1, K0, K1, K2
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, HOM1CR, HOM2CR 
      DOUBLE PRECISION :: LXR, MXR, LXR_CR, LXR_CI, MXR_CR, MXR_CI

!  Upwelling field
!  ---------------

      IF ( DO_UPWELLING ) THEN
        DO L = -1, 0
          NL = NCHOICE + L ; L1 = L + 2 ; N = NL + 1

!  Lowest level

          IF ( NL .EQ. NLAYERS ) THEN
            KO1 = HSOL%K_REAL(NL) + 1
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              DO O1 = 1, NSTOKES
                SHOM_R = ZERO
                DO K = 1, HSOL%K_REAL(NL)
                  LXR = PBVP%LCON(K,NL)*HSOL%SOLA_XPOS(I1,O1,K,NL)
                  MXR = PBVP%MCON(K,NL)*HSOL%SOLB_XNEG(I1,O1,K,NL)
                  SHOM_R = SHOM_R + LXR * HSOL%T_DELT_EIGEN(K,NL) + MXR
                ENDDO
                SHOM_CR = ZERO
                DO K = 1, HSOL%K_COMPLEX(NL)
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                  LXR_CR =  PBVP%LCON(K1,NL) * HSOL%SOLA_XPOS(I1,O1,K1,NL) - PBVP%LCON(K2,NL) * HSOL%SOLA_XPOS(I1,O1,K2,NL)
                  LXR_CI =  PBVP%LCON(K1,NL) * HSOL%SOLA_XPOS(I1,O1,K2,NL) + PBVP%LCON(K2,NL) * HSOL%SOLA_XPOS(I1,O1,K1,NL)
                  MXR_CR =  PBVP%MCON(K1,NL) * HSOL%SOLB_XNEG(I1,O1,K1,NL) - PBVP%MCON(K2,NL) * HSOL%SOLB_XNEG(I1,O1,K2,NL)
                  HOM1CR =  LXR_CR*HSOL%T_DELT_EIGEN(K1,NL) - LXR_CI*HSOL%T_DELT_EIGEN(K2,NL)
                  SHOM_CR = SHOM_CR + HOM1CR + MXR_CR
                ENDDO
                QSTOKES_F(I,O1,L1) = SHOM_R + SHOM_CR + PSOL%WLOWER(I1,O1,NL)
              ENDDO
            ENDDO

!  other levels
!mick fix 6/6/2023 - changed index of HSOL%K_COMPLEX for K loop from NL to N

          ELSE
            KO1 = HSOL%K_REAL(N) + 1
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              DO O1 = 1, NSTOKES
                SHOM_R = ZERO
                DO K = 1, HSOL%K_REAL(N)
                  LXR = PBVP%LCON(K,N)*HSOL%SOLA_XPOS(I1,O1,K,N)
                  MXR = PBVP%MCON(K,N)*HSOL%SOLB_XNEG(I1,O1,K,N)
                  SHOM_R = SHOM_R + LXR +  MXR * HSOL%T_DELT_EIGEN(K,N)
                ENDDO
                SHOM_CR = ZERO
                DO K = 1, HSOL%K_COMPLEX(N)
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                  LXR_CR =  PBVP%LCON(K1,N) * HSOL%SOLA_XPOS(I1,O1,K1,N) - PBVP%LCON(K2,N) * HSOL%SOLA_XPOS(I1,O1,K2,N)
                  MXR_CR =  PBVP%MCON(K1,N) * HSOL%SOLB_XNEG(I1,O1,K1,N) - PBVP%MCON(K2,N) * HSOL%SOLB_XNEG(I1,O1,K2,N)
                  MXR_CI =  PBVP%MCON(K1,N) * HSOL%SOLB_XNEG(I1,O1,K2,N) + PBVP%MCON(K2,N) * HSOL%SOLB_XNEG(I1,O1,K1,N)
                  HOM2CR =  MXR_CR * HSOL%T_DELT_EIGEN(K1,N) - MXR_CI * HSOL%T_DELT_EIGEN(K2,N)
                  SHOM_CR = SHOM_CR + LXR_CR + HOM2CR
                ENDDO
                QSTOKES_F(I,O1,L1) = SHOM_R + SHOM_CR + PSOL%WUPPER(I1,O1,N)
              ENDDO
            ENDDO
          ENDIF

        ENDDO
      ENDIF

!  Downwelling field
!  -----------------

      IF ( .not.DO_UPWELLING ) THEN
        DO L = -1, 0
          NL = NCHOICE + L ; L1 = L + 2 ; N = NL

!  Highest level

          IF ( NL .EQ. 0 ) THEN
            KO1 = HSOL%K_REAL(NL) + 1
            DO I = 1, NSTREAMS
              DO O1 = 1, NSTOKES
                QSTOKES_F(I,O1,L1) = ZERO
              ENDDO
            ENDDO

!  other levels

          ELSE
            KO1 = HSOL%K_REAL(N) + 1
            DO I = 1, NSTREAMS
              DO O1 = 1, NSTOKES
                SHOM_R = ZERO
                DO K = 1, HSOL%K_REAL(N)
                  LXR = PBVP%LCON(K,N)*HSOL%SOLA_XPOS(I,O1,K,N)
                  MXR = PBVP%MCON(K,N)*HSOL%SOLB_XNEG(I,O1,K,N)
                  SHOM_R = SHOM_R + LXR * HSOL%T_DELT_EIGEN(K,N) + MXR 
                ENDDO
                SHOM_CR = ZERO
                DO K = 1, HSOL%K_COMPLEX(NL)
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                  LXR_CR =  PBVP%LCON(K1,N) * HSOL%SOLA_XPOS(I,O1,K1,N) - PBVP%LCON(K2,N) * HSOL%SOLA_XPOS(I,O1,K2,N)
                  MXR_CR =  PBVP%MCON(K1,N) * HSOL%SOLB_XNEG(I,O1,K1,N) - PBVP%MCON(K2,N) * HSOL%SOLB_XNEG(I,O1,K2,N)
                  MXR_CI =  PBVP%MCON(K1,N) * HSOL%SOLB_XNEG(I,O1,K2,N) + PBVP%MCON(K2,N) * HSOL%SOLB_XNEG(I,O1,K1,N)
                  HOM1CR =  LXR_CR * HSOL%T_DELT_EIGEN(K1,N) - LXR_CI * HSOL%T_DELT_EIGEN(K2,N)
                  SHOM_CR = SHOM_CR + MXR_CR + HOM1CR
                ENDDO
                QSTOKES_F(I,O1,L1) = SHOM_R + SHOM_CR + PSOL%WLOWER(I,O1,N)
              ENDDO
            ENDDO
          ENDIF

        ENDDO
      ENDIF

!  done

      RETURN
      END SUBROUTINE GET_QuadIntens_Levels

!

      END MODULE MSSeg_Partic_PostProcess_m

