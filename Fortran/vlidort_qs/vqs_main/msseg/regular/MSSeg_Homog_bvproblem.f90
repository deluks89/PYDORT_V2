
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

      MODULE MSSeg_Homog_BVProblem_m

!  05 April 2023. Revision

      PUBLIC  :: MSSEG_MATRIXSETUP_MASTER
      PRIVATE :: BVP_SURFACE_SETUP_HOM, BVP_MATRIX_INIT, BVP_MATRIX_SETUP, BVP_MATRIX_SVD

      CONTAINS

      SUBROUTINE MSSEG_MATRIXSETUP_MASTER ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, SURFACE_FACTOR,    &
        FOURIER_COMPONENT, QUAD_STRMWTS, LAMBERTIAN_ALBEDO, BRDF_F,   &
        NSTOKES, NSTREAMS, NLAYERS, NSTKS_NSTRMS, MUELLER_INDEX,      &
        NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS_2,     &
        HSOL, HBVP, STATUS, MESSAGE, TRACE )

!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension, Only use quantities you need

      USE MSSeg_PARS_m, Only : MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAXSTREAMS_2, MAXSTRMSTKS_2, &
                               MAXEVALUES, MAXSTOKES_SQ, MAXBANDTOTAL, MAXTOTAL, MAXMOMENTS,  &
                               VQS_SUCCESS, VQS_SERIOUS

      use MSSeg_HPTypes_Def_m

      IMPLICIT NONE

!  Inputs
!  ------

      LOGICAL, INTENT (IN) ::            DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::            DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) ::   SURFACE_FACTOR

      INTEGER, INTENT (IN) ::            FOURIER_COMPONENT
      DOUBLE PRECISION, INTENT (IN) ::   QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::   LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::   BRDF_F ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )

      INTEGER, INTENT (IN) ::            NSTOKES
      INTEGER, INTENT (IN) ::            NSTREAMS
      INTEGER, INTENT (IN) ::            NLAYERS
      INTEGER, INTENT (IN) ::            NSTKS_NSTRMS

      INTEGER, INTENT (IN) ::            NSTREAMS_2
      INTEGER, INTENT (IN) ::            NTOTAL
      INTEGER, INTENT (IN) ::            N_SUBDIAG
      INTEGER, INTENT (IN) ::            N_SUPDIAG
      INTEGER, INTENT (IN) ::            NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::            MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

!  input solution type structure
!  -----------------------------

      type(Homog_Solutions), intent(inout) :: HSOL

!  outputs
!  -------

!  output type structure, should be intent(inout)

      type(Homog_BVPSetup), intent(inout) :: HBVP

!  Exception handling

      INTEGER          , INTENT (OUT)   :: STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  local variables
!  ---------------

      INTEGER :: STATUS_SUB

      INTEGER :: NMIN ( MAXTOTAL )
      INTEGER :: NMAX ( MAXTOTAL )
      INTEGER :: KALL

!  Initialize exception handling

      STATUS  = 0
      MESSAGE = ' '
      TRACE   = ' '

!  Additional setups for the albedo layer

      CALL BVP_SURFACE_SETUP_HOM ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, FOURIER_COMPONENT, &
        SURFACE_FACTOR, NSTOKES, NSTREAMS, NLAYERS,                   &
        LAMBERTIAN_ALBEDO, BRDF_F, QUAD_STRMWTS, MUELLER_INDEX,       &
        HSOL%K_REAL, HSOL%K_COMPLEX, HSOL%SOLA_XPOS, HSOL%SOLB_XNEG,  &
        HBVP%R2_HOMP, HBVP%R2_HOMM, HBVP%AXBID_F )

!  initialize compression matrix (Do this for every Fourier component)

      CALL BVP_MATRIX_INIT ( &
        NLAYERS, NSTREAMS_2, NTOTAL, &
        N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        NMIN, NMAX, KALL, HBVP%BANDMAT2, HBVP%SMAT2 )

!  set up boundary values matrix in compressed form (the "A" as in AX=B)

      CALL BVP_MATRIX_SETUP ( &
        DO_INCLUDE_SURFACE, FOURIER_COMPONENT, &
        NSTOKES, NSTREAMS, NLAYERS, NSTREAMS_2, NTOTAL, NSTKS_NSTRMS, NSTKS_NSTRMS_2,   &
        HSOL%K_REAL, HSOL%K_COMPLEX, HSOL%SOLA_XPOS, HSOL%SOLB_XNEG, HSOL%T_DELT_EIGEN, &
        HBVP%R2_HOMP, HBVP%R2_HOMM, NMIN, NMAX, KALL, HBVP%BANDMAT2, HBVP%SMAT2 )

!  SVD decomposition of compressed boundary values matrix

      CALL BVP_MATRIX_SVD ( &
        NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG, &
        HBVP%BANDMAT2, HBVP%SMAT2, HBVP%IPIVOT, HBVP%SIPIVOT, &
        STATUS_SUB, MESSAGE, TRACE )

!  error tracing

      IF ( STATUS_SUB .NE. VQS_SUCCESS ) THEN
        STATUS = VQS_SERIOUS ; RETURN
      ENDIF

!  finish

      RETURN
      END SUBROUTINE MSSEG_MATRIXSETUP_MASTER

!

      SUBROUTINE BVP_SURFACE_SETUP_HOM ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, FOURIER_COMPONENT, &
        SURFACE_FACTOR, NSTOKES, NSTREAMS, NLAYERS,                   &
        LAMBERTIAN_ALBEDO, BRDF_F, QUAD_STRMWTS, MUELLER_INDEX,       &
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,                      &
        R2_HOMP, R2_HOMM, AXBID_F )

!  This is the Lambertian or BRDF surface routine
!     Reflected homogeneous solutions

!  merged July 26, 2010. Using new BRDF material
!   RT SOLUTIONS Inc,

!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension, Only use quantities you need

!  Mk3, 8/27/21. .
!   -- BRDF input array has lost its "MAXMOMENTS" dimension. Now locally defined for each Fourier.

      USE MSSeg_PARS_m, only : MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAXSTREAMS_2, MAXSTRMSTKS_2, &
                                 MAXMOMENTS, MAXBANDTOTAL, MAXTOTAL, MAXEVALUES, MAXSTOKES_SQ, zero

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE
      INTEGER, INTENT (IN) ::           FOURIER_COMPONENT
      DOUBLE PRECISION, INTENT (IN) ::  SURFACE_FACTOR

      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS

      DOUBLE PRECISION, INTENT (IN) ::  LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  BRDF_F ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )

      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

      INTEGER, INTENT (IN) ::           K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (OUT) :: R2_HOMP &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (OUT) :: R2_HOMM &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (OUT) :: AXBID_F &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )

!  local variables
!  ---------------

      INTEGER          :: I,J,O1,O2,OM,K,KO1,K0,K1,K2,M,NL
      DOUBLE PRECISION :: AXJ, KMULT, FACTOR_KERNEL
      DOUBLE PRECISION :: H_1,   H_2,   H_1_S,   H_2_S
      DOUBLE PRECISION :: H_1_CR,   H_2_CR,   H_1_CI,   H_2_CI
      DOUBLE PRECISION :: H_1_S_CR, H_2_S_CR, H_1_S_CI, H_2_S_CI

!  Initialization
!  ==============

!  Zero total reflected contributions

     R2_HOMP = zero
     R2_HOMM = zero
     AXBID_F = zero

!  Return with Zeroed values if albedo flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Fourier component

      M = FOURIER_COMPONENT

!  Last layer

      NL = NLAYERS

!  Offset

      KO1 = K_REAL(NL) + 1

!  Return if Fourier component not zero for Lambertian surface

      IF ( DO_LAMBERTIAN_SURFACE .and. FOURIER_COMPONENT .GT. 0 ) RETURN
      IF ( FOURIER_COMPONENT .GT. 0 ) RETURN

!  Go to total BRDF case for non-Lambertian surface

      IF ( .not.  DO_LAMBERTIAN_SURFACE ) goto 877

!  For Lambertian reflectance, all streams are the same
!  ----------------------------------------------------

!  Integrate Downward streams of particular solutions

      KMULT = SURFACE_FACTOR * LAMBERTIAN_ALBEDO

!  Homogeneous real solutions

      DO K = 1, K_REAL(NLAYERS)
        DO I = 1, NSTREAMS
          H_1 = ZERO
          H_2 = ZERO
          DO J = 1, NSTREAMS
            AXJ = QUAD_STRMWTS(J)
            H_1 = H_1 + AXJ*SOLA_XPOS(J,1,K,NLAYERS)
            H_2 = H_2 + AXJ*SOLB_XNEG(J,1,K,NLAYERS)
          ENDDO
          R2_HOMP(I,1,K) = H_1 * KMULT
          R2_HOMM(I,1,K) = H_2 * KMULT
        ENDDO
      ENDDO

!  Homogeneous complex solutions

      KO1 = K_REAL(NLAYERS) + 1
      DO K = 1, K_COMPLEX(NLAYERS)
        K0 = 2*K - 2
        K1 = KO1 + K0
        K2 = K1  + 1
        DO I = 1, NSTREAMS
          H_1_CR = ZERO
          H_2_CR = ZERO
          H_1_CI = ZERO
          H_2_CI = ZERO
          DO J = 1, NSTREAMS
            AXJ = QUAD_STRMWTS(J)
            H_1_CR = H_1_CR + AXJ*SOLA_XPOS(J,1,K1,NLAYERS)
            H_2_CR = H_2_CR + AXJ*SOLB_XNEG(J,1,K1,NLAYERS)
            H_1_CI = H_1_CI + AXJ*SOLA_XPOS(J,1,K2,NLAYERS)
            H_2_CI = H_2_CI + AXJ*SOLB_XNEG(J,1,K2,NLAYERS)
          ENDDO
          R2_HOMP(I,1,K1) = H_1_CR * KMULT
          R2_HOMM(I,1,K1) = H_2_CR * KMULT
          R2_HOMP(I,1,K2) = H_1_CI * KMULT
          R2_HOMM(I,1,K2) = H_2_CI * KMULT
        ENDDO
      ENDDO

!  debug Lambertian

!      DO I = 1, NSTREAMS
!        DO O1 = 1, NSTOKES
!          DO K = 1, NSTKS_NSTRMS
!           write(*,'(3i3,1p2e14.6)')
!     &             I,O1,K,R2_HOMP(I,O1,K),R2_HOMM(I,O1,K)
!          ENDDO
!        ENDDO
!      ENDDO
!      pause'r2homp'

!  return from Lambertian

      RETURN

!  Total BRDF case
!  ===============

!  Continuation point

877   continue

!  Integrate Downward streams of particular solutions

      FACTOR_KERNEL = SURFACE_FACTOR

!  help variable
!   8/27/21. Mk3. BRDF_F array has no Fourier dimension. Define locally now.

      DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES
          DO J = 1, NSTREAMS
            DO O2 = 1, NSTOKES
              OM = MUELLER_INDEX(O1,O2)
!mick fix - 2/15/2011 - reactivated 1st statement; commented out 2nd statement
              AXBID_F(I,J,OM) = QUAD_STRMWTS(J) * BRDF_F(OM,I,J)
!              AXBID_F(I,J,OM) = QUAD_STRMWTS(J) * BRDF_F(M,OM,I,J)
!              AXBID_F(I,J,OM) = QUAD_STRMWTS(J) * BRDF_F(M,OM,J,I)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  homogeneous real solutions

      DO K = 1, K_REAL(NL)
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            H_1 = ZERO
            H_2 = ZERO
            DO J = 1, NSTREAMS
              H_1_S = ZERO
              H_2_S = ZERO
              DO O2 = 1, NSTOKES
                OM = MUELLER_INDEX(O1,O2)
                H_1_S = H_1_S + AXBID_F(I,J,OM) * SOLA_XPOS(J,O2,K,NL)
                H_2_S = H_2_S + AXBID_F(I,J,OM) * SOLB_XNEG(J,O2,K,NL)
              ENDDO
              H_1 = H_1 + H_1_S
              H_2 = H_2 + H_2_S
            ENDDO
            R2_HOMP(I,O1,K) = FACTOR_KERNEL * H_1
            R2_HOMM(I,O1,K) = FACTOR_KERNEL * H_2
          ENDDO
        ENDDO
      ENDDO

!  homogeneous complex solutions

      DO K = 1, K_COMPLEX(NL)
        K0 = 2*K - 2
        K1 = KO1 + K0
        K2 = K1  + 1
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            H_1_CR = ZERO
            H_2_CR = ZERO
            H_1_CI = ZERO
            H_2_CI = ZERO
            DO J = 1, NSTREAMS
              H_1_S_CR = ZERO
              H_2_S_CR = ZERO
              H_1_S_CI = ZERO
              H_2_S_CI = ZERO
              DO O2 = 1, NSTOKES
                OM = MUELLER_INDEX(O1,O2)
                H_1_S_CR = H_1_S_CR + &
                    AXBID_F(I,J,OM) * SOLA_XPOS(J,O2,K1,NL)
                H_2_S_CR = H_2_S_CR + &
                    AXBID_F(I,J,OM) * SOLB_XNEG(J,O2,K1,NL)
                H_1_S_CI = H_1_S_CI + &
                    AXBID_F(I,J,OM) * SOLA_XPOS(J,O2,K2,NL)
                H_2_S_CI = H_2_S_CI + &
                    AXBID_F(I,J,OM) * SOLB_XNEG(J,O2,K2,NL)
              ENDDO
              H_1_CR = H_1_CR + H_1_S_CR
              H_2_CR = H_2_CR + H_2_S_CR
              H_1_CI = H_1_CI + H_1_S_CI
              H_2_CI = H_2_CI + H_2_S_CI
            ENDDO
            R2_HOMP(I,O1,K1) = FACTOR_KERNEL * H_1_CR
            R2_HOMM(I,O1,K1) = FACTOR_KERNEL * H_2_CR
            R2_HOMP(I,O1,K2) = FACTOR_KERNEL * H_1_CI
            R2_HOMM(I,O1,K2) = FACTOR_KERNEL * H_2_CI
          ENDDO
        ENDDO
      ENDDO

!  debug total  BRDF case

!      DO I = 1, NSTREAMS
!        DO O1 = 1, NSTOKES
!          DO K = 1, NSTKS_NSTRMS
!           write(*,'(3i3,1p2e14.6)')
!     &             2,1,K,R2_HOMP(2,1,K),R2_HOMM(2,1,K)
!          ENDDO
!        ENDDO
!      ENDDO
!      pause'r2homp'

!  Finish

      RETURN
      END SUBROUTINE BVP_SURFACE_SETUP_HOM

!

      SUBROUTINE BVP_MATRIX_INIT ( &
        NLAYERS, NSTREAMS_2, NTOTAL, &
        N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        NMIN, NMAX, KALL, &
        BANDMAT2, SMAT2 )

!  Initialise the compressed matrix

!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension, Only use quantities you need

      USE MSSeg_PARS_m, Only: MAXBANDTOTAL, MAXTOTAL, MAXSTRMSTKS_2, ZERO

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2

      INTEGER, INTENT (OUT) ::          NMIN ( MAXTOTAL )
      INTEGER, INTENT (OUT) ::          NMAX ( MAXTOTAL )
      INTEGER, INTENT (OUT) ::          KALL
      DOUBLE PRECISION, INTENT (OUT) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      DOUBLE PRECISION, INTENT (OUT) :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )

!  Local variables
!  ---------------

!  compression function

!      INTEGER  BMAT_ROWMASK
!      EXTERNAL BMAT_ROWMASK

!  help variables

      INTEGER :: I, J, N3, JS, JF, IS, LF, L, I1, NMJ, NXJ

!  special case

      IF ( NLAYERS .EQ. 1 ) THEN
        DO I = 1, NTOTAL
          DO J = 1, NTOTAL
            SMAT2(I,J)     = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  compression row indices, NMIN, NMAX and KALL
!    Results are now stored in commons

      DO J = 1, N_SUPDIAG + 1
        NMIN(J) = 1
      ENDDO
      DO J = N_SUPDIAG + 2, NTOTAL
        NMIN(J) = J - N_SUPDIAG
      ENDDO
      DO J = 1, NTOTAL - N_SUBDIAG
        NMAX(J) = J + N_SUBDIAG
      ENDDO
      DO J = NTOTAL - N_SUBDIAG + 1, NTOTAL
        NMAX(J) = NTOTAL
      ENDDO
      KALL = N_SUBDIAG + N_SUPDIAG + 1

!  Former code to assign Array - now replaced by external function
!    Jukaa Kujanpaa, FMI, August 2005
!      DO I = 1, NTOTAL
!        DO J = 1, NTOTAL
!          IF ( (I.GE.NMIN(J)) .AND. (I.LE.NMAX(J)) ) THEN
!            BMAT_ROWMASK(I,J) = KALL + I - J
!          ENDIF
!        ENDDO
!      ENDDO

!  compression matrix zeroing

      N3 = NSTKS_NSTRMS_2 + NSTKS_NSTRMS
      LF = NLAYERS - 2

!  upper band top

      JS = NSTKS_NSTRMS_2 + 1
      JF = N3 - 1
      DO I = 1, NSTKS_NSTRMS
        DO J = JS, JF + I
          NMJ = NMIN(J); NXJ = NMAX(J)
          BANDMAT2(BMAT_ROWMASK(I,J,NMJ,NXJ,KALL),J) = ZERO
        ENDDO
      ENDDO

!  upper band

      DO L = 1, LF
        IS = L*NSTKS_NSTRMS_2 - NSTKS_NSTRMS + 1
        JS = IS + N3
        JF = JS - 1
        DO I = 1, NSTKS_NSTRMS_2-1
          I1 = I + IS
          DO J = JS, JF + I
            NMJ = NMIN(J); NXJ = NMAX(J)
            BANDMAT2(BMAT_ROWMASK(I1,J,NMJ,NXJ,KALL),J) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  lower band

      DO L = 1, LF
        IS = L*NSTKS_NSTRMS_2 + NSTKS_NSTRMS
        JS = IS - N3 + 1
        JF = IS - NSTKS_NSTRMS
        DO I = 1, NSTREAMS_2-1
          I1 = I + IS
          DO J = JS + I, JF
            NMJ = NMIN(J); NXJ = NMAX(J)
            BANDMAT2(BMAT_ROWMASK(I1,J,NMJ,NXJ,KALL),J) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  lower band bottom

      JS = LF * NSTKS_NSTRMS_2 + 1
      IS = JS + N3 - 1
      JF = IS - NSTKS_NSTRMS
      DO I = 1, NSTKS_NSTRMS
        I1 = I + IS
        DO J = JS + I, JF
          NMJ = NMIN(J); NXJ = NMAX(J)
          BANDMAT2(BMAT_ROWMASK(I1,J,NMJ,NXJ,KALL),J) = ZERO
        ENDDO
      ENDDO

!  finish

      RETURN
      END SUBROUTINE BVP_MATRIX_INIT

!

      SUBROUTINE BVP_MATRIX_SETUP ( &
        DO_INCLUDE_SURFACE, FOURIER_COMPONENT, &
        NSTOKES, NSTREAMS, &
        NLAYERS, &
        NSTREAMS_2, NTOTAL, &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        T_DELT_EIGEN, &
        R2_HOMP, R2_HOMM, &
        NMIN, NMAX, KALL, &
        BANDMAT2, SMAT2 )

!  Fills up the compressed matrix directly

!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension, Only use quantities you need

      USE MSSeg_PARS_m, only : MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAXSTREAMS_2, &
                                 MAXSTRMSTKS_2, MAXEVALUES, MAXBANDTOTAL, MAXTOTAL, ZERO

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      INTEGER, INTENT (IN) ::           FOURIER_COMPONENT

      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NSTREAMS_2
      INTEGER, INTENT (IN) ::           NTOTAL
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::           K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  R2_HOMP &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (IN) ::  R2_HOMM &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      INTEGER, INTENT (IN) ::           NMIN ( MAXTOTAL )
      INTEGER, INTENT (IN) ::           NMAX ( MAXTOTAL )
      INTEGER, INTENT (IN) ::           KALL

      DOUBLE PRECISION, INTENT (INOUT) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      DOUBLE PRECISION, INTENT (INOUT) :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )

!  compression function
!  --------------------

!      INTEGER  BMAT_ROWMASK
!      EXTERNAL BMAT_ROWMASK

!  Local variables
!  ---------------

      INTEGER ::          I, I1, N, N1, C0, CE_OFFSET, O1, M
      INTEGER ::          EP, EM, EPC, EPS, EMS
      INTEGER ::          CP, CM, CEP, CEM, CEPS, CEMS, CEM1, CEP1
      INTEGER ::          ROWEP, ROWEM, ROWEPS, ROWEMS, ROWCEP, ROWCEM
      INTEGER ::          ROWCEPS, ROWCEMS, ROWCEP1, ROWCEM1
      INTEGER ::          IR, IROW, KO11, KO1, K0, K1, K2
      DOUBLE PRECISION :: XPNET, XMNET
      DOUBLE PRECISION :: X1R, X2R, X1I, X2I, X1R_H, X1I_H

!  Initialization

      EP = 0

!  Fourier component

      M = FOURIER_COMPONENT

!  If Nlayers = 1, go to special case

      IF ( NLAYERS .EQ. 1 ) GO TO 345

!mick fix 2/17/11 - turn back on initialise (fix under test!)
!  Initialise - makes no difference
      DO EP = 1, NTOTAL
        DO EM = 1, MAXBANDTOTAL
          BANDMAT2(EM,EP) = ZERO
        ENDDO
      ENDDO

!  top BC for layer 1: no downward diffuse radiation
!  -------------------------------------------------

      N = 1
      KO1 = K_REAL(N) + 1
      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1

!  real solutions

          DO EP = 1, K_REAL(N)
            EM = EP + NSTKS_NSTRMS
            ROWEP = BMAT_ROWMASK(IROW,EP,NMIN(EP),NMAX(EP),KALL)
            ROWEM = BMAT_ROWMASK(IROW,EM,NMIN(EM),NMAX(EM),KALL)
            BANDMAT2(ROWEP,EP) = SOLA_XPOS(I,O1,EP,N)
            BANDMAT2(ROWEM,EM) = SOLB_XNEG(I,O1,EP,N)*T_DELT_EIGEN(EP,N)
          ENDDO

!  complex solutions (REWORKED)

          DO EPC = 1, K_COMPLEX(N)
            K0 = 2*EPC - 2
            K1 = KO1 + K0
            K2 = K1  + 1

            EP  = K_REAL(N) + EPC
            EPS = K_REAL(N) + K_COMPLEX(N) + EPC
            EM  = EP  + NSTKS_NSTRMS
            EMS = EPS + NSTKS_NSTRMS

            X1R = SOLA_XPOS(I,O1,K1,N)
            X1I = SOLA_XPOS(I,O1,K2,N)
            X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N) &
                 -SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
            X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N) &
                 +SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)

            ROWEP  = BMAT_ROWMASK(IROW,EP,NMIN(EP),NMAX(EP),KALL)
            ROWEM  = BMAT_ROWMASK(IROW,EM,NMIN(EM),NMAX(EM),KALL)
            ROWEPS = BMAT_ROWMASK(IROW,EPS,NMIN(EPS),NMAX(EPS),KALL)
            ROWEMS = BMAT_ROWMASK(IROW,EMS,NMIN(EMS),NMAX(EMS),KALL)

            BANDMAT2(ROWEP,EP)   =   X1R
            BANDMAT2(ROWEPS,EPS) = - X1I
            BANDMAT2(ROWEM,EM)   =   X2R
            BANDMAT2(ROWEMS,EMS) = - X2I
          ENDDO

        ENDDO
      ENDDO

!  intermediate layer boundaries
!  -----------------------------

      C0 = - NSTKS_NSTRMS
      DO N = 2, NLAYERS
        N1 = N - 1
        KO1  = K_REAL(N) + 1
        KO11 = K_REAL(N1) + 1
        C0   = C0 + NSTKS_NSTRMS_2
        CE_OFFSET = C0 - NSTKS_NSTRMS
        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW

!  real solutions for layer above the boundary

            DO EP = 1, K_REAL(N1)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS

              ROWCEP = BMAT_ROWMASK(CM,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM = BMAT_ROWMASK(CM,CEM,NMIN(CEM),NMAX(CEM),KALL)

              BANDMAT2(ROWCEP,CEP)   = &
                       T_DELT_EIGEN(EP,N1)*SOLA_XPOS(I,O1,EP,N1)
              BANDMAT2(ROWCEM,CEM)   = &
                         SOLB_XNEG(I,O1,EP,N1)
            ENDDO

!  real solutions for layer below the boundary
!   ( Note the change of sign !!! )

            DO EP = 1, K_REAL(N)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS
              CEP1 = CEP + NSTKS_NSTRMS_2
              CEM1 = CEM + NSTKS_NSTRMS_2

              ROWCEP1 = BMAT_ROWMASK(CM,CEP1,NMIN(CEP1),NMAX(CEP1),KALL)
              ROWCEM1 = BMAT_ROWMASK(CM,CEM1,NMIN(CEM1),NMAX(CEM1),KALL)

              BANDMAT2(ROWCEP1,CEP1) = &
                         - SOLA_XPOS(I,O1,EP,N)
              BANDMAT2(ROWCEM1,CEM1) = &
                         - T_DELT_EIGEN(EP,N)*SOLB_XNEG(I,O1,EP,N)
            ENDDO

!  complex solutions for layer above boundary

            DO EPC = 1, K_COMPLEX(N1)
              K0 = 2*EPC - 2
              K1 = KO11 + K0
              K2 = K1  + 1

              EP  = K_REAL(N1) + EPC
              EPS = K_REAL(N1) + K_COMPLEX(N1) + EPC
              CEP  = CE_OFFSET + EP
              CEM  = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS
              CEMS = CEPS + NSTKS_NSTRMS

              X1R = SOLA_XPOS(I,O1,K1,N1)*T_DELT_EIGEN(K1,N1) &
                   -SOLA_XPOS(I,O1,K2,N1)*T_DELT_EIGEN(K2,N1)
              X1I = SOLA_XPOS(I,O1,K1,N1)*T_DELT_EIGEN(K2,N1) &
                   +SOLA_XPOS(I,O1,K2,N1)*T_DELT_EIGEN(K1,N1)
              X2R = SOLB_XNEG(I,O1,K1,N1)
              X2I = SOLB_XNEG(I,O1,K2,N1)

              ROWCEP  = BMAT_ROWMASK(CM,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM  = BMAT_ROWMASK(CM,CEM,NMIN(CEM),NMAX(CEM),KALL)
              ROWCEPS = BMAT_ROWMASK(CM,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
              ROWCEMS = BMAT_ROWMASK(CM,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

              BANDMAT2(ROWCEP ,CEP)  =   X1R
              BANDMAT2(ROWCEPS,CEPS) = - X1I
              BANDMAT2(ROWCEM ,CEM)  =   X2R
              BANDMAT2(ROWCEMS,CEMS) = - X2I
            ENDDO

!  complex solutions for layer below boundary
!   ( Note the change of sign !!! )

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              CEP  = CE_OFFSET + EP + NSTKS_NSTRMS_2
              CEM  = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS + NSTKS_NSTRMS_2
              CEMS = CEPS + NSTKS_NSTRMS

              X1R = SOLA_XPOS(I,O1,K1,N)
              X1I = SOLA_XPOS(I,O1,K2,N)
              X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N) &
                   -SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
              X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N) &
                   +SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)

              ROWCEP  = BMAT_ROWMASK(CM,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM  = BMAT_ROWMASK(CM,CEM,NMIN(CEM),NMAX(CEM),KALL)
              ROWCEPS = BMAT_ROWMASK(CM,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
              ROWCEMS = BMAT_ROWMASK(CM,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

              BANDMAT2(ROWCEP,CEP)   = - X1R
              BANDMAT2(ROWCEPS,CEPS) =   X1I
              BANDMAT2(ROWCEM,CEM)   = - X2R
              BANDMAT2(ROWCEMS,CEMS) =   X2I
            ENDDO

          ENDDO
        ENDDO
      ENDDO

!  bottom BC (with albedo additions if flagged)

      N   = NLAYERS
      KO1 = K_REAL(N) + 1
      C0  = C0 + NSTKS_NSTRMS_2
      CE_OFFSET = C0 - NSTKS_NSTRMS

!  With albedo
!  -----------

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = C0 + IROW

!  real solutions

            DO EP = 1, K_REAL(N)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS

              XPNET = SOLA_XPOS(I1,O1,EP,N) - R2_HOMP(I,O1,EP)
              XMNET = SOLB_XNEG(I1,O1,EP,N) - R2_HOMM(I,O1,EP)

              ROWCEP  = BMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM  = BMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)

              BANDMAT2(ROWCEP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
              BANDMAT2(ROWCEM,CEM) = XMNET
            ENDDO

!  Complex solutions

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS
              CEMS = CEPS + NSTKS_NSTRMS

              X1R_H = SOLA_XPOS(I1,O1,K1,N) - R2_HOMP(I,O1,K1)
              X1I_H = SOLA_XPOS(I1,O1,K2,N) - R2_HOMP(I,O1,K2)
              X2R = SOLB_XNEG(I1,O1,K1,N) - R2_HOMM(I,O1,K1)
              X2I = SOLB_XNEG(I1,O1,K2,N) - R2_HOMM(I,O1,K2)

              X1R = X1R_H * T_DELT_EIGEN(K1,N) - &
                    X1I_H * T_DELT_EIGEN(K2,N)
              X1I = X1R_H * T_DELT_EIGEN(K2,N) + &
                    X1I_H * T_DELT_EIGEN(K1,N)

              ROWCEP  = BMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM  = BMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)
              ROWCEPS = BMAT_ROWMASK(CP,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
              ROWCEMS = BMAT_ROWMASK(CP,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

              BANDMAT2(ROWCEP,CEP)   =   X1R
              BANDMAT2(ROWCEPS,CEPS) = - X1I
              BANDMAT2(ROWCEM,CEM)   =   X2R
              BANDMAT2(ROWCEMS,CEMS) = - X2I

!         iunit = 20
!        IF (do_fdtest)iunit=19
!         write(iunit,4567)I,O1,EPC,
!     &         R2_HOMP(I,O1,EPC)*T_DELT_EIGEN(EPC,N)+
!     &         R2_HOMP(I,O1,EPC)*T_DELT_EIGEN(EPC,N),
!     &                   R2_HOMM(I,O1,EPC)
! 4567  format   (3i4,1p2e20.10)

            ENDDO

          ENDDO
        ENDDO

!  No albedo

      ELSE

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = C0 + IROW

!  real solutions

            DO EP = 1, K_REAL(N)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS

              XPNET = SOLA_XPOS(I1,O1,EP,N)
              XMNET = SOLB_XNEG(I1,O1,EP,N)

              ROWCEP  = BMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM  = BMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)

              BANDMAT2(ROWCEP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
              BANDMAT2(ROWCEM,CEM) = XMNET
            ENDDO

!  Complex solutions
!    ----- Deep bug found for X1I, 17 December 2005.

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              CEP  = CE_OFFSET + EP
              CEM  = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS
              CEMS = CEPS + NSTKS_NSTRMS

              X1R = SOLA_XPOS(I1,O1,K1,N) * T_DELT_EIGEN(K1,N) &
                   -SOLA_XPOS(I1,O1,K2,N) * T_DELT_EIGEN(K2,N)
              X1I = SOLA_XPOS(I1,O1,K1,N) * T_DELT_EIGEN(K2,N) &
                   +SOLA_XPOS(I1,O1,K2,N) * T_DELT_EIGEN(K1,N)
              X2R = SOLB_XNEG(I1,O1,K1,N)
              X2I = SOLB_XNEG(I1,O1,K2,N)

              ROWCEP  = BMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
              ROWCEM  = BMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)
              ROWCEPS = BMAT_ROWMASK(CP,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
              ROWCEMS = BMAT_ROWMASK(CP,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

              BANDMAT2(ROWCEP,CEP)   =   X1R
              BANDMAT2(ROWCEPS,CEPS) = - X1I
              BANDMAT2(ROWCEM,CEM)   =   X2R
              BANDMAT2(ROWCEMS,CEMS) = - X2I
           ENDDO

          ENDDO
        ENDDO

      ENDIF

!  normal completion

      RETURN

!  special case for 1 layer only
!  =============================

345   CONTINUE

!  top BC for layer 1: no downward diffuse radiation
!  -------------------------------------------------

      N = 1
      KO1 = K_REAL(N) + 1
      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1

!  real solutions

          DO EP = 1, K_REAL(N)
            EM = EP + NSTKS_NSTRMS
            SMAT2(IROW,EP) = SOLA_XPOS(I,O1,EP,N)
            SMAT2(IROW,EM) = SOLB_XNEG(I,O1,EP,N)*T_DELT_EIGEN(EP,N)
          ENDDO

!  REWORKED complex solutions

          DO EPC = 1, K_COMPLEX(N)
            K0 = 2*EPC - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            EP  = K_REAL(N) + EPC
            EPS = K_REAL(N) + K_COMPLEX(N) + EPC
            EM  = EP  + NSTKS_NSTRMS
            EMS = EPS + NSTKS_NSTRMS
            X1R = SOLA_XPOS(I,O1,K1,N)
            X1I = SOLA_XPOS(I,O1,K2,N)
            X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N) &
                 -SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
            X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N) &
                 +SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)
            SMAT2(IROW,EP)  =   X1R
            SMAT2(IROW,EPS) = - X1I
            SMAT2(IROW,EM)  =   X2R
            SMAT2(IROW,EMS) = - X2I
          ENDDO

        ENDDO
      ENDDO

!  bottom BC (with albedo additions)
!  ---------------------------------

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = NSTKS_NSTRMS + IROW

!  real solutions

            DO EP = 1, K_REAL(N)
              CEP = EP
              CEM = CEP + NSTKS_NSTRMS
              XPNET = SOLA_XPOS(I1,O1,EP,N) - R2_HOMP(I,O1,EP)
              XMNET = SOLB_XNEG(I1,O1,EP,N) - R2_HOMM(I,O1,EP)
              SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
              SMAT2(CP,CEM) = XMNET
            ENDDO

!  Complex solutions

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              EM  = EP  + NSTKS_NSTRMS
              CEP = EP
              CEM = EM
              CEPS = EPS
              CEMS = CEPS + NSTKS_NSTRMS
              X1R_H = SOLA_XPOS(I1,O1,K1,N) - R2_HOMP(I,O1,K1)
              X1I_H = SOLA_XPOS(I1,O1,K2,N) - R2_HOMP(I,O1,K2)
              X2R = SOLB_XNEG(I1,O1,K1,N) - R2_HOMM(I,O1,K1)
              X2I = SOLB_XNEG(I1,O1,K2,N) - R2_HOMM(I,O1,K2)
              X1R = X1R_H * T_DELT_EIGEN(K1,N) - &
                    X1I_H * T_DELT_EIGEN(K2,N)
              X1I = X1R_H * T_DELT_EIGEN(K2,N) + &
                    X1I_H * T_DELT_EIGEN(K1,N)
              SMAT2(CP,CEP)  =   X1R
              SMAT2(CP,CEPS) = - X1I
              SMAT2(CP,CEM)  =   X2R
              SMAT2(CP,CEMS) = - X2I
            ENDDO

          ENDDO
        ENDDO

!  bottom BC (with No albedo)
!  --------------------------

      ELSE

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = NSTKS_NSTRMS + IROW

!  real solutions

            DO EP = 1, K_REAL(N)
              CEP = EP
              CEM = CEP + NSTKS_NSTRMS
              XPNET = SOLA_XPOS(I1,O1,EP,N)
              XMNET = SOLB_XNEG(I1,O1,EP,N)
              SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
              SMAT2(CP,CEM) = XMNET
            ENDDO

!  REWORKED Complex solutions

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              EM  = EP  + NSTKS_NSTRMS
              CEP = EP
              CEM = EM
              CEPS = EPS
              CEMS = CEPS + NSTKS_NSTRMS
              X1R = SOLA_XPOS(I1,O1,K1,N)*T_DELT_EIGEN(K1,N) &
                   -SOLA_XPOS(I1,O1,K2,N)*T_DELT_EIGEN(K2,N)
              X1I = SOLA_XPOS(I1,O1,K1,N)*T_DELT_EIGEN(K2,N) &
                   +SOLA_XPOS(I1,O1,K2,N)*T_DELT_EIGEN(K1,N)
              X2R = SOLB_XNEG(I1,O1,K1,N)
              X2I = SOLB_XNEG(I1,O1,K2,N)
              SMAT2(CP,CEP)  =   X1R
              SMAT2(CP,CEPS) = - X1I
              SMAT2(CP,CEM)  =   X2R
              SMAT2(CP,CEMS) = - X2I
            ENDDO

          ENDDO
        ENDDO

      ENDIF

!  normal return and finish

      RETURN
      END SUBROUTINE BVP_MATRIX_SETUP

!

      SUBROUTINE BVP_MATRIX_SVD ( &
        NLAYERS, NTOTAL, &
        N_SUBDIAG, N_SUPDIAG, &
        BANDMAT2, SMAT2, &
        IPIVOT, SIPIVOT, &
        STATUS, MESSAGE, TRACE )

!  Solves the boundary value problem.

!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension, Only use quantities you need

      USE MSSeg_PARS_m, only : MAXBANDTOTAL, MAXTOTAL, MAXSTRMSTKS_2, VQS_SUCCESS, VQS_SERIOUS
      USE LAPACK_TOOLS_m, only : DGBTRF, DGETRF

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG

      DOUBLE PRECISION, INTENT (INOUT) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      DOUBLE PRECISION, INTENT (INOUT) :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )

      INTEGER, INTENT (OUT) ::            IPIVOT ( MAXTOTAL )
      INTEGER, INTENT (OUT) ::            SIPIVOT ( MAXSTRMSTKS_2 )
      INTEGER, INTENT (OUT) ::            STATUS
      CHARACTER (LEN=*), INTENT (INOUT) ::  MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) ::  TRACE

!  local variables
!  ---------------

      INTEGER ::           INFO
      CHARACTER (LEN=3) :: CI

!  Intialize Exception handling

      STATUS = VQS_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  SVD the BVP matrix: With compression (multilayers)
!  --------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK LU-decomposition for band matrix

        CALL DGBTRF &
           ( NTOTAL, NTOTAL, N_SUBDIAG, N_SUPDIAG, &
             BANDMAT2, MAXBANDTOTAL, IPIVOT, INFO )

!  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE   = 'DGBTRF call (nlayers>1) in BVP_MATRIX_SVD'
          STATUS  = VQS_SERIOUS
          RETURN
        ELSE IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRF call (nlayers>1) in BVP_MATRIX_SVD'
          STATUS  = VQS_SERIOUS
          RETURN
        ENDIF

!  SVD the BVP matrix: No compression, Single Layer only
!  -----------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK LU-decomposition for  matrix

        CALL DGETRF &
           ( NTOTAL, NTOTAL, SMAT2, MAXSTRMSTKS_2, SIPIVOT, INFO )

!  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRF call (nlayers=1) in BVP_MATRIX_SVD'
          STATUS  = VQS_SERIOUS
          RETURN
        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BVP_MATRIX_SVD


      INTEGER FUNCTION BMAT_ROWMASK ( &
        I, J, NMINJ, NMAXJ, KALL )

!  Integer function for Band-matrix compression
!   Replaces Array BMAT_ROWMASK(I,J) which was too large.
!   Replacement by J. Kujanpaa, FMI, August 2005.

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: I
      INTEGER, INTENT (IN) :: J
      INTEGER, INTENT (IN) :: NMINJ
      INTEGER, INTENT (IN) :: NMAXJ
      INTEGER, INTENT (IN) :: KALL

!  Assign the function (Keep the hard-wired Stop for now)
!    KALL and NMIN,NMAX are assigned in BVP_MATRIX_INIT

      IF ( (I.GE.NMINJ) .AND. (I.LE.NMAXJ) ) THEN
         BMAT_ROWMASK = KALL + I - J
      ELSE
         BMAT_ROWMASK = 0
         WRITE(*,*) 'ERROR=out of bounds, BMAT_ROWMASK:',I,J
         STOP
      ENDIF

!  Finish

      RETURN
      END FUNCTION BMAT_ROWMASK

      END MODULE MSSeg_Homog_BVProblem_m
