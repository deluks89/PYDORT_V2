
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
! # Subroutines in this Module                                  #
! #                                                             #
! #  For linearizations involving atmospheric parameters        #
! #                                                             #
! #            MSSEG_LC_BVP_SOLUTION_MASTER  (calling)          #
! #                                                             #
! #                MSSEG_LC_BVP_COLUMN_SETUP                    #
! #                LC_BEAMSOLUTION_NEQK                         #
! #                                                             #
! #                L1_BVP_BACKSUB         (generic linear)      #
! #                L1_BVP_SURFACE_SETUP   (generic linear)      #
! #                                                             #
! ###############################################################

!  Mk2a. 9/21/21. LC routines created by R. Spurr

!  Mk2, 6/21/21. Several Changes.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- ASYMTX_TOLERANCE now introduced and passed into the QHOM call to ASYMTX
!   -- Flux output (LOCAL_FLUXES) controlled by DO_LOCAL_FLUXES input flag
!   -- DO_FOURIER0_NSTOKES2 performance flag introduced.
!   -- Call subroutine specially written to get the discrete-ordinate field for Fluxes
!   -- Use public routine GETQUAD2 instead of GAULEG

!  Mk3. 4/11/23. (1/31/22 created). Modified for MK3 use

      MODULE MSSeg_LC_BVProblem_m

      use MSSEG_HPTypes_Def_m
      use MSSEG_L_HPTypes_Def_m

      PRIVATE
      PUBLIC ::  MSSEG_LC_BVPSOLUTION_MASTER

      CONTAINS

      SUBROUTINE MSSEG_LC_BVPSOLUTION_MASTER ( &
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, DO_LAMBERTIAN_SURFACE,       &
        IBEAM, FOURIER, N_COLUMN_WFS, NCHOICE, NSTOKES, NSTREAMS, NLAYERS,      &
        NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        SURFACE_FACTOR, LAMBERTIAN_ALBEDO, BRDF_F, QUAD_STRMWTS, DIRECT_BEAM,   &
        DO_LAYER_SCATTERING, MUELLER_INDEX, DELTAU_SLANT, L_DELTAU_VERT,        &
        BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, LC_INITIAL_TRANS, LC_T_DELT_MUBAR, &
        HSOL, HBVP, PSOL, PBVP, L_HSOL, L_PSOL, L_PBVP, STATUS, MESSAGE, TRACE )

!  Mk2.  6/21/21. Use module names now have the "_m" extension, Only use quantities you need
!  Mk2a. 9/21/21. Created by  R. Spurr
!  Mk3.  1/31/22. Argument list re-organized for clarity, output status changed
!  Mk3.  1/31/22. Use type structure inputs and outputs directly

      USE MSSeg_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAX_ATMOSWFS, &
                               MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, MAXSTRMSTKS_2, MAXSTOKES_SQ,   &
                               MAXBANDTOTAL, MAXTOTAL, VQS_SUCCESS, VQS_SERIOUS

      IMPLICIT NONE

!  Input Arguments
!  ---------------

!  Flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE

!  Control integers

      INTEGER, INTENT (IN) ::          FOURIER, IBEAM




      INTEGER, INTENT (IN) ::          N_COLUMN_WFS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NCHOICE

!  Bookkeeping

      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG

!  Bookkeeping and optical

      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )
      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT  ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  surface stuff. Mk3 4/11/23. Drop Fourier index from BRDF_F

      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )

!  Average secant stuff

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF      ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR     ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS    ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Type structure arguments
!  ------------------------

!  input solution type structures

      type(Homog_Solutions) , intent(in) :: HSOL
      type(Partic_Solutions), intent(in) :: PSOL

!  input BVProblem type structures

      type(Homog_BVPSetup)  , intent(in) :: HBVP
      type(Partic_BVPSolve) , intent(in) :: PBVP

!  linearized solutions (mostly input)

      type(L_Homog_Solutions) , intent(in)    :: L_HSOL
      type(L_Partic_Solutions), intent(inout) :: L_PSOL

!  output linearized BVP stuff

      type(L_Partic_BVPSolve) , intent(inout) :: L_PBVP

!  exception handling

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  error tracing variables

      INTEGER :: STATUS_SUB

!  module status and message initialization

      STATUS  = VQS_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Linearization of the regular BVP case
!  =====================================

      CALL MSSEG_LC_BVPCOLUMN_SETUP ( &
          DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,     & ! Input flags
          FOURIER, IBEAM, NCHOICE, NSTOKES, NSTREAMS, NLAYERS, N_COLUMN_WFS,    & ! Input control integers
          NSTREAMS_2, NTOTAL, NSTKS_NSTRMS, NSTKS_NSTRMS_2, MUELLER_INDEX,      & ! input bookkeeping numbers
          SURFACE_FACTOR, LAMBERTIAN_ALBEDO, BRDF_F, QUAD_STRMWTS, DIRECT_BEAM, & ! Input surface quantities
          DO_LAYER_SCATTERING, DELTAU_SLANT, L_DELTAU_VERT, BEAM_CUTOFF,        & ! Input bookkeeping/Optical
          T_DELT_MUBAR, INITIAL_TRANS, LC_T_DELT_MUBAR, LC_INITIAL_TRANS,       & ! Input average-scant stuff
          HSOL%K_REAL, HSOL%K_COMPLEX, HSOL%SOLA_XPOS, HSOL%SOLB_XNEG,          & ! TS Input  Homog solutions
          HSOL%T_DELT_EIGEN, PSOL%BVEC, PBVP%LCON, PBVP%MCON, L_PSOL%LC_BVEC,   & ! TS Input  PI Sols/BVP
          L_HSOL%L_T_DELT_EIGEN, L_HSOL%L_SOLA_XPOS, L_HSOL%L_SOLB_XNEG,        & ! TS Input  L_Homog solutions
          L_PSOL%L_WLOWER, L_PSOL%L_WUPPER, L_PBVP%COL2_WF, L_PBVP%SCOL2_WF )     ! TS Output linearized BVP

!  Back-substitution

      CALL L1_BVP_BACKSUB ( &
          N_COLUMN_WFS, NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, & ! Input numbers
          HSOL%K_REAL, HSOL%K_COMPLEX, HBVP%BANDMAT2, HBVP%IPIVOT, HBVP%SMAT2, HBVP%SIPIVOT, & ! TS Input Homog Sol/BVP
          L_PBVP%COL2_WF, L_PBVP%SCOL2_WF, L_PBVP%NCON, L_PBVP%PCON,                         & ! TS Input/Output linearized BVP
          STATUS_SUB, MESSAGE, TRACE )

!  error tracing

      IF ( STATUS_SUB .NE. VQS_SUCCESS ) THEN
        STATUS = VQS_SERIOUS
        RETURN
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE MSSEG_LC_BVPSOLUTION_MASTER

!

      SUBROUTINE MSSEG_LC_BVPCOLUMN_SETUP ( &
          DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,     & ! Input  flags
          FOURIER, IBEAM, NCHOICE, NSTOKES, NSTREAMS, NLAYERS, N_COLUMN_WFS,    & ! Input  control integers
          NSTREAMS_2, NTOTAL, NSTKS_NSTRMS, NSTKS_NSTRMS_2, MUELLER_INDEX,      & ! input  bookkeeping numbers
          SURFACE_FACTOR, LAMBERTIAN_ALBEDO, BRDF_F, QUAD_STRMWTS, DIRECT_BEAM, & ! Input  surface quantities
          DO_LAYER_SCATTERING, DELTAU_SLANT, L_DELTAU_VERT, BEAM_CUTOFF,        & ! Input  bookkeeping/Optical
          T_DELT_MUBAR, INITIAL_TRANS, LC_T_DELT_MUBAR, LC_INITIAL_TRANS,       & ! Input  average-scant stuff
          K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, BVEC,          & ! Input  Hom/PI Solutions
          LCON, MCON, LC_BVEC, L_T_DELT_EIGEN, L_SOLA_XPOS, L_SOLB_XNEG,        & ! Input  BVP, Lin Hom/PI Solutions
          L_WLOWER, L_WUPPER, COL2_WF, SCOL2_WF )                                 ! Output linearized BVP

!  Mk2.  6/21/21. Use module names now have the "_m" extension, Only use quantities you need
!  Mk3.  1/31/22. Argument list re-organized for clarity, output status changed

      USE MSSeg_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAX_ATMOSWFS, &
                               MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, MAXSTRMSTKS_2, MAXSTOKES_SQ,   &
                               MAXBANDTOTAL, MAXTOTAL, ZERO

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE

!  Control integers

      INTEGER, INTENT (IN) ::          FOURIER, IBEAM
      INTEGER, INTENT (IN) ::          NCHOICE     ! debug only
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_COLUMN_WFS

!  Bookkeeping numbers

      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

!  Surface inputs. Mk3 4/11/23. Drop Fourier index from BRDF_F

      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  optical and bookkeeping

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: DIRECT_BEAM  ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  Average secant and linearization

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR    ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS    ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Type structure inputs
!  ---------------------

!  Homog solutions

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  BVP constants

      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Linearized Homog solutions

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  particular integral and linearization

      DOUBLE PRECISION, INTENT (IN) :: BVEC    ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LC_BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  output BVP vector and Linearized PI solutions at the boundaries
!     Quantities must be Intent(in/out) now, as they are in type structures.

      DOUBLE PRECISION, INTENT (INOUT) :: COL2_WF  ( MAXTOTAL, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: SCOL2_WF ( MAXSTRMSTKS_2, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: L_WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  local variables
!  ---------------

!  Help variables

      INTEGER ::          Q, N, N1, I, I1, IR, IROW, CM, C0, O1
      INTEGER ::          K, KO1, K0, K1, K2
      DOUBLE PRECISION :: CPOS, CNEG, L_HOM_R, L_HOM_CR, L_BEAM, FAC
      DOUBLE PRECISION :: T1, T2, T1R, T1I, T2R, T2I
      DOUBLE PRECISION :: L_HOM_U_R  , L_HOM_D_R
      DOUBLE PRECISION :: L_HOM_U_CR , L_HOM_D_CR
      LOGICAL ::          MODIFIED_BOUNDARY

!  Linearized surface reflectance solutions

      DOUBLE PRECISION :: R2_L_BEAM ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: R2_L_HOMP ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: R2_L_HOMM ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

!  initialise
!  ----------

!  zero the results vectors

      DO I = 1, NTOTAL
        DO Q = 1, MAX_ATMOSWFS
          COL2_WF(I,Q) = ZERO
        ENDDO
      ENDDO

!    This is a very important zeroing.................!!!!!

      DO I = 1, NSTREAMS_2
        DO O1 = 1, NSTOKES
          DO Q = 1, N_COLUMN_WFS
            DO N = 1, NLAYERS
              L_WUPPER(I,O1,N,Q) = ZERO
              L_WLOWER(I,O1,N,Q) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  BCL1 or BCL3M - top of first layer (TOA), UPPER boundary condition
!  ------------------------------------------------------------------

      N = 1

!  Get the linearized beam solution for the varying layer

      CALL LC_BEAMSOLUTION_NEQK ( &
             FOURIER, IBEAM, N, N_COLUMN_WFS, NSTOKES, NSTREAMS_2,          & ! Input numbers
             DO_LAYER_SCATTERING, BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, & ! Input av. secant
             BVEC, LC_T_DELT_MUBAR, LC_INITIAL_TRANS, LC_BVEC,              & ! Input linearized stuff
             L_WUPPER, L_WLOWER )                                             ! Output at level boundaries

!  Start  loops

      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          DO Q = 1, N_COLUMN_WFS

!  beam solution linearization at top of layer

            L_BEAM = - L_WUPPER(I,O1,N,Q)

!  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
              CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + &
                   L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) - &
                   L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
              T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                   - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                   - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
              T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                   + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N) &
                   + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

!  Final contribution

            COL2_WF(IROW,Q) = L_BEAM - L_HOM_R - L_HOM_CR

!  end loops

          ENDDO
        ENDDO
      ENDDO

!  Intermediate boundary conditions
!  --------------------------------

      DO N = 1, NLAYERS - 1

!  N1 is the layer below, C0 is the offset

        N1 = N + 1
        C0 = N*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

!  Get the linearized beam solution for the next layer
!    -- 1/31/21. Version 2.8.3. Additional arguments to handle Green's function treatment :--
!         DO_CLASSICAL_SOLUTION, TAYLOR_ORDER, AVERAGE_SECANT, LC_AVERAGE_SECANT, NSTKS_NSTRMS, DELTAU_VERT, L_DELTAU_VERT
!         CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE

        CALL LC_BEAMSOLUTION_NEQK ( &
             FOURIER, IBEAM, N1, N_COLUMN_WFS, NSTOKES, NSTREAMS_2,         & ! Input numbers
             DO_LAYER_SCATTERING, BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, & ! Input av. secant
             BVEC, LC_T_DELT_MUBAR, LC_INITIAL_TRANS, LC_BVEC,              & ! Input linearized stuff
             L_WUPPER, L_WLOWER )                                             ! Output at level boundaries

!  .. 2 contributions to L_BEAM, from variations L_WUPPER L_WLOWER
!  .. 2 contributions to L_HOM,  from variations above and below

        DO I = 1, NSTREAMS_2
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          CM = C0 + IROW
          DO Q = 1, N_COLUMN_WFS

!  Beam contributions

            L_BEAM  = + L_WUPPER(I,O1,N1,Q) - L_WLOWER(I,O1,N,Q)

!  Linearized Real homogeneous solution contributions above

!  Enhancement # 4, 6/27/16. NOT IMPLEMENTED, LACKS CLARITY
!            L_HOM_R =sum( ( LCON(1:K_REAL(N),N) *  L_SOLA_XPOS(I,O1,1:K_REAL(N),N,Q) ) + &
!                ( MCON(1:K_REAL(N),N) * ( T_DELT_EIGEN(1:K_REAL(N),N)   * L_SOLB_XNEG(I,O1,1:K_REAL(N),N,Q) + &
!                 L_T_DELT_EIGEN(1:K_REAL(N),N,Q) *   SOLB_XNEG(I,O1,1:K_REAL(N),N) ) ) )
            L_HOM_U_R = ZERO
            DO K = 1, K_REAL(N1)
              CPOS = L_SOLA_XPOS(I,O1,K,N1,Q)
              CNEG = T_DELT_EIGEN(K,N1)   * L_SOLB_XNEG(I,O1,K,N1,Q) + &
                   L_T_DELT_EIGEN(K,N1,Q) *   SOLB_XNEG(I,O1,K,N1)
              T1 = LCON(K,N1) * CPOS
              T2 = MCON(K,N1) * CNEG
              L_HOM_U_R = L_HOM_U_R + T1 + T2
            ENDDO

!  Linearized Real homogeneous solution contributions below

            L_HOM_D_R = ZERO
            DO K = 1, K_REAL(N)
              CNEG = L_SOLB_XNEG(I,O1,K,N,Q)
              CPOS = T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q) + &
                   L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_D_R = L_HOM_D_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions above

            L_HOM_U_CR  = ZERO
            KO1 = K_REAL(N1) + 1
            DO K = 1, K_COMPLEX(N1)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N1,Q) * LCON(K1,N1) - &
                   L_SOLA_XPOS(I,O1,K2,N1,Q) * LCON(K2,N1)
              T2R =  T_DELT_EIGEN(K1,N1)   * L_SOLB_XNEG(I,O1,K1,N1,Q) &
                   - T_DELT_EIGEN(K2,N1)   * L_SOLB_XNEG(I,O1,K2,N1,Q) &
                   + L_T_DELT_EIGEN(K1,N1,Q) *   SOLB_XNEG(I,O1,K1,N1) &
                   - L_T_DELT_EIGEN(K2,N1,Q) *   SOLB_XNEG(I,O1,K2,N1)
              T2I =  T_DELT_EIGEN(K1,N1)   * L_SOLB_XNEG(I,O1,K2,N1,Q) &
                   + T_DELT_EIGEN(K2,N1)   * L_SOLB_XNEG(I,O1,K1,N1,Q) &
                   + L_T_DELT_EIGEN(K1,N1,Q) *   SOLB_XNEG(I,O1,K2,N1) &
                   + L_T_DELT_EIGEN(K2,N1,Q) *   SOLB_XNEG(I,O1,K1,N1)
              T2 =  T2R * MCON(K1,N1) - T2I * MCON(K2,N1)
              L_HOM_U_CR = L_HOM_U_CR + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions below

            L_HOM_D_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T2 = L_SOLB_XNEG(I,O1,K1,N,Q) * MCON(K1,N) - &
                   L_SOLB_XNEG(I,O1,K2,N,Q) * MCON(K2,N)
              T1R =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) &
                   - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K2,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N) &
                   - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
              T1I =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K2,N,Q) &
                   + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N) &
                   + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              L_HOM_D_CR = L_HOM_D_CR + T1 + T2
            ENDDO

!  Final contribution

            L_HOM_R  =  L_HOM_U_R  - L_HOM_D_R
            L_HOM_CR =  L_HOM_U_CR - L_HOM_D_CR
            COL2_WF(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR

          ENDDO
         ENDDO
        ENDDO

!  End layer

      ENDDO

!  FINAL layer
!  -----------

      N = NLAYERS
      MODIFIED_BOUNDARY = .TRUE.

!  get the linearized downward-reflected term
!    -- 1/31/21. Version 2.8.3. BRDF_F array defined locally for each Fourier.

      CALL L1_BVP_SURFACE_SETUP ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, MODIFIED_BOUNDARY, & ! flag    control
        FOURIER, IBEAM, N_COLUMN_WFS, NSTOKES, NSTREAMS, NLAYERS,     & ! integer Control
        QUAD_STRMWTS, NSTKS_NSTRMS, MUELLER_INDEX,          & ! Bookkeeping
        SURFACE_FACTOR, LAMBERTIAN_ALBEDO, BRDF_F,          & ! Surface
        T_DELT_EIGEN, K_REAL, K_COMPLEX, SOLA_XPOS,         & ! Input Homog
        L_T_DELT_EIGEN, L_SOLA_XPOS, L_SOLB_XNEG, L_WLOWER, & ! Input linearized Hom/PI
        R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )                     ! Output reflected solutions

!  offsets

      C0  = (N-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS

!  start loops

      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        I1 = I + NSTREAMS
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          CM = C0 + IROW
          DO Q = 1, N_COLUMN_WFS

!  Beam contributions

            L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)

!  Linearized Real homogeneous solution contributions

!  Enhancement # 5, 6/27/16. NOT IMPLEMENTED, LACKS CLARITY
!            L_HOM_R = sum( ( LCON(1:K_REAL(N),N) * ( ( T_DELT_EIGEN(1:K_REAL(N),N)   * L_SOLA_XPOS(I1,O1,1:K_REAL(N),N,Q) &
!                + L_T_DELT_EIGEN(1:K_REAL(N),N,Q) *   SOLA_XPOS(I1,O1,1:K_REAL(N),N) ) - R2_L_HOMP(I,O1,1:K_REAL(N),Q) ) ) + &
!                ( MCON(1:K_REAL(N),N) * ( L_SOLB_XNEG(I1,O1,1:K_REAL(N),N,Q) - R2_L_HOMM(I,O1,1:K_REAL(N),Q) ) ) )
            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) &
                  + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
              CPOS = CPOS - R2_L_HOMP(I,O1,K,Q)
              CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
              CNEG = CNEG - R2_L_HOMM(I,O1,K,Q)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions
!    Bug Fixed 16 December 2005.

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                  - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                  + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) &
                  - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
              T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                  + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                  + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N) &
                  + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
              T1R = T1R - R2_L_HOMP(I,O1,K1,Q)
              T1I = T1I - R2_L_HOMP(I,O1,K2,Q)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              T2R = L_SOLB_XNEG(I1,O1,K1,N,Q) - R2_L_HOMM(I,O1,K1,Q)
              T2I = L_SOLB_XNEG(I1,O1,K2,N,Q) - R2_L_HOMM(I,O1,K2,Q)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

!  Final contributions

            COL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

          ENDDO
        ENDDO
      ENDDO


!  Add direct beam variation to Final boundary
!  -------------------------------------------

!  Must use the reflected beam here
      
      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
         DO Q = 1, N_COLUMN_WFS
            DO I = 1, NSTREAMS
               IR = NSTOKES*(I-1) ; I1 = I + NSTREAMS
               DO O1 = 1, NSTOKES
                  IROW = IR + O1 ;  CM = C0 + IROW
                  FAC = - DIRECT_BEAM(I,IBEAM,O1)
                  L_BEAM = DOT_PRODUCT(L_DELTAU_VERT(Q,1:NLAYERS),DELTAU_SLANT(N,1:NLAYERS,IBEAM))
!                  COL2_WF(CM,Q) = COL2_WF(CM,Q) - FAC * L_BEAM   Bug 8/12/19
                  COL2_WF(CM,Q) = COL2_WF(CM,Q) + FAC * L_BEAM
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!  copy the single layer vector

      IF ( NLAYERS .EQ. 1 ) THEN
!  Enhancement # 6, 6/27/16
        SCOL2_WF(1:NTOTAL,1:N_COLUMN_WFS) = COL2_WF(1:NTOTAL,1:N_COLUMN_WFS)
      ENDIF

!  finish

      RETURN
      END SUBROUTINE MSSEG_LC_BVPCOLUMN_SETUP

!

      SUBROUTINE LC_BEAMSOLUTION_NEQK ( &
        FOURIER, IB, N, N_COLUMN_WFS, NSTOKES, NSTREAMS_2,             & ! Input numbers
        DO_LAYER_SCATTERING, BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, & ! Input av. secant
        BVEC, LC_T_DELT_MUBAR, LC_INITIAL_TRANS, LC_BVEC,              & ! Input linearized stuff
        L_WUPPER, L_WLOWER )                                             ! Output at level boundaries

!  Linearization of beam particular integral in one layer.

!  In this module, this is the Layer that contains the variation.
!                   N = LAYER_TO_VARY

!  Mk2.  6/21/21. Use module names now have the "_m" extension, Only use quantities you need
!  Mk2a. 9/21/21. New subroutine  by R. Spurr
!  Mk3.  1/31/22. Argument list re-organized for clarity, output status changed

      USE MSSeg_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXBEAMS, MAX_ATMOSWFS, &
                               MAXSTREAMS_2, MAXSTRMSTKS, MAXEVALUES, TAYLOR_SMALL, ONE

      IMPLICIT NONE

!  Input arguments
!  ---------------

!   Control

      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          IB
      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          N_COLUMN_WFS
      INTEGER, INTENT (IN) ::          NSTOKES

!   bookkeeping


      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  Average secant and linearization

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  PI solution and linearization

      DOUBLE PRECISION, INTENT (IN) :: BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LC_BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  outputs 

      DOUBLE PRECISION, INTENT (INOUT) :: L_WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      INTEGER          :: I, O1, Q
      DOUBLE PRECISION :: CONST, WDEL, VAR1, VAR2, VAR_U, TRANS2

!  Local linearized particular integrals at the layer boundaries

      DOUBLE PRECISION :: LBSOL_LOWER(MAXSTREAMS_2,MAXSTOKES,MAX_ATMOSWFS)
      DOUBLE PRECISION :: LBSOL_UPPER(MAXSTREAMS_2,MAXSTOKES,MAX_ATMOSWFS)

!  No linearized particular solution beyond the cutoff layer. ALSO--
!  Nothing if the solution saving mode is on and layer is inactive
!    Exit (Solutions have already been zeroed).

      IF (.NOT.DO_LAYER_SCATTERING(FOURIER,N) .OR. (N .GT.BEAM_CUTOFF(IB))) RETURN

!  Classical solution
!  ==================

!  Very simple, same code for all situations

      CONST   = INITIAL_TRANS(N,IB)
      WDEL    = T_DELT_MUBAR(N,IB)
      TRANS2  = CONST * WDEL
      DO Q = 1, N_COLUMN_WFS
         VAR1 = LC_T_DELT_MUBAR (N,IB,Q) * CONST
         VAR2 = LC_INITIAL_TRANS(N,IB,Q)
         DO I = 1, NSTREAMS_2
            DO O1 = 1, NSTOKES
               VAR_U = VAR2 * BVEC(I,O1,N) + LC_BVEC(I,O1,N,Q)
               LBSOL_UPPER(I,O1,Q) = CONST  * VAR_U
               LBSOL_LOWER(I,O1,Q) = TRANS2 * VAR_U + VAR1 * BVEC(I,O1,N)
            ENDDO
         ENDDO
      ENDDO

! Add to existing solution

      DO Q = 1, N_COLUMN_WFS
        DO I = 1, NSTREAMS_2
          DO O1 = 1, NSTOKES
            L_WUPPER(I,O1,N,Q) = L_WUPPER(I,O1,N,Q) + LBSOL_UPPER(I,O1,Q)
            L_WLOWER(I,O1,N,Q) = L_WLOWER(I,O1,N,Q) + LBSOL_LOWER(I,O1,Q)
          ENDDO
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LC_BEAMSOLUTION_NEQK

!

      SUBROUTINE L1_BVP_BACKSUB ( &
        N_COLUMN_WFS, NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        K_REAL, K_COMPLEX, BANDMAT2, IPIVOT, SMAT2, SIPIVOT, COL2_WF, SCOL2_WF, &
        NCON, PCON, STATUS, MESSAGE, TRACE )

!  Solves the linearized boundary value problem.

!  Mk2.  6/21/21. Use module names now have the "_m" extension, Only use quantities you need
!  Mk2a. 9/21/21. Subroutine copied by R. Spurr
!  Mk3.  1/31/22. Argument list re-organized for clarity, output status changed

      USE MSSeg_PARS_m, Only : MAXLAYERS, MAX_ATMOSWFS, MAXSTRMSTKS, MAXSTRMSTKS_2, &
                               MAXBANDTOTAL, MAXTOTAL, VQS_SUCCESS, VQS_SERIOUS

      USE LAPACK_TOOLS_m, Only : DGBTRS, DGETRS

      IMPLICIT NONE

!  Inputs
!  ------

! Basic control

      INTEGER, INTENT (IN) ::          N_COLUMN_WFS
      INTEGER, INTENT (IN) ::          NLAYERS

!  Bookkeeping control

      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2


      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

!  BVP problem (LU decomp)

      DOUBLE PRECISION, INTENT (IN) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::          IPIVOT ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::          SIPIVOT ( MAXSTRMSTKS_2 )

!  BVP linearized vector
  
      DOUBLE PRECISION, INTENT (INOUT) :: COL2_WF  ( MAXTOTAL, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: SCOL2_WF ( MAXSTRMSTKS_2, MAX_ATMOSWFS )

!  output (must be Intent(inout))

      DOUBLE PRECISION, INTENT (INOUT) ::  NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) ::  PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Exception handling

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  local variables
!  ---------------

      INTEGER ::           C0, LAY, K, Q, KO1, K0, K1, K2
      INTEGER ::           IROW, IROW1, IROW_S, IROW1_S
      INTEGER ::           INFO
      CHARACTER (LEN=3) :: CI

!  Intialize status

      STATUS = VQS_SUCCESS

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS &
           ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_COLUMN_WFS, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, &
              COL2_WF, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'Column Wfs:  DGBTRS call in L1_BVP_BACKSUB'
          STATUS  = VQS_SERIOUS
          RETURN
        ENDIF

!  Set Linearized integration constants NCON and PCON, all layers

        DO LAY = 1, NLAYERS
          C0 = (LAY-1)*NSTKS_NSTRMS_2
          KO1 = K_REAL(LAY) + 1
          DO K = 1, K_REAL(LAY)
            IROW = K
            IROW1 = IROW + NSTKS_NSTRMS
            DO Q = 1, N_COLUMN_WFS
              NCON(K,LAY,Q) = COL2_WF(C0+IROW,Q)
              PCON(K,LAY,Q) = COL2_WF(C0+IROW1,Q)
            ENDDO
          ENDDO
          DO K = 1, K_COMPLEX(LAY)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(LAY)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(LAY)
            IROW1_S = IROW_S + NSTKS_NSTRMS
            DO Q = 1, N_COLUMN_WFS
              NCON(K1,LAY,Q) = COL2_WF(C0+IROW,   Q)
              NCON(K2,LAY,Q) = COL2_WF(C0+IROW_S, Q)
              PCON(K1,LAY,Q) = COL2_WF(C0+IROW1,  Q)
              PCON(K2,LAY,Q) = COL2_WF(C0+IROW1_S,Q)
            ENDDO
          ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2

        CALL DGETRS &
           ( 'N', NTOTAL, N_COLUMN_WFS, SMAT2, MAXSTRMSTKS_2, SIPIVOT, &
              SCOL2_WF, MAXSTRMSTKS_2, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE = 'Atmos_Wfs for 1-layer: DGETRS call in L1_BVP_BACKSUB'
          STATUS = VQS_SERIOUS
          RETURN
        ENDIF

!  Set linearized integration constants NCON and PCON, 1 layer

        LAY = 1
        KO1 = K_REAL(LAY) + 1
        DO K = 1, K_REAL(LAY)
          IROW = K
          IROW1 = IROW + NSTKS_NSTRMS
          DO Q = 1, N_COLUMN_WFS
            NCON(K,LAY,Q) = SCOL2_WF(IROW, Q)
            PCON(K,LAY,Q) = SCOL2_WF(IROW1,Q)
          ENDDO
        ENDDO
        DO K = 1, K_COMPLEX(LAY)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW    = K + K_REAL(LAY)
          IROW1   = IROW + NSTKS_NSTRMS
          IROW_S  = K + K_REAL(LAY) + K_COMPLEX(LAY)
          IROW1_S = IROW_S + NSTKS_NSTRMS
          DO Q = 1, N_COLUMN_WFS
            NCON(K1,LAY,Q) = SCOL2_WF(IROW,    Q)
            NCON(K2,LAY,Q) = SCOL2_WF(IROW_S,  Q)
            PCON(K1,LAY,Q) = SCOL2_WF(IROW1,   Q)
            PCON(K2,LAY,Q) = SCOL2_WF(IROW1_S, Q)
          ENDDO
        ENDDO

      ENDIF

!  finish

      RETURN
      END SUBROUTINE L1_BVP_BACKSUB

!

      SUBROUTINE L1_BVP_SURFACE_SETUP ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, MODIFIED_BCL4, & ! flag    control
        FOURIER, IBEAM, N_COLUMN_WFS, NSTOKES, NSTREAMS, NLAYERS, & ! integer Control
        QUAD_STRMWTS, NSTKS_NSTRMS, MUELLER_INDEX,          & ! Bookkeeping
        SURFACE_FACTOR, LAMBERTIAN_ALBEDO, BRDF_F,          & ! Surface
        T_DELT_EIGEN, K_REAL, K_COMPLEX, SOLA_XPOS,         & ! Input Homog
        L_T_DELT_EIGEN, L_SOLA_XPOS, L_SOLB_XNEG, L_WLOWER, & ! Input linearized Hom/PI
        R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )                     ! Output reflected solutions

!  Mk2.  6/21/21. Use module names now have the "_m" extension, Only use quantities you need
!  Mk2a. 9/21/21. Subroutine copied by R. Spurr
!  Mk3.  1/31/22. Argument list re-organized for clarity.

      USE MSSeg_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTREAMS, MAX_ATMOSWFS, &
                               MAXSTREAMS_2, MAXEVALUES, MAXSTOKES_SQ, ZERO

      IMPLICIT NONE

!  inputs
!  ------

!  Flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::          MODIFIED_BCL4

!  Control numbers

      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          IBEAM
      INTEGER, INTENT (IN) ::          N_COLUMN_WFS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS

!  bookkeeping

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

!  surface. Mk3 4/11/23. Drop Fourier index from BRDF_F

      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )

!  homogeneous solutions

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Linearized homogeneous solutions

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized solution at lower boundary

      DOUBLE PRECISION, INTENT (INOUT) :: L_WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  output surface contributions linearized

      DOUBLE PRECISION, INTENT (OUT) :: R2_L_BEAM ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: R2_L_HOMP ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: R2_L_HOMM ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

!  Local variables
!  ---------------

!  storage arrays

      DOUBLE PRECISION :: PV_W ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: HV_P ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: HV_M ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

!  Help

      INTEGER ::          I, J, O1, O2, OM, IB, Q, N, M, NELEMENTS
      INTEGER ::          K, KO1, K0, K1, K2

      DOUBLE PRECISION :: H1R, H1I, REFL_P, REFL_B, REFL_M, KMULT
      DOUBLE PRECISION :: H1, H2, H1_CR,   H2_CR,   H1_CI,   H2_CI
      DOUBLE PRECISION :: H1_S_CR, H2_S_CR, H1_S_CI, H2_S_CI, HM, HP, BIREF

!  Initial section
!  ---------------

!  Beam index

      IB = IBEAM

!  Fourier component

      M = FOURIER

!  Always zero the result to start

      DO I = 1, NSTREAMS
       DO O1 = 1, NSTOKES
        DO Q = 1, N_COLUMN_WFS
          R2_L_BEAM(I,O1,Q)    = ZERO
          DO K = 1, NSTKS_NSTRMS
            R2_L_HOMP(I,O1,K,Q) = ZERO
            R2_L_HOMM(I,O1,K,Q) = ZERO
          ENDDO
        ENDDO
       ENDDO
      ENDDO

!  Return if no albedo

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Return if Fourier component > 0 (Lambertian)

      IF ( DO_LAMBERTIAN_SURFACE .and. M.gt.0 ) RETURN

!  Set up Auxiliary arrays
!  -----------------------

!  Only require (1,1) component

      O1 = 1

!  Last layer

      N = NLAYERS

!  number of elements

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        NELEMENTS = 1
      ELSE
        NELEMENTS = NSTOKES
      ENDIF

!  Need Beam arrays regardless
!   Always something from layers above the PBL and also from PBL

      DO J = 1, NSTREAMS
        DO O1 = 1, NSTOKES
          DO Q = 1, N_COLUMN_WFS
            PV_W(J,O1,Q) = L_WLOWER(J,O1,N,Q) * QUAD_STRMWTS(J)
          ENDDO
        ENDDO
      ENDDO

!   Modified boundary condition
!     This applies when PBL is varying layer.
!       Require homogeneous solutions linearizations

      IF ( MODIFIED_BCL4 ) THEN

!  start loops

        DO J = 1, NSTREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, N_COLUMN_WFS

!  real homogeneous solution contributions

           DO K = 1, K_REAL(N)
            H1 = L_SOLA_XPOS(J,O1,K,N,Q) *   T_DELT_EIGEN(K,N) + &
                   SOLA_XPOS(J,O1,K,N)   * L_T_DELT_EIGEN(K,N,Q)
            H2 = L_SOLB_XNEG(J,O1,K,N,Q)
            HV_P(J,O1,K,Q) = QUAD_STRMWTS(J)*H1
            HV_M(J,O1,K,Q) = QUAD_STRMWTS(J)*H2
           ENDDO

!  Complex homogeneous solution contributions

           KO1 = K_REAL(N) + 1
           DO K = 1, K_COMPLEX(N)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            H1R = L_SOLA_XPOS(J,O1,K1,N,Q) *   T_DELT_EIGEN(K1,N) &
                - L_SOLA_XPOS(J,O1,K2,N,Q) *   T_DELT_EIGEN(K2,N) &
                +   SOLA_XPOS(J,O1,K1,N)   * L_T_DELT_EIGEN(K1,N,Q) &
                -   SOLA_XPOS(J,O1,K2,N)   * L_T_DELT_EIGEN(K2,N,Q)
            H1I = L_SOLA_XPOS(J,O1,K1,N,Q) *   T_DELT_EIGEN(K2,N) &
                + L_SOLA_XPOS(J,O1,K2,N,Q) *   T_DELT_EIGEN(K1,N) &
                +   SOLA_XPOS(J,O1,K1,N)   * L_T_DELT_EIGEN(K2,N,Q) &
                +   SOLA_XPOS(J,O1,K2,N)   * L_T_DELT_EIGEN(K1,N,Q)
            HV_P(J,O1,K1,Q) = QUAD_STRMWTS(J)* H1R
            HV_P(J,O1,K2,Q) = QUAD_STRMWTS(J)* H1I
            HV_M(J,O1,K1,Q) = QUAD_STRMWTS(J)* L_SOLB_XNEG(J,O1,K1,N,Q)
            HV_M(J,O1,K2,Q) = QUAD_STRMWTS(J)* L_SOLB_XNEG(J,O1,K2,N,Q)
           ENDDO

!  End loops

          ENDDO
         ENDDO
        ENDDO

!  End modified BCL4 condition

      ENDIF

!  Lambertian condition
!  ====================

!  Skip if BRDF surface

      if ( .not. DO_LAMBERTIAN_SURFACE ) go to 988

!  reflection

      KMULT = SURFACE_FACTOR * LAMBERTIAN_ALBEDO

!  only 1 component

      O1 = 1

!  Integrated Downward reflection (Calculation)
!  --------------------------------------------

      DO Q = 1, N_COLUMN_WFS

!  Particular solution (only for the first Stokes component)

       REFL_B = ZERO
       DO J = 1, NSTREAMS
         REFL_B = REFL_B + PV_W(J,O1,Q)
       ENDDO
       REFL_B = REFL_B * KMULT
       DO I = 1, NSTREAMS
         R2_L_BEAM(I,O1,Q) = REFL_B
       ENDDO

!  Homogeneous solutions for the modified condition

       IF ( MODIFIED_BCL4 ) THEN

!  Homogeneous real solutions

         DO K = 1, K_REAL(NLAYERS)
          REFL_P = ZERO
          REFL_M = ZERO
          DO J = 1, NSTREAMS
           REFL_P = REFL_P + HV_P(J,O1,K,Q)
           REFL_M = REFL_M + HV_M(J,O1,K,Q)
          ENDDO
          REFL_P = REFL_P * KMULT
          REFL_M = REFL_M * KMULT
          DO I = 1, NSTREAMS
           R2_L_HOMP(I,O1,K,Q) = REFL_P
           R2_L_HOMM(I,O1,K,Q) = REFL_M
          ENDDO
         ENDDO

!  Homogeneous complex solutions

         KO1 = K_REAL(NLAYERS) + 1
         DO K = 1, K_COMPLEX(NLAYERS)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          H1_CR = ZERO
          H1_CI = ZERO
          H2_CR = ZERO
          H2_CI = ZERO
          DO J = 1, NSTREAMS
           H1_CR = H1_CR + HV_P(J,O1,K1,Q)
           H1_CI = H1_CI + HV_P(J,O1,K2,Q)
           H2_CR = H2_CR + HV_M(J,O1,K1,Q)
           H2_CI = H2_CI + HV_M(J,O1,K2,Q)
          ENDDO
          H1_CR = H1_CR * KMULT
          H1_CI = H1_CI * KMULT
          H2_CR = H2_CR * KMULT
          H2_CI = H2_CI * KMULT
          DO I = 1, NSTREAMS
            R2_L_HOMP(I,O1,K1,Q) = H1_CR
            R2_L_HOMP(I,O1,K2,Q) = H1_CI
            R2_L_HOMM(I,O1,K1,Q) = H2_CR
            R2_L_HOMM(I,O1,K2,Q) = H2_CI
          ENDDO
         ENDDO

!  End modified boundary condition clause

        ENDIF

!  end parameter loop

      ENDDO

!  Complete Lambertian case

      RETURN

!  BRDF surface condition
!  ======================

!  Continuation point

988   continue

!  start loops

      DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES
          DO Q = 1, N_COLUMN_WFS

!  Particular solution

            REFL_B = ZERO
            DO J = 1, NSTREAMS
              H1 = ZERO
              DO O2 = 1, NSTOKES
               OM = MUELLER_INDEX(O1,O2)
               H1 = H1 + PV_W(J,O2,Q)*BRDF_F(OM,J,I)
              ENDDO
              REFL_B = REFL_B + H1
            ENDDO
            REFL_B = REFL_B * SURFACE_FACTOR
            R2_L_BEAM(I,O1,Q) = R2_L_BEAM(I,O1,Q) + REFL_B

!  Homogeneous solutions for the modified condition

            IF ( MODIFIED_BCL4 ) THEN

!  Homogeneous real solutions

              DO K = 1, K_REAL(NLAYERS)
               REFL_P = ZERO
               REFL_M = ZERO
               DO J = 1, NSTREAMS
                HP = ZERO
                HM = ZERO
                DO O2 = 1, NSTOKES
                 OM = MUELLER_INDEX(O1,O2)
                 HP = HP + HV_P(J,O2,K,Q)*BRDF_F(OM,J,I)
                 HM = HM + HV_M(J,O2,K,Q)*BRDF_F(OM,J,I)
                ENDDO
                REFL_P = REFL_P + HP
                REFL_M = REFL_M + HM
               ENDDO
               R2_L_HOMP(I,O1,K,Q) =  REFL_P * SURFACE_FACTOR
               R2_L_HOMM(I,O1,K,Q) =  REFL_M * SURFACE_FACTOR
              ENDDO

!  homogeneous complex solutions

              KO1 = K_REAL(NLAYERS) + 1
              DO K = 1, K_COMPLEX(NLAYERS)
               K0 = 2*K - 2
               K1 = KO1 + K0
               K2 = K1  + 1
               H1_CR = ZERO
               H2_CR = ZERO
               H1_CI = ZERO
               H2_CI = ZERO
               DO J = 1, NSTREAMS
                H1_S_CR = ZERO
                H2_S_CR = ZERO
                H1_S_CI = ZERO
                H2_S_CI = ZERO
                DO O2 = 1, NSTOKES
                 OM    = MUELLER_INDEX(O1,O2)
                 BIREF = BRDF_F(OM,J,I)
                 H1_S_CR = H1_S_CR + HV_P(J,O2,K1,Q)*BIREF
                 H2_S_CR = H2_S_CR + HV_M(J,O2,K1,Q)*BIREF
                 H1_S_CI = H1_S_CI + HV_P(J,O2,K2,Q)*BIREF
                 H2_S_CI = H2_S_CI + HV_M(J,O2,K2,Q)*BIREF
                ENDDO
                H1_CR = H1_CR + H1_S_CR
                H2_CR = H2_CR + H2_S_CR
                H1_CI = H1_CI + H1_S_CI
                H2_CI = H2_CI + H2_S_CI
               ENDDO
               R2_L_HOMP(I,O1,K1,Q) = H1_CR * SURFACE_FACTOR
               R2_L_HOMM(I,O1,K1,Q) = H2_CR * SURFACE_FACTOR
               R2_L_HOMP(I,O1,K2,Q) = H1_CI * SURFACE_FACTOR
               R2_L_HOMM(I,O1,K2,Q) = H2_CI * SURFACE_FACTOR
              ENDDO

!  End modified condition

            ENDIF

!  end parameter and stream and Stokes loops

          ENDDO
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE L1_BVP_SURFACE_SETUP

!  End module

end module MSSeg_LC_BVProblem_m

