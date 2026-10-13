
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

!  Mk2a. 9/21/21. LC Subroutines created
!  Mk3.  1/31/22. Homogeneous setups only, LC stuff hived off elsewhere.
!  Mk3.  4/11/23. Modified for MK3 use

      MODULE MSSeg_LA_MiscSetups_m

      PRIVATE :: MSSEG_LA_DELTAMSCALE, MSSEG_LA_SSALBINIT
      PUBLIC  :: MSSEG_LA_MISCSETUPS

      CONTAINS

!

      SUBROUTINE MSSEG_LA_MISCSETUPS ( &
        DO_DELTAM_SCALING, NSTOKES, NLAYERS, NMOMENTS,        & ! Input
        MUELLER_INDEX, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,    & ! Input
        OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, OMEGA_GREEK, & ! Input Optical/Deltam
        TRUNC_FACTOR, FAC1, L_DELTAU_VERT_INPUT,              & ! Input Optical/Deltam
        L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,          & ! Input Linearized Optical
        L_OMEGA_TOTAL, L_DELTAU_VERT, L_GREEKMAT_TOTAL,       & ! Output
        DO_SCATMAT_VARIATION, L_TRUNC_FACTOR, L_OMEGA_GREEK )   ! Output

!  Mk2a, 9/21/21. LLAC top-level subroutine created by R. Spurr
!   -- Modeled after LAP subroutine of the same name

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- Remove SFI_OPDEPS array from argument list.

      USE MSSeg_PARS_m, Only : MAXSTOKES, MAXLAYERS, MAXSTREAMS, MAX_ATMOSWFS, MAX_GEOMETRIES, &
                               MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, ZERO

      IMPLICIT NONE

!  Inputs
!  ------

!  control

      LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NMOMENTS

!  bookkeeping

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG   ( MAXLAYERS )

!  optical

      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

!  Linearized optical

      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_TOTAL_INPUT    ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT_INPUT    ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_GREEKMAT_TOTAL_INPUT ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

!  Input from original DELTAMSCALE routine

      DOUBLE PRECISION, INTENT (IN) :: TRUNC_FACTOR ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: FAC1 ( MAXLAYERS )

!  Input from origianl SSALBINIT

      DOUBLE PRECISION, INTENT (IN) :: OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  output
!  ------

!  from LA_DELTAMSCALE

      DOUBLE PRECISION, INTENT (OUT) :: L_OMEGA_TOTAL    ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: L_DELTAU_VERT    ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: L_GREEKMAT_TOTAL ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      LOGICAL, INTENT (OUT)          ::  DO_SCATMAT_VARIATION ( MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_TRUNC_FACTOR ( MAX_ATMOSWFS, MAXLAYERS )

!  from MSSEG_LA_SSALBINIT

      DOUBLE PRECISION, INTENT (OUT) :: L_OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  miscellaneous setup operations, Master routine

!  Delta-m scaling of input quantities
!   Mk2, 6/21/21. remove SFI_OPDEPS array from argument list

      CALL MSSEG_LA_DELTAMSCALE ( &
        DO_DELTAM_SCALING, NSTOKES, NLAYERS, NMOMENTS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, TRUNC_FACTOR, FAC1,                       &
        L_OMEGA_TOTAL_INPUT, L_DELTAU_VERT_INPUT, L_GREEKMAT_TOTAL_INPUT,                  &
        L_OMEGA_TOTAL, L_DELTAU_VERT, L_GREEKMAT_TOTAL, DO_SCATMAT_VARIATION, L_TRUNC_FACTOR )

!  initialise single scatter albedo terms
!    GREEKMAT x OMEGA

      CALL MSSEG_LA_SSALBINIT ( &
        NSTOKES, NLAYERS, NMOMENTS, MUELLER_INDEX,       &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, OMEGA_GREEK, &
        L_OMEGA_TOTAL, L_GREEKMAT_TOTAL, L_OMEGA_GREEK )

!  Finish

      RETURN
      END SUBROUTINE MSSEG_LA_MISCSETUPS

!

      SUBROUTINE MSSEG_LA_DELTAMSCALE ( &
        DO_DELTAM_SCALING, NSTOKES, NLAYERS, NMOMENTS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, TRUNC_FACTOR, FAC1,                       &
        L_OMEGA_TOTAL_INPUT, L_DELTAU_VERT_INPUT, L_GREEKMAT_TOTAL_INPUT,                  &
        L_OMEGA_TOTAL, L_DELTAU_VERT, L_GREEKMAT_TOTAL, DO_SCATMAT_VARIATION, L_TRUNC_FACTOR )

!  module, dimensions and numbers

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- remove SFI_OPDEPS array from argument list

      USE MSSeg_PARS_m, Only : MAX_ATMOSWFS, MAXMOMENTS, MAXLAYERS, &
                               MAX_GEOMETRIES, MAXBEAMS, MAXSTOKES_SQ, zero, one, SMALLNUM

      IMPLICIT NONE

!  Inputs
!  ------

!  control

      LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NMOMENTS

!  bookkeeping

      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG   ( MAXLAYERS )

!  optical

      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

!  Linearized optical

      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_TOTAL_INPUT    ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT_INPUT    ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_GREEKMAT_TOTAL_INPUT ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

!  Input from original DELTAMSCALE routine

      DOUBLE PRECISION, INTENT (IN) :: TRUNC_FACTOR ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: FAC1 ( MAXLAYERS )

!  output
!  ------

      DOUBLE PRECISION, INTENT (OUT) :: L_OMEGA_TOTAL    ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: L_DELTAU_VERT    ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: L_GREEKMAT_TOTAL ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      LOGICAL, INTENT (OUT)          ::  DO_SCATMAT_VARIATION ( MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_TRUNC_FACTOR ( MAX_ATMOSWFS, MAXLAYERS )

!  local variables

      DOUBLE PRECISION :: BLD, DL, OF1, F, F1, UQ, EQ, T1
      DOUBLE PRECISION :: FZM, ZMQ, ZLQ, UZQ_SCALE, ZQ_SCALE, UQ_SCALE
      INTEGER          :: N, Q, NM1, L, K, K1, K2
      LOGICAL          :: LOOP

!  Indexing

      INTEGER :: KTYPE1(4), KTYPE2(4)

      KTYPE1 = (/ 1, 6, 11, 16 /)
      KTYPE2 = (/ 2, 5, 12, 15 /)

!  Set Scattering matrix linearization flag
!   Only examine the (1,1) entry (the phase function coefficient)
!    Dimensioning bug, 21 March 2007. Care with using NM1

      K1 = 1
      IF ( DO_DELTAM_SCALING ) THEN
        NM1 = NMOMENTS+1
      ELSE
        NM1 = NMOMENTS
      ENDIF

      DO N = 1, NLAYERS
         IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LOOP = .TRUE. ; L = 0
              DO WHILE (LOOP.AND.L.LT.NM1)
                L = L + 1
                LOOP = ( DABS(L_GREEKMAT_TOTAL_INPUT(Q,L,N,K1)) .LT.1000.0*SMALLNUM )
              ENDDO
              DO_SCATMAT_VARIATION(N,Q) = .NOT.LOOP
            ENDDO
         ENDIF
      ENDDO

!  DELTAM SCALING
!  ==============

!mick fix 2/17/11 - initialize
      L_GREEKMAT_TOTAL = ZERO

!  New section added 21 December 2005 by R. Spurr
!   Linearization of the scaled delta-M inputs

      IF ( DO_DELTAM_SCALING ) THEN

         NM1 = NMOMENTS+1

!  start loop over layers

         DO N = 1, NLAYERS

!  If the layer is varying

            IF ( LAYER_VARY_FLAG(N) ) THEN

              OF1 = ( ONE - FAC1(N) ) / FAC1(N)
              F   = TRUNC_FACTOR(N)
              F1  = ONE - F

!  Start loop over varying parameters in this layer

              DO Q = 1, LAYER_VARY_NUMBER(N)

!  scale scattering matrix linearization additionally

                IF ( DO_SCATMAT_VARIATION(N,Q) ) THEN

!  set bulk property values

                  UQ  = L_OMEGA_TOTAL_INPUT(Q,N)
                  EQ  = L_DELTAU_VERT_INPUT(Q,N)
                  ZMQ = L_GREEKMAT_TOTAL_INPUT(Q,NM1,N,1)
                  FZM = F * ZMQ
                  UZQ_SCALE = ( UQ + ZMQ ) * OF1
                  ZQ_SCALE  = FZM / F1
                  L_TRUNC_FACTOR(Q,N) = FZM
                  L_OMEGA_TOTAL(Q,N)      = UQ + UZQ_SCALE - ZQ_SCALE
                  L_DELTAU_VERT(Q,N)      = EQ - UZQ_SCALE

!  do the phase function first
!Rob Fix. 1/17/14. Cannot have BLD anf F both zero.

                  K1 = 1
                  L_GREEKMAT_TOTAL(Q,0,N,K1) = ZERO
                  FZM = F * L_GREEKMAT_TOTAL_INPUT(Q,NM1,N,1)
                  DO L = 1, NMOMENTS
                    DL  = DBLE(2*L+1)
                    ZLQ = L_GREEKMAT_TOTAL_INPUT(Q,L,N,K1)
                    BLD = GREEKMAT_TOTAL_INPUT(L,N,K1) / DL
                    IF ( BLD.ne.ZERO.or.F.ne.zero ) THEN
                      T1  = ( BLD*ZLQ - FZM ) / ( BLD - F )
                      L_GREEKMAT_TOTAL(Q,L,N,K1) = T1 + ZQ_SCALE
                    ENDIF
                  ENDDO

!  If there is polarization.........

                  IF ( NSTOKES .GT. 1 ) THEN

!  Now do the other Type 1 elements
!Rob Fix. 1/17/14. Cannot have BLD anf F both zero.

                    DO L = 0, NMOMENTS
                      DL  = DBLE(2*L+1)
                      DO K = 2, 4
                        K1 = KTYPE1(K)
                        ZLQ = L_GREEKMAT_TOTAL_INPUT(Q,L,N,K1)
                        BLD = GREEKMAT_TOTAL_INPUT(L,N,K1) / DL
                        IF ( BLD.ne.ZERO.or.F.ne.zero ) THEN
                          T1  = ( BLD*ZLQ - FZM ) / ( BLD - F )
                          L_GREEKMAT_TOTAL(Q,L,N,K1) = T1 + ZQ_SCALE
                        ENDIF
                      ENDDO
                    ENDDO

!  Now do the Type 2 elements
!    Bug discovered 24 September 2007
!    Formerly, we were overwriting the Type 1 elements

                    DO L = 0, NMOMENTS
                      DL  = DBLE(2*L+1)
                      DO K = 1, 4
                        K2 = KTYPE2(K)
                        ZLQ = L_GREEKMAT_TOTAL_INPUT(Q,L,N,K2)
!  !!!! WRONG           L_GREEKMAT_TOTAL(Q,L,N,K1) = ZLQ + ZQ_SCALE
                        L_GREEKMAT_TOTAL(Q,L,N,K2) = ZLQ + ZQ_SCALE
                      ENDDO
                    ENDDO

!  Finish scaled Greek matrix linearization

                  ENDIF

!  No scattering matrix linearization

                ELSE

!  Bulk property linearization

                  UQ = L_OMEGA_TOTAL_INPUT(Q,N)
                  EQ = L_DELTAU_VERT_INPUT(Q,N)
                  L_TRUNC_FACTOR(Q,N) = ZERO
                  UQ_SCALE = UQ * OF1
                  L_OMEGA_TOTAL(Q,N) = UQ + UQ_SCALE
                  L_DELTAU_VERT(Q,N) = EQ - UQ_SCALE

!   Zero all linearized scattering matrix quantities now;
!     Not necessary if pre-zeroed
!                  DO L = 0, NMOMENTS
!                   DO O1 = 1, NSTOKES_SQ
!                    L_GREEKMAT_TOTAL(Q,L,N,O1) = ZERO
!                   ENDDO
!                  ENDDO

                ENDIF

              ENDDO

!  end layer loop

            ENDIF
         ENDDO

!  NO DELTAM SCALING
!  =================

!  move input geophysical variables to Workspace quantities

      ELSE

         DO N = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(N) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                L_TRUNC_FACTOR(Q,N) = ZERO
                L_DELTAU_VERT(Q,N)  = L_DELTAU_VERT_INPUT(Q,N)
                L_OMEGA_TOTAL(Q,N)  = L_OMEGA_TOTAL_INPUT(Q,N)
                IF ( DO_SCATMAT_VARIATION(N,Q) ) THEN
                  DO L = 0, MAXMOMENTS
                    L_GREEKMAT_TOTAL(Q,L,N,:) = L_GREEKMAT_TOTAL_INPUT(Q,L,N,:)
                  ENDDO
                ELSE
                  DO L = 0, MAXMOMENTS
                    L_GREEKMAT_TOTAL(Q,L,N,:) = ZERO
                  ENDDO
                ENDIF
              ENDDO
            ENDIF
         ENDDO

      ENDIF

!  Finish module

      RETURN
      END SUBROUTINE MSSEG_LA_DELTAMSCALE

!

      SUBROUTINE MSSEG_LA_SSALBINIT ( &
        NSTOKES, NLAYERS, NMOMENTS, MUELLER_INDEX,       &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, OMEGA_GREEK, &
        L_OMEGA_TOTAL, L_GREEKMAT_TOTAL, L_OMEGA_GREEK )

!  Mk2, 6/21/21. Use module names now have the "_m" extension, Only use quantities you need
!  Mk3, 1/31/22. Not changed

      USE MSSeg_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES_SQ, MAX_ATMOSWFS

      IMPLICIT NONE

! inputs
! ------

!  Control

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NMOMENTS

!  bookkeeping

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )

!  optical

      DOUBLE PRECISION, INTENT (IN) :: OMEGA_GREEK      ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_TOTAL    ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_GREEKMAT_TOTAL ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

!  output
!  ------

      DOUBLE PRECISION, INTENT (OUT) :: L_OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          N, L, Q, O1, O2, OM
      DOUBLE PRECISION :: VL

!  phase moment-weighted OMEGA and linearizations
!  Including phase function linearization

      DO N = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(N) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(N)
           DO L = 0, NMOMENTS
            DO O1 = 1, NSTOKES
             DO O2 = 1, NSTOKES
              OM = MUELLER_INDEX(O1,O2)
              VL = L_OMEGA_TOTAL   (Q,N)   + &
                   L_GREEKMAT_TOTAL(Q,L,N,OM)
              L_OMEGA_GREEK(L,N,O1,O2,Q) = OMEGA_GREEK(L,N,O1,O2) * VL
             ENDDO
            ENDDO
           ENDDO
          ENDDO
        ENDIF
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE MSSEG_LA_SSALBINIT

!  Finish Module

END MODULE MSSeg_LA_MiscSetups_m

